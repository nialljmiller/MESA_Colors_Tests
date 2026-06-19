#!/usr/bin/env python3
"""
Fine-tune a MESA RSP RR Lyrae model to the Kepler light curve of FN Lyr
(KIC 6936115) with a Nelder-Mead simplex.

This is stage 2 of a two-stage fit. Stage 1 was a coarse grid search that
located the basin (T~6700, L~47.7, kick=40, alfam~0.15, M=0.65) matching
FN Lyr's period (~0.527 d) and amplitude (~1.05 mag) but leaving a
shock-split maximum. Stage 2 refines four parameters automatically:

    RSP_Teff   strip position
    RSP_L      luminosity
    RSP_alfam  eddy viscosity (amplitude lever)
    RSP_mass   the lever that decouples amplitude from the surface shock

against the OBSERVED Kepler waveform, which is the quantity actually being
judged. Period and amplitude are held near their observed values by soft
penalties; the parameters are held in physical RRab ranges by hard bounds.

Objective for each evaluation:
    1. write inlist, run ./rn
    2. if the run crashed / wrote no history  -> large finite penalty
    3. if the limit cycle has NOT settled       -> graded penalty
       (amp_last10_range above threshold; scoring a transient is invalid)
    4. else: fold the best converged cycle, fit phase shift + vertical
       offset against the binned FN Lyr curve, return RMS residual
       + soft period penalty + soft amplitude guard.

Every evaluation is appended to optimize_log.csv -- that log is the paper
artifact showing Colors driving convergence onto a real dataset. The best
run's history.data and inlist are copied to optimize_best/.

Run from:
    ~/MESA/MESA_Colors_Tests/rsp_RR_Lyrae/FN_Lyr

Usage:
    # fast / rough (write the section today, rerun longer afterwards)
    python optimize_fn_lyr.py --mode fast | tee optimize_fn_lyr_fast.log

    # production (overnight; final numbers)
    python optimize_fn_lyr.py --mode production | tee optimize_fn_lyr_prod.log

    # resume: re-seed the simplex from the best rows already logged
    python optimize_fn_lyr.py --mode production --resume | tee -a optimize_fn_lyr_prod.log
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize


# ----------------------------------------------------------------------------
# Observed FN Lyr targets
# ----------------------------------------------------------------------------
OBS_PERIOD = 0.527398471          # days, Nemec et al. 2011
OBS_AMP = 1.0464                  # P99-P01 Kepler amplitude (mag)
OBS_FILE = "FN_Lyr_KIC6936115_kepler_binned.csv"

INLIST = Path("inlist_rsp_RR_Lyrae")
RUN_CMD = ["./rn"]

LOG_CSV = Path("optimize_log.csv")
BEST_DIR = Path("optimize_best")
BEST_DIR.mkdir(exist_ok=True)


# ----------------------------------------------------------------------------
# Fixed (non-optimized) parameters
# ----------------------------------------------------------------------------
# kick is a seed, not a fit parameter -- fixed at the value that settled
# fastest in stage 1. alfat/gammar were swept and rejected (alfat does nothing
# to the waveform; gammar amplifies and destabilizes), so set_inlist forces
# BOTH to 0 every eval. This matters: the working inlist may carry a leftover
# RSP_gammar (e.g. 2d0 from the last gammar trial), and pinning it here
# overwrites that so every eval runs the radiative-set baseline, not the
# worst-shock configuration.
FIXED_KICK = 40.0
FIXED_ALFAT = 0.0
FIXED_GAMMAR = 0.0
FIXED_Z = 0.0004


# ----------------------------------------------------------------------------
# Framing B: the genuinely-independent observables are PERIOD (a dynamical
# quantity, ~exact) and metallicity (from the Fourier phi31 relation, largely
# independent of the amplitude/shape we fit). These anchor the fit. The
# published M, L, Teff (Nemec et al. 2011, Tables 4 & 6) are derived from the
# same light curve via empirical Fourier calibrations, so we do NOT fix them
# as inputs -- we let them float within their measured ranges and report where
# the converged fit lands relative to the Nemec central values. The claim is
# then "fitting the Kepler waveform recovers M, L, Teff consistent with the
# independent Fourier-decomposition estimates", which avoids circularity.
#
# Nemec et al. (2011) FN Lyr central values (for distance-from-literature
# reporting; NOT used as constraints):
NEMEC = {
    "teff": 6520.0,   # KW01/S06, Kp & V span ~6482-6559; midpoint
    "lum":  48.0,     # L(puls) ~40-43, L(evol) ~51-56; midpoint of the union
    "mass": 0.64,     # M(puls)=0.595, M(evol)=0.69; midpoint of the disagreement
}


@dataclass
class Param:
    name: str
    lo: float
    hi: float
    seed: float


# Bounds are the Nemec measured ranges (mass spans puls->evol, Teff spans the
# KW01/S06 Kp&V estimates, L spans puls->evol), each padded by a hair so the
# simplex can sit fully inside the measured box rather than on its edge. alfam
# is the one genuinely free convective parameter (not measured) and carries a
# physical RRab range. The seed is placed INSIDE the measured box -- NOT at the
# old grid basin (T6700), which was hotter and heavier than FN Lyr actually is.
PARAMS = [
    Param("teff",  6470.0, 6570.0, 6520.0),   # measured 6482-6559
    Param("lum",     40.0,   56.0,   48.0),    # L(puls)->L(evol)
    Param("alfam",   0.08,    0.30,   0.15),   # free convective lever
    Param("mass",    0.59,    0.69,   0.62),    # M(puls)=0.595 -> M(evol)=0.69
]


# ----------------------------------------------------------------------------
# Objective weights / gates
# ----------------------------------------------------------------------------
# RMS waveform residual is ~0.06 for a good fit. Weights are chosen so the
# soft penalties are comparable to RMS only when a target is meaningfully
# violated, not for small deviations.
PERIOD_TOL = 0.015        # fractional period error tolerated with no penalty
PERIOD_W = 1.0            # penalty per unit fractional period error beyond tol
AMP_TOL = 0.10            # fractional amplitude deviation tolerated (RMS already
                          #   captures amplitude; this is just a guard rail)
AMP_W = 0.5

# Convergence gate. Stage-1 settled models had amp_last10_range ~0.001;
# de-converged runs reached ~0.01+. Reject anything not settled to here.
CONV_RANGE_MAX = 0.02

# Penalties (finite, so Nelder-Mead keeps moving rather than seeing inf).
PENALTY_CRASH = 5.0       # no history / nonzero exit
PENALTY_UNCONV_BASE = 2.0 # unsettled run; graded by how unsettled


# ----------------------------------------------------------------------------
# Mode configuration
# ----------------------------------------------------------------------------
MODE_CFG = {
    # FAST: ~120 periods still mostly settles near this basin; capped evals
    # and a tight initial simplex keep wall time to roughly half a day.
    # The amplitude will read slightly low from truncation -- expected.
    "fast": dict(periods=120, max_evals=30, simplex_frac=0.06, max_model_number=240000),
    # PRODUCTION: full settling, larger eval budget. Overnight+.
    "production": dict(periods=200, max_evals=80, simplex_frac=0.10, max_model_number=320000),
}


# ============================================================================
# Inlist editing (same discipline as the grid driver: full-line anchors only)
# ============================================================================
def replace_or_insert(text: str, pattern: str, replacement: str, anchor: str | None = None) -> str:
    new, n = re.subn(pattern, replacement, text, flags=re.MULTILINE)
    if n > 0:
        return new
    if anchor is not None and anchor in text:
        return text.replace(anchor, anchor + "\n   " + replacement)
    return text + "\n   " + replacement + "\n"


def set_inlist(teff: float, lum: float, alfam: float, mass: float, periods: int,
               max_model_number: int) -> None:
    s = INLIST.read_text()

    s = replace_or_insert(s, r"^\s*create_RSP_model\s*=.*$", "   create_RSP_model = .true.")
    s = replace_or_insert(s, r"^\s*load_saved_model\s*=.*$", "   load_saved_model = .false.")
    s = re.sub(r"^\s*load_model_filename\s*=.*\n", "", s, flags=re.MULTILINE)
    s = re.sub(r"^\s*saved_model_name\s*=.*\n", "", s, flags=re.MULTILINE)

    # Stellar parameters.
    s = replace_or_insert(s, r"^\s*RSP_mass\s*=.*$", f"   RSP_mass = {mass:.5g}d0")
    s = replace_or_insert(s, r"^\s*RSP_Teff\s*=.*$", f"   RSP_Teff = {teff:.0f}d0")
    s = replace_or_insert(s, r"^\s*RSP_L\s*=.*$", f"   RSP_L = {lum:.5g}d0")

    # Kick (fixed seed).
    s = replace_or_insert(
        s, r"^\s*RSP_kick_vsurf_km_per_sec\s*=.*$",
        f"   RSP_kick_vsurf_km_per_sec = {FIXED_KICK:.4g}d0",
        anchor="RSP_Z",
    )

    # Eddy viscosity -- anchor on the full RSP_Z assignment line (safe: not a
    # prefix of any keyword, stays inside &controls).
    z_anchor = f"RSP_Z = {FIXED_Z}d0"
    if "RSP_alfam" in s:
        s = replace_or_insert(s, r"^\s*RSP_alfam\s*=.*$", f"   RSP_alfam = {alfam:.5g}d0")
    else:
        s = replace_or_insert(s, r"^\s*RSP_alfam\s*=.*$", f"   RSP_alfam = {alfam:.5g}d0",
                              anchor=z_anchor)

    # Fixed convective terms (swept and rejected in stage 1).
    alfam_anchor = f"RSP_alfam = {alfam:.5g}d0"
    if "RSP_alfat" in s:
        s = replace_or_insert(s, r"^\s*RSP_alfat\s*=.*$", f"   RSP_alfat = {FIXED_ALFAT:.4g}d0")
    else:
        s = replace_or_insert(s, r"^\s*RSP_alfat\s*=.*$", f"   RSP_alfat = {FIXED_ALFAT:.4g}d0",
                              anchor=alfam_anchor)
    alfat_anchor = f"RSP_alfat = {FIXED_ALFAT:.4g}d0"
    if "RSP_gammar" in s:
        s = replace_or_insert(s, r"^\s*RSP_gammar\s*=.*$", f"   RSP_gammar = {FIXED_GAMMAR:.4g}d0")
    else:
        s = replace_or_insert(s, r"^\s*RSP_gammar\s*=.*$", f"   RSP_gammar = {FIXED_GAMMAR:.4g}d0",
                              anchor=alfat_anchor)

    # Period target / number of periods / model cap.
    s = replace_or_insert(s, r"^\s*x_ctrl\(1\)\s*=.*$",
                          f"   x_ctrl(1) = {OBS_PERIOD:.9f}d0    ! observed FN Lyr period")
    s = replace_or_insert(s, r"^\s*x_integer_ctrl\(1\)\s*=.*$",
                          f"   x_integer_ctrl(1) = {periods}   ! which period to check")
    s = replace_or_insert(s, r"^\s*max_model_number\s*=.*$",
                          f"   max_model_number = {max_model_number}")

    # Keep output light.
    s = replace_or_insert(s, r"^\s*make_csv\s*=.*$", "   make_csv = .false.")
    s = replace_or_insert(s, r"^\s*sed_per_model\s*=.*$", "   sed_per_model = .false.")
    s = replace_or_insert(s, r"^\s*save_model_filename\s*=.*$",
                          "      save_model_filename = 'final.mod'")

    INLIST.write_text(s)


# ============================================================================
# History reading + scoring
# ============================================================================
def read_history(path: Path) -> pd.DataFrame:
    lines = path.read_text().splitlines()
    header_idx = next(i for i, line in enumerate(lines) if "model_number" in line.split())
    return pd.read_csv(path, sep=r"\s+", skiprows=header_idx)


def get_obs() -> pd.DataFrame:
    obs = pd.read_csv(OBS_FILE)
    obs = obs[np.isfinite(obs["phase"]) & np.isfinite(obs["mag"])].copy()
    return obs.sort_values("phase")


def periodic_interp(phase, mag, phase_new):
    phase = np.asarray(phase, float) % 1.0
    mag = np.asarray(mag, float)
    phase_new = np.asarray(phase_new, float) % 1.0
    good = np.isfinite(phase) & np.isfinite(mag)
    phase, mag = phase[good], mag[good]
    order = np.argsort(phase)
    phase, mag = phase[order], mag[order]
    _, uniq = np.unique(np.round(phase, 8), return_index=True)
    uniq = np.sort(uniq)
    phase, mag = phase[uniq], mag[uniq]
    x = np.concatenate([phase - 1.0, phase, phase + 1.0])
    y = np.concatenate([mag, mag, mag])
    return np.interp(phase_new, x, y)


def cycle_amp(h, mag_col, n):
    m = h[h["rsp_num_periods"].astype(int) == n]
    if len(m) < 200:
        return np.nan
    return np.nanmax(m[mag_col]) - np.nanmin(m[mag_col])


def convergence_range(h, mag_col, nmax, n_last=10):
    amps = [cycle_amp(h, mag_col, n) for n in range(max(0, nmax - n_last + 1), nmax + 1)]
    amps = [a for a in amps if np.isfinite(a)]
    if not amps:
        return np.nan, np.nan
    return float(np.nanmedian(amps)), float(np.nanmax(amps) - np.nanmin(amps))


def best_cycle_rms(h, mag_col, obs_phase, obs_mag, nmax, n_scan=10):
    """Over the last n_scan cycles, fit phase shift + offset and take min RMS."""
    shifts = np.linspace(0, 1, 1500, endpoint=False)
    best = None
    for n in range(max(0, nmax - n_scan + 1), nmax + 1):
        m = h[h["rsp_num_periods"].astype(int) == n]
        if len(m) < 200:
            continue
        ph = np.asarray(m["rsp_phase"], float) % 1.0
        mg = np.asarray(m[mag_col], float)
        good = np.isfinite(ph) & np.isfinite(mg)
        ph, mg = ph[good], mg[good]
        if len(ph) < 200 or (np.nanmax(ph) - np.nanmin(ph)) < 0.80:
            continue
        for shift in shifts:
            m_eval = periodic_interp(ph, mg, (obs_phase - shift) % 1.0)
            offset = np.nanmedian(obs_mag - m_eval)
            resid = obs_mag - (m_eval + offset)
            rms = float(np.sqrt(np.nanmean(resid ** 2)))
            if best is None or rms < best[0]:
                best = (rms, n, float(shift), float(offset))
    return best  # (rms, cycle, shift, offset) or None


# ============================================================================
# Parameter <-> normalized space (Nelder-Mead behaves badly across scales)
# ============================================================================
def to_norm(phys):
    return np.array([(v - p.lo) / (p.hi - p.lo) for v, p in zip(phys, PARAMS)])


def to_phys(norm):
    return np.array([p.lo + np.clip(x, 0.0, 1.0) * (p.hi - p.lo) for x, p in zip(norm, PARAMS)])


# ============================================================================
# Objective
# ============================================================================
_eval_counter = {"n": 0}
_best = {"obj": np.inf}


def clean_run_products():
    for name in ["LOGS", "photos", "png", "SED"]:
        p = Path(name)
        if p.exists():
            shutil.rmtree(p)
    Path("SED").mkdir(exist_ok=True)
    Path("png").mkdir(exist_ok=True)
    if Path("final.mod").exists():
        Path("final.mod").unlink()


def log_row(row: dict):
    new = not LOG_CSV.exists()
    # Stable column order.
    fields = ["eval", "timestamp", "teff", "lum", "alfam", "mass",
              "status", "period", "period_err_pct", "amp", "amp_frac",
              "amp_last10_median", "amp_last10_range", "vmax_cs_last",
              "best_cycle", "rms", "objective",
              "d_teff_nemec", "d_lum_nemec", "d_mass_nemec", "runtime_s"]
    with LOG_CSV.open("a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        if new:
            w.writeheader()
        w.writerow({k: row.get(k, "") for k in fields})


def objective(norm, cfg, obs_phase, obs_mag):
    _eval_counter["n"] += 1
    ev = _eval_counter["n"]
    teff, lum, alfam, mass = to_phys(norm)
    t0 = time.time()

    print(f"\n--- eval {ev}: Teff={teff:.1f} L={lum:.3f} "
          f"alfam={alfam:.4f} M={mass:.4f} ---", flush=True)

    set_inlist(teff, lum, alfam, mass, cfg["periods"], cfg["max_model_number"])
    clean_run_products()

    proc = subprocess.run(RUN_CMD, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          text=True, check=False)

    row = dict(eval=ev, timestamp=time.strftime("%Y-%m-%d %H:%M:%S"),
               teff=teff, lum=lum, alfam=alfam, mass=mass,
               d_teff_nemec=round(teff - NEMEC["teff"], 1),
               d_lum_nemec=round(lum - NEMEC["lum"], 3),
               d_mass_nemec=round(mass - NEMEC["mass"], 4))

    hist = Path("LOGS/history.data")
    if proc.returncode != 0 or not hist.exists():
        row.update(status="crash", objective=PENALTY_CRASH,
                   runtime_s=round(time.time() - t0, 1))
        log_row(row)
        print(f"    CRASH (rc={proc.returncode}) -> penalty {PENALTY_CRASH}", flush=True)
        return PENALTY_CRASH

    try:
        h = read_history(hist)
        mag_col = "Kepler" if "Kepler" in h.columns else "K"
        nmax = int(np.nanmax(h["rsp_num_periods"]))
        amp_med, amp_range = convergence_range(h, mag_col, nmax)
        vmax = (np.nanmax(np.abs(h[h["rsp_num_periods"].astype(int) == nmax]["max_abs_v_div_cs"]))
                if "max_abs_v_div_cs" in h.columns else np.nan)
    except Exception as exc:
        row.update(status=f"read_error:{exc!r}", objective=PENALTY_CRASH,
                   runtime_s=round(time.time() - t0, 1))
        log_row(row)
        print(f"    READ ERROR {exc!r} -> penalty {PENALTY_CRASH}", flush=True)
        return PENALTY_CRASH

    # Period from the last cycle.
    m_last = h[h["rsp_num_periods"].astype(int) == nmax]
    period = float(np.nanmax(m_last["star_age_day"]) - np.nanmin(m_last["star_age_day"]))
    period_err = abs(period - OBS_PERIOD) / OBS_PERIOD

    row.update(period=round(period, 6), period_err_pct=round(100 * period_err, 4),
               amp=round(float(amp_med), 5) if np.isfinite(amp_med) else "",
               amp_frac=round(float(amp_med) / OBS_AMP, 4) if np.isfinite(amp_med) else "",
               amp_last10_median=round(float(amp_med), 5) if np.isfinite(amp_med) else "",
               amp_last10_range=round(float(amp_range), 6) if np.isfinite(amp_range) else "",
               vmax_cs_last=round(float(vmax), 4) if np.isfinite(vmax) else "")

    # Convergence gate: scoring a transient is invalid.
    if not np.isfinite(amp_range) or amp_range > CONV_RANGE_MAX:
        graded = PENALTY_UNCONV_BASE + 10.0 * (amp_range if np.isfinite(amp_range) else 1.0)
        row.update(status="unconverged", rms="", objective=round(graded, 4),
                   runtime_s=round(time.time() - t0, 1))
        log_row(row)
        print(f"    UNCONVERGED (amp_last10_range={amp_range}) -> penalty {graded:.3f}",
              flush=True)
        return graded

    # Waveform residual on the best converged cycle.
    bc = best_cycle_rms(h, mag_col, obs_phase, obs_mag, nmax)
    if bc is None:
        row.update(status="no_scorable_cycle", rms="", objective=PENALTY_UNCONV_BASE,
                   runtime_s=round(time.time() - t0, 1))
        log_row(row)
        print("    NO SCORABLE CYCLE -> penalty", flush=True)
        return PENALTY_UNCONV_BASE

    rms, cyc, _, _ = bc
    amp_frac_dev = abs(float(amp_med) - OBS_AMP) / OBS_AMP

    obj = (rms
           + PERIOD_W * max(0.0, period_err - PERIOD_TOL)
           + AMP_W * max(0.0, amp_frac_dev - AMP_TOL))

    row.update(status="ok", best_cycle=cyc, rms=round(rms, 5),
               objective=round(obj, 5), runtime_s=round(time.time() - t0, 1))
    log_row(row)
    print(f"    ok  cycle={cyc} period_err={100*period_err:.3f}% "
          f"amp={amp_med:.4f} vmax/cs={vmax:.2f} RMS={rms:.5f} -> obj={obj:.5f}",
          flush=True)

    # Save best run's products.
    if obj < _best["obj"]:
        _best["obj"] = obj
        shutil.copy2(hist, BEST_DIR / "history.data")
        shutil.copy2(INLIST, BEST_DIR / "inlist_rsp_RR_Lyrae")
        (BEST_DIR / "best_params.txt").write_text(
            f"eval={ev}\nTeff={teff:.2f}\nL={lum:.4f}\nalfam={alfam:.5f}\nmass={mass:.5f}\n"
            f"period={period:.6f} d (err {100*period_err:.3f}%)\n"
            f"amp_last10_median={amp_med:.5f} (frac {amp_med/OBS_AMP:.4f})\n"
            f"vmax_cs_last={vmax:.4f}\nbest_cycle={cyc}\nRMS={rms:.5f}\nobjective={obj:.5f}\n"
        )
        print(f"    *** new best: obj={obj:.5f} (saved to {BEST_DIR}/) ***", flush=True)

    return obj


# ============================================================================
# Initial simplex
# ============================================================================
def build_initial_simplex(frac, resume):
    n = len(PARAMS)
    seed = to_norm([p.seed for p in PARAMS])

    if resume and LOG_CSV.exists():
        df = pd.read_csv(LOG_CSV)
        df = df[df["status"] == "ok"].copy()
        if len(df) >= n + 1:
            df = df.sort_values("objective").head(n + 1)
            verts = [to_norm([r.teff, r.lum, r.alfam, r.mass]) for r in df.itertuples()]
            print(f"Resuming: seeded simplex from {n+1} best logged evals.", flush=True)
            return np.array(verts)
        print("Resume requested but too few 'ok' evals logged; using fresh simplex.",
              flush=True)

    simplex = [seed]
    for i in range(n):
        v = seed.copy()
        v[i] = np.clip(v[i] + frac, 0.0, 1.0)
        simplex.append(v)
    return np.array(simplex)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=list(MODE_CFG), default="fast")
    ap.add_argument("--resume", action="store_true")
    args = ap.parse_args()

    cfg = MODE_CFG[args.mode]
    print(f"Mode: {args.mode}  periods={cfg['periods']}  max_evals={cfg['max_evals']}",
          flush=True)
    print(f"Free params: {[p.name for p in PARAMS]}", flush=True)
    print(f"Bounds: {[(p.lo, p.hi) for p in PARAMS]}", flush=True)

    obs = get_obs()
    obs_phase = obs["phase"].to_numpy(float)
    obs_mag = obs["mag"].to_numpy(float)

    # Snapshot the inlist before touching it.
    (BEST_DIR / "inlist_original_before_opt").write_text(INLIST.read_text())

    init = build_initial_simplex(cfg["simplex_frac"], args.resume)

    res = minimize(
        objective,
        x0=init[0],
        args=(cfg, obs_phase, obs_mag),
        method="Nelder-Mead",
        options=dict(
            initial_simplex=init,
            maxfev=cfg["max_evals"],
            xatol=0.01,     # ~1% of a normalized bound range
            fatol=0.002,    # ~RMS resolution we care about
            disp=True,
        ),
    )

    best_phys = to_phys(res.x)
    print("\n" + "=" * 60)
    print("OPTIMIZATION FINISHED")
    print("=" * 60)
    print(f"evals run:   {_eval_counter['n']}")
    print(f"best obj:    {_best['obj']:.5f}")
    print(f"best params: Teff={best_phys[0]:.1f}  L={best_phys[1]:.4f}  "
          f"alfam={best_phys[2]:.5f}  M={best_phys[3]:.5f}")
    print(f"\nFull trace:  {LOG_CSV}")
    print(f"Best run:    {BEST_DIR}/ (history.data, inlist, best_params.txt)")
    print("\nNext: run compare_fn_lyr_kepler_absolute.py against "
          f"{BEST_DIR}/history.data to inspect the waveform.")


if __name__ == "__main__":
    main()
