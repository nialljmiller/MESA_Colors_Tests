#!/usr/bin/env python3
"""
Stage 2 morphology-targeted Nelder-Mead simplex for FN Lyr RSP+Colors fit.

Stage 1 (optimize_fn_lyr_morphology.py, RSP_cq free) showed RSP_cq has no
effect on the shock: vmax/cs stayed pinned at 7.7-8.1 across cq=4.98-6.48,
and all morphology improvement came from Teff/L drifting upward. cq is
therefore dropped here and replaced with RSP_zsh, the shock onset threshold
(default 0.1) -- this gates WHEN artificial viscosity switches on, which is
the more likely lever if the shock is triggering at the wrong compression
phase rather than being under-damped once triggered.

Free parameters
---------------
  RSP_zsh    shock onset threshold       [0.02, 0.40]   seed 0.10
  RSP_alfam  eddy viscosity               [0.10, 0.20]   seed 0.1524
  RSP_Teff   strip position               [6650, 6800]   seed 6728.5
  RSP_L      luminosity                   [46.0, 52.0]   seed 49.10

Seeded at the stage-1 eval-23 basin (best objective: 0.05576).

Fixed
-----
  RSP_mass = 0.65, kick = 40, alfat = 0, gammar = 0, RSP_cq = 4.0 (default),
  Z = 0.0004

Objective is identical to stage 1: amplitude-normalized shape RMS
(survey_fn_lyr.py-style _unit_norm comparison) plus soft period/amplitude
guards.

Usage
-----
  python optimize_fn_lyr_morphology_zsh.py --mode fast | tee opt_morph_zsh_fast.log
  python optimize_fn_lyr_morphology_zsh.py --mode production | tee opt_morph_zsh_prod.log
  python optimize_fn_lyr_morphology_zsh.py --mode production --resume | tee -a opt_morph_zsh_prod.log

Run from:
  ~/MESA/MESA_Colors_Tests/rsp_RR_Lyrae/FN_Lyr
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


# ---------------------------------------------------------------------------
# Observed targets
# ---------------------------------------------------------------------------
OBS_PERIOD = 0.527398471      # days, Nemec et al. 2011
OBS_AMP    = 1.0464           # Kepler amplitude (mag)
OBS_FILE   = "FN_Lyr_KIC6936115_kepler_binned.csv"

INLIST  = Path("inlist_rsp_RR_Lyrae")
RUN_CMD = ["./rn"]

LOG_CSV  = Path("optimize_morph_zsh_log.csv")
BEST_DIR = Path("optimize_morph_zsh_best")
BEST_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Fixed parameters
# ---------------------------------------------------------------------------
FIXED_MASS   = 0.65
FIXED_KICK   = 40.0
FIXED_ALFAT  = 0.0
FIXED_GAMMAR = 0.0
FIXED_CQ     = 4.0     # back to RSP default -- stage 1 showed cq has no effect
FIXED_Z      = 0.0004

# ---------------------------------------------------------------------------
# Free parameters -- seeded at stage-1 eval-23 basin
# ---------------------------------------------------------------------------
@dataclass
class Param:
    name:  str
    lo:    float
    hi:    float
    seed:  float

PARAMS = [
    # RSP_zsh: shock onset threshold. Default 0.1. Lower values trigger
    # artificial viscosity earlier/more readily (in zsh, larger -> later
    # onset, more restrictive). Range spans well below and above default.
    Param("zsh",   0.02,  0.40,  0.10),
    Param("alfam", 0.10,  0.20,  0.1524),
    Param("teff",  6650., 6800., 6728.5),
    Param("lum",   46.0,  52.0,  49.0954),
]

# ---------------------------------------------------------------------------
# Objective weights / gates  (unchanged from stage 1)
# ---------------------------------------------------------------------------
PERIOD_TOL = 0.015
PERIOD_W   = 1.5
AMP_TOL    = 0.08
AMP_W      = 0.5

CONV_RANGE_MAX = 0.02

PENALTY_CRASH   = 5.0
PENALTY_UNCONV_BASE = 2.0

MIN_POINTS_PER_CYCLE = 250
MIN_PHASE_SPAN       = 0.80
SEARCH_LAST_N_CYCLES = 20
N_SHIFTS             = 1500

# ---------------------------------------------------------------------------
# Mode config
# ---------------------------------------------------------------------------
MODE_CFG = {
    "fast":       dict(periods=150, max_evals=25,  simplex_frac=0.08, max_model_number=260000),
    "production": dict(periods=200, max_evals=70,  simplex_frac=0.10, max_model_number=320000),
}

# ---------------------------------------------------------------------------
# Global state
# ---------------------------------------------------------------------------
_eval_counter = {"n": 0}
_best = {"obj": np.inf}


# ===========================================================================
# Normalisation helpers
# ===========================================================================
def to_norm(phys: list[float]) -> np.ndarray:
    return np.array([(v - p.lo) / (p.hi - p.lo) for v, p in zip(phys, PARAMS)],
                    dtype=float)

def to_phys(norm: np.ndarray) -> list[float]:
    return [float(np.clip(n, 0.0, 1.0) * (p.hi - p.lo) + p.lo)
            for n, p in zip(norm, PARAMS)]


# ===========================================================================
# Inlist editing
# ===========================================================================
def replace_or_insert(text: str, pattern: str, replacement: str,
                      anchor: str | None = None) -> str:
    new, n = re.subn(pattern, replacement, text, flags=re.MULTILINE)
    if n > 0:
        return new
    if anchor is not None and anchor in text:
        return text.replace(anchor, anchor + "\n   " + replacement)
    return text + "\n   " + replacement + "\n"


def set_inlist(zsh: float, alfam: float, teff: float, lum: float, cfg: dict) -> None:
    s = INLIST.read_text()

    # Fresh model every time.
    s = replace_or_insert(s, r"^\s*create_RSP_model\s*=.*$",
                          "   create_RSP_model = .true.")
    s = replace_or_insert(s, r"^\s*load_saved_model\s*=.*$",
                          "   load_saved_model = .false.")
    s = re.sub(r"^\s*load_model_filename\s*=.*\n", "", s, flags=re.MULTILINE)

    # Stellar parameters.
    s = replace_or_insert(s, r"^\s*RSP_Teff\s*=.*$",  f"   RSP_Teff = {teff:.2f}d0")
    s = replace_or_insert(s, r"^\s*RSP_L\s*=.*$",     f"   RSP_L = {lum:.4f}d0")
    s = replace_or_insert(s, r"^\s*RSP_mass\s*=.*$",  f"   RSP_mass = {FIXED_MASS:.4f}d0")

    # Kick.
    s = replace_or_insert(s, r"^\s*RSP_kick_vsurf_km_per_sec\s*=.*$",
                          f"   RSP_kick_vsurf_km_per_sec = {FIXED_KICK:.1f}d0",
                          anchor="RSP_Z = 0.0004d0")

    # Convective parameters, chained anchors.
    z_anchor = "RSP_Z = 0.0004d0"
    s = replace_or_insert(s, r"^\s*RSP_alfam\s*=.*$",
                          f"   RSP_alfam = {alfam:.5f}d0", anchor=z_anchor)

    alfam_anchor = f"RSP_alfam = {alfam:.5f}d0"
    s = replace_or_insert(s, r"^\s*RSP_alfat\s*=.*$",
                          f"   RSP_alfat = {FIXED_ALFAT:.1f}d0", anchor=alfam_anchor)

    alfat_anchor = f"RSP_alfat = {FIXED_ALFAT:.1f}d0"
    s = replace_or_insert(s, r"^\s*RSP_gammar\s*=.*$",
                          f"   RSP_gammar = {FIXED_GAMMAR:.1f}d0", anchor=alfat_anchor)

    # cq pinned back to RSP default (stage 1 showed it does nothing here).
    gammar_anchor = f"RSP_gammar = {FIXED_GAMMAR:.1f}d0"
    s = replace_or_insert(s, r"^\s*RSP_cq\s*=.*$",
                          f"   RSP_cq = {FIXED_CQ:.4f}d0", anchor=gammar_anchor)

    # RSP_zsh: the new morphology lever, stage 2.
    cq_anchor = f"RSP_cq = {FIXED_CQ:.4f}d0"
    s = replace_or_insert(s, r"^\s*RSP_zsh\s*=.*$",
                          f"   RSP_zsh = {zsh:.5f}d0", anchor=cq_anchor)

    # Run length.
    s = replace_or_insert(s, r"^\s*x_integer_ctrl\(1\)\s*=.*$",
                          f"   x_integer_ctrl(1) = {cfg['periods']}")
    s = replace_or_insert(s, r"^\s*x_ctrl\(1\)\s*=.*$",
                          f"   x_ctrl(1) = {OBS_PERIOD:.9f}d0")
    s = replace_or_insert(s, r"^\s*max_model_number\s*=.*$",
                          f"   max_model_number = {cfg['max_model_number']}")

    INLIST.write_text(s)


# ===========================================================================
# History reading helpers (identical to stage 1 / survey_fn_lyr.py)
# ===========================================================================
def read_history(path: Path) -> pd.DataFrame:
    lines = path.read_text(errors="ignore").splitlines()
    hi = next(i for i, l in enumerate(lines) if "model_number" in l.split())
    return pd.read_csv(path, sep=r"\s+", skiprows=hi)


def periodic_interp(phase, mag, pnew):
    phase = np.asarray(phase, float) % 1.0
    mag   = np.asarray(mag,   float)
    pnew  = np.asarray(pnew,  float) % 1.0
    g = np.isfinite(phase) & np.isfinite(mag)
    phase, mag = phase[g], mag[g]
    o = np.argsort(phase); phase, mag = phase[o], mag[o]
    _, u = np.unique(np.round(phase, 8), return_index=True)
    phase, mag = phase[np.sort(u)], mag[np.sort(u)]
    x = np.concatenate([phase - 1, phase, phase + 1])
    y = np.concatenate([mag, mag, mag])
    return np.interp(pnew, x, y)


def _unit_norm(mag):
    mag = np.asarray(mag, float)
    rng = np.nanmax(mag) - np.nanmin(mag)
    if not np.isfinite(rng) or rng <= 0:
        return None
    return (mag - np.nanmin(mag)) / rng


def last_full_cycle(h):
    mc = None
    for c in ("Kepler", "K"):
        if c in h.columns:
            mc = c; break
    if mc is None or "rsp_num_periods" not in h.columns:
        return None
    nmax = int(np.nanmax(h["rsp_num_periods"]))
    for n in range(nmax, max(0, nmax - SEARCH_LAST_N_CYCLES) - 1, -1):
        cyc = h[h["rsp_num_periods"].astype(int) == n]
        if len(cyc) < MIN_POINTS_PER_CYCLE:
            continue
        ph = np.asarray(cyc["rsp_phase"], float) % 1.0
        mg = np.asarray(cyc[mc], float)
        g  = np.isfinite(ph) & np.isfinite(mg)
        ph, mg = ph[g], mg[g]
        if (ph.max() - ph.min()) < MIN_PHASE_SPAN:
            continue
        o = np.argsort(ph); ph, mg = ph[o], mg[o]
        if "rsp_period_in_days" in h.columns:
            pv = pd.to_numeric(cyc["rsp_period_in_days"], errors="coerce")
            pv = pv[(pv > 0.1) & (pv < 2.0)]
            period = float(np.nanmedian(pv)) if len(pv) else np.nan
        else:
            period = float(cyc["star_age_day"].max() - cyc["star_age_day"].min())
        vmax = (float(np.nanmax(np.abs(cyc["max_abs_v_div_cs"])))
                if "max_abs_v_div_cs" in h.columns else np.nan)
        return ph, mg, period, vmax, n
    return None


def convergence_range(h, mc, nmax):
    amps = []
    for n in range(max(0, nmax - 10), nmax + 1):
        cyc = h[h["rsp_num_periods"].astype(int) == n]
        if len(cyc) < MIN_POINTS_PER_CYCLE:
            continue
        mg = cyc[mc].to_numpy(float)
        amps.append(np.nanmax(mg) - np.nanmin(mg))
    if len(amps) < 3:
        return np.nan, np.nan
    return float(np.nanmedian(amps)), float(np.nanmax(amps) - np.nanmin(amps))


def shape_residual_score(mph, mmg, obs_phase, obs_mag):
    shifts = np.linspace(0, 1, N_SHIFTS, endpoint=False)
    best = None
    om = _unit_norm(obs_mag)
    if om is None:
        return None

    for s in shifts:
        m_on_obs = periodic_interp(mph, mmg, (obs_phase - s) % 1.0)
        mm = _unit_norm(m_on_obs)
        if mm is None:
            continue
        resid = om - mm
        rms = float(np.sqrt(np.nanmean(resid ** 2)))
        off = float(np.nanmedian(obs_mag - m_on_obs))
        if best is None or rms < best["shape_rms"]:
            best = dict(shape_rms=rms, shift=float(s), offset=off)

    return best


# ===========================================================================
# Logging
# ===========================================================================
LOG_FIELDS = [
    "eval", "zsh", "alfam", "teff", "lum",
    "status", "period", "period_err_pct",
    "amp", "amp_frac", "amp_last10_range",
    "shape_rms", "objective", "vmax_cs_last", "runtime_s",
]

def log_row(row: dict) -> None:
    write_header = not LOG_CSV.exists()
    with LOG_CSV.open("a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=LOG_FIELDS, extrasaction="ignore")
        if write_header:
            w.writeheader()
        w.writerow(row)


# ===========================================================================
# Objective
# ===========================================================================
def get_obs():
    obs = pd.read_csv(OBS_FILE)
    obs = obs[np.isfinite(obs["phase"]) & np.isfinite(obs["mag"])].sort_values("phase")
    return obs


def objective(norm_x: np.ndarray, cfg: dict,
              obs_phase: np.ndarray, obs_mag: np.ndarray) -> float:
    _eval_counter["n"] += 1
    ev = _eval_counter["n"]

    zsh, alfam, teff, lum = to_phys(norm_x)

    print(f"\n[eval {ev}]  zsh={zsh:.4f}  alfam={alfam:.4f}  "
          f"Teff={teff:.1f}  L={lum:.3f}", flush=True)

    row = dict(eval=ev, zsh=round(zsh, 5), alfam=round(alfam, 5),
               teff=round(teff, 2), lum=round(lum, 4))

    set_inlist(zsh, alfam, teff, lum, cfg)

    for name in ("LOGS", "photos", "png"):
        p = Path(name)
        if p.exists():
            shutil.rmtree(p)
    if Path("final.mod").exists():
        Path("final.mod").unlink()

    t0 = time.time()
    proc = subprocess.run(RUN_CMD, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True, check=False)

    hist = Path("LOGS/history.data")
    if proc.returncode != 0 or not hist.exists():
        row.update(status="crash", objective=PENALTY_CRASH,
                   runtime_s=round(time.time() - t0, 1))
        log_row(row)
        print(f"    CRASH -> penalty {PENALTY_CRASH}", flush=True)
        return PENALTY_CRASH

    try:
        h  = read_history(hist)
        mc = "Kepler" if "Kepler" in h.columns else "K"
        nmax = int(np.nanmax(h["rsp_num_periods"]))
        amp_med, amp_range = convergence_range(h, mc, nmax)
    except Exception as exc:
        row.update(status=f"read_error:{exc!r}", objective=PENALTY_CRASH,
                   runtime_s=round(time.time() - t0, 1))
        log_row(row)
        print(f"    READ ERROR {exc!r} -> penalty", flush=True)
        return PENALTY_CRASH

    res = last_full_cycle(h)
    if res is None:
        row.update(status="no_cycle", objective=PENALTY_UNCONV_BASE,
                   runtime_s=round(time.time() - t0, 1))
        log_row(row); return PENALTY_UNCONV_BASE

    mph, mmg, period, vmax, cyc = res

    row.update(
        period=round(period, 6) if np.isfinite(period) else "",
        period_err_pct=round(100 * abs(period - OBS_PERIOD) / OBS_PERIOD, 4)
                       if np.isfinite(period) else "",
        amp=round(float(amp_med), 5) if np.isfinite(amp_med) else "",
        amp_frac=round(float(amp_med) / OBS_AMP, 4) if np.isfinite(amp_med) else "",
        amp_last10_range=round(float(amp_range), 6) if np.isfinite(amp_range) else "",
        vmax_cs_last=round(float(vmax), 3) if np.isfinite(vmax) else "",
    )

    if not np.isfinite(amp_range) or amp_range > CONV_RANGE_MAX:
        graded = PENALTY_UNCONV_BASE + 10.0 * (amp_range if np.isfinite(amp_range) else 1.0)
        row.update(status="unconverged", shape_rms="", objective=round(graded, 4),
                   runtime_s=round(time.time() - t0, 1))
        log_row(row)
        print(f"    UNCONVERGED (range={amp_range:.4f}) -> penalty {graded:.3f}", flush=True)
        return graded

    period_err = abs(period - OBS_PERIOD) / OBS_PERIOD if np.isfinite(period) else 1.0
    period_pen = PERIOD_W * max(0.0, period_err - PERIOD_TOL)

    amp_frac_dev = abs(float(amp_med) - OBS_AMP) / OBS_AMP if np.isfinite(amp_med) else 1.0
    amp_pen = AMP_W * max(0.0, amp_frac_dev - AMP_TOL)

    sc = shape_residual_score(mph, mmg, obs_phase, obs_mag)
    if sc is None:
        row.update(status="no_shape", objective=PENALTY_UNCONV_BASE,
                   runtime_s=round(time.time() - t0, 1))
        log_row(row); return PENALTY_UNCONV_BASE

    shape_rms = sc["shape_rms"]
    obj = shape_rms + period_pen + amp_pen

    row.update(status="ok", shape_rms=round(shape_rms, 5),
               objective=round(obj, 5), runtime_s=round(time.time() - t0, 1))
    log_row(row)
    print(f"    ok  cycle={cyc}  period_err={100*period_err:.2f}%  "
          f"amp={amp_med:.4f}  vmax/cs={vmax:.2f}  "
          f"shape_rms={shape_rms:.5f}  obj={obj:.5f}", flush=True)

    if obj < _best["obj"]:
        _best["obj"] = obj
        shutil.copy2(hist, BEST_DIR / "history.data")
        shutil.copy2(INLIST, BEST_DIR / "inlist_rsp_RR_Lyrae")
        (BEST_DIR / "best_params.txt").write_text(
            f"eval={ev}\nzsh={zsh:.5f}\nalfam={alfam:.5f}\n"
            f"teff={teff:.2f}\nlum={lum:.4f}\n"
            f"period={period:.6f} d  (err {100*period_err:.3f}%)\n"
            f"amp_last10_median={amp_med:.5f}\n"
            f"vmax_cs_last={vmax:.4f}\n"
            f"shape_rms={shape_rms:.5f}\nobjective={obj:.5f}\n"
        )
        print(f"    *** new best: obj={obj:.5f} -> {BEST_DIR}/ ***", flush=True)

    return obj


# ===========================================================================
# Initial simplex
# ===========================================================================
def build_initial_simplex(frac: float, resume: bool) -> np.ndarray:
    n    = len(PARAMS)
    seed = to_norm([p.seed for p in PARAMS])

    if resume and LOG_CSV.exists():
        df = pd.read_csv(LOG_CSV)
        df = df[df["status"] == "ok"].copy()
        if len(df) >= n + 1:
            df = df.sort_values("objective").head(n + 1)
            verts = [to_norm([r.zsh, r.alfam, r.teff, r.lum])
                     for r in df.itertuples()]
            print(f"Resuming: seeded from {n+1} best logged evals.", flush=True)
            return np.array(verts)
        print("Too few ok evals to resume; using fresh simplex.", flush=True)

    simplex = [seed]
    for i in range(n):
        v = seed.copy()
        v[i] = np.clip(v[i] + frac, 0.0, 1.0)
        simplex.append(v)
    return np.array(simplex)


# ===========================================================================
# Main
# ===========================================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=list(MODE_CFG), default="fast")
    ap.add_argument("--resume", action="store_true")
    args = ap.parse_args()

    cfg = MODE_CFG[args.mode]
    print(f"Mode: {args.mode}  periods={cfg['periods']}  max_evals={cfg['max_evals']}",
          flush=True)
    print(f"Free params: {[p.name for p in PARAMS]}", flush=True)
    print(f"Seeds: zsh={PARAMS[0].seed}  alfam={PARAMS[1].seed}  "
          f"teff={PARAMS[2].seed}  lum={PARAMS[3].seed}", flush=True)
    print(f"Fixed: cq={FIXED_CQ} (RSP default; stage 1 found no effect)",
          flush=True)

    obs       = get_obs()
    obs_phase = obs["phase"].to_numpy(float)
    obs_mag   = obs["mag"].to_numpy(float)

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
            xatol=0.01,
            fatol=0.002,
            disp=True,
        ),
    )

    best_phys = to_phys(res.x)
    print("\n" + "=" * 60)
    print("OPTIMIZATION FINISHED")
    print("=" * 60)
    print(f"evals:       {_eval_counter['n']}")
    print(f"best obj:    {_best['obj']:.5f}")
    print(f"best params: zsh={best_phys[0]:.4f}  alfam={best_phys[1]:.4f}  "
          f"Teff={best_phys[2]:.1f}  L={best_phys[3]:.4f}")
    print(f"\nLog:      {LOG_CSV}")
    print(f"Best run: {BEST_DIR}/")
    print("\nNext: python plot_fn_lyr_fit.py "
          f"--history {BEST_DIR}/history.data --out figures/fig_fn_lyr_fit.pdf")


if __name__ == "__main__":
    main()
