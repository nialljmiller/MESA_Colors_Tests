#!/usr/bin/env python3
"""
Take stock of the entire FN Lyr fitting campaign -- every grid run, every
parameter sweep, and the simplex -- scored on the THREE axes that define a
good fit: period, amplitude, and morphology.

Morphology is scored to capture the specific defects seen by eye that a plain
RMS misses:
  - a spurious secondary maximum on the descending branch (the "phase ~0.25
    bump") that the real FN Lyr light curve does not have, and
  - a split / early-peaking maximum (the shock shoulder).
We therefore report, per model:
    period_err   |P - P_obs| / P_obs
    amp_err      |A - A_obs| / A_obs
    r            Pearson corr of the aligned curve (overall shape)
    bump         amplitude of the largest spurious secondary maximum on the
                 descending branch, in mag (0 = none; bigger = worse)
    rise_err     |risetime_model - risetime_obs| (phase units)
    vmax_cs      max |v|/c_s on the scored cycle (shock proxy)
A combined, equally-weighted morphology score and an overall score are formed
so the campaign can be ranked, and a period-vs-amplitude trade plot (colored by
morphology) shows whether all three can be minimized at once.

Outputs:
    fn_lyr_campaign_scores.csv     one row per model, all axes + params
    fn_lyr_campaign_tradeoff.pdf   period_err vs amp_err, colored by morphology

Usage:
    python survey_fn_lyr.py
    python survey_fn_lyr.py --min-bytes 1000000 --obs FN_Lyr_KIC6936115_kepler_binned.csv
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ---- observed reference (Nemec et al. 2011) -------------------------------
P_OBS = 0.527398471
A_OBS = 1.0464
# Observed risetime (fraction of period from min to max light). Nemec Table 1
# gives RT(Kp) ~ 0.118 for FN Lyr; we compute the model's the same way.
RT_OBS = 0.118

MIN_POINTS_PER_CYCLE = 250
MIN_PHASE_SPAN = 0.80
SEARCH_LAST_N_CYCLES = 20
N_SHIFTS = 1500


# ---------------------------------------------------------------------------
def read_history(path: Path) -> pd.DataFrame:
    lines = path.read_text(errors="ignore").splitlines()
    hi = next(i for i, l in enumerate(lines) if "model_number" in l.split())
    return pd.read_csv(path, sep=r"\s+", skiprows=hi)


def mag_col(h):
    for c in ("Kepler", "K"):
        if c in h.columns:
            return c
    cand = [c for c in h.columns if "kepler" in c.lower()]
    return cand[0] if cand else None


def periodic_interp(phase, mag, pnew):
    phase = np.asarray(phase, float) % 1.0
    mag = np.asarray(mag, float)
    pnew = np.asarray(pnew, float) % 1.0
    g = np.isfinite(phase) & np.isfinite(mag)
    phase, mag = phase[g], mag[g]
    o = np.argsort(phase)
    phase, mag = phase[o], mag[o]
    _, u = np.unique(np.round(phase, 8), return_index=True)
    u = np.sort(u)
    phase, mag = phase[u], mag[u]
    x = np.concatenate([phase - 1, phase, phase + 1])
    y = np.concatenate([mag, mag, mag])
    return np.interp(pnew, x, y)


def last_full_cycle(h):
    """Return (phase, mag) of the last cycle with a real number of points,
    plus its period from rsp_period_in_days. Skips the degenerate final cycle."""
    mc = mag_col(h)
    if mc is None or "rsp_num_periods" not in h.columns:
        return None
    nmax = int(np.nanmax(h["rsp_num_periods"]))
    for n in range(nmax, max(0, nmax - SEARCH_LAST_N_CYCLES) - 1, -1):
        cyc = h[h["rsp_num_periods"].astype(int) == n]
        if len(cyc) < MIN_POINTS_PER_CYCLE:
            continue
        ph = np.asarray(cyc["rsp_phase"], float) % 1.0
        mg = np.asarray(cyc[mc], float)
        good = np.isfinite(ph) & np.isfinite(mg)
        ph, mg = ph[good], mg[good]
        if (ph.max() - ph.min()) < MIN_PHASE_SPAN:
            continue
        o = np.argsort(ph)
        ph, mg = ph[o], mg[o]
        # period
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


def _unit_norm(mag):
    """Normalize a magnitude curve to unit peak-to-peak range, max light = 0.
    Removes amplitude so SHAPE can be compared independently of amplitude."""
    mag = np.asarray(mag, float)
    rng = np.nanmax(mag) - np.nanmin(mag)
    if not np.isfinite(rng) or rng <= 0:
        return None
    return (mag - np.nanmin(mag)) / rng   # 0 at brightest, 1 at faintest


def shape_residual(model_phase, model_mag, obs_phase, obs_mag, shift):
    """RMS residual between the AMPLITUDE-NORMALIZED model and observed curves,
    on the observed phase grid, using the already-chosen phase shift.

    This is the honest morphology metric: both curves are scaled to unit range
    first, so amplitude (scored separately) does not leak in. The residual is
    large only where the SHAPE differs -- e.g. where the model over-produces
    the descending-branch bump or sharpens the rise relative to the real star.
    Returns (total, rise_part, desc_part), all in normalized-flux units."""
    om = _unit_norm(obs_mag)
    if om is None:
        return np.nan, np.nan, np.nan
    # model on obs grid, then normalize the SAME way
    m_on_obs = periodic_interp(model_phase, model_mag, (obs_phase - shift) % 1.0)
    mm = _unit_norm(m_on_obs)
    if mm is None:
        return np.nan, np.nan, np.nan
    resid = om - mm
    total = float(np.sqrt(np.nanmean(resid ** 2)))
    # split at max light (phase of brightest observed point) into rise vs descent
    grid = obs_phase % 1.0
    imax_obs = int(np.nanargmin(om))           # brightest obs point
    ph_max = grid[imax_obs]
    # "rise" = the ~20% of phase just before max light; "descent" = the rest
    dphase = (grid - ph_max) % 1.0
    rise_mask = (dphase > 0.80) | (dphase < 0.02)   # just before & at max
    desc_mask = ~rise_mask
    rise_part = float(np.sqrt(np.nanmean(resid[rise_mask] ** 2))) if rise_mask.any() else np.nan
    desc_part = float(np.sqrt(np.nanmean(resid[desc_mask] ** 2))) if desc_mask.any() else np.nan
    return total, rise_part, desc_part


def risetime(phase, mag):
    """Fraction of cycle from min light (faintest) to max light (brightest)."""
    grid = np.linspace(0, 1, 400, endpoint=False)
    mg = periodic_interp(phase, mag, grid)
    imax = int(np.argmin(mg))   # brightest
    imin = int(np.argmax(mg))   # faintest
    rt = (imax - imin) / 400.0 % 1.0
    return float(rt)


def align_score(mph, mmg, obs_phase, obs_mag):
    """Best phase shift + offset; return r, rms, and aligned model on obs grid."""
    shifts = np.linspace(0, 1, N_SHIFTS, endpoint=False)
    best = None
    for s in shifts:
        m = periodic_interp(mph, mmg, (obs_phase - s) % 1.0)
        off = np.nanmedian(obs_mag - m)
        resid = obs_mag - (m + off)
        rms = float(np.sqrt(np.nanmean(resid ** 2)))
        if best is None or rms < best[0]:
            best = (rms, s, off, m + off)
    rms, s, off, pred = best
    r = float(np.corrcoef(obs_mag - obs_mag.mean(), pred - pred.mean())[0, 1])
    return r, rms, s, off


PARAM_RE = re.compile(r"T(\d+)_L(\d+)p(\d+)(?:_kick(\d+))?"
                      r"(?:_alfam(\d+)p(\d+))?(?:_alfat(\d+)p(\d+))?"
                      r"(?:_gammar(\d+)p(\d+))?")


def parse_params(path: Path):
    """Pull T/L/kick/alfam/alfat/gammar out of the directory name if present."""
    name = path.parent.name
    m = PARAM_RE.search(name)
    d = {}
    if m:
        g = m.groups()
        d["teff"] = float(g[0]) if g[0] else np.nan
        d["lum"] = float(f"{g[1]}.{g[2]}") if g[1] else np.nan
        d["kick"] = float(g[3]) if g[3] else np.nan
        d["alfam"] = float(f"{g[4]}.{g[5]}") if g[4] else np.nan
        d["alfat"] = float(f"0.{g[7]}") if g[7] else np.nan
        d["gammar"] = float(f"{g[8]}.{g[9]}") if g[8] else np.nan
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--obs", default="FN_Lyr_KIC6936115_kepler_binned.csv")
    ap.add_argument("--min-bytes", type=int, default=1_000_000,
                    help="ignore history.data smaller than this (stubs)")
    ap.add_argument("--root", default=".")
    ap.add_argument("--out-csv", default="fn_lyr_campaign_scores.csv")
    ap.add_argument("--out-fig", default="fn_lyr_campaign_tradeoff.pdf")
    ap.add_argument("--simplex-log", default="optimize_log.csv")
    args = ap.parse_args()

    obs = pd.read_csv(args.obs)
    obs = obs[np.isfinite(obs["phase"]) & np.isfinite(obs["mag"])].sort_values("phase")
    obs_phase = obs["phase"].to_numpy(float)
    obs_mag = obs["mag"].to_numpy(float)

    rows = []
    histories = [p for p in Path(args.root).rglob("history.data")
                 if p.stat().st_size >= args.min_bytes]
    print(f"Found {len(histories)} real history.data files (>= {args.min_bytes} bytes)\n")

    for p in sorted(histories):
        tag = str(p.parent.relative_to(args.root))
        try:
            h = read_history(p)
            res = last_full_cycle(h)
            if res is None:
                print(f"  skip (no scorable cycle): {tag}")
                continue
            mph, mmg, period, vmax, cyc = res
            amp = float(np.nanmax(mmg) - np.nanmin(mmg))
            r, rms, shift, off = align_score(mph, mmg, obs_phase, obs_mag)
            shp, shp_rise, shp_desc = shape_residual(mph, mmg, obs_phase, obs_mag, shift)
            rt = risetime(mph, mmg)
            row = dict(source="grid/sweep", tag=tag, cycle=cyc,
                       period=period, period_err=abs(period - P_OBS) / P_OBS,
                       amp=amp, amp_err=abs(amp - A_OBS) / A_OBS,
                       r=r, rms=rms,
                       shape_resid=shp, shape_rise=shp_rise, shape_desc=shp_desc,
                       risetime=rt, rise_err=abs(rt - RT_OBS),
                       vmax_cs=vmax)
            row.update(parse_params(p))
            rows.append(row)
            print(f"  scored {tag:55s} P={period:.4f} "
                  f"({100*row['period_err']:.1f}%) A={amp:.3f} "
                  f"shape={shp:.3f} (rise {shp_rise:.3f}/desc {shp_desc:.3f}) "
                  f"vmax/cs={vmax:.2f}")
        except Exception as exc:
            print(f"  ERROR {tag}: {exc!r}")

    df = pd.DataFrame(rows)

    # ---- fold in simplex evals from optimize_log.csv (summary only) --------
    slog = Path(args.simplex_log)
    if slog.exists():
        s = pd.read_csv(slog)
        s = s[s["status"] == "ok"].copy()
        for col in ("period", "amp"):
            s[col] = pd.to_numeric(s[col], errors="coerce")
        s_rows = pd.DataFrame(dict(
            source="simplex",
            tag="simplex_eval_" + s["eval"].astype(str),
            period=s["period"],
            period_err=(s["period"] - P_OBS).abs() / P_OBS,
            amp=s["amp"], amp_err=(s["amp"] - A_OBS).abs() / A_OBS,
            r=np.nan, rms=pd.to_numeric(s["rms"], errors="coerce"),
            shape_resid=np.nan, shape_rise=np.nan, shape_desc=np.nan,
            risetime=np.nan, rise_err=np.nan,
            vmax_cs=pd.to_numeric(s["vmax_cs_last"], errors="coerce"),
            teff=s["teff"], lum=s["lum"], alfam=s["alfam"], mass=s["mass"]))
        df = pd.concat([df, s_rows], ignore_index=True)
        print(f"\nFolded in {len(s_rows)} simplex evals from {slog}")

    if df.empty:
        print("No models scored. Check --min-bytes / paths.")
        return

    # ---- combined scores ---------------------------------------------------
    # Morphology = amplitude-normalized shape residual against the REAL FN Lyr
    # curve (lower = better). This compares shape independent of amplitude, so
    # flat/dead models no longer score as "good morphology", and the model is
    # only penalized where its shape DIFFERS from the observed star (e.g. an
    # over-produced descending-branch bump or an over-sharp rise).
    df["morph_score"] = df["shape_resid"]
    # Overall: equal weight period / amplitude / morphology (per decision).
    # Normalizations chosen so a "good" model scores ~1 on each axis:
    #   2% period, 5% amplitude, 0.06 normalized-shape-residual.
    pe = df["period_err"] / 0.02
    ae = df["amp_err"] / 0.05
    me = df["morph_score"] / 0.06
    df["overall"] = pe + ae + me

    df = df.sort_values("overall").reset_index(drop=True)
    df.to_csv(args.out_csv, index=False)

    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 30)
    show = ["source", "tag", "period_err", "amp_err",
            "shape_resid", "shape_rise", "shape_desc", "vmax_cs", "overall"]
    print("\n=== TOP 15 (equal-weight period / amplitude / morphology) ===")
    print(df[show].head(15).to_string(index=False))
    print(f"\nFull table -> {args.out_csv}")

    # ---- trade-off plot ----------------------------------------------------
    g = df[df["source"] != "simplex"].copy()  # morphology only exists here
    fig, ax = plt.subplots(figsize=(7, 5.5))
    sc = ax.scatter(100 * g["period_err"], 100 * g["amp_err"],
                    c=g["shape_resid"], s=70, cmap="viridis_r",
                    edgecolor="k", linewidth=0.4)
    # simplex evals (no morphology) as small grey points for context
    si = df[df["source"] == "simplex"]
    ax.scatter(100 * si["period_err"], 100 * si["amp_err"],
               s=18, color="0.6", marker="x", label="simplex evals (period/amp only)")
    ax.axvline(2.0, ls=":", color="0.4", lw=0.8)
    ax.axhline(5.0, ls=":", color="0.4", lw=0.8)
    ax.set_xlabel("period error (\\%)")
    ax.set_ylabel("amplitude error (\\%)")
    cb = fig.colorbar(sc)
    cb.set_label("morphology score (lower = better)")
    ax.legend(loc="upper right", fontsize=8, frameon=False)
    ax.set_title("FN Lyr campaign: period / amplitude / morphology trade-off")
    fig.tight_layout()
    fig.savefig(args.out_fig, bbox_inches="tight")
    fig.savefig(str(Path(args.out_fig).with_suffix(".png")), dpi=200, bbox_inches="tight")
    print(f"Trade-off plot -> {args.out_fig}")


if __name__ == "__main__":
    main()
