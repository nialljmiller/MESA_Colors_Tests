#!/usr/bin/env python3
"""
ApJ-quality FN Lyr light-curve comparison figure.

Two stacked panels sharing the x-axis with no vertical gap:
  - top (2 height units): observed binned Kepler light curve + best-fit
    MESA RSP + MESA Custom Colors model, both over two pulsation phases
  - bottom (1 height unit): residuals (observed - model) over two phases

The model is aligned to the data by the only two physically-arbitrary
degrees of freedom: a phase shift (epoch is arbitrary) and a vertical
magnitude offset (the observed curve is relative Kepler magnitude). No
amplitude rescaling or shape warping is applied -- the shape is set by the
physics of the converged RSP limit cycle.

Usage:
    python plot_fn_lyr_fit.py \
        --history optimize_best/history.data \
        --obs FN_Lyr_KIC6936115_kepler_binned.csv \
        --out figures/fig_fn_lyr_fit.pdf

If --history is omitted it defaults to LOGS/history.data.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


# ---- observed targets (Nemec et al. 2011) ---------------------------------
P_OBS = 0.527398471
OBS_AMP = 1.0464

# Cycle-selection guards (match the comparison script).
MIN_POINTS_PER_CYCLE = 250
MIN_PHASE_SPAN = 0.80
SEARCH_LAST_N_CYCLES = 10
N_SHIFTS = 2500


# ---------------------------------------------------------------------------
# IO + alignment helpers (self-contained; mirror compare_fn_lyr_kepler_absolute)
# ---------------------------------------------------------------------------
def read_history(path: Path) -> pd.DataFrame:
    path = Path(path)
    lines = path.read_text().splitlines()
    header_idx = next(i for i, l in enumerate(lines) if "model_number" in l.split())
    return pd.read_csv(path, sep=r"\s+", skiprows=header_idx)


def get_mag_col(h: pd.DataFrame) -> str:
    for c in ("Kepler", "K"):
        if c in h.columns:
            print(h.columns)
            return c
    cand = [c for c in h.columns if "kepler" in c.lower() or c.lower() in {"kp", "kplr"}]
    if cand:
        return cand[0]
    raise RuntimeError(f"No Kepler magnitude column. Columns: {list(h.columns)}")


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


def get_cycle(h, mag_col, n):
    cyc = h[h["rsp_num_periods"].astype(int) == n]
    ph = np.asarray(cyc["rsp_phase"], float) % 1.0
    mg = np.asarray(cyc[mag_col], float)
    good = np.isfinite(ph) & np.isfinite(mg)
    ph, mg = ph[good], mg[good]
    order = np.argsort(ph)
    return ph[order], mg[order]


def best_alignment(h, mag_col, obs_phase, obs_mag):
    """Pick the late cycle and (phase shift, vertical offset) with lowest RMS."""
    nmax = int(np.nanmax(h["rsp_num_periods"]))
    shifts = np.linspace(0, 1, N_SHIFTS, endpoint=False)
    best = None
    for n in range(max(0, nmax - SEARCH_LAST_N_CYCLES), nmax + 1):
        cyc = h[h["rsp_num_periods"].astype(int) == n]
        if len(cyc) < MIN_POINTS_PER_CYCLE:
            continue
        ph, mg = get_cycle(h, mag_col, n)
        if len(ph) < MIN_POINTS_PER_CYCLE or (np.nanmax(ph) - np.nanmin(ph)) < MIN_PHASE_SPAN:
            continue
        for shift in shifts:
            m_eval = periodic_interp(ph, mg, (obs_phase - shift) % 1.0)
            offset = np.nanmedian(obs_mag - m_eval)
            resid = obs_mag - (m_eval + offset)
            rms = float(np.sqrt(np.nanmean(resid ** 2)))
            if best is None or rms < best["rms"]:
                best = dict(cycle=n, shift=float(shift), offset=float(offset), rms=rms)
    if best is None:
        raise RuntimeError("No scorable cycle found in history.")
    return best


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--history", default="LOGS/history.data")
    ap.add_argument("--obs", default="FN_Lyr_KIC6936115_kepler_binned.csv")
    ap.add_argument("--out", default="figures/fig_fn_lyr_fit.pdf")
    ap.add_argument("--title", default=None, help="optional suptitle")
    args = ap.parse_args()

    # ApJ-ish style: serif, modest line weights, no chartjunk.
    plt.rcParams.update({
        "font.family": "serif",
        "mathtext.fontset": "dejavuserif",
        "font.size": 10,
        "axes.linewidth": 0.8,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "legend.frameon": False,
        "figure.dpi": 150,
        # ApJ requires embedded Type 42 (TrueType) fonts, not Type 3 bitmaps.
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    obs = pd.read_csv(args.obs)
    obs = obs[np.isfinite(obs["phase"]) & np.isfinite(obs["mag"])].sort_values("phase")
    obs_phase = obs["phase"].to_numpy(float)
    obs_mag = obs["mag"].to_numpy(float)

    h = read_history(args.history)
    mag_col = get_mag_col(h)
    al = best_alignment(h, mag_col, obs_phase, obs_mag)

    model_phase, model_mag = get_cycle(h, mag_col, al["cycle"])

    # Model curve over two phases, shifted + offset.
    mp_shift = (model_phase + al["shift"]) % 1.0
    order = np.argsort(mp_shift)
    mp1 = mp_shift[order]
    mm1 = model_mag[order] + al["offset"]
    mp2 = np.concatenate([mp1, mp1 + 1.0])
    mm2 = np.concatenate([mm1, mm1])

    # Observed over two phases.
    op2 = np.concatenate([obs_phase, obs_phase + 1.0])
    om2 = np.concatenate([obs_mag, obs_mag])

    # Residuals (obs - model) at observed phases, one phase, then duplicated.
    m_at_obs = periodic_interp(model_phase, model_mag, (obs_phase - al["shift"]) % 1.0)
    resid = obs_mag - (m_at_obs + al["offset"])
    rp2 = np.concatenate([obs_phase, obs_phase + 1.0])
    rr2 = np.concatenate([resid, resid])

    model_amp = float(np.nanmax(model_mag) - np.nanmin(model_mag))
    model_cycle = h[h["rsp_num_periods"].astype(int) == al["cycle"]]
    model_period = float(np.nanmedian(model_cycle["rsp_period_in_days"]))
    mad = float(np.nanmedian(np.abs(resid)))
    # shape correlation after mean subtraction (no amplitude rescale)
    pred = m_at_obs + al["offset"]
    corr = float(np.corrcoef(obs_mag - np.nanmean(obs_mag),
                             pred - np.nanmean(pred))[0, 1])

    # ---- figure: 2:1 stacked, shared x, zero gap ----
    fig = plt.figure(figsize=(7.0, 5.0))
    gs = GridSpec(2, 1, height_ratios=[2, 1], hspace=0.0, figure=fig)
    ax_lc = fig.add_subplot(gs[0])
    ax_re = fig.add_subplot(gs[1], sharex=ax_lc)

    # Light-curve panel.
    ax_lc.plot(op2, om2, ".", ms=3.0, color="#3b6ea5", alpha=0.85,
               label="FN Lyr / KIC 6936115 (Kepler)", rasterized=True)
    ax_lc.plot(mp2, mm2, "-", lw=1.6, color="#d1622b",
               label="MESA RSP + MESA Custom Colors")
    ax_lc.invert_yaxis()
    ax_lc.set_ylabel("Relative Kepler magnitude")
    ax_lc.legend(loc="upper center", fontsize=8.5, handlelength=1.6)
    plt.setp(ax_lc.get_xticklabels(), visible=False)

    # Annotation block with the headline numbers.
    txt = (f"$P_{{\\rm obs}} = {P_OBS:.6f}$, $P_{{\\rm mod}} = {model_period:.6f}$ d\n"
           f"$A_{{\\rm obs}} = {OBS_AMP:.3f}$, $A_{{\\rm mod}} = {model_amp:.3f}$ mag")
    ax_lc.text(0.015, 0.04, txt, transform=ax_lc.transAxes, fontsize=7.8,
               va="bottom", ha="left",
               bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.7", lw=0.6))

    # Residual panel.
    ax_re.axhline(0.0, lw=0.8, color="0.4")
    ax_re.plot(rp2, rr2, ".", ms=3.0, color="#3b6ea5", alpha=0.85, rasterized=True)
    ax_re.set_xlabel("Pulsation phase")
    ax_re.set_ylabel("O $-$ C")
    rmax = np.nanmax(np.abs(rr2)) * 1.15 if np.isfinite(np.nanmax(np.abs(rr2))) else 0.1
    ax_re.set_ylim(rmax, -rmax)  # inverted to match magnitude convention
    ax_re.set_xlim(0, 2)

    if args.title:
        fig.suptitle(args.title, y=0.98, fontsize=11)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    # also a PNG companion for quick viewing
    png = str(Path(args.out).with_suffix(".png"))
    fig.savefig(png, bbox_inches="tight", dpi=250)

    print(f"cycle={al['cycle']}  shift={al['shift']:.4f}  offset={al['offset']:.4f}")
    print(f"model_period={model_period:.6f}  model_amp={model_amp:.4f}  r={corr:.4f}  MAD={mad:.4f}  RMS={al['rms']:.4f}")
    print(f"wrote {args.out} and {png}")


if __name__ == "__main__":
    main()
