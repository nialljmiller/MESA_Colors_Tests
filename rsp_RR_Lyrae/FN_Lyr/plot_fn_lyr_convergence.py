#!/usr/bin/env python3
"""
Appendix figure: convergence of the Nelder-Mead fit of the MESA RSP+Colors
model to the FN Lyr Kepler light curve.

Reads optimize_log.csv (written by optimize_fn_lyr.py) and shows, as a
function of evaluation number, how the optimizer drives the synthetic light
curve onto the observed one. This is the figure that demonstrates Colors
closes the loop between a nonlinear pulsation model and real photometry.

Layout (5 stacked panels, shared x = evaluation number):
  (a) objective and running-best objective
  (b) waveform RMS residual (the thing actually minimized), with converged
      vs rejected (crash/unconverged) evals marked
  (c) period error (%), with the +/-1.5% tolerance band
  (d) amplitude fraction A_model / A_obs, with the observed value at 1.0
  (e) the four free parameters, each normalized within its bound so they
      share one axis; the Nemec measured central values are marked where
      defined (Teff, L, mass) so you can see the fit settle inside the
      measured box.

Usage:
    python plot_fn_lyr_convergence.py \
        --log optimize_log.csv --out figures/fig_fn_lyr_convergence.pdf
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


# Bounds + Nemec centres must match optimize_fn_lyr.py.
BOUNDS = {
    "teff": (6470.0, 6570.0),
    "lum":  (40.0, 56.0),
    "alfam": (0.08, 0.30),
    "mass": (0.59, 0.69),
}
NEMEC = {"teff": 6520.0, "lum": 48.0, "mass": 0.64}  # alfam not measured
PERIOD_TOL = 1.5   # percent
OBS_AMP = 1.0464


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--log", default="optimize_log.csv")
    ap.add_argument("--out", default="figures/fig_fn_lyr_convergence.pdf")
    args = ap.parse_args()

    plt.rcParams.update({
        "font.family": "serif",
        "mathtext.fontset": "dejavuserif",
        "font.size": 9,
        "axes.linewidth": 0.8,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "legend.frameon": False,
        "figure.dpi": 150,
    })

    df = pd.read_csv(args.log)
    df = df.sort_values("eval")
    ev = df["eval"].to_numpy(float)

    ok = df["status"] == "ok"
    bad = ~ok

    # numeric coercion (rejected rows have blank rms/amp)
    def num(col):
        return pd.to_numeric(df[col], errors="coerce").to_numpy(float)

    obj = num("objective")
    rms = num("rms")
    perr = num("period_err_pct")
    ampf = num("amp_frac")

    running_best = np.fmin.accumulate(np.where(np.isfinite(obj), obj, np.inf))

    fig = plt.figure(figsize=(7.0, 9.2))
    gs = GridSpec(4, 1, height_ratios=[1, 1, 1, 1.3], hspace=0.0, figure=fig)
    axes = [fig.add_subplot(gs[i]) for i in range(4)]
    for a in axes[:-1]:
        plt.setp(a.get_xticklabels(), visible=False)

    blue, orange, green, grey = "#3b6ea5", "#d1622b", "#3a8a5f", "0.6"

    # (a) objective + running best
    ax = axes[0]
    ax.plot(ev, obj, "o", ms=3.2, color=blue, alpha=0.7, label="evaluation")
    ax.step(ev, running_best, where="post", color=orange, lw=1.5, label="running best")
    ax.set_ylabel("objective")

    leg = ax.legend(
        loc="upper right",
        fontsize=7.5,
        ncol=1,
        frameon=True,
        fancybox=False,
        framealpha=0.95,
        facecolor="white",
        edgecolor="0.35",
        borderpad=0.35,
        handlelength=1.6,
        handletextpad=0.5,
    )
    
    leg.get_frame().set_linewidth(0.5)
    leg.set_zorder(10)

    ax.text(0.012, 0.9, "(a)", transform=ax.transAxes, va="top", fontsize=9)

    # (b) RMS, converged vs rejected
    ax = axes[1]
    ax.plot(ev[ok.values], rms[ok.values], "o", ms=3.4, color=blue, label="converged")
    # rejected evals: show at the top as markers (no RMS)
    if bad.any():
        ymark = np.nanmax(rms[np.isfinite(rms)]) if np.isfinite(rms).any() else 1.0
        ax.plot(ev[bad.values], np.full(bad.sum(), ymark), "x", ms=4,
                color=grey, label="rejected (crash/unconverged)")
    ax.set_ylabel("waveform RMS\n(mag)")
    ax.text(0.012, 0.9, "(b)", transform=ax.transAxes, va="top", fontsize=9)


    # (d) amplitude fraction
    ax = axes[2]
    ax.axhline(1.0, lw=0.8, color=orange, label="$A_{\\rm obs}$")
    ax.plot(ev[ok.values], ampf[ok.values], "o", ms=3.2, color=blue)
    ax.set_ylabel("$A_{\\rm mod}/A_{\\rm obs}$")
    ax.legend(loc="lower right", fontsize=7.5)
    ax.text(0.012, 0.9, "(d)", transform=ax.transAxes, va="top", fontsize=9)

    # (e) normalized parameters
    ax = axes[3]
    pcolors = {"teff": "#b5482a", "lum": "#3a8a5f", "alfam": "#7a4fa3", "mass": "#2a6fb5"}
    plabels = {"teff": "$T_{\\rm eff}$", "lum": "$L$",
               "alfam": r"$\alpha_{\rm m}$", "mass": "$M$"}
    for name, (lo, hi) in BOUNDS.items():
        vals = (num(name) - lo) / (hi - lo)
        ax.plot(ev, vals, "-o", ms=2.6, lw=1.2, color=pcolors[name],
                label=plabels[name], alpha=0.9)
        if name in NEMEC:
            nv = (NEMEC[name] - lo) / (hi - lo)
            ax.axhline(nv, color=pcolors[name], lw=1, ls=":", alpha=1)
    ax.set_ylim(-0.05, 1.05)
    ax.set_ylabel("parameter\n(normalized in bounds)")
    ax.set_xlabel("evaluation number")
    ax.legend(loc="upper right", fontsize=7.5, ncol=4, handlelength=1.4)
    ax.text(0.012, 0.9, "(e)", transform=ax.transAxes, va="top", fontsize=9)
    ax.text(0.99, 0.02,
            "dotted lines: Nemec et al. (2011) central values",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.8, color="0.4")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    png = str(Path(args.out).with_suffix(".png"))
    fig.savefig(png, bbox_inches="tight", dpi=200)

    nbest = int(df.loc[df["objective"].idxmin(), "eval"]) if np.isfinite(obj).any() else -1
    print(f"evals: {len(df)}  converged: {int(ok.sum())}  rejected: {int(bad.sum())}")
    print(f"best eval: {nbest}  best objective: {np.nanmin(obj):.5f}")
    print(f"wrote {args.out} and {png}")


if __name__ == "__main__":
    main()
