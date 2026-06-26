#!/usr/bin/env python3

from pathlib import Path
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
from matplotlib.gridspec import GridSpec


def read_gaia_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".parquet":
        df = pd.read_parquet(path)
    elif path.suffix.lower() in {".csv", ".txt"}:
        df = pd.read_csv(path)
    else:
        raise ValueError(f"Unsupported Gaia file format: {path}")

    df = df.copy()

    if "bp_rp" not in df.columns:
        raise KeyError("Gaia table must contain 'bp_rp'.")

    if "abs_g" not in df.columns:
        if {"phot_g_mean_mag", "parallax"}.issubset(df.columns):
            df["abs_g"] = df["phot_g_mean_mag"] + 5.0 * np.log10(df["parallax"]) - 10.0
        else:
            raise KeyError("Gaia table must contain 'abs_g', or both 'phot_g_mean_mag' and 'parallax'.")

    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=["bp_rp", "abs_g"])

    return df


def read_mesa_history(path: Path) -> pd.DataFrame:
    lines = path.read_text().splitlines()

    header_idx = None
    for i, line in enumerate(lines):
        cols = line.split()
        if "model_number" in cols:
            header_idx = i

    if header_idx is None:
        raise RuntimeError(f"Could not find history header in {path}")

    cols = lines[header_idx].split()
    rows = []

    for line in lines[header_idx + 1:]:
        parts = line.split()
        if len(parts) != len(cols):
            continue
        try:
            rows.append([float(x) for x in parts])
        except ValueError:
            continue

    return pd.DataFrame(rows, columns=cols)


def prepare_mesa_track(path: Path) -> pd.DataFrame:
    df = read_mesa_history(path)

    missing = {"G", "Gbp", "Grp"} - set(df.columns)
    if missing:
        raise KeyError(f"MESA history missing columns: {missing}")

    df = df.replace([np.inf, -np.inf], np.nan)
    df["bp_rp"] = df["Gbp"] - df["Grp"]
    df["M_G"] = df["G"]
    df = df.dropna(subset=["bp_rp", "M_G"])

    if "model_number" in df.columns:
        df = df.sort_values("model_number")
    elif "star_age" in df.columns:
        df = df.sort_values("star_age")

    return df


def density_hist(x, y, xlim, ylim, dx, dy):
    good = (
        np.isfinite(x)
        & np.isfinite(y)
        & (x >= xlim[0])
        & (x <= xlim[1])
        & (y >= ylim[0])
        & (y <= ylim[1])
    )

    x = x[good]
    y = y[good]

    xedges = np.arange(xlim[0], xlim[1] + dx, dx)
    yedges = np.arange(ylim[0], ylim[1] + dy, dy)

    H, _, _ = np.histogram2d(x, y, bins=[xedges, yedges])
    H = H.T
    H = np.ma.masked_where(H <= 0, H)

    return H, xedges, yedges


def draw_density(ax, H, xedges, yedges, norm):
    return ax.pcolormesh(
        xedges,
        yedges,
        H,
        cmap="Greys",
        norm=norm,
        shading="auto",
        rasterized=True,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gaia", default="gaia_cmd_all.parquet")
    parser.add_argument("--history", default="LOGS_wd_cool/history.data")
    parser.add_argument("--out", default="figures/gaia_wd_cmd_mesa_overlay_apj")
    args = parser.parse_args()

    gaia = read_gaia_table(Path(args.gaia))
    mesa = prepare_mesa_track(Path(args.history))

    outbase = Path(args.out)
    outbase.parent.mkdir(parents=True, exist_ok=True)

    full_xlim = (-1.0, 5.0)
    full_ylim = (-5.4, 18.0)

    wd_xlim = (-0.8, 1.3)
    wd_ylim = (7.4, 16.0)

    x = gaia["bp_rp"].to_numpy()
    y = gaia["abs_g"].to_numpy()

    H_full, xe_full, ye_full = density_hist(
        x, y, full_xlim, full_ylim, dx=0.020, dy=0.060
    )

    H_wd, xe_wd, ye_wd = density_hist(
        x, y, wd_xlim, wd_ylim, dx=0.015, dy=0.04
    )

    counts = np.concatenate(
        [
            np.asarray(H_full.compressed(), dtype=float),
            np.asarray(H_wd.compressed(), dtype=float),
        ]
    )

    if len(counts) == 0:
        raise RuntimeError("No Gaia points survived plotting cuts.")

    vmax = np.nanpercentile(counts, 96.8)
    vmax = max(vmax, 2.0)
    norm = LogNorm(vmin=0.1, vmax=vmax)

    plt.rcParams.update(
        {
            "font.size": 16,
            "axes.linewidth": 1.0,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.major.size": 5,
            "ytick.major.size": 5,
            "xtick.minor.size": 2.5,
            "ytick.minor.size": 2.5,
        }
    )

    fig = plt.figure(figsize=(13.8, 6.4))

    outer = GridSpec(
        1,
        2,
        width_ratios=[1.0, 1.12],
        wspace=0.15,
        left=0.065,
        right=0.93,
        bottom=0.12,
        top=0.955,
        figure=fig,
    )

    ax1 = fig.add_subplot(outer[0, 0])

    right = outer[0, 1].subgridspec(
        1,
        2,
        width_ratios=[1.0, 0.055],
        wspace=0.0,
    )

    ax2 = fig.add_subplot(right[0, 0])
    cax = fig.add_subplot(right[0, 1])

    m1 = draw_density(ax1, H_full, xe_full, ye_full, norm)
    m2 = draw_density(ax2, H_wd, xe_wd, ye_wd, norm)

    track_color = "#0072B2"


    ax1.plot(
        mesa["bp_rp"],
        mesa["M_G"],
        color=track_color,
        lw=2,
        zorder=10,
        #solid_capstyle="round",
    )

    ax2.plot(
        mesa["bp_rp"],
        mesa["M_G"],
        color=track_color,
        lw=2.4,
        zorder=10,
        #solid_capstyle="round",
        label="MESA Colors WD cooling track",
    )

    ax2.scatter(
        [mesa["bp_rp"].iloc[0], mesa["bp_rp"].iloc[-1]],
        [mesa["M_G"].iloc[0], mesa["M_G"].iloc[-1]],
        s=30,
        facecolor=track_color,
        linewidth=0.9,
        zorder=11,
    )

    rect = Rectangle(
        (wd_xlim[0], wd_ylim[0]),
        wd_xlim[1] - wd_xlim[0],
        wd_ylim[1] - wd_ylim[0],
        fill=False,
        edgecolor="black",
        linewidth=1.35,
        zorder=12,
    )
    ax1.add_patch(rect)

    ax1.set_xlim(*full_xlim)
    ax1.set_ylim(full_ylim[1], full_ylim[0])

    ax2.set_xlim(*wd_xlim)
    ax2.set_ylim(wd_ylim[1], wd_ylim[0])

    ax1.set_xlabel(r"$G_{\rm BP} - G_{\rm RP}$", fontsize=18)
    ax1.set_ylabel(r"$M_G$", fontsize=18)

    ax2.set_xlabel(r"$G_{\rm BP} - G_{\rm RP}$", fontsize=18)
    ax2.set_ylabel(r"$M_G$", fontsize=18)

    ax1.minorticks_on()
    ax2.minorticks_on()

    ax1.tick_params(axis="both", which="major", labelsize=16, pad=5)
    ax2.tick_params(axis="both", which="major", labelsize=16, pad=5)

    for ax in (ax1, ax2):
        for spine in ax.spines.values():
            spine.set_linewidth(1.0)

    cb = fig.colorbar(m2, cax=cax)
    cb.set_label("Source density", fontsize=18)
    cb.ax.tick_params(labelsize=15)
    cb.outline.set_linewidth(1.0)

    leg = ax2.legend(
        loc="lower left",
        fontsize=14,
        frameon=True,
        framealpha=0.92,
        borderpad=0.35,
        handlelength=3.0,
    )
    leg.get_frame().set_linewidth(0.9)

    fig.savefig(f"{outbase}.png", dpi=500)
    fig.savefig(f"{outbase}.pdf")

    print(f"Gaia points plotted: {len(gaia):,}")
    print(f"wrote {outbase}.png")
    print(f"wrote {outbase}.pdf")


if __name__ == "__main__":
    main()
