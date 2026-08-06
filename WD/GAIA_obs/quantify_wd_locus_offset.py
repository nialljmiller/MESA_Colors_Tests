#!/usr/bin/env python3

from pathlib import Path
import argparse

import numpy as np
import pandas as pd


def read_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".parquet":
        df = pd.read_parquet(path)
    elif path.suffix.lower() in {".csv", ".txt"}:
        df = pd.read_csv(path)
    else:
        raise ValueError(f"Unsupported file format: {path}")

    df = df.copy()
    df = df.replace([np.inf, -np.inf], np.nan)
    return df


def read_gaia_table(path: Path) -> pd.DataFrame:
    df = read_table(path)

    if "bp_rp" not in df.columns:
        raise KeyError("Gaia table must contain 'bp_rp'.")

    if "abs_g" not in df.columns:
        if {"phot_g_mean_mag", "parallax"}.issubset(df.columns):
            # Gaia parallax is in mas: M_G = G + 5 log10(parallax_mas) - 10
            df["abs_g"] = df["phot_g_mean_mag"] + 5.0 * np.log10(df["parallax"]) - 10.0
        else:
            raise KeyError(
                "Gaia table must contain 'abs_g', or both 'phot_g_mean_mag' and 'parallax'."
            )

    df = df.dropna(subset=["bp_rp", "abs_g"])
    return df


def read_mesa_history(path: Path) -> pd.DataFrame:
    lines = path.read_text(errors="ignore").splitlines()

    header_idx = None
    for i, line in enumerate(lines):
        cols = line.split()
        if "model_number" in cols:
            header_idx = i

    if header_idx is None:
        raise RuntimeError(f"Could not find MESA history header in {path}")

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

    if not rows:
        raise RuntimeError(f"No numeric history rows read from {path}")

    return pd.DataFrame(rows, columns=cols)


def prepare_mesa_track(path: Path) -> pd.DataFrame:
    df = read_mesa_history(path)

    missing = {"G", "Gbp", "Grp"} - set(df.columns)
    if missing:
        raise KeyError(f"MESA history missing required Gaia columns: {sorted(missing)}")

    df = df.replace([np.inf, -np.inf], np.nan)
    df["bp_rp"] = df["Gbp"] - df["Grp"]
    df["M_G"] = df["G"]
    df = df.dropna(subset=["bp_rp", "M_G"])

    if len(df) < 2:
        raise RuntimeError("Need at least two finite MESA WD track points for interpolation.")

    if "model_number" in df.columns:
        df = df.sort_values("model_number")
    elif "star_age" in df.columns:
        df = df.sort_values("star_age")

    return df


def select_gaia_wds(
    gaia: pd.DataFrame,
    wd_prob_min: float,
    color_min: float,
    color_max: float,
    mg_min: float,
    mg_max: float,
) -> pd.DataFrame:
    df = gaia.copy()

    base = (
        np.isfinite(df["bp_rp"])
        & np.isfinite(df["abs_g"])
        & (df["bp_rp"] >= color_min)
        & (df["bp_rp"] <= color_max)
        & (df["abs_g"] >= mg_min)
        & (df["abs_g"] <= mg_max)
    )

    # Preferred path for your gaia_cmd_all_plus_dsc_wds.parquet table:
    # use the Gaia DR3 DSC white-dwarf probability where available.
    if "wd_prob" in df.columns and df["wd_prob"].notna().any():
        good = base & (df["wd_prob"] >= wd_prob_min)
        selected = df.loc[good].copy()
        selection = f"wd_prob >= {wd_prob_min:g} plus CMD window"
    else:
        # Fallback for a plain CMD table. This is not as clean as the DSC selection,
        # but it lets the script run on older gaia_cmd_all.parquet products.
        selected = df.loc[base].copy()
        selection = "CMD window only; no wd_prob column present"

    if selected.empty:
        raise RuntimeError(
            "No Gaia WD points survived the selection. "
            "Check --gaia, --wd-prob-min, and the CMD bounds."
        )

    selected.attrs["selection"] = selection
    return selected


def robust_bin_median(values: np.ndarray, nsigma: float | None) -> tuple[float, float, int]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]

    if len(values) == 0:
        return np.nan, np.nan, 0

    if nsigma is not None and len(values) >= 5:
        med0 = np.nanmedian(values)
        mad0 = 1.4826 * np.nanmedian(np.abs(values - med0))
        if np.isfinite(mad0) and mad0 > 0:
            keep = np.abs(values - med0) <= nsigma * mad0
            values = values[keep]

    med = np.nanmedian(values)
    mad = 1.4826 * np.nanmedian(np.abs(values - med)) if len(values) > 0 else np.nan
    return med, mad, len(values)


def make_model_interpolator(mesa: pd.DataFrame):
    # np.interp needs a single-valued function with ascending x.
    # WD cooling should be effectively single-valued in M_G versus color, but
    # we average duplicate M_G values after rounding to avoid repeated x entries.
    tmp = mesa[["M_G", "bp_rp"]].copy()
    tmp = tmp.dropna().sort_values("M_G")
    tmp["M_G_round"] = tmp["M_G"].round(6)
    tmp = tmp.groupby("M_G_round", as_index=False).agg({"M_G": "mean", "bp_rp": "mean"})
    tmp = tmp.sort_values("M_G")

    model_mg = tmp["M_G"].to_numpy(dtype=float)
    model_color = tmp["bp_rp"].to_numpy(dtype=float)

    if len(model_mg) < 2:
        raise RuntimeError("MESA track has too few unique M_G values for interpolation.")

    def interp_color(mg_values: np.ndarray) -> np.ndarray:
        mg_values = np.asarray(mg_values, dtype=float)
        out = np.full_like(mg_values, np.nan, dtype=float)
        inside = (mg_values >= model_mg.min()) & (mg_values <= model_mg.max())
        out[inside] = np.interp(mg_values[inside], model_mg, model_color)
        return out

    return interp_color, model_mg.min(), model_mg.max()


def quantify_offsets(
    gaia_wd: pd.DataFrame,
    mesa: pd.DataFrame,
    mg_min: float,
    mg_max: float,
    bin_width: float,
    min_per_bin: int,
    model_at: str,
    sigma_clip: float | None,
) -> tuple[pd.DataFrame, dict]:
    interp_color, model_mg_min, model_mg_max = make_model_interpolator(mesa)

    edges = np.arange(mg_min, mg_max + 0.5 * bin_width, bin_width)
    if edges[-1] < mg_max:
        edges = np.append(edges, mg_max)

    rows = []

    for lo, hi in zip(edges[:-1], edges[1:]):
        in_bin = (gaia_wd["abs_g"] >= lo) & (gaia_wd["abs_g"] < hi)
        sub = gaia_wd.loc[in_bin]

        if len(sub) < min_per_bin:
            continue

        color_med, color_mad, n_used = robust_bin_median(
            sub["bp_rp"].to_numpy(dtype=float), sigma_clip
        )

        if n_used < min_per_bin or not np.isfinite(color_med):
            continue

        mg_mid = 0.5 * (lo + hi)
        mg_med = float(np.nanmedian(sub["abs_g"].to_numpy(dtype=float)))
        mg_ref = mg_mid if model_at == "bin_center" else mg_med

        model_col = float(interp_color(np.array([mg_ref]))[0])
        if not np.isfinite(model_col):
            continue

        delta = model_col - color_med

        rows.append(
            {
                "M_G_lo": lo,
                "M_G_hi": hi,
                "M_G_mid": mg_mid,
                "M_G_median_gaia": mg_med,
                "M_G_reference": mg_ref,
                "N_gaia": int(len(sub)),
                "N_gaia_used": int(n_used),
                "gaia_bp_rp_median": color_med,
                "gaia_bp_rp_mad": color_mad,
                "model_bp_rp": model_col,
                "delta_bp_rp_model_minus_gaia": delta,
                "abs_delta_bp_rp": abs(delta),
            }
        )

    table = pd.DataFrame(rows)

    if table.empty:
        raise RuntimeError(
            "No bins could be compared. "
            f"MESA model M_G range is {model_mg_min:.3f}--{model_mg_max:.3f}; "
            "try changing --mg-min/--mg-max or --min-per-bin."
        )

    delta = table["delta_bp_rp_model_minus_gaia"].to_numpy(dtype=float)
    abs_delta = np.abs(delta)
    n_weights = table["N_gaia_used"].to_numpy(dtype=float)

    summary = {
        "n_bins": int(len(table)),
        "n_gaia_total_selected": int(len(gaia_wd)),
        "n_gaia_used_in_bins": int(table["N_gaia_used"].sum()),
        "mg_min_compared": float(table["M_G_lo"].min()),
        "mg_max_compared": float(table["M_G_hi"].max()),
        "model_mg_min": float(model_mg_min),
        "model_mg_max": float(model_mg_max),
        "median_abs_delta_bp_rp": float(np.nanmedian(abs_delta)),
        "mean_abs_delta_bp_rp": float(np.nanmean(abs_delta)),
        "rms_delta_bp_rp": float(np.sqrt(np.nanmean(delta**2))),
        "max_abs_delta_bp_rp": float(np.nanmax(abs_delta)),
        "weighted_mean_abs_delta_bp_rp": float(np.average(abs_delta, weights=n_weights)),
        "weighted_rms_delta_bp_rp": float(np.sqrt(np.average(delta**2, weights=n_weights))),
    }

    return table, summary


def write_summary(path: Path, args, gaia_wd: pd.DataFrame, table: pd.DataFrame, summary: dict):
    path.parent.mkdir(parents=True, exist_ok=True)

    selection = gaia_wd.attrs.get("selection", "unknown")

    lines = []
    lines.append("White dwarf Gaia-locus comparison")
    lines.append("===================================")
    lines.append("")
    lines.append(f"Gaia input:       {args.gaia}")
    lines.append(f"MESA history:     {args.history}")
    lines.append(f"Gaia selection:   {selection}")
    lines.append(f"Model reference:  {args.model_at.replace('_', ' ')}")
    lines.append(f"Magnitude bins:   {args.mg_min:g} <= M_G <= {args.mg_max:g}, width {args.bin_width:g} mag")
    lines.append("")
    lines.append(f"Bins compared:                    {summary['n_bins']}")
    lines.append(f"Gaia WD points selected:          {summary['n_gaia_total_selected']:,}")
    lines.append(f"Gaia WD points used in bins:      {summary['n_gaia_used_in_bins']:,}")
    lines.append(f"Compared M_G range:               {summary['mg_min_compared']:.2f}--{summary['mg_max_compared']:.2f}")
    lines.append(f"Model M_G range:                  {summary['model_mg_min']:.3f}--{summary['model_mg_max']:.3f}")
    lines.append("")
    lines.append(f"median |Delta(G_BP-G_RP)|:        {summary['median_abs_delta_bp_rp']:.4f} mag")
    lines.append(f"mean |Delta(G_BP-G_RP)|:          {summary['mean_abs_delta_bp_rp']:.4f} mag")
    lines.append(f"RMS Delta(G_BP-G_RP):             {summary['rms_delta_bp_rp']:.4f} mag")
    lines.append(f"max |Delta(G_BP-G_RP)|:           {summary['max_abs_delta_bp_rp']:.4f} mag")
    lines.append(f"weighted mean |Delta|:            {summary['weighted_mean_abs_delta_bp_rp']:.4f} mag")
    lines.append(f"weighted RMS Delta:               {summary['weighted_rms_delta_bp_rp']:.4f} mag")
    lines.append("")
    lines.append("Paper-ready sentence:")
    lines.append(
        "Across $M_G = "
        f"{summary['mg_min_compared']:.0f}$--${summary['mg_max_compared']:.0f}$, "
        "the binned Gaia white-dwarf locus and the $0.6\\,M_\\odot$ DA Colors track "
        "differ by a median absolute color offset of "
        f"${summary['median_abs_delta_bp_rp']:.3f}$ mag, "
        "with an RMS offset of "
        f"${summary['rms_delta_bp_rp']:.3f}$ mag and a maximum bin-wise offset of "
        f"${summary['max_abs_delta_bp_rp']:.3f}$ mag."
    )
    lines.append("")
    lines.append("Bin table:")
    lines.append(table.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    lines.append("")

    path.write_text("\n".join(lines))


def make_plot(path: Path, table: pd.DataFrame, mesa: pd.DataFrame):
    import matplotlib.pyplot as plt

    path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(6.2, 5.2))

    ax.plot(mesa["bp_rp"], mesa["M_G"], lw=2.0, label="MESA Colors 0.6 $M_\\odot$ DA")
    ax.errorbar(
        table["gaia_bp_rp_median"],
        table["M_G_reference"],
        xerr=table["gaia_bp_rp_mad"],
        fmt="o",
        ms=4,
        capsize=2,
        label="Binned Gaia WD locus",
    )

    ax.set_xlabel(r"$G_{\rm BP} - G_{\rm RP}$")
    ax.set_ylabel(r"$M_G$")
    ax.invert_yaxis()
    ax.legend(frameon=True, fontsize=9)
    ax.minorticks_on()
    fig.tight_layout()
    fig.savefig(path.with_suffix(".pdf"))
    fig.savefig(path.with_suffix(".png"), dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Quantify the color offset between a MESA Colors white-dwarf cooling track "
            "and the empirical Gaia white-dwarf locus in bins of M_G."
        )
    )

    parser.add_argument("--gaia", default="gaia_cmd_all_plus_dsc_wds.parquet")
    parser.add_argument("--history", default="LOGS_wd_cool/history.data")
    parser.add_argument("--out", default="wd_locus_offsets.csv")
    parser.add_argument("--summary", default="wd_locus_offsets_summary.txt")
    parser.add_argument("--plot", default=None, help="Optional output base for a diagnostic plot, e.g. figures/wd_locus_offset")

    parser.add_argument("--mg-min", type=float, default=9.0)
    parser.add_argument("--mg-max", type=float, default=15.0)
    parser.add_argument("--bin-width", type=float, default=1.0)
    parser.add_argument("--min-per-bin", type=int, default=50)

    parser.add_argument("--wd-prob-min", type=float, default=0.5)
    parser.add_argument("--color-min", type=float, default=-0.8)
    parser.add_argument("--color-max", type=float, default=1.3)
    parser.add_argument(
        "--model-at",
        choices=["bin_center", "median_mg"],
        default="bin_center",
        help="M_G value at which to interpolate the model in each bin.",
    )
    parser.add_argument(
        "--sigma-clip",
        type=float,
        default=None,
        help="Optional robust sigma clipping of Gaia colors within each M_G bin, e.g. 3.",
    )

    args = parser.parse_args()

    gaia = read_gaia_table(Path(args.gaia))
    mesa = prepare_mesa_track(Path(args.history))
    gaia_wd = select_gaia_wds(
        gaia,
        wd_prob_min=args.wd_prob_min,
        color_min=args.color_min,
        color_max=args.color_max,
        mg_min=args.mg_min,
        mg_max=args.mg_max,
    )

    table, summary = quantify_offsets(
        gaia_wd=gaia_wd,
        mesa=mesa,
        mg_min=args.mg_min,
        mg_max=args.mg_max,
        bin_width=args.bin_width,
        min_per_bin=args.min_per_bin,
        model_at=args.model_at,
        sigma_clip=args.sigma_clip,
    )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out, index=False)

    summary_path = Path(args.summary)
    write_summary(summary_path, args, gaia_wd, table, summary)

    if args.plot is not None:
        make_plot(Path(args.plot), table, mesa)

    print()
    print(f"Wrote bin table: {out}")
    print(f"Wrote summary:   {summary_path}")
    if args.plot is not None:
        print(f"Wrote plot:      {Path(args.plot).with_suffix('.pdf')}")
        print(f"Wrote plot:      {Path(args.plot).with_suffix('.png')}")

    print()
    print("Summary:")
    print(f"  bins compared:             {summary['n_bins']}")
    print(f"  Gaia WD points selected:   {summary['n_gaia_total_selected']:,}")
    print(f"  median |Delta color|:      {summary['median_abs_delta_bp_rp']:.4f} mag")
    print(f"  RMS Delta color:           {summary['rms_delta_bp_rp']:.4f} mag")
    print(f"  max |Delta color|:         {summary['max_abs_delta_bp_rp']:.4f} mag")
    print()
    print("Paper-ready sentence:")
    print(
        "Across $M_G = "
        f"{summary['mg_min_compared']:.0f}$--${summary['mg_max_compared']:.0f}$, "
        "the binned Gaia white-dwarf locus and the $0.6\\,M_\\odot$ DA Colors track "
        "differ by a median absolute color offset of "
        f"${summary['median_abs_delta_bp_rp']:.3f}$ mag, "
        "with an RMS offset of "
        f"${summary['rms_delta_bp_rp']:.3f}$ mag and a maximum bin-wise offset of "
        f"${summary['max_abs_delta_bp_rp']:.3f}$ mag."
    )


if __name__ == "__main__":
    main()
