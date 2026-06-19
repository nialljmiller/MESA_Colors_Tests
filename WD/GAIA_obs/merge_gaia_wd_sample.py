#!/usr/bin/env python3

from pathlib import Path
import argparse

import numpy as np
import pandas as pd
from matplotlib.path import Path as MplPath


# Hand-tuned Gaia WD-locus polygon in CMD space:
# x = G_BP - G_RP
# y = M_G
#
# This deliberately excludes most of the rectangular WD-query box and keeps
# only the diagonal WD sequence / WD binary region.
WD_POLYGON = np.array(
    [
        [-0.75, 7.8],
        [-0.55, 9.3],
        [-0.35, 10.1],
        [-0.10, 10.8],
        [0.25, 11.6],
        [0.60, 12.5],
        [0.95, 13.5],
        [1.25, 14.5],
        [1.45, 15.7],

        [0.95, 16.0],
        [0.45, 15.3],
        [0.00, 14.3],
        [-0.35, 13.2],
        [-0.65, 11.8],
        [-0.85, 10.0],
    ],
    dtype=float,
)


def read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)

    if path.suffix.lower() == ".parquet":
        df = pd.read_parquet(path)
    elif path.suffix.lower() in {".csv", ".txt"}:
        df = pd.read_csv(path)
    else:
        raise ValueError(f"Unsupported file type: {path}")

    df = df.copy()
    df = df.replace([np.inf, -np.inf], np.nan)

    if "bp_rp" not in df.columns:
        raise KeyError(f"{path} has no bp_rp column.")

    if "abs_g" not in df.columns:
        if {"phot_g_mean_mag", "parallax"}.issubset(df.columns):
            df["abs_g"] = df["phot_g_mean_mag"] + 5.0 * np.log10(df["parallax"]) - 10.0
        else:
            raise KeyError(
                f"{path} has no abs_g column and cannot compute it "
                "because phot_g_mean_mag/parallax are missing."
            )

    df = df.dropna(subset=["bp_rp", "abs_g"])

    return df


def select_wd_locus(df: pd.DataFrame, polygon: np.ndarray) -> pd.DataFrame:
    x = df["bp_rp"].to_numpy()
    y = df["abs_g"].to_numpy()

    finite = np.isfinite(x) & np.isfinite(y)

    points = np.column_stack([x[finite], y[finite]])

    path = MplPath(polygon)
    inside = np.zeros(len(df), dtype=bool)
    inside[np.where(finite)[0]] = path.contains_points(points)

    return df.loc[inside].copy()


def main():
    parser = argparse.ArgumentParser(
        description="Merge broad Gaia CMD sample with only the WD-locus part of a WD-region sample."
    )

    parser.add_argument(
        "--base",
        default="gaia_cmd_all.parquet",
        help="Broad Gaia CMD file.",
    )

    parser.add_argument(
        "--wd",
        default="gaia_cmd_wd_region.parquet",
        help="Dense WD-region Gaia file.",
    )

    parser.add_argument(
        "--out",
        default="gaia_cmd_all_plus_wd_locus.parquet",
        help="Output merged Gaia CMD file.",
    )

    parser.add_argument(
        "--selected-wd-out",
        default="gaia_cmd_wd_locus_selected.parquet",
        help="Optional output file containing only the selected WD-locus supplement.",
    )

    parser.add_argument(
        "--write-csv",
        action="store_true",
        help="Also write CSV versions.",
    )

    args = parser.parse_args()

    base_path = Path(args.base)
    wd_path = Path(args.wd)
    out_path = Path(args.out)
    selected_wd_out = Path(args.selected_wd_out)

    print(f"Reading base Gaia sample: {base_path}")
    base = read_table(base_path)

    print(f"Reading WD-region Gaia sample: {wd_path}")
    wd_region = read_table(wd_path)

    print()
    print(f"Base rows before merge:      {len(base):,}")
    print(f"WD-region rows before cut:   {len(wd_region):,}")

    wd_selected = select_wd_locus(wd_region, WD_POLYGON)

    print(f"WD-locus rows after cut:     {len(wd_selected):,}")

    wd_selected.to_parquet(selected_wd_out, index=False)
    print(f"Wrote selected WD locus:     {selected_wd_out}")

    if args.write_csv:
        selected_csv = selected_wd_out.with_suffix(".csv")
        wd_selected.to_csv(selected_csv, index=False)
        print(f"Wrote selected WD CSV:       {selected_csv}")

    # Make sure both tables have the same columns.
    all_columns = sorted(set(base.columns) | set(wd_selected.columns))
    base = base.reindex(columns=all_columns)
    wd_selected = wd_selected.reindex(columns=all_columns)

    merged = pd.concat([base, wd_selected], ignore_index=True)

    before_dedup = len(merged)

    if "source_id" in merged.columns:
        merged = merged.drop_duplicates(subset="source_id", keep="first")
    else:
        merged = merged.drop_duplicates(subset=["bp_rp", "abs_g"], keep="first")

    after_dedup = len(merged)

    print()
    print(f"Rows after concat:           {before_dedup:,}")
    print(f"Rows after dedup:            {after_dedup:,}")
    print(f"New unique WD rows added:    {after_dedup - len(base):,}")

    merged.to_parquet(out_path, index=False)
    print()
    print(f"Wrote merged file:           {out_path}")

    if args.write_csv:
        csv_path = out_path.with_suffix(".csv")
        merged.to_csv(csv_path, index=False)
        print(f"Wrote merged CSV:            {csv_path}")

    print()
    print("Merged CMD ranges:")
    print(f"  bp_rp: {merged['bp_rp'].min():.3f} to {merged['bp_rp'].max():.3f}")
    print(f"  abs_g: {merged['abs_g'].min():.3f} to {merged['abs_g'].max():.3f}")


if __name__ == "__main__":
    main()
