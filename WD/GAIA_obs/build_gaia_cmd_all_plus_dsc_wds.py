#!/usr/bin/env python3

from pathlib import Path
import argparse
import time

import numpy as np
import pandas as pd
from astroquery.gaia import Gaia


Gaia.ROW_LIMIT = -1


QUERY_TEMPLATE = """
SELECT
    gs.source_id,
    gs.phot_g_mean_mag,
    gs.bp_rp,
    gs.parallax,
    gs.parallax_over_error,
    gs.ruwe,
    ap.classprob_dsc_combmod_whitedwarf AS wd_prob
FROM gaiadr3.gaia_source AS gs
JOIN gaiadr3.astrophysical_parameters AS ap
    ON gs.source_id = ap.source_id
WHERE
    gs.random_index >= {lo}
    AND gs.random_index < {hi}
    AND ap.classprob_dsc_combmod_whitedwarf >= {prob_min}
    AND gs.parallax > 0
    AND gs.parallax_over_error > {parallax_over_error}
    AND gs.bp_rp IS NOT NULL
    AND gs.phot_g_mean_mag IS NOT NULL
    AND gs.ruwe < {ruwe_max}
    AND gs.phot_g_mean_mag < {g_max}
"""


def read_cmd_table(path: Path) -> pd.DataFrame:
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


def fetch_wd_chunk(
    lo,
    hi,
    prob_min,
    parallax_over_error,
    ruwe_max,
    g_max,
    retries,
    sleep,
):
    query = QUERY_TEMPLATE.format(
        lo=lo,
        hi=hi,
        prob_min=prob_min,
        parallax_over_error=parallax_over_error,
        ruwe_max=ruwe_max,
        g_max=g_max,
    )

    for attempt in range(1, retries + 1):
        try:
            print(f"Querying random_index {lo:,} -> {hi:,}  attempt {attempt}/{retries}")

            # Sync TAP avoids the async missing-result problem you hit earlier.
            job = Gaia.launch_job(query)
            table = job.get_results()
            df = table.to_pandas()

            df = df.replace([np.inf, -np.inf], np.nan)

            if len(df) > 0:
                df["abs_g"] = df["phot_g_mean_mag"] + 5.0 * np.log10(df["parallax"]) - 10.0
                df = df.dropna(subset=["bp_rp", "abs_g"])

            return df

        except Exception as e:
            print(f"  failed: {e}")
            if attempt < retries:
                print(f"  sleeping {sleep} s then retrying")
                time.sleep(sleep)

    print(f"  giving up on chunk {lo:,} -> {hi:,}")
    return None


def combine_wd_chunks(chunk_dir: Path, wd_out: Path) -> pd.DataFrame:
    chunks = sorted(chunk_dir.glob("wd_chunk_*.parquet"))

    if not chunks:
        raise RuntimeError(f"No successful WD chunks found in {chunk_dir}")

    print()
    print(f"Combining {len(chunks)} WD chunks")

    dfs = []
    for p in chunks:
        df = pd.read_parquet(p)
        if len(df) > 0:
            dfs.append(df)

    if not dfs:
        raise RuntimeError("WD chunks exist, but all are empty.")

    wd = pd.concat(dfs, ignore_index=True)

    if "source_id" in wd.columns:
        wd = wd.drop_duplicates(subset="source_id", keep="first")
    else:
        wd = wd.drop_duplicates(subset=["bp_rp", "abs_g"], keep="first")

    wd = wd.replace([np.inf, -np.inf], np.nan)
    wd = wd.dropna(subset=["bp_rp", "abs_g"])

    wd.to_parquet(wd_out, index=False)

    print(f"WD rows after dedup: {len(wd):,}")
    print(f"Wrote WD catalogue:  {wd_out}")

    return wd


def merge_base_and_wds(base_path: Path, wd_path: Path, out_path: Path, write_csv: bool):
    print()
    print(f"Reading base CMD:    {base_path}")
    base = read_cmd_table(base_path)

    print(f"Reading WD catalog:  {wd_path}")
    wd = read_cmd_table(wd_path)

    print()
    print(f"Base rows:           {len(base):,}")
    print(f"WD rows:             {len(wd):,}")

    all_columns = sorted(set(base.columns) | set(wd.columns))
    base = base.reindex(columns=all_columns)
    wd = wd.reindex(columns=all_columns)

    merged = pd.concat([base, wd], ignore_index=True)

    before = len(merged)

    if "source_id" in merged.columns:
        merged = merged.drop_duplicates(subset="source_id", keep="first")
    else:
        merged = merged.drop_duplicates(subset=["bp_rp", "abs_g"], keep="first")

    after = len(merged)

    print(f"Rows after concat:   {before:,}")
    print(f"Rows after dedup:    {after:,}")
    print(f"New rows added:      {after - len(base):,}")

    merged.to_parquet(out_path, index=False)
    print()
    print(f"Wrote merged CMD:    {out_path}")

    if write_csv:
        csv_path = out_path.with_suffix(".csv")
        merged.to_csv(csv_path, index=False)
        print(f"Wrote merged CSV:    {csv_path}")

    print()
    print("Merged CMD ranges:")
    print(f"  bp_rp: {merged['bp_rp'].min():.3f} to {merged['bp_rp'].max():.3f}")
    print(f"  abs_g: {merged['abs_g'].min():.3f} to {merged['abs_g'].max():.3f}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build a broad Gaia CMD supplemented with Gaia DR3 DSC-classified "
            "white dwarfs, without CMD polygon cuts."
        )
    )

    parser.add_argument("--base", default="gaia_cmd_all.parquet")
    parser.add_argument("--out", default="gaia_cmd_all_plus_dsc_wds.parquet")
    parser.add_argument("--wd-out", default="gaia_cmd_dsc_wds.parquet")
    parser.add_argument("--chunk-dir", default="gaia_dsc_wd_chunks")

    parser.add_argument("--prob-min", type=float, default=0.5)
    parser.add_argument("--max-random-index", type=int, default=2_000_000_000)
    parser.add_argument("--chunk-size", type=int, default=10_000_000)

    parser.add_argument("--parallax-over-error", type=float, default=10.0)
    parser.add_argument("--ruwe-max", type=float, default=1.4)
    parser.add_argument("--g-max", type=float, default=21.0)

    parser.add_argument("--retries", type=int, default=5)
    parser.add_argument("--sleep", type=float, default=10.0)

    parser.add_argument("--combine-only", action="store_true")
    parser.add_argument("--merge-only", action="store_true")
    parser.add_argument("--write-csv", action="store_true")

    args = parser.parse_args()

    base_path = Path(args.base)
    out_path = Path(args.out)
    wd_out = Path(args.wd_out)
    chunk_dir = Path(args.chunk_dir)
    chunk_dir.mkdir(parents=True, exist_ok=True)

    if not args.merge_only:
        if not args.combine_only:
            for lo in range(0, args.max_random_index, args.chunk_size):
                hi = min(lo + args.chunk_size, args.max_random_index)
                chunk_path = chunk_dir / f"wd_chunk_{lo:012d}_{hi:012d}.parquet"

                if chunk_path.exists():
                    print(f"Skipping existing {chunk_path}")
                    continue

                df = fetch_wd_chunk(
                    lo=lo,
                    hi=hi,
                    prob_min=args.prob_min,
                    parallax_over_error=args.parallax_over_error,
                    ruwe_max=args.ruwe_max,
                    g_max=args.g_max,
                    retries=args.retries,
                    sleep=args.sleep,
                )

                if df is None:
                    continue

                print(f"  got {len(df):,} WD rows")
                df.to_parquet(chunk_path, index=False)
                print(f"  wrote {chunk_path}")

        combine_wd_chunks(chunk_dir, wd_out)

    merge_base_and_wds(
        base_path=base_path,
        wd_path=wd_out,
        out_path=out_path,
        write_csv=args.write_csv,
    )


if __name__ == "__main__":
    main()
