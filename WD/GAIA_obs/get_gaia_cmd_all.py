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
    source_id,
    phot_g_mean_mag,
    bp_rp,
    parallax,
    parallax_over_error,
    ruwe
FROM gaiadr3.gaia_source
WHERE
    random_index >= {lo}
    AND random_index < {hi}
    AND parallax > 0
    AND parallax_over_error > {parallax_over_error}
    AND bp_rp IS NOT NULL
    AND phot_g_mean_mag IS NOT NULL
    AND ruwe < {ruwe_max}
    AND phot_g_mean_mag < {g_max}
"""


def fetch_chunk(lo, hi, parallax_over_error, ruwe_max, g_max, retries, sleep):
    query = QUERY_TEMPLATE.format(
        lo=lo,
        hi=hi,
        parallax_over_error=parallax_over_error,
        ruwe_max=ruwe_max,
        g_max=g_max,
    )

    for attempt in range(1, retries + 1):
        try:
            print(f"Querying random_index {lo:,} -> {hi:,}  attempt {attempt}/{retries}")

            # Use synchronous TAP. The Gaia async result store is what is failing.
            job = Gaia.launch_job(query)
            table = job.get_results()
            df = table.to_pandas()

            df["abs_g"] = df["phot_g_mean_mag"] + 5.0 * np.log10(df["parallax"]) - 10.0
            df = df.replace([np.inf, -np.inf], np.nan)
            df = df.dropna(subset=["bp_rp", "abs_g"])

            return df

        except Exception as e:
            print(f"  failed: {e}")
            if attempt < retries:
                print(f"  sleeping {sleep} s then retrying")
                time.sleep(sleep)

    print(f"  giving up on chunk {lo:,} -> {hi:,}")
    return None


def combine_chunks(chunk_dir, out):
    chunks = sorted(chunk_dir.glob("chunk_*.parquet"))

    if not chunks:
        raise RuntimeError(f"No successful chunks found in {chunk_dir}")

    print(f"Combining {len(chunks)} chunks")

    dfs = []
    for p in chunks:
        print(f"  reading {p}")
        dfs.append(pd.read_parquet(p))

    df = pd.concat(dfs, ignore_index=True)

    if "source_id" in df.columns:
        df = df.drop_duplicates(subset="source_id")
    else:
        df = df.drop_duplicates(subset=["bp_rp", "abs_g"])

    df.to_parquet(out, index=False)
    df.to_csv(out.with_suffix(".csv"), index=False)

    print(f"Total rows: {len(df):,}")
    print(f"wrote {out}")
    print(f"wrote {out.with_suffix('.csv')}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="gaia_cmd_all.parquet")
    parser.add_argument("--max-random-index", type=int, default=80_000_000)
    parser.add_argument("--chunk-size", type=int, default=50_000)
    parser.add_argument("--parallax-over-error", type=float, default=10.0)
    parser.add_argument("--ruwe-max", type=float, default=1.4)
    parser.add_argument("--g-max", type=float, default=21.0)
    parser.add_argument("--retries", type=int, default=5)
    parser.add_argument("--sleep", type=float, default=10.0)
    parser.add_argument("--chunk-dir", default="gaia_chunks")
    parser.add_argument("--combine-only", action="store_true")
    args = parser.parse_args()

    out = Path(args.out)
    chunk_dir = Path(args.chunk_dir)
    chunk_dir.mkdir(parents=True, exist_ok=True)

    if not args.combine_only:
        for lo in range(0, args.max_random_index, args.chunk_size):
            hi = min(lo + args.chunk_size, args.max_random_index)
            chunk_path = chunk_dir / f"chunk_{lo:012d}_{hi:012d}.parquet"

            if chunk_path.exists():
                print(f"Skipping existing {chunk_path}")
                continue

            df = fetch_chunk(
                lo=lo,
                hi=hi,
                parallax_over_error=args.parallax_over_error,
                ruwe_max=args.ruwe_max,
                g_max=args.g_max,
                retries=args.retries,
                sleep=args.sleep,
            )

            if df is None:
                continue

            print(f"  got {len(df):,} rows")
            df.to_parquet(chunk_path, index=False)
            print(f"  wrote {chunk_path}")

    combine_chunks(chunk_dir, out)


if __name__ == "__main__":
    main()
