#!/usr/bin/env python3
"""
Download Gaia DR3 CMD data for WD comparison.

Output:
  gaia_cmd_wd_region.parquet
  gaia_cmd_background.parquet
"""

from pathlib import Path

import pandas as pd
from astroquery.gaia import Gaia


OUT = Path(".")
Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
Gaia.ROW_LIMIT = -1


def run_query(adql: str, outfile: Path):
    if outfile.exists():
        print(f"exists: {outfile}")
        return

    print(f"running query -> {outfile}")
    job = Gaia.launch_job_async(adql, dump_to_file=False)
    tab = job.get_results()
    df = tab.to_pandas()
    df.to_parquet(outfile, index=False)
    print(f"wrote {outfile} with {len(df):,} rows")


# Background CMD: high-quality, downsampled, all-sky.
# random_index is useful because it gives a stable random-ish subset.
background_query = """
SELECT
    source_id,
    phot_g_mean_mag,
    bp_rp,
    parallax,
    parallax_over_error,
    ruwe,
    phot_bp_rp_excess_factor,
    phot_g_mean_mag + 5 * LOG10(parallax) - 10 AS abs_g
FROM gaiadr3.gaia_source
WHERE
    parallax > 0
    AND parallax_over_error > 10
    AND phot_g_mean_mag IS NOT NULL
    AND bp_rp IS NOT NULL
    AND ruwe < 1.4
    AND random_index < 2000000
"""

# WD-region query: not downsampled as aggressively.
# This is the useful science background for the WD cooling track.
wd_region_query = """
SELECT
    source_id,
    phot_g_mean_mag,
    bp_rp,
    parallax,
    parallax_over_error,
    ruwe,
    phot_bp_rp_excess_factor,
    phot_g_mean_mag + 5 * LOG10(parallax) - 10 AS abs_g
FROM gaiadr3.gaia_source
WHERE
    parallax > 0
    AND parallax_over_error > 10
    AND phot_g_mean_mag IS NOT NULL
    AND bp_rp IS NOT NULL
    AND ruwe < 1.4
    AND bp_rp BETWEEN -0.8 AND 2.0
    AND phot_g_mean_mag + 5 * LOG10(parallax) - 10 BETWEEN 8.0 AND 17.5
"""

if __name__ == "__main__":
    run_query(background_query, OUT / "gaia_cmd_background.parquet")
    run_query(wd_region_query, OUT / "gaia_cmd_wd_region.parquet")
