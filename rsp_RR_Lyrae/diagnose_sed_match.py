#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd

from sed_model import load_grid, load_filters_from_instrument_dir, run_forward

HISTORY = "LOGS_overnight/history_last100.data"
DATA_ROOT = Path("/home/njm/SED_Tools/data")

FILTER_DIR = DATA_ROOT / "filters/Roman/WFI"
VEGA_SED = DATA_ROOT / "stellar_models/vega_flam.csv"

GRID_DIRS = [
    DATA_ROOT / "stellar_models/Kurucz2003all",
    DATA_ROOT / "stellar_models/Kurucz2003",
    DATA_ROOT / "stellar_models/Kurucz2003_dense",
    DATA_ROOT / "stellar_models/combined_models",
    DATA_ROOT / "stellar_models/my_combined_grid",
]

METAS = [-2.5, -2.0, -1.53, -1.5, -1.0, -0.5, 0.0, 0.0004]
INTERPS = ["hermite", "linear"]
MAG_SYSTEMS = ["AB", "Vega"]

DISTANCE_CM = 3.0857e19
RSUN_CM = 6.957e10
FCOLS = ["F062", "F087", "F106", "F129", "F146", "F158", "F184", "F213"]

def read_history(path):
    with open(path) as f:
        for i, l in enumerate(f):
            if "model_number" in l.split():
                return pd.read_csv(path, sep=r"\s+", skiprows=i)
    raise RuntimeError("Could not find MESA history header")

h = read_history(HISTORY)

nmax = int(h["rsp_num_periods"].max())
if (h["rsp_num_periods"].astype(int) == nmax).sum() < 200:
    nmax -= 1

h = h[h["rsp_num_periods"].astype(int).between(nmax - 2, nmax)].copy()
h = h.iloc[::max(1, len(h)//80)].copy()

fcols = [c for c in FCOLS if c in h.columns]

filters = load_filters_from_instrument_dir(
    str(FILTER_DIR),
    vega_sed_path=str(VEGA_SED) if VEGA_SED.exists() else None,
)

if "photosphere_r" in h.columns:
    h["R_cm_for_sed"] = h["photosphere_r"].astype(float) * RSUN_CM
elif "log_R" in h.columns:
    h["R_cm_for_sed"] = (10.0 ** h["log_R"].astype(float)) * RSUN_CM
else:
    raise RuntimeError("No radius column found")

if "effective_T" in h.columns:
    h["Teff_for_sed"] = h["effective_T"].astype(float)
else:
    h["Teff_for_sed"] = 10.0 ** h["log_Teff"].astype(float)

print(f"sample rows: {len(h)}")
print(f"filters compared: {fcols}")

scores = []

for grid_dir in GRID_DIRS:
    if not grid_dir.exists():
        continue

    print(f"\nloading grid: {grid_dir}")
    grid = load_grid(str(grid_dir))

    for meta in METAS:
        for interp in INTERPS:
            for mag_system in MAG_SYSTEMS:
                diffs = []

                try:
                    for row in h.itertuples(index=False):
                        res = run_forward(
                            teff=float(row.Teff_for_sed),
                            logg=float(row.log_g),
                            meta=float(meta),
                            R=float(row.R_cm_for_sed),
                            d=DISTANCE_CM,
                            grid=grid,
                            filters=filters,
                            mag_system=mag_system,
                            interp_method=interp,
                        )

                        row_dict = row._asdict()
                        for fn in fcols:
                            diffs.append(res.magnitudes[fn] - float(row_dict[fn]))

                    diffs = np.asarray(diffs, dtype=float)
                    score = np.nanmedian(np.abs(diffs))
                    signed = np.nanmedian(diffs)

                    scores.append((score, signed, grid_dir.name, meta, interp, mag_system))
                    print(
                        f"{score:9.5f}  signed={signed:+9.5f}  "
                        f"grid={grid_dir.name:22s} meta={meta:7g} "
                        f"interp={interp:7s} mag={mag_system}"
                    )

                except Exception as e:
                    print(
                        f"FAILED grid={grid_dir.name} meta={meta} "
                        f"interp={interp} mag={mag_system}: {e}"
                    )

print("\n=== BEST MATCHES ===")
for score, signed, grid, meta, interp, mag_system in sorted(scores)[:20]:
    print(
        f"{score:9.5f}  signed={signed:+9.5f}  "
        f"grid={grid:22s} meta={meta:7g} interp={interp:7s} mag={mag_system}"
    )
