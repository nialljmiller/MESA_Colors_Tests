"""
yrec_sed_interface.py
=====================
Reads YREC SPOTS evolutionary tracks and runs every timestep through
SED_Model, producing per-track output directories structured as:

    <output_root>/
        m070_f034/
            m070_f034_colors.track      # input columns + SED_Model colors appended
            SEDs/
                step_00001_hot.csv      # surface (T_hot) SED
                step_00001_cool.csv     # spot (T_cool) SED
                step_00002_hot.csv
                step_00002_cool.csv
                ...

The hot SED is the unspotted photosphere (T_hot).
The cool SED is the spot component (T_cool).
Both are full ForwardResults — wavelengths, surface flux, observed flux —
so downstream code can combine them however it likes (e.g. rotational
modulation as a function of visible spot fraction).

Usage
-----
    python yrec_sed_interface.py \\
        --tracks  /path/to/YREC/ \\
        --grid    /path/to/atmosphere/grid/ \\
        --filters /path/to/filters/ \\
        --output  /path/to/output/ \\
        --mag-system AB \\
        --distance-pc 10.0
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

from sed_model import (
    load_grid,
    load_filters_from_instrument_dir,
    run_forward,
    RSUN_TO_CM,
    PC_TO_CM,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SOLAR_ZX = 0.0229   # Grevesse & Sauval 1998 — same reference as SPOTS

# ---------------------------------------------------------------------------
# Track column names (SPOTS format)
# ---------------------------------------------------------------------------

_COLNAMES = [
    "StepNum", "Mass", "Fspot", "Xspot", "Age_Gyr",
    "log_L", "log_R", "log_g", "log_Teff",
    "log_T_hot", "log_T_cool",
    "TauCZ", "total_I", "CZ_I",
    "X_cen", "X_surf", "ZX_surf",
    "Li_Li0", "H2_H2_0",
    "B_mag", "V_mag", "Rc_mag", "Ic_mag",
    "J_mag", "H_mag", "Ks_mag", "W1_mag",
    "G_mag", "BP_mag", "RP_mag",
]

# ---------------------------------------------------------------------------
# Track reader
# ---------------------------------------------------------------------------

def read_track(path: str | Path) -> pd.DataFrame:
    """Read one SPOTS .track file into a DataFrame with derived columns."""
    df = pd.read_csv(
        Path(path),
        sep=r"\s+",
        comment="#",
        header=None,
        names=_COLNAMES,
    )
    df["Teff"]   = 10.0 ** df["log_Teff"]
    df["T_hot"]  = 10.0 ** df["log_T_hot"]
    df["T_cool"] = 10.0 ** df["log_T_cool"]
    df["R_cm"]   = (10.0 ** df["log_R"]) * RSUN_TO_CM
    df["meta"]   = np.log10((df["ZX_surf"] / SOLAR_ZX).clip(lower=1e-10))
    return df


def read_all_tracks(directory: str | Path) -> dict[str, pd.DataFrame]:
    """Read every .track file in a directory, keyed by stem."""
    return {
        p.stem: read_track(p)
        for p in sorted(Path(directory).glob("*.track"))
    }


def parse_track_name(stem: str) -> tuple[float, float]:
    """Parse (mass_msun, fspot) from a stem like 'm070_f034'."""
    m = re.match(r"m(\d+)_f(\d+)", stem)
    if not m:
        raise ValueError(f"Cannot parse track name: {stem!r}")
    return int(m.group(1)) / 100.0, int(m.group(2)) / 1000.0


# ---------------------------------------------------------------------------
# SED writer
# ---------------------------------------------------------------------------

def _write_sed(path: Path, result, label: str) -> None:
    """Write a ForwardResult SED to CSV with a label in the header."""
    header = (
        f"# SED_Model SED output [{label}]\n"
        f"# Teff={result.teff:.1f} K  logg={result.logg:.3f}  "
        f"[M/H]={result.meta:.3f}\n"
        f"# R={result.R:.4e} cm  d={result.d:.4e} cm\n"
        f"# wavelength_AA,surface_flux_erg_s_cm2_AA,observed_flux_erg_s_cm2_AA"
    )
    data = np.column_stack([
        result.wavelengths,
        result.surface_flux,
        result.observed_flux,
    ])
    np.savetxt(path, data, delimiter=",", header=header, comments="")


# ---------------------------------------------------------------------------
# Per-track runner
# ---------------------------------------------------------------------------

def run_track(
    stem:          str,
    df:            pd.DataFrame,
    grid,
    filters:       list,
    output_dir:    Path,
    d:             float,
    mag_system:    str,
    interp_method: str,
    verbose:       bool,
) -> None:
    """
    Run every timestep of one track through SED_Model.

    Produces:
        <output_dir>/<stem>/SEDs/step_NNNNN_hot.csv
        <output_dir>/<stem>/SEDs/step_NNNNN_cool.csv
        <output_dir>/<stem>/<stem>_colors.track
    """
    track_dir = output_dir / stem
    sed_dir   = track_dir / "SEDs"
    sed_dir.mkdir(parents=True, exist_ok=True)

    filter_names = [f.name for f in filters]

    hot_rows  = []
    cool_rows = []

    for _, row in df.iterrows():
        step = int(row["StepNum"])
        logg = float(row["log_g"])
        meta = float(row["meta"])
        R    = float(row["R_cm"])

        # Hot component: unspotted photosphere
        res_hot = run_forward(
            teff=float(row["T_hot"]),
            logg=logg, meta=meta, R=R, d=d,
            grid=grid, filters=filters,
            mag_system=mag_system, interp_method=interp_method,
        )

        # Cool component: spot
        res_cool = run_forward(
            teff=float(row["T_cool"]),
            logg=logg, meta=meta, R=R, d=d,
            grid=grid, filters=filters,
            mag_system=mag_system, interp_method=interp_method,
        )

        # Save both SEDs
        step_tag = f"step_{step:05d}"
        _write_sed(sed_dir / f"{step_tag}_hot.csv",  res_hot,  label="hot/surface")
        _write_sed(sed_dir / f"{step_tag}_cool.csv", res_cool, label="cool/spot")

        # Collect magnitudes for the augmented track file
        hot_row  = {"bol_mag_hot":  res_hot.bol_mag}
        cool_row = {"bol_mag_cool": res_cool.bol_mag}
        for fname in filter_names:
            hot_row[f"{fname}_hot"]   = res_hot.magnitudes.get(fname,  float("nan"))
            cool_row[f"{fname}_cool"] = res_cool.magnitudes.get(fname, float("nan"))
        hot_rows.append(hot_row)
        cool_rows.append(cool_row)

        if verbose:
            print(f"    step {step:5d}  T_hot={res_hot.teff:.0f} K  "
                  f"T_cool={res_cool.teff:.0f} K")

    # Augmented track file: original columns + SED_Model columns
    out_df = pd.concat(
        [df.reset_index(drop=True),
         pd.DataFrame(hot_rows),
         pd.DataFrame(cool_rows)],
        axis=1,
    )

    out_path = track_dir / f"{stem}_colors.track"
    _write_colors_track(out_path, out_df, stem, mag_system)


def _write_colors_track(
    path: Path, df: pd.DataFrame, stem: str, mag_system: str
) -> None:
    """Write the augmented track in the original SPOTS column style."""
    lines = [
        "# SPOTS track augmented with SED_Model synthetic photometry",
        f"# Track: {stem}",
        f"# Magnitude system: {mag_system}",
        "# Original SPOTS columns preserved; SED_Model columns appended.",
        "# _hot = photosphere (T_hot)  |  _cool = spot (T_cool)",
        "# " + "  ".join(df.columns),
    ]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
        df.to_csv(fh, sep=" ", index=False, header=False,
                  float_format="%.8e", na_rep="-99.0")


# ---------------------------------------------------------------------------
# Top-level runner
# ---------------------------------------------------------------------------

def run_all_tracks(
    track_dir:     str | Path,
    grid_dir:      str | Path,
    filter_dir:    str | Path,
    output_dir:    str | Path,
    distance_pc:   float = 10.0,
    mag_system:    str   = "AB",
    interp_method: str   = "hermite",
    verbose:       bool  = True,
) -> None:
    """Run every .track file in track_dir through SED_Model."""

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if verbose:
        print("Loading atmosphere grid ...", grid_dir, Path(grid_dir), Path(grid_dir).resolve)
    grid = load_grid(grid_dir)

    if verbose:
        print("Loading filters ...")
    filters = load_filters_from_instrument_dir(filter_dir)
    if verbose:
        print(f"  {len(filters)} filters: {[f.name for f in filters]}")

    d = distance_pc * PC_TO_CM

    tracks = read_all_tracks(track_dir)
    if verbose:
        print(f"\nFound {len(tracks)} tracks in {track_dir}")

    for stem, df in tracks.items():
        mass, fspot = parse_track_name(stem)
        if verbose:
            print(f"\n{stem}  (M={mass:.2f} Msun  Fspot={fspot:.3f}  "
                  f"{len(df)} steps)")
        run_track(
            stem=stem, df=df,
            grid=grid, filters=filters,
            output_dir=output_dir,
            d=d, mag_system=mag_system,
            interp_method=interp_method,
            verbose=verbose,
        )

    if verbose:
        print(f"\nDone. Output in {output_dir}")




if __name__ == "__main__":

    tracks = '/home/njm/MESA/MESA_Colors_Tests/YREC/tracks/'
    grid   = '/home/njm/SED_Tools/data/stellar_models/Kurucz2003all__alpha_00/'
    filters = '/home/njm/SED_Tools/data/filters/Generic/Johnson/'
    output = '/home/njm/MESA/MESA_Colors_Tests/YREC/'
    distance_pc = 10.0
    mag_system = "AB"
    interp_method = 'hermite'

    run_all_tracks(
        track_dir=tracks,
        grid_dir=grid,
        filter_dir=filters,
        output_dir=output,
        distance_pc=distance_pc,
        mag_system=mag_system,
        interp_method=interp_method,
        verbose=True
    )
