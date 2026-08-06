#!/usr/bin/env python3
"""
diagnose_koester_12k.py

Targeted diagnostics for the discontinuity seen in the MESA Colors WD run
near Teff ~ 12,000 K.

What this checks
----------------
1. Whether the stellar structure is smooth while Colors magnitudes jump.
2. Whether Colors Mag_bol agrees with the MESA luminosity-derived Mbol trend.
3. Whether the Koester raw SED files or the precomputed flux_cube.bin contain
   a normalization discontinuity near the relevant Teff/logg grid cell.
4. Whether the cube values agree with the raw lookup-table SED files.

Run from your WD/python_helpers directory, e.g.

python diagnose_koester_12k.py \
  --history ../LOGS/history.data \
  --grid-dir /home/njm/SED_Tools/data/stellar_models/koester2 \
  --out-dir figures/koester_12k_diagnostics
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

SIGMA_SB_CGS = 5.670374419e-5  # erg s^-1 cm^-2 K^-4
MBOL_SUN_DEFAULT = 4.74

BANDS = ["u", "g", "r", "i", "z", "y"]
COLOR_PAIRS = [("u", "g"), ("g", "r"), ("r", "i"), ("i", "z"), ("z", "y")]

mpl.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "legend.fontsize": 8,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "figure.dpi": 150,
})


def read_history(path: Path) -> pd.DataFrame:
    lines = path.read_text(errors="ignore").splitlines()
    header_idx = None
    for i, line in enumerate(lines):
        cols = line.strip().split()
        if "model_number" in cols and "star_age" in cols:
            header_idx = i
            break
    if header_idx is None:
        raise RuntimeError(f"Could not find history header in {path}")

    df = pd.read_csv(path, sep=r"\s+", skiprows=header_idx, comment="#")
    if "model_number" in df.columns:
        n0 = len(df)
        df = df.drop_duplicates("model_number", keep="last").reset_index(drop=True)
        if len(df) != n0:
            print(f"Removed {n0-len(df)} duplicate model_number rows; kept last occurrence.")
    return df


def normalize_lookup_columns(df: pd.DataFrame) -> Tuple[pd.DataFrame, str, str, str, str]:
    original = list(df.columns)
    mapping = {c: c.strip().lower().replace(" ", "_") for c in df.columns}
    df = df.rename(columns=mapping)

    # filename column: prefer explicit names; otherwise first column.
    filename_candidates = ["filename", "file", "spectrum", "sed", "model", "path"]
    filename_col = next((c for c in filename_candidates if c in df.columns), df.columns[0])

    teff_col = next((c for c in ["teff", "t_eff", "temperature", "temp"] if c in df.columns), None)
    logg_col = next((c for c in ["logg", "log_g", "gravity", "loggrav"] if c in df.columns), None)
    meta_col = next((c for c in ["meta", "feh", "metallicity", "m_h", "mh"] if c in df.columns), None)

    missing = []
    if teff_col is None:
        missing.append("teff")
    if logg_col is None:
        missing.append("logg/log_g")
    if meta_col is None:
        missing.append("meta/feh/metallicity")
    if missing:
        raise RuntimeError(
            "Could not identify required lookup columns: "
            + ", ".join(missing)
            + f"\nOriginal columns: {original}\nNormalized columns: {list(df.columns)}"
        )

    for c in [teff_col, logg_col, meta_col]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=[teff_col, logg_col, meta_col]).copy()
    return df, filename_col, teff_col, logg_col, meta_col


def read_lookup(grid_dir: Path) -> Tuple[pd.DataFrame, str, str, str, str]:
    path = grid_dir / "lookup_table.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path)
    return normalize_lookup_columns(df)


def parse_numeric_two_col(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    wave: List[float] = []
    flux: List[float] = []
    with path.open("r", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            # Accept whitespace or comma separated files.
            s = s.replace(",", " ")
            parts = s.split()
            if len(parts) < 2:
                continue
            try:
                x = float(parts[0])
                y = float(parts[1])
            except ValueError:
                # Header line.
                continue
            if np.isfinite(x) and np.isfinite(y):
                wave.append(x)
                flux.append(y)
    if len(wave) < 4:
        raise RuntimeError(f"Could not read numeric SED columns from {path}")
    w = np.asarray(wave, dtype=float)
    f = np.asarray(flux, dtype=float)
    order = np.argsort(w)
    return w[order], f[order]


def integrate_flux(wave: np.ndarray, flux: np.ndarray) -> float:
    good = np.isfinite(wave) & np.isfinite(flux) & (wave > 0)
    if good.sum() < 4:
        return np.nan
    w = wave[good]
    f = flux[good]
    order = np.argsort(w)
    return float(np.trapezoid(f[order], w[order]))


def read_cube(grid_dir: Path) -> Dict[str, object]:
    path = grid_dir / "flux_cube.bin"
    if not path.exists():
        raise FileNotFoundError(path)

    with path.open("rb") as f:
        dims = np.fromfile(f, dtype=np.int32, count=4)
        if len(dims) != 4:
            raise RuntimeError("Could not read 4 int32 cube dimensions")
        nteff, nlogg, nmeta, nlam = [int(x) for x in dims]
        if min(nteff, nlogg, nmeta, nlam) <= 0 or max(nteff, nlogg, nmeta, nlam) > 2_000_000:
            raise RuntimeError(f"Implausible cube dimensions: {tuple(dims)}")
        teff_grid = np.fromfile(f, dtype=np.float64, count=nteff)
        logg_grid = np.fromfile(f, dtype=np.float64, count=nlogg)
        meta_grid = np.fromfile(f, dtype=np.float64, count=nmeta)
        wave_grid = np.fromfile(f, dtype=np.float64, count=nlam)
        flat = np.fromfile(f, dtype=np.float64)

    expected = nteff * nlogg * nmeta * nlam
    if flat.size != expected:
        raise RuntimeError(f"Cube payload size mismatch: got {flat.size}, expected {expected}")

    return {
        "path": path,
        "dims": (nteff, nlogg, nmeta, nlam),
        "teff": teff_grid,
        "logg": logg_grid,
        "meta": meta_grid,
        "wave": wave_grid,
        "flat": flat,
    }


def nearest_index(arr: np.ndarray, value: float) -> int:
    return int(np.argmin(np.abs(arr - value)))


def bracket(arr: np.ndarray, value: float) -> Tuple[float, float]:
    arr = np.asarray(arr, dtype=float)
    if value <= arr.min():
        return float(arr.min()), float(arr.min())
    if value >= arr.max():
        return float(arr.max()), float(arr.max())
    hi = int(np.searchsorted(arr, value, side="left"))
    lo = hi - 1
    return float(arr[lo]), float(arr[hi])


def lookup_sed_file(
    lookup: pd.DataFrame,
    filename_col: str,
    teff_col: str,
    logg_col: str,
    meta_col: str,
    teff: float,
    logg: float,
    meta: float,
    grid_dir: Path,
) -> Optional[Path]:
    x = lookup.copy()
    sel = (
        np.isclose(x[teff_col].to_numpy(float), teff, rtol=0, atol=1e-8)
        & np.isclose(x[logg_col].to_numpy(float), logg, rtol=0, atol=1e-8)
        & np.isclose(x[meta_col].to_numpy(float), meta, rtol=0, atol=1e-8)
    )
    rows = x.loc[sel]
    if rows.empty:
        return None
    fn = str(rows.iloc[0][filename_col])
    p = Path(fn)
    if not p.is_absolute():
        p = grid_dir / p
    return p


def choose_cube_order(
    cube_info: Dict[str, object],
    lookup: pd.DataFrame,
    filename_col: str,
    teff_col: str,
    logg_col: str,
    meta_col: str,
    grid_dir: Path,
) -> Tuple[str, np.ndarray, List[str]]:
    nteff, nlogg, nmeta, nlam = cube_info["dims"]
    flat = cube_info["flat"]
    wave = cube_info["wave"]
    teff_grid = cube_info["teff"]
    logg_grid = cube_info["logg"]
    meta_grid = cube_info["meta"]

    candidates = {
        "Fortran": flat.reshape((nteff, nlogg, nmeta, nlam), order="F"),
        "C": flat.reshape((nteff, nlogg, nmeta, nlam), order="C"),
    }

    sample_points = []
    for T in [10000, 12000, 12250, 15000, 20000, 40000]:
        if teff_grid.min() <= T <= teff_grid.max():
            sample_points.append((float(teff_grid[nearest_index(teff_grid, T)]), 8.0, float(meta_grid[0])))
    sample_points = list(dict.fromkeys(sample_points))

    notes: List[str] = []
    scores: Dict[str, float] = {}
    for name, cube in candidates.items():
        rels = []
        for T, g, m in sample_points:
            if not (logg_grid.min() <= g <= logg_grid.max()):
                continue
            p = lookup_sed_file(lookup, filename_col, teff_col, logg_col, meta_col, T, g, m, grid_dir)
            if p is None or not p.exists():
                continue
            try:
                rw, rf = parse_numeric_two_col(p)
            except Exception:
                continue
            iT = nearest_index(teff_grid, T)
            ig = nearest_index(logg_grid, g)
            im = nearest_index(meta_grid, m)
            raw_int = integrate_flux(rw, rf)
            cube_int = integrate_flux(wave, cube[iT, ig, im, :])
            if np.isfinite(raw_int) and np.isfinite(cube_int) and raw_int != 0:
                rels.append(abs(cube_int - raw_int) / abs(raw_int))
        score = float(np.nanmedian(rels)) if rels else np.inf
        scores[name] = score
        notes.append(f"cube order candidate {name}: median integrated-flux rel. diff vs raw = {score:.3e}")

    order = min(scores, key=scores.get)
    notes.append(f"selected cube order: {order}")
    return order, candidates[order], notes


def find_largest_jump(df: pd.DataFrame, cols: Iterable[str]) -> Tuple[str, int, float]:
    best_col = ""
    best_idx = -1
    best_abs = -np.inf
    for c in cols:
        if c not in df.columns:
            continue
        d = df[c].diff().to_numpy(float)
        if len(d) <= 1:
            continue
        idx = int(np.nanargmax(np.abs(d[1:])) + 1)
        val = abs(d[idx])
        if val > best_abs:
            best_col = c
            best_idx = idx
            best_abs = val
    return best_col, best_idx, best_abs


def write_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    pd.DataFrame(rows).to_csv(path, index=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--history", type=Path, default=Path("../LOGS/history.data"))
    ap.add_argument("--grid-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, default=Path("figures/koester_12k_diagnostics"))
    ap.add_argument("--mbol-sun", type=float, default=MBOL_SUN_DEFAULT)
    ap.add_argument("--teff-min", type=float, default=9000.0)
    ap.add_argument("--teff-max", type=float, default=15000.0)
    ap.add_argument("--logg-values", type=float, nargs="*", default=[8.0, 8.25])
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading history: {args.history}")
    h = read_history(args.history)
    print(f"Rows after cleaning: {len(h)}")

    required = ["model_number", "star_age", "log_Teff", "log_L"]
    missing = [c for c in required if c not in h.columns]
    if missing:
        raise RuntimeError(f"Missing required history columns: {missing}")

    h = h.copy()
    h["cooling_age"] = h["star_age"] - h["star_age"].iloc[0]
    h["Teff"] = 10.0 ** h["log_Teff"]
    h["Mbol_from_logL"] = args.mbol_sun - 2.5 * h["log_L"]
    if "Mag_bol" in h.columns:
        h["Mag_bol_minus_logL_Mbol"] = h["Mag_bol"] - h["Mbol_from_logL"]

    jump_cols = ["Mag_bol"] + BANDS
    best_col, best_idx, best_abs = find_largest_jump(h, jump_cols)
    if best_idx < 1:
        raise RuntimeError("Could not identify a jump")

    print(f"Largest jump among {jump_cols}: {best_col} at row {best_idx-1}->{best_idx}, |delta|={best_abs:.6g}")
    row_before = h.iloc[best_idx - 1]
    row_after = h.iloc[best_idx]
    suspect_teff = float(row_after["Teff"])
    suspect_logg = float(row_after["log_g"] if "log_g" in h.columns else np.nan)
    print(f"Suspect after-row model={int(row_after['model_number'])}, age={row_after['cooling_age']/1e9:.6f} Gyr, Teff={suspect_teff:.1f} K, logg={suspect_logg:.5f}")

    print(f"Reading lookup/cube from: {args.grid_dir}")
    lookup, filename_col, teff_col, logg_col, meta_col = read_lookup(args.grid_dir)
    cube_info = read_cube(args.grid_dir)
    print(f"cube dims: {cube_info['dims']} = (n_teff, n_logg, n_meta, n_lambda)")

    order_name, cube, order_notes = choose_cube_order(
        cube_info, lookup, filename_col, teff_col, logg_col, meta_col, args.grid_dir
    )
    for note in order_notes:
        print(note)

    teff_grid = np.asarray(cube_info["teff"], dtype=float)
    logg_grid = np.asarray(cube_info["logg"], dtype=float)
    meta_grid = np.asarray(cube_info["meta"], dtype=float)
    wave_grid = np.asarray(cube_info["wave"], dtype=float)
    meta0 = float(meta_grid[0])

    print(f"Teff bracket for suspect: {bracket(teff_grid, suspect_teff)}")
    if np.isfinite(suspect_logg):
        print(f"logg bracket for suspect: {bracket(logg_grid, suspect_logg)}")

    # History local table.
    local = h.iloc[max(0, best_idx - 6): min(len(h), best_idx + 7)].copy()
    local_cols = [
        "model_number", "cooling_age", "Teff", "log_Teff", "log_L", "Mbol_from_logL",
        "Mag_bol", "Mag_bol_minus_logL_Mbol", "log_R", "log_g", "Interp_rad",
    ] + BANDS
    local_cols = [c for c in local_cols if c in local.columns]
    local.to_csv(args.out_dir / "history_local_around_jump.csv", index=False, columns=local_cols)

    # Grid flux integrals around 12k.
    grid_rows: List[Dict[str, object]] = []
    for gval in args.logg_values:
        g = float(logg_grid[nearest_index(logg_grid, gval)])
        ig = nearest_index(logg_grid, g)
        im = nearest_index(meta_grid, meta0)
        for iT, T in enumerate(teff_grid):
            if not (args.teff_min <= T <= args.teff_max):
                continue
            cube_flux = cube[iT, ig, im, :]
            cube_int = integrate_flux(wave_grid, cube_flux)
            cube_ratio = cube_int / (SIGMA_SB_CGS * T**4) if np.isfinite(cube_int) else np.nan
            raw_int = np.nan
            raw_ratio = np.nan
            raw_path = lookup_sed_file(lookup, filename_col, teff_col, logg_col, meta_col, float(T), g, meta0, args.grid_dir)
            raw_vs_cube_rel_int = np.nan
            raw_vs_cube_max_rel = np.nan
            if raw_path is not None and raw_path.exists():
                rw, rf = parse_numeric_two_col(raw_path)
                raw_int = integrate_flux(rw, rf)
                raw_ratio = raw_int / (SIGMA_SB_CGS * T**4) if np.isfinite(raw_int) else np.nan
                if np.isfinite(raw_int) and raw_int != 0:
                    raw_vs_cube_rel_int = (cube_int - raw_int) / raw_int
                # Compare shapes on cube wavelength grid, avoiding deep zeros.
                rf_i = np.interp(wave_grid, rw, rf, left=np.nan, right=np.nan)
                good = np.isfinite(rf_i) & np.isfinite(cube_flux) & (np.abs(rf_i) > 0)
                if good.sum() > 10:
                    rel = np.abs((cube_flux[good] - rf_i[good]) / rf_i[good])
                    raw_vs_cube_max_rel = float(np.nanpercentile(rel, 99.0))
            grid_rows.append({
                "Teff": float(T),
                "logg": g,
                "meta": meta0,
                "cube_fint": cube_int,
                "raw_fint": raw_int,
                "cube_fint_over_sigmaT4": cube_ratio,
                "raw_fint_over_sigmaT4": raw_ratio,
                "cube_minus_raw_rel_integrated": raw_vs_cube_rel_int,
                "cube_vs_raw_p99_abs_rel_flux": raw_vs_cube_max_rel,
                "raw_file": str(raw_path) if raw_path is not None else "",
            })

    grid_df = pd.DataFrame(grid_rows)
    # Adjacent jumps in integrated flux ratios.
    jump_rows = []
    for g, sub in grid_df.groupby("logg"):
        sub = sub.sort_values("Teff").reset_index(drop=True)
        for ratio_col in ["raw_fint_over_sigmaT4", "cube_fint_over_sigmaT4"]:
            vals = sub[ratio_col].to_numpy(float)
            Ts = sub["Teff"].to_numpy(float)
            # Use log-ratio jumps; robust to arbitrary constant normalization.
            logvals = np.log10(np.abs(vals))
            d = np.diff(logvals)
            for idx in np.argsort(np.abs(d))[::-1][:8]:
                jump_rows.append({
                    "logg": g,
                    "quantity": ratio_col,
                    "Teff_lo": Ts[idx],
                    "Teff_hi": Ts[idx + 1],
                    "delta_log10_ratio": d[idx],
                    "factor_hi_over_lo": 10.0 ** d[idx],
                    "ratio_lo": vals[idx],
                    "ratio_hi": vals[idx + 1],
                })
    jump_df = pd.DataFrame(jump_rows).sort_values("delta_log10_ratio", key=lambda s: np.abs(s), ascending=False)

    grid_df.to_csv(args.out_dir / "koester_grid_integrals_around_12k.csv", index=False)
    jump_df.to_csv(args.out_dir / "koester_grid_integral_jumps.csv", index=False)

    # Print the rows closest to the culprit cell.
    print("\nGrid integral ratios near suspect temperature:")
    near = grid_df[(grid_df["Teff"] >= 11500) & (grid_df["Teff"] <= 12750)].copy()
    print(near[[
        "Teff", "logg", "raw_fint_over_sigmaT4", "cube_fint_over_sigmaT4",
        "cube_minus_raw_rel_integrated", "raw_file",
    ]].to_string(index=False, max_colwidth=70))

    print("\nLargest grid-integral log jumps:")
    print(jump_df.head(12).to_string(index=False))

    # Figure 1: history structure and Mbol residual.
    fig, axes = plt.subplots(5, 1, figsize=(5.2, 7.0), sharex=True)
    x = h["cooling_age"] / 1e9
    axes[0].plot(x, h["Teff"], lw=1.1)
    axes[0].set_ylabel(r"$T_{\rm eff}$")
    if "log_g" in h.columns:
        axes[1].plot(x, h["log_g"], lw=1.1)
        axes[1].set_ylabel(r"$\log g$")
    axes[2].plot(x, h["log_L"], lw=1.1, label=r"$\log L$")
    axes[2].set_ylabel(r"$\log L$")
    axes[3].plot(x, h["Mbol_from_logL"], lw=1.1, label="from log_L")
    if "Mag_bol" in h.columns:
        axes[3].plot(x, h["Mag_bol"], lw=1.1, label="Colors Mag_bol")
    axes[3].invert_yaxis()
    axes[3].set_ylabel(r"$M_{\rm bol}$")
    axes[3].legend(framealpha=0.8)
    if "Mag_bol_minus_logL_Mbol" in h.columns:
        axes[4].plot(x, h["Mag_bol_minus_logL_Mbol"], lw=1.1)
    axes[4].set_ylabel("Colors - logL")
    axes[4].set_xlabel("Cooling age (Gyr)")
    for ax in axes:
        ax.axvline(row_after["cooling_age"] / 1e9, ls="--", lw=0.8, color="0.6")
    fig.tight_layout()
    fig.savefig(args.out_dir / "history_structure_and_mbol_check.pdf", bbox_inches="tight")
    plt.close(fig)

    # Figure 2: local zoom around jump.
    fig, axes = plt.subplots(3, 1, figsize=(5.4, 5.6), sharex=True)
    lx = local["cooling_age"] / 1e9
    axes[0].plot(lx, local["Mbol_from_logL"], "o-", ms=3, label="Mbol from log_L")
    if "Mag_bol" in local.columns:
        axes[0].plot(lx, local["Mag_bol"], "o-", ms=3, label="Colors Mag_bol")
    axes[0].invert_yaxis()
    axes[0].legend(framealpha=0.8)
    axes[0].set_ylabel(r"$M_{\rm bol}$")
    for b in BANDS:
        if b in local.columns:
            axes[1].plot(lx, local[b], "o-", ms=3, label=b)
    axes[1].invert_yaxis()
    axes[1].set_ylabel("AB mag")
    axes[1].legend(ncol=3, framealpha=0.8)
    for b1, b2 in COLOR_PAIRS:
        if b1 in local.columns and b2 in local.columns:
            axes[2].plot(lx, local[b1] - local[b2], "o-", ms=3, label=f"{b1}-{b2}")
    axes[2].set_ylabel("color")
    axes[2].set_xlabel("Cooling age (Gyr)")
    axes[2].legend(ncol=3, framealpha=0.8)
    for ax in axes:
        ax.axvline(row_after["cooling_age"] / 1e9, ls="--", lw=0.8, color="0.6")
    fig.tight_layout()
    fig.savefig(args.out_dir / "history_local_zoom_around_jump.pdf", bbox_inches="tight")
    plt.close(fig)

    # Figure 3: grid integral ratios vs Teff.
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    for g, sub in grid_df.groupby("logg"):
        sub = sub.sort_values("Teff")
        ax.plot(sub["Teff"], sub["raw_fint_over_sigmaT4"], "o-", ms=3, lw=1.0, label=f"raw logg={g:g}")
        ax.plot(sub["Teff"], sub["cube_fint_over_sigmaT4"], "x--", ms=3, lw=1.0, label=f"cube logg={g:g}")
    ax.axvline(suspect_teff, ls="--", lw=0.8, color="0.5")
    ax.set_xlabel(r"$T_{\rm eff}$ (K)")
    ax.set_ylabel(r"$\int F_\lambda d\lambda / \sigma T_{\rm eff}^4$")
    ax.set_yscale("log")
    ax.legend(ncol=2, framealpha=0.8)
    fig.tight_layout()
    fig.savefig(args.out_dir / "koester_integrated_flux_ratio_vs_teff.pdf", bbox_inches="tight")
    plt.close(fig)

    # Figure 4: raw SEDs bracketing the suspect cell.
    Tlo, Thi = bracket(teff_grid, suspect_teff)
    if np.isfinite(suspect_logg):
        glo, ghi = bracket(logg_grid, suspect_logg)
    else:
        glo, ghi = args.logg_values[0], args.logg_values[-1]
    bracket_points = [(Tlo, glo), (Thi, glo), (Tlo, ghi), (Thi, ghi)]
    fig, ax = plt.subplots(figsize=(5.4, 3.8))
    for T, g in bracket_points:
        p = lookup_sed_file(lookup, filename_col, teff_col, logg_col, meta_col, T, g, meta0, args.grid_dir)
        if p is None or not p.exists():
            continue
        rw, rf = parse_numeric_two_col(p)
        ax.plot(rw, rf, lw=1.0, label=f"raw T={T:g}, logg={g:g}")
    ax.set_xlim(2500, 11000)
    ax.set_yscale("log")
    ax.set_xlabel(r"Wavelength ($\AA$)")
    ax.set_ylabel(r"Raw $F_\lambda$")
    ax.legend(framealpha=0.8)
    fig.tight_layout()
    fig.savefig(args.out_dir / "koester_raw_seds_bracketing_jump.pdf", bbox_inches="tight")
    plt.close(fig)

    # Figure 5: cube vs raw relative integrated difference.
    fig, ax = plt.subplots(figsize=(5.2, 3.5))
    for g, sub in grid_df.groupby("logg"):
        sub = sub.sort_values("Teff")
        ax.plot(sub["Teff"], sub["cube_minus_raw_rel_integrated"], "o-", ms=3, lw=1.0, label=f"logg={g:g}")
    ax.axhline(0, lw=0.8, color="0.5")
    ax.axvline(suspect_teff, ls="--", lw=0.8, color="0.5")
    ax.set_xlabel(r"$T_{\rm eff}$ (K)")
    ax.set_ylabel(r"$(F_{\rm cube}-F_{\rm raw})/F_{\rm raw}$")
    ax.legend(framealpha=0.8)
    fig.tight_layout()
    fig.savefig(args.out_dir / "koester_cube_minus_raw_integrated_flux.pdf", bbox_inches="tight")
    plt.close(fig)

    # Text report.
    report_lines: List[str] = []
    report_lines.append("Koester 12k discontinuity diagnostic")
    report_lines.append("====================================")
    report_lines.append(f"history: {args.history}")
    report_lines.append(f"grid_dir: {args.grid_dir}")
    report_lines.append(f"rows used: {len(h)}")
    report_lines.append("")
    report_lines.append("Largest Colors jump")
    report_lines.append("-------------------")
    report_lines.append(f"column: {best_col}")
    report_lines.append(f"rows/models: {best_idx-1}->{best_idx} / {int(row_before['model_number'])}->{int(row_after['model_number'])}")
    report_lines.append(f"cooling age after row: {row_after['cooling_age']/1e9:.9f} Gyr")
    report_lines.append(f"Teff after row: {suspect_teff:.3f} K")
    if np.isfinite(suspect_logg):
        report_lines.append(f"logg after row: {suspect_logg:.6f}")
    report_lines.append(f"delta {best_col}: {row_after[best_col] - row_before[best_col]:+.9f}")
    report_lines.append(f"delta log_L: {row_after['log_L'] - row_before['log_L']:+.9f}")
    report_lines.append(f"expected delta Mbol from log_L: {-2.5*(row_after['log_L'] - row_before['log_L']):+.9f}")
    if "Mag_bol" in h.columns:
        report_lines.append(f"actual delta Colors Mag_bol: {row_after['Mag_bol'] - row_before['Mag_bol']:+.9f}")
    report_lines.append("")
    report_lines.append("Cube")
    report_lines.append("----")
    report_lines.append(f"dims: {cube_info['dims']} = (n_teff, n_logg, n_meta, n_lambda)")
    report_lines.append(f"selected memory order: {order_name}")
    report_lines.extend(order_notes)
    report_lines.append(f"Teff grid range: {teff_grid.min():.3f}->{teff_grid.max():.3f}")
    report_lines.append(f"logg grid range: {logg_grid.min():.3f}->{logg_grid.max():.3f}")
    report_lines.append(f"meta grid range: {meta_grid.min():.3f}->{meta_grid.max():.3f}")
    report_lines.append(f"suspect Teff bracket: {bracket(teff_grid, suspect_teff)}")
    if np.isfinite(suspect_logg):
        report_lines.append(f"suspect logg bracket: {bracket(logg_grid, suspect_logg)}")
    report_lines.append("")
    report_lines.append("Largest grid-integral jumps")
    report_lines.append("---------------------------")
    for _, r in jump_df.head(12).iterrows():
        report_lines.append(
            f"logg={r['logg']:.3f} {r['quantity']} "
            f"Teff {r['Teff_lo']:.0f}->{r['Teff_hi']:.0f}: "
            f"delta_log10={r['delta_log10_ratio']:+.6f}, factor={r['factor_hi_over_lo']:.6g}"
        )
    report_lines.append("")
    report_lines.append("Files written")
    report_lines.append("-------------")
    for name in [
        "history_local_around_jump.csv",
        "koester_grid_integrals_around_12k.csv",
        "koester_grid_integral_jumps.csv",
        "history_structure_and_mbol_check.pdf",
        "history_local_zoom_around_jump.pdf",
        "koester_integrated_flux_ratio_vs_teff.pdf",
        "koester_raw_seds_bracketing_jump.pdf",
        "koester_cube_minus_raw_integrated_flux.pdf",
    ]:
        report_lines.append(str(args.out_dir / name))

    report_path = args.out_dir / "koester_12k_diagnostic_report.txt"
    report_path.write_text("\n".join(report_lines) + "\n")
    print(f"\nSaved report: {report_path}")
    print(f"Saved local history table: {args.out_dir / 'history_local_around_jump.csv'}")
    print(f"Saved grid-integral table: {args.out_dir / 'koester_grid_integrals_around_12k.csv'}")
    print("Done.")


if __name__ == "__main__":
    main()
