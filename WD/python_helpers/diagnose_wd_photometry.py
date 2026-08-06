#!/usr/bin/env python3
"""
diagnose_wd_photometry.py

Diagnostic plots/tables for MESA Colors white-dwarf cooling tracks.

This is deliberately more diagnostic than publication-quality. It checks whether
features in the WD photometry are coming from the MESA cooling model, the
atmosphere interpolation, the filter magnitudes, or history.data ordering.

Typical use from the same directory as plot_wd.py:

    python diagnose_wd_photometry.py \
        --history ../LOGS/history.data \
        --grid-dir /home/njm/SED_Tools/data/stellar_models/koester2 \
        --out-dir figures/wd_diagnostics

Outputs:
    wd_diagnostic_report.txt
    wd_diagnostic_top_jumps.csv
    wd_diagnostic_structure_vs_age.pdf
    wd_diagnostic_mags_colors_vs_age.pdf
    wd_diagnostic_color_teff_raw.pdf
    wd_diagnostic_cmd_colored_by_logg.pdf
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


DEFAULT_HISTORY = Path("../LOGS/history.data")
DEFAULT_GRID_DIR = Path("/home/njm/SED_Tools/data/stellar_models/koester2")
DEFAULT_OUT_DIR = Path("figures/wd_diagnostics")

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


@dataclass
class GridInfo:
    teff_grid: np.ndarray | None = None
    logg_grid: np.ndarray | None = None
    meta_grid: np.ndarray | None = None
    cube_dims: tuple[int, int, int, int] | None = None
    lookup_path: Path | None = None
    cube_path: Path | None = None


def normalize_name(name: str) -> str:
    return name.strip().lower().replace("_", "").replace("-", "")


def read_history(path: Path) -> dict[str, np.ndarray]:
    """Read a MESA history.data file into a dict of numpy arrays."""
    if not path.exists():
        raise FileNotFoundError(f"history.data not found: {path}")

    lines = path.read_text(errors="replace").splitlines()
    header_idx = None

    for i, line in enumerate(lines):
        fields = line.split()
        if "model_number" in fields and "star_age" in fields:
            header_idx = i
            break

    if header_idx is None:
        # fallback matching the original plotting script logic
        for i, line in enumerate(lines):
            stripped = line.strip()
            if stripped.startswith("model_number") or "star_age" in stripped.split():
                header_idx = i
                break

    if header_idx is None:
        raise ValueError(f"Could not find MESA history header in {path}")

    cols = lines[header_idx].split()
    rows: list[list[float]] = []

    for line in lines[header_idx + 1:]:
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        fields = s.split()
        if len(fields) != len(cols):
            # history.data may include repeated header blocks after restarts.
            if "model_number" in fields and "star_age" in fields:
                continue
            print(f"WARNING: skipping malformed line with {len(fields)} fields, expected {len(cols)}")
            continue
        try:
            rows.append([float(x) for x in fields])
        except ValueError:
            continue

    if not rows:
        raise ValueError(f"No numeric history rows found in {path}")

    data = np.asarray(rows, dtype=float)
    return {c: data[:, i] for i, c in enumerate(cols)}


def find_col(h: dict[str, np.ndarray], aliases: Iterable[str], required: bool = False) -> str | None:
    wanted = {normalize_name(a) for a in aliases}
    for c in h:
        if normalize_name(c) in wanted:
            return c
    if required:
        raise KeyError(f"Could not find required column. Tried aliases: {list(aliases)}")
    return None


def clean_history(h: dict[str, np.ndarray]) -> tuple[dict[str, np.ndarray], list[str]]:
    """Keep the last occurrence of duplicate model_number rows and sort by model_number.

    MESA history files can contain stale rows after restarts/retries. Those can create
    fake jumps in photometry. This cleaning is diagnostic-safe: the report states
    whether anything was removed.
    """
    notes: list[str] = []
    n0 = len(next(iter(h.values())))

    model_col = find_col(h, ["model_number"])
    age_col = find_col(h, ["star_age", "age"], required=True)

    keep_idx: np.ndarray
    if model_col is not None:
        model = h[model_col].astype(int)
        last_for_model: dict[int, int] = {}
        for i, m in enumerate(model):
            last_for_model[int(m)] = i
        keep_idx = np.array(sorted(last_for_model.values(), key=lambda i: model[i]), dtype=int)
        n_removed = n0 - len(keep_idx)
        if n_removed:
            notes.append(f"Removed {n_removed} duplicate model_number rows by keeping the last occurrence.")
        else:
            notes.append("No duplicate model_number rows found.")
    else:
        keep_idx = np.arange(n0, dtype=int)
        notes.append("No model_number column found; duplicate-model cleaning skipped.")

    # If ages are still not monotonic, sort by star_age for derivative diagnostics.
    age_after = h[age_col][keep_idx]
    if np.any(np.diff(age_after) < 0):
        notes.append("WARNING: star_age is non-monotonic after duplicate cleaning; sorting by star_age for diagnostics.")
        keep_idx = keep_idx[np.argsort(age_after)]
    else:
        notes.append("star_age is monotonic after cleaning.")

    return {c: v[keep_idx] for c, v in h.items()}, notes


def finite_common_mask(arrays: Iterable[np.ndarray]) -> np.ndarray:
    arrays = list(arrays)
    mask = np.ones(len(arrays[0]), dtype=bool)
    for a in arrays:
        mask &= np.isfinite(a)
    return mask


def read_lookup_grid(grid_dir: Path) -> GridInfo:
    info = GridInfo()
    lookup = grid_dir / "lookup_table.csv"
    info.lookup_path = lookup if lookup.exists() else None

    if not lookup.exists():
        return info

    with lookup.open(newline="", errors="replace") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            return info

        norm_to_name = {normalize_name(name): name for name in reader.fieldnames}
        teff_name = None
        logg_name = None
        meta_name = None

        for candidate in ["teff", "temperature", "t"]:
            if candidate in norm_to_name:
                teff_name = norm_to_name[candidate]
                break
        for candidate in ["logg", "gravity", "g"]:
            if candidate in norm_to_name:
                logg_name = norm_to_name[candidate]
                break
        for candidate in ["meta", "feh", "mh", "metallicity", "z"]:
            if candidate in norm_to_name:
                meta_name = norm_to_name[candidate]
                break

        teff_vals: list[float] = []
        logg_vals: list[float] = []
        meta_vals: list[float] = []

        for row in reader:
            if teff_name is not None and row.get(teff_name, "") != "":
                teff_vals.append(float(row[teff_name]))
            if logg_name is not None and row.get(logg_name, "") != "":
                logg_vals.append(float(row[logg_name]))
            if meta_name is not None and row.get(meta_name, "") != "":
                meta_vals.append(float(row[meta_name]))

    if teff_vals:
        info.teff_grid = np.unique(np.asarray(teff_vals, dtype=float))
    if logg_vals:
        info.logg_grid = np.unique(np.asarray(logg_vals, dtype=float))
    if meta_vals:
        info.meta_grid = np.unique(np.asarray(meta_vals, dtype=float))

    return info


def read_cube_header(grid_dir: Path, info: GridInfo) -> GridInfo:
    cube = grid_dir / "flux_cube.bin"
    info.cube_path = cube if cube.exists() else None
    if not cube.exists():
        return info

    try:
        dims = np.fromfile(cube, dtype=np.int32, count=4)
    except Exception:
        return info

    if len(dims) != 4:
        return info

    dims_tuple = tuple(int(x) for x in dims)
    plausible = (
        1 <= dims_tuple[0] <= 10000
        and 1 <= dims_tuple[1] <= 10000
        and 1 <= dims_tuple[2] <= 10000
        and 10 <= dims_tuple[3] <= 10_000_000
    )

    if not plausible:
        return info

    info.cube_dims = dims_tuple

    # Try reading the axes assuming float64 immediately after the int32 dims.
    n_teff, n_logg, n_meta, n_lambda = dims_tuple
    try:
        with cube.open("rb") as f:
            f.seek(4 * 4)
            teff = np.fromfile(f, dtype=np.float64, count=n_teff)
            logg = np.fromfile(f, dtype=np.float64, count=n_logg)
            meta = np.fromfile(f, dtype=np.float64, count=n_meta)
            # wavelengths follow, but not needed for this diagnostic
        if len(teff) == n_teff and np.all(np.isfinite(teff)):
            info.teff_grid = teff
        if len(logg) == n_logg and np.all(np.isfinite(logg)):
            info.logg_grid = logg
        if len(meta) == n_meta and np.all(np.isfinite(meta)):
            info.meta_grid = meta
    except Exception:
        pass

    return info


def bracket_value(grid: np.ndarray | None, value: float) -> tuple[float | None, float | None, float | None]:
    if grid is None or len(grid) == 0 or not np.isfinite(value):
        return None, None, None
    grid = np.asarray(grid, dtype=float)
    nearest = float(grid[np.argmin(np.abs(grid - value))])
    lower_vals = grid[grid <= value]
    upper_vals = grid[grid >= value]
    lower = float(lower_vals[-1]) if len(lower_vals) else None
    upper = float(upper_vals[0]) if len(upper_vals) else None
    return lower, upper, nearest


def collect_series(h: dict[str, np.ndarray]) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    bands = {b: h[b] for b in BANDS if b in h}
    colors = {}
    for b1, b2 in COLOR_PAIRS:
        if b1 in bands and b2 in bands:
            colors[f"{b1}-{b2}"] = bands[b1] - bands[b2]
    return bands, colors


def top_jumps(
    h: dict[str, np.ndarray],
    cooling_age: np.ndarray,
    teff: np.ndarray,
    logg: np.ndarray | None,
    series: dict[str, np.ndarray],
    grid: GridInfo,
    top_n: int,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    model_col = find_col(h, ["model_number"])
    logl_col = find_col(h, ["log_L", "logL", "log_l"])
    magbol_col = find_col(h, ["Mag_bol", "mag_bol", "Mbol", "M_bol"])
    interp_col = find_col(h, ["Interp_rad", "interp_rad", "interp_radius"])

    model = h[model_col].astype(int) if model_col else np.arange(len(cooling_age))
    logl = h[logl_col] if logl_col else np.full_like(cooling_age, np.nan)
    magbol = h[magbol_col] if magbol_col else np.full_like(cooling_age, np.nan)
    interp = h[interp_col] if interp_col else np.full_like(cooling_age, np.nan)
    logg_arr = logg if logg is not None else np.full_like(cooling_age, np.nan)

    for name, y in series.items():
        if len(y) < 2:
            continue
        dy = np.diff(y)
        dt_gyr = np.diff(cooling_age) / 1e9
        with np.errstate(divide="ignore", invalid="ignore"):
            slope = dy / dt_gyr
        score = np.abs(dy)
        finite = np.isfinite(score) & np.isfinite(slope)
        idxs = np.where(finite)[0]
        idxs = idxs[np.argsort(score[idxs])[::-1]][:top_n]

        for i in idxs:
            j = i + 1
            tlo, thi, tnear = bracket_value(grid.teff_grid, float(teff[j]))
            glo, ghi, gnear = bracket_value(grid.logg_grid, float(logg_arr[j]))
            rows.append({
                "series": name,
                "i_before": int(i),
                "i_after": int(j),
                "model_before": int(model[i]),
                "model_after": int(model[j]),
                "cooling_age_before_yr": float(cooling_age[i]),
                "cooling_age_after_yr": float(cooling_age[j]),
                "cooling_age_after_gyr": float(cooling_age[j] / 1e9),
                "Teff_after_K": float(teff[j]),
                "logg_after": float(logg_arr[j]),
                "logL_after": float(logl[j]),
                "Mag_bol_after": float(magbol[j]),
                "Interp_rad_after": float(interp[j]),
                "value_before": float(y[i]),
                "value_after": float(y[j]),
                "delta_value": float(dy[i]),
                "abs_delta_value": float(abs(dy[i])),
                "slope_mag_per_Gyr": float(slope[i]),
                "teff_grid_lower": tlo,
                "teff_grid_upper": thi,
                "teff_grid_nearest": tnear,
                "logg_grid_lower": glo,
                "logg_grid_upper": ghi,
                "logg_grid_nearest": gnear,
            })

    rows.sort(key=lambda r: float(r["abs_delta_value"]), reverse=True)
    return rows[:top_n]


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("No jumps found.\n")
        return
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_report(
    path: Path,
    history_path: Path,
    h_raw: dict[str, np.ndarray],
    h: dict[str, np.ndarray],
    cleaning_notes: list[str],
    grid: GridInfo,
    age: np.ndarray,
    cooling_age: np.ndarray,
    teff: np.ndarray,
    logg: np.ndarray | None,
    bands: dict[str, np.ndarray],
    colors: dict[str, np.ndarray],
    jumps: list[dict[str, object]],
) -> None:
    lines: list[str] = []
    lines.append("WD photometry diagnostic report")
    lines.append("=" * 34)
    lines.append(f"history: {history_path}")
    lines.append(f"raw rows: {len(next(iter(h_raw.values())))}")
    lines.append(f"clean rows used for diagnostics: {len(age)}")
    lines.append("")
    lines.append("History cleaning")
    lines.append("----------------")
    lines.extend(f"- {note}" for note in cleaning_notes)
    lines.append("")

    lines.append("Available diagnostics")
    lines.append("---------------------")
    lines.append(f"bands found: {', '.join(bands) if bands else 'none'}")
    lines.append(f"colors built: {', '.join(colors) if colors else 'none'}")
    lines.append(f"cooling_age range: {np.nanmin(cooling_age):.6e} -> {np.nanmax(cooling_age):.6e} yr")
    lines.append(f"Teff range: {np.nanmin(teff):.3f} -> {np.nanmax(teff):.3f} K")
    if logg is not None:
        lines.append(f"model logg range: {np.nanmin(logg):.6f} -> {np.nanmax(logg):.6f}")
        lines.append(f"model logg span: {np.nanmax(logg) - np.nanmin(logg):.6f} dex")
    else:
        lines.append("model logg column: NOT FOUND")
    lines.append("")

    lines.append("Atmosphere grid")
    lines.append("---------------")
    if grid.lookup_path:
        lines.append(f"lookup_table.csv: {grid.lookup_path}")
    else:
        lines.append("lookup_table.csv: not found/read")
    if grid.cube_path:
        lines.append(f"flux_cube.bin: {grid.cube_path}")
    else:
        lines.append("flux_cube.bin: not found/read")
    if grid.cube_dims:
        lines.append(f"flux cube dims: {grid.cube_dims} = (n_teff, n_logg, n_meta, n_lambda)")
    if grid.teff_grid is not None:
        lines.append(f"Teff grid: n={len(grid.teff_grid)}, range={np.nanmin(grid.teff_grid):.3f}->{np.nanmax(grid.teff_grid):.3f}")
    else:
        lines.append("Teff grid: unavailable")
    if grid.logg_grid is not None:
        lines.append(f"logg grid: n={len(grid.logg_grid)}, range={np.nanmin(grid.logg_grid):.3f}->{np.nanmax(grid.logg_grid):.3f}")
    else:
        lines.append("logg grid: unavailable")
    if grid.meta_grid is not None:
        lines.append(f"meta grid: n={len(grid.meta_grid)}, range={np.nanmin(grid.meta_grid):.3f}->{np.nanmax(grid.meta_grid):.3f}")
    else:
        lines.append("meta grid: unavailable")
    lines.append("")

    lines.append("Largest adjacent jumps")
    lines.append("----------------------")
    if not jumps:
        lines.append("No jumps found.")
    else:
        for r in jumps[:10]:
            lines.append(
                f"{r['series']:>6s}: model {r['model_before']} -> {r['model_after']}, "
                f"age={float(r['cooling_age_after_gyr']):.6f} Gyr, "
                f"Teff={float(r['Teff_after_K']):.1f} K, "
                f"logg={float(r['logg_after']):.5f}, "
                f"delta={float(r['delta_value']):+.6f} mag"
            )
    lines.append("")
    lines.append("Interpretation hints")
    lines.append("--------------------")
    lines.append("- If Mag_bol/log_L are smooth but LSST bands/colors jump, the feature is atmospheric/filter-related.")
    lines.append("- If Teff or log_L jumps at the same place, the feature is in the cooling model/history, not Colors.")
    lines.append("- If the jump occurs exactly at a Teff or logg grid node, inspect the adjacent Koester SEDs/cube interpolation.")
    lines.append("- A small logg span is expected for a single-mass WD cooling sequence; that can make logg visually invisible in the CMD.")
    lines.append("")

    path.write_text("\n".join(lines))


def mark_ages(ax, ages_gyr: list[float], **kwargs) -> None:
    for x in ages_gyr:
        ax.axvline(x, **kwargs)


def plot_structure_vs_age(
    out: Path,
    cooling_age: np.ndarray,
    teff: np.ndarray,
    logg: np.ndarray | None,
    h: dict[str, np.ndarray],
    feature_ages_gyr: list[float],
) -> None:
    logl_col = find_col(h, ["log_L", "logL", "log_l"])
    magbol_col = find_col(h, ["Mag_bol", "mag_bol", "Mbol", "M_bol"])
    interp_col = find_col(h, ["Interp_rad", "interp_rad", "interp_radius"])

    panels = [(r"$T_{\rm eff}$ (K)", teff)]
    if logg is not None:
        panels.append((r"$\log g$", logg))
    if logl_col is not None:
        panels.append((r"$\log L$", h[logl_col]))
    if magbol_col is not None:
        panels.append((r"$M_{\rm bol}$", h[magbol_col]))
    if interp_col is not None:
        panels.append(("Interp_rad", h[interp_col]))

    fig, axes = plt.subplots(len(panels), 1, figsize=(6.5, 1.65 * len(panels)), sharex=True)
    if len(panels) == 1:
        axes = [axes]

    x = cooling_age / 1e9
    for ax, (ylabel, y) in zip(axes, panels):
        ax.plot(x, y, lw=1.0)
        mark_ages(ax, feature_ages_gyr, color="0.4", ls="--", lw=0.7, alpha=0.6)
        ax.set_ylabel(ylabel)
        if "M_{\\rm bol}" in ylabel:
            ax.invert_yaxis()
    axes[-1].set_xlabel("Cooling age (Gyr)")
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def plot_mags_colors_vs_age(
    out: Path,
    cooling_age: np.ndarray,
    bands: dict[str, np.ndarray],
    colors: dict[str, np.ndarray],
    feature_ages_gyr: list[float],
) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(6.5, 6.0), sharex=True)
    x = cooling_age / 1e9

    for b, y in bands.items():
        axes[0].plot(x, y, lw=1.0, label=b)
    axes[0].invert_yaxis()
    axes[0].set_ylabel("AB mag")
    axes[0].legend(ncol=3)
    mark_ages(axes[0], feature_ages_gyr, color="0.4", ls="--", lw=0.7, alpha=0.6)

    for name, y in colors.items():
        axes[1].plot(x, y, lw=1.0, label=name)
    axes[1].set_ylabel("color (mag)")
    axes[1].set_xlabel("Cooling age (Gyr)")
    axes[1].legend(ncol=3)
    mark_ages(axes[1], feature_ages_gyr, color="0.4", ls="--", lw=0.7, alpha=0.6)

    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def plot_color_teff_raw(
    out: Path,
    teff: np.ndarray,
    colors: dict[str, np.ndarray],
    jumps: list[dict[str, object]],
) -> None:
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    for name, y in colors.items():
        ax.plot(teff, y, lw=1.0, label=name)

    # Mark the top jumps by color series.
    for r in jumps[:10]:
        if str(r["series"]) in colors:
            j = int(r["i_after"])
            name = str(r["series"])
            ax.plot(teff[j], colors[name][j], marker="o", ms=3.5, mec="k", mfc="none", mew=0.6)

    ax.set_xscale("log")
    ax.invert_xaxis()
    ax.set_xticks([5000, 10000, 20000, 40000, 80000])
    ax.get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    ax.set_xlabel(r"$T_{\rm eff}$ (K)")
    ax.set_ylabel("raw color (AB mag)")
    ax.legend(ncol=2)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def plot_cmd_colored_by_logg(
    out: Path,
    bands: dict[str, np.ndarray],
    logg: np.ndarray | None,
    cooling_age: np.ndarray,
) -> None:
    if "g" not in bands or "r" not in bands:
        return

    x = bands["g"] - bands["r"]
    y = bands["g"]
    c = logg if logg is not None else np.log10(np.clip(cooling_age, 1e3, None))
    label = r"model $\log g$" if logg is not None else r"$\log_{10}$(cooling age / yr)"

    fig, ax = plt.subplots(figsize=(4.3, 5.0))
    sc = ax.scatter(x, y, c=c, s=8, linewidths=0, cmap="viridis")
    ax.invert_yaxis()
    ax.set_xlabel(r"$(g-r)_{\rm AB}$")
    ax.set_ylabel(r"$M_g$ (AB mag)")
    cb = fig.colorbar(sc, ax=ax, pad=0.02)
    cb.set_label(label)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Diagnose WD photometry features in MESA Colors output.")
    parser.add_argument("--history", type=Path, default=DEFAULT_HISTORY, help="Path to MESA LOGS/history.data")
    parser.add_argument("--grid-dir", type=Path, default=DEFAULT_GRID_DIR, help="Atmosphere grid directory, e.g. koester2")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR, help="Output directory for diagnostics")
    parser.add_argument("--top", type=int, default=30, help="Number of largest jumps to report")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading history: {args.history}")
    h_raw = read_history(args.history)
    h, cleaning_notes = clean_history(h_raw)
    for note in cleaning_notes:
        print(note)

    age_col = find_col(h, ["star_age", "age"], required=True)
    logteff_col = find_col(h, ["log_Teff", "logTeff", "log_teff"], required=True)
    logg_col = find_col(h, ["log_g", "logg", "surface_log_g"])

    age = h[age_col]
    log_teff = h[logteff_col]
    teff = 10.0 ** log_teff
    cooling_age = age - age[0]
    logg = h[logg_col] if logg_col is not None else None

    bands, colors = collect_series(h)
    if len(bands) < 2:
        raise RuntimeError(f"Found too few photometry bands. Available matched bands: {list(bands)}")

    # Remove rows with bad common values. Keep this simple and explicit.
    common_arrays = [age, teff] + list(bands.values()) + list(colors.values())
    if logg is not None:
        common_arrays.append(logg)
    mask = finite_common_mask(common_arrays)
    if not np.all(mask):
        print(f"Removing {np.count_nonzero(~mask)} rows with non-finite diagnostic values.")
        h = {c: v[mask] for c, v in h.items()}
        age = age[mask]
        teff = teff[mask]
        cooling_age = cooling_age[mask]
        if logg is not None:
            logg = logg[mask]
        bands, colors = collect_series(h)

    print(f"Rows used: {len(age)}")
    print(f"Cooling-age range: {cooling_age.min():.4e} -> {cooling_age.max():.4e} yr")
    print(f"Teff range: {teff.min():.1f} -> {teff.max():.1f} K")
    if logg is not None:
        print(f"model logg range: {logg.min():.5f} -> {logg.max():.5f}  span={logg.max() - logg.min():.5f} dex")
    else:
        print("WARNING: no model logg column found in history.data")

    print(f"Reading grid diagnostics from: {args.grid_dir}")
    grid = read_lookup_grid(args.grid_dir)
    grid = read_cube_header(args.grid_dir, grid)
    if grid.cube_dims:
        print(f"flux_cube.bin dims: {grid.cube_dims} = (n_teff, n_logg, n_meta, n_lambda)")
    if grid.teff_grid is not None:
        print(f"Teff grid n={len(grid.teff_grid)} range={grid.teff_grid.min():.1f}->{grid.teff_grid.max():.1f}")
    if grid.logg_grid is not None:
        print(f"logg grid n={len(grid.logg_grid)} range={grid.logg_grid.min():.2f}->{grid.logg_grid.max():.2f}")

    all_series = {**{f"mag_{b}": v for b, v in bands.items()}, **colors}
    magbol_col = find_col(h, ["Mag_bol", "mag_bol", "Mbol", "M_bol"])
    if magbol_col is not None:
        all_series["Mag_bol"] = h[magbol_col]

    jumps = top_jumps(h, cooling_age, teff, logg, all_series, grid, args.top)
    feature_ages_gyr = []
    for r in jumps[:5]:
        x = float(r["cooling_age_after_gyr"])
        if all(abs(x - y) > 0.02 for y in feature_ages_gyr):
            feature_ages_gyr.append(x)

    csv_path = args.out_dir / "wd_diagnostic_top_jumps.csv"
    report_path = args.out_dir / "wd_diagnostic_report.txt"
    write_csv(csv_path, jumps)
    write_report(report_path, args.history, h_raw, h, cleaning_notes, grid, age, cooling_age, teff, logg, bands, colors, jumps)

    print(f"Writing plots to: {args.out_dir}")
    plot_structure_vs_age(args.out_dir / "wd_diagnostic_structure_vs_age.pdf", cooling_age, teff, logg, h, feature_ages_gyr)
    plot_mags_colors_vs_age(args.out_dir / "wd_diagnostic_mags_colors_vs_age.pdf", cooling_age, bands, colors, feature_ages_gyr)
    plot_color_teff_raw(args.out_dir / "wd_diagnostic_color_teff_raw.pdf", teff, colors, jumps)
    plot_cmd_colored_by_logg(args.out_dir / "wd_diagnostic_cmd_colored_by_logg.pdf", bands, logg, cooling_age)

    print("\nLargest adjacent jumps:")
    for r in jumps[:10]:
        print(
            f"  {str(r['series']):>7s}  "
            f"model {int(r['model_before'])}->{int(r['model_after'])}  "
            f"age={float(r['cooling_age_after_gyr']):.5f} Gyr  "
            f"Teff={float(r['Teff_after_K']):.0f} K  "
            f"logg={float(r['logg_after']):.4f}  "
            f"delta={float(r['delta_value']):+.5f} mag  "
            f"Tgrid=({r['teff_grid_lower']},{r['teff_grid_upper']})  "
            f"ggrid=({r['logg_grid_lower']},{r['logg_grid_upper']})"
        )

    print(f"\nSaved report: {report_path}")
    print(f"Saved jump table: {csv_path}")
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
