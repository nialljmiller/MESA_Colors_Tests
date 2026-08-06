#!/usr/bin/env python3
"""
diagnose_runtime_vs_independent_cube.py

Compare the MESA Colors history columns against an independent Python
interpolation/convolution through the same Koester flux_cube.bin.

Purpose
-------
If the raw Koester SEDs and flux_cube.bin are smooth, but MESA Colors history
magnitudes jump, this script checks whether a direct Python evaluation of the
cube at the same history rows is smooth. That localizes the problem to the
runtime Colors path rather than the atmosphere data.

Example
-------
python diagnose_runtime_vs_independent_cube.py \
  --history ../LOGS/history.data \
  --grid-dir /home/njm/SED_Tools/data/stellar_models/koester2 \
  --instrument /home/njm/SED_Tools/data/filters/LSST/LSST \
  --out-dir figures/runtime_vs_independent_cube
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

R_SUN_CGS = 6.957e10
DEFAULT_DISTANCE = 3.0857e19
BANDS_DEFAULT = ["u", "g", "r", "i", "z", "y"]
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


def read_cube(grid_dir: Path, order: str = "F") -> Dict[str, np.ndarray]:
    path = grid_dir / "flux_cube.bin"
    with path.open("rb") as f:
        dims = np.fromfile(f, dtype=np.int32, count=4)
        if len(dims) != 4:
            raise RuntimeError("Could not read cube dimensions")
        nteff, nlogg, nmeta, nlam = map(int, dims)
        teff = np.fromfile(f, dtype=np.float64, count=nteff)
        logg = np.fromfile(f, dtype=np.float64, count=nlogg)
        meta = np.fromfile(f, dtype=np.float64, count=nmeta)
        wave = np.fromfile(f, dtype=np.float64, count=nlam)
        flat = np.fromfile(f, dtype=np.float64)
    expected = nteff * nlogg * nmeta * nlam
    if flat.size != expected:
        raise RuntimeError(f"Cube payload mismatch: got {flat.size}, expected {expected}")
    cube = flat.reshape((nteff, nlogg, nmeta, nlam), order=order)
    return {"dims": dims, "teff": teff, "logg": logg, "meta": meta, "wave": wave, "cube": cube}


def parse_two_col(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    xs: List[float] = []
    ys: List[float] = []
    with path.open("r", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.replace(",", " ").split()
            if len(parts) < 2:
                continue
            try:
                x = float(parts[0])
                y = float(parts[1])
            except ValueError:
                continue
            if np.isfinite(x) and np.isfinite(y):
                xs.append(x)
                ys.append(y)
    if len(xs) < 4:
        raise RuntimeError(f"Could not parse numeric two-column file: {path}")
    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)
    idx = np.argsort(x)
    return x[idx], y[idx]


def discover_filters(instrument: Path) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    index = instrument / instrument.name
    files: List[Path] = []
    if index.exists():
        for line in index.read_text(errors="ignore").splitlines():
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            p = Path(s)
            if not p.is_absolute():
                p = instrument / p
            files.append(p)
    else:
        files = sorted(instrument.glob("*.dat"))
    filters: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    for p in files:
        if not p.exists():
            print(f"WARNING: filter file listed but missing: {p}")
            continue
        name = p.name.split(".")[0]
        filters[name] = parse_two_col(p)
    if not filters:
        raise RuntimeError(f"No filters found in {instrument}")
    return filters


def lower_bracket_index(grid: np.ndarray, x: float) -> Tuple[int, int, float]:
    """Return lo, hi, frac for ascending grid, clamped to boundaries."""
    if x <= grid[0]:
        return 0, 0, 0.0
    if x >= grid[-1]:
        n = len(grid) - 1
        return n, n, 0.0
    hi = int(np.searchsorted(grid, x, side="left"))
    lo = hi - 1
    dx = grid[hi] - grid[lo]
    frac = 0.0 if dx == 0 else (x - grid[lo]) / dx
    return lo, hi, float(frac)


def interp_cube_sed(cube_info: Dict[str, np.ndarray], teff: float, logg: float, meta: float = 0.0) -> np.ndarray:
    Tg = cube_info["teff"]
    gg = cube_info["logg"]
    mg = cube_info["meta"]
    C = cube_info["cube"]

    i0, i1, ft = lower_bracket_index(Tg, teff)
    j0, j1, fg = lower_bracket_index(gg, logg)
    k0, k1, fm = lower_bracket_index(mg, meta)

    sed = np.zeros(C.shape[-1], dtype=float)
    for ii, wt in [(i0, 1.0 - ft), (i1, ft)]:
        if wt == 0.0 and i0 != i1:
            continue
        for jj, wg in [(j0, 1.0 - fg), (j1, fg)]:
            if wg == 0.0 and j0 != j1:
                continue
            for kk, wm in [(k0, 1.0 - fm), (k1, fm)]:
                if wm == 0.0 and k0 != k1:
                    continue
                sed += wt * wg * wm * C[ii, jj, kk, :]
    return sed


def trapz_positive(wave: np.ndarray, flux: np.ndarray) -> float:
    good = np.isfinite(wave) & np.isfinite(flux) & (wave > 0) & (flux > 0)
    if good.sum() < 4:
        return np.nan
    return float(np.trapezoid(flux[good], wave[good]))


def weighted_filter_flux(wave: np.ndarray, flux: np.ndarray, fw: np.ndarray, ft: np.ndarray) -> float:
    S = np.interp(wave, fw, ft, left=0.0, right=0.0)
    good = np.isfinite(wave) & np.isfinite(flux) & np.isfinite(S) & (S > 0)
    if good.sum() < 4:
        return np.nan
    num = np.trapezoid(flux[good] * S[good], wave[good])
    den = np.trapezoid(S[good], wave[good])
    if den <= 0 or num <= 0:
        return np.nan
    return float(num / den)


def mag_rel(flux: float) -> float:
    if not np.isfinite(flux) or flux <= 0:
        return np.nan
    return float(-2.5 * np.log10(flux))


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


def align_to_history(history: pd.Series, independent: pd.Series, exclude_idx: int, window: int = 6) -> Tuple[pd.Series, float]:
    mask = np.isfinite(history.to_numpy(float)) & np.isfinite(independent.to_numpy(float))
    lo = max(0, exclude_idx - window)
    hi = min(len(history), exclude_idx + window + 1)
    mask[lo:hi] = False
    if mask.sum() < 5:
        mask = np.isfinite(history.to_numpy(float)) & np.isfinite(independent.to_numpy(float))
    offset = float(np.nanmedian(history.to_numpy(float)[mask] - independent.to_numpy(float)[mask]))
    return independent + offset, offset


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--history", type=Path, default=Path("../LOGS/history.data"))
    ap.add_argument("--grid-dir", type=Path, required=True)
    ap.add_argument("--instrument", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, default=Path("figures/runtime_vs_independent_cube"))
    ap.add_argument("--distance", type=float, default=DEFAULT_DISTANCE)
    ap.add_argument("--cube-order", choices=["F", "C"], default="F", help="Fortran order is expected for flux_cube.bin")
    ap.add_argument("--meta", type=float, default=0.0)
    ap.add_argument("--bands", nargs="*", default=BANDS_DEFAULT)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading history: {args.history}")
    h = read_history(args.history)
    h = h.copy()
    h["cooling_age"] = h["star_age"] - h["star_age"].iloc[0]
    h["Teff"] = 10.0 ** h["log_Teff"]
    if "log_g" not in h.columns:
        raise RuntimeError("history.data needs log_g for this diagnostic")
    if "log_R" not in h.columns:
        raise RuntimeError("history.data needs log_R for this diagnostic")

    bands = [b for b in args.bands if b in h.columns]
    if not bands:
        raise RuntimeError(f"None of requested bands found in history: {args.bands}")

    best_col, best_idx, best_abs = find_largest_jump(h, ["Mag_bol"] + bands)
    print(f"Largest MESA Colors jump: {best_col} row {best_idx-1}->{best_idx}, |delta|={best_abs:.6g}")
    print(
        f"Suspect model {int(h.loc[best_idx, 'model_number'])}: "
        f"age={h.loc[best_idx, 'cooling_age']/1e9:.9f} Gyr, "
        f"Teff={h.loc[best_idx, 'Teff']:.2f} K, logg={h.loc[best_idx, 'log_g']:.6f}"
    )

    print(f"Reading cube: {args.grid_dir / 'flux_cube.bin'}")
    cube_info = read_cube(args.grid_dir, order=args.cube_order)
    print(f"cube dims: {tuple(map(int, cube_info['dims']))} = (n_teff, n_logg, n_meta, n_lambda)")
    print(f"Teff grid: {cube_info['teff'].min():.1f}->{cube_info['teff'].max():.1f}; logg grid: {cube_info['logg'].min():.2f}->{cube_info['logg'].max():.2f}")

    print(f"Reading filters: {args.instrument}")
    filters = discover_filters(args.instrument)
    filters = {k: v for k, v in filters.items() if k in bands}
    if not filters:
        raise RuntimeError("No matching filters between instrument and history bands")
    print("filters used:", ", ".join(filters))

    wave = cube_info["wave"]
    independent: Dict[str, List[float]] = {"ind_Mbol_rel": []}
    for b in filters:
        independent[f"ind_{b}_rel"] = []

    print("Evaluating independent cube photometry for every history row...")
    for _, row in h.iterrows():
        sed_surface = interp_cube_sed(cube_info, float(row["Teff"]), float(row["log_g"]), args.meta)
        radius = (10.0 ** float(row["log_R"])) * R_SUN_CGS
        dilution = (radius / args.distance) ** 2
        sed_obs = sed_surface * dilution
        fbol = trapz_positive(wave, sed_obs)
        independent["ind_Mbol_rel"].append(mag_rel(fbol))
        for b, (fw, ft) in filters.items():
            fb = weighted_filter_flux(wave, sed_obs, fw, ft)
            independent[f"ind_{b}_rel"].append(mag_rel(fb))

    for k, vals in independent.items():
        h[k] = np.asarray(vals, dtype=float)

    # Align arbitrary relative magnitudes to MESA history columns so residuals are easy to inspect.
    align_cols = []
    if "Mag_bol" in h.columns:
        aligned, off = align_to_history(h["Mag_bol"], h["ind_Mbol_rel"], best_idx)
        h["ind_Mbol_aligned"] = aligned
        h["resid_Mbol"] = h["Mag_bol"] - h["ind_Mbol_aligned"]
        align_cols.append(("Mag_bol", "ind_Mbol_aligned", "resid_Mbol", off))
    for b in filters:
        aligned, off = align_to_history(h[b], h[f"ind_{b}_rel"], best_idx)
        h[f"ind_{b}_aligned"] = aligned
        h[f"resid_{b}"] = h[b] - h[f"ind_{b}_aligned"]
        align_cols.append((b, f"ind_{b}_aligned", f"resid_{b}", off))

    # Local comparison table.
    lo = max(0, best_idx - 6)
    hi = min(len(h), best_idx + 7)
    local = h.iloc[lo:hi].copy()
    local_cols = ["model_number", "cooling_age", "Teff", "log_Teff", "log_L", "log_R", "log_g"]
    for c, a, r, _ in align_cols:
        local_cols.extend([c, a, r])
    local_cols = [c for c in local_cols if c in local.columns]
    local.to_csv(args.out_dir / "runtime_vs_independent_local.csv", index=False, columns=local_cols)

    # Delta table for the exact suspect step.
    delta_rows = []
    for c, a, r, off in align_cols:
        d_hist = float(h.loc[best_idx, c] - h.loc[best_idx - 1, c])
        d_ind = float(h.loc[best_idx, a] - h.loc[best_idx - 1, a])
        d_res = float(h.loc[best_idx, r] - h.loc[best_idx - 1, r])
        delta_rows.append({
            "quantity": c,
            "history_delta": d_hist,
            "independent_cube_delta": d_ind,
            "residual_delta_history_minus_independent": d_res,
            "alignment_offset": off,
        })
    delta_df = pd.DataFrame(delta_rows)
    delta_df.to_csv(args.out_dir / "runtime_vs_independent_suspect_step_deltas.csv", index=False)

    print("\nSuspect-step comparison: history delta vs independent-cube delta")
    print(delta_df.to_string(index=False))

    # Figure 1: local aligned comparison.
    x = local["cooling_age"] / 1e9
    n_panels = 3
    fig, axes = plt.subplots(n_panels, 1, figsize=(6.0, 6.2), sharex=True)
    if "Mag_bol" in local.columns:
        axes[0].plot(x, local["Mag_bol"], "o-", ms=3, lw=1.0, label="MESA Colors Mag_bol")
        axes[0].plot(x, local["ind_Mbol_aligned"], "s--", ms=3, lw=1.0, label="independent cube, aligned")
        axes[0].invert_yaxis()
        axes[0].set_ylabel(r"$M_{bol}$")
        axes[0].legend(framealpha=0.8)
    for b in filters:
        axes[1].plot(x, local[b], "o-", ms=3, lw=1.0, label=f"MESA {b}")
        axes[1].plot(x, local[f"ind_{b}_aligned"], "--", lw=0.9, alpha=0.8, label=f"ind {b}")
    axes[1].invert_yaxis()
    axes[1].set_ylabel("AB mag")
    axes[1].legend(ncol=3, framealpha=0.75)
    for c, a, r, _ in align_cols:
        axes[2].plot(x, local[r], "o-", ms=3, lw=1.0, label=c)
    axes[2].axhline(0, color="0.5", lw=0.8)
    axes[2].set_ylabel("history - independent")
    axes[2].set_xlabel("Cooling age (Gyr)")
    axes[2].legend(ncol=4, framealpha=0.75)
    for ax in axes:
        ax.axvline(h.loc[best_idx, "cooling_age"] / 1e9, ls="--", color="0.55", lw=0.8)
    fig.tight_layout()
    fig.savefig(args.out_dir / "runtime_vs_independent_local_comparison.pdf", bbox_inches="tight")
    plt.close(fig)

    # Figure 2: global residuals.
    fig, axes = plt.subplots(2, 1, figsize=(6.0, 5.2), sharex=True)
    gx = h["cooling_age"] / 1e9
    if "resid_Mbol" in h.columns:
        axes[0].plot(gx, h["resid_Mbol"], lw=1.0, label="Mbol")
    for b in filters:
        axes[0].plot(gx, h[f"resid_{b}"], lw=1.0, label=b)
    axes[0].axhline(0, color="0.5", lw=0.8)
    axes[0].set_ylabel("history - independent")
    axes[0].legend(ncol=4, framealpha=0.8)
    for b in filters:
        axes[1].plot(gx, h[f"ind_{b}_aligned"], lw=1.0, label=f"ind {b}")
    axes[1].invert_yaxis()
    axes[1].set_ylabel("independent aligned mag")
    axes[1].set_xlabel("Cooling age (Gyr)")
    axes[1].legend(ncol=3, framealpha=0.8)
    for ax in axes:
        ax.axvline(h.loc[best_idx, "cooling_age"] / 1e9, ls="--", color="0.55", lw=0.8)
    fig.tight_layout()
    fig.savefig(args.out_dir / "runtime_vs_independent_global_residuals.pdf", bbox_inches="tight")
    plt.close(fig)

    report = []
    report.append("Runtime MESA Colors vs independent cube diagnostic")
    report.append("===================================================")
    report.append(f"history: {args.history}")
    report.append(f"grid_dir: {args.grid_dir}")
    report.append(f"instrument: {args.instrument}")
    report.append(f"rows used: {len(h)}")
    report.append(f"cube dims: {tuple(map(int, cube_info['dims']))} = (n_teff, n_logg, n_meta, n_lambda)")
    report.append(f"largest MESA Colors jump: {best_col} row {best_idx-1}->{best_idx}, |delta|={best_abs:.9f}")
    report.append(f"suspect model after step: {int(h.loc[best_idx, 'model_number'])}")
    report.append(f"suspect age Gyr: {h.loc[best_idx, 'cooling_age']/1e9:.9f}")
    report.append(f"suspect Teff: {h.loc[best_idx, 'Teff']:.6f}")
    report.append(f"suspect logg: {h.loc[best_idx, 'log_g']:.6f}")
    report.append("")
    report.append("Suspect-step deltas")
    report.append("-------------------")
    report.append(delta_df.to_string(index=False))
    report.append("")
    report.append("Interpretation")
    report.append("--------------")
    report.append("If independent_cube_delta is smooth/small while history_delta is large, the atmosphere grid and cube are not the source of the jump.")
    report.append("That points to the runtime MESA Colors interpolation/convolution/scaling path.")
    report.append("If independent_cube_delta is also large, the Python interpolation is reproducing the issue and the next step is to inspect the cube interpolation formula or filter convolution.")
    report.append("")
    report.append("Files written")
    report.append("-------------")
    for name in [
        "runtime_vs_independent_local.csv",
        "runtime_vs_independent_suspect_step_deltas.csv",
        "runtime_vs_independent_local_comparison.pdf",
        "runtime_vs_independent_global_residuals.pdf",
    ]:
        report.append(str(args.out_dir / name))
    (args.out_dir / "runtime_vs_independent_report.txt").write_text("\n".join(report) + "\n")

    print(f"\nSaved report: {args.out_dir / 'runtime_vs_independent_report.txt'}")
    print(f"Saved local table: {args.out_dir / 'runtime_vs_independent_local.csv'}")
    print("Done.")


if __name__ == "__main__":
    main()
