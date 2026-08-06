#!/usr/bin/env python3
"""
verify_rrlyrae_physical.py

Quick physical/sanity verification for a MESA RSP RR Lyrae + Colors run.

Run from the MESA work directory, e.g.

    cd ~/MESA/MESA_Colors_Tests/rsp_RR_Lyrae
    python verify_rrlyrae_physical.py

It checks whether the run looks like a real settled RR Lyrae simulation rather
than a bookkeeping/plotting artefact. It is intentionally diagnostic, not a
paper-quality validation suite.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd

try:
    import matplotlib.pyplot as plt
except Exception:  # pragma: no cover
    plt = None

TSUN = 5772.0

# Effective wavelengths in Angstrom, only for ordering/amplitude sanity checks.
FILTER_LAMBDA_A = {
    # Roman/WFI-ish
    "F062": 6200.0,
    "F087": 8700.0,
    "F106": 10600.0,
    "F129": 12900.0,
    "F146": 14600.0,
    "F158": 15800.0,
    "F184": 18400.0,
    "F213": 21300.0,
    # LSST-ish
    "u": 3671.0,
    "g": 4827.0,
    "r": 6223.0,
    "i": 7546.0,
    "z": 8691.0,
    "y": 9711.0,
    # Gaia-ish
    "Gbp": 5320.0,
    "G": 6730.0,
    "Grp": 7970.0,
    # Other common broad bands
    "U": 3650.0,
    "B": 4450.0,
    "V": 5510.0,
    "R": 6580.0,
    "I": 8060.0,
    "Kepler": 6400.0,
    "TESS": 7865.0,
}

BAD_SENTINELS = {-99.0, -999.0, 99.0, 999.0}


@dataclass
class CheckResult:
    name: str
    status: str  # PASS / WARN / FAIL / INFO
    message: str


def status_symbol(status: str) -> str:
    return {"PASS": "✓", "WARN": "!", "FAIL": "✗", "INFO": "·"}.get(status, "?")


def add(results: list[CheckResult], name: str, status: str, message: str) -> None:
    results.append(CheckResult(name=name, status=status, message=message))


def clean_numeric_series(s: pd.Series) -> pd.Series:
    out = pd.to_numeric(s, errors="coerce")
    for bad in BAD_SENTINELS:
        out = out.mask(np.isclose(out, bad, rtol=0, atol=1e-12))
    return out


def read_mesa_history(path: Path) -> pd.DataFrame:
    text = path.read_text(errors="ignore").splitlines()
    header_idx = None
    for i, line in enumerate(text):
        toks = line.split()
        if "model_number" in toks:
            header_idx = i
            break
    if header_idx is None:
        raise ValueError(f"Could not find history header containing model_number in {path}")
    df = pd.read_csv(path, sep=r"\s+", skiprows=header_idx, comment="!", engine="python")
    # MESA may repeat headers after restarts; coerce model_number and drop repeats.
    if "model_number" in df.columns:
        df["model_number"] = pd.to_numeric(df["model_number"], errors="coerce")
        df = df.dropna(subset=["model_number"]).copy()
        df["model_number"] = df["model_number"].astype(int)
        df = df.drop_duplicates(subset=["model_number"], keep="last")
        df = df.sort_values("model_number").reset_index(drop=True)
    for c in df.columns:
        if c != "model_number":
            df[c] = pd.to_numeric(df[c], errors="ignore")
    return df


def find_history(explicit: Optional[str]) -> Path:
    if explicit:
        p = Path(explicit)
        if not p.exists():
            raise FileNotFoundError(p)
        return p
    candidates = [
        Path("LOGS_colors/history.data"),
        Path("LOGS/history.data"),
        Path("LOGS_settle/history.data"),
        Path("LOGS_overnight/history.data"),
        Path("LOGS_overnight/history_last100.data"),
    ]
    existing = [p for p in candidates if p.exists()]
    if not existing:
        raise FileNotFoundError("No history.data found in LOGS_colors, LOGS, LOGS_settle, or LOGS_overnight")
    return max(existing, key=lambda p: p.stat().st_mtime)


def finite_values(df: pd.DataFrame, col: str) -> np.ndarray:
    if col not in df.columns:
        return np.array([], dtype=float)
    x = clean_numeric_series(df[col]).to_numpy(dtype=float)
    return x[np.isfinite(x)]


def get_teff(df: pd.DataFrame) -> Optional[pd.Series]:
    if "Teff" in df.columns:
        return clean_numeric_series(df["Teff"])
    if "effective_T" in df.columns:
        return clean_numeric_series(df["effective_T"])
    if "log_Teff" in df.columns:
        return 10.0 ** clean_numeric_series(df["log_Teff"])
    return None


def infer_radius_rsun(df: pd.DataFrame, teff: pd.Series, log_l: pd.Series) -> tuple[Optional[pd.Series], str, Optional[float]]:
    """Infer radius in Rsun. Test candidate radius columns against Stefan-Boltzmann."""
    candidates: list[tuple[str, pd.Series]] = []
    for col in ["radius", "photosphere_r", "R", "star_radius"]:
        if col in df.columns:
            raw = clean_numeric_series(df[col])
            if raw.notna().sum() > 5:
                candidates.append((col, raw))
    for col in ["log_R", "log_R_phot", "lg_R"]:
        if col in df.columns:
            raw = 10.0 ** clean_numeric_series(df[col])
            if raw.notna().sum() > 5:
                candidates.append((col, raw))

    if not candidates:
        return None, "none", None

    best_name = "none"
    best_r = None
    best_mad = np.inf
    mask_common = np.isfinite(teff) & np.isfinite(log_l) & (teff > 0)
    for name, r_raw in candidates:
        # Try as linear radius and as log10 radius; choose what best satisfies L=4*piR^2sigmaT^4.
        options = [(name, r_raw)]
        if np.nanmedian(r_raw.to_numpy(dtype=float)) < 3.0:
            options.append((name + " interpreted as log10(R/Rsun)", 10.0 ** r_raw))
        for opt_name, r in options:
            mask = mask_common & np.isfinite(r) & (r > 0)
            if mask.sum() < 5:
                continue
            pred = 2.0 * np.log10(r[mask].to_numpy(dtype=float)) + 4.0 * np.log10(teff[mask].to_numpy(dtype=float) / TSUN)
            mad = float(np.nanmedian(np.abs(pred - log_l[mask].to_numpy(dtype=float))))
            if mad < best_mad:
                best_mad = mad
                best_name = opt_name
                best_r = r
    if best_r is None:
        return None, "none", None
    return best_r, best_name, best_mad


def select_last_periods(df: pd.DataFrame, n_periods: int) -> pd.DataFrame:
    if "rsp_num_periods" not in df.columns:
        return df.tail(min(len(df), 5000)).copy()
    periods = clean_numeric_series(df["rsp_num_periods"])
    if periods.notna().sum() == 0:
        return df.tail(min(len(df), 5000)).copy()
    pmax = int(np.nanmax(periods.to_numpy(dtype=float)))
    sub = df[periods >= (pmax - n_periods + 1)].copy()
    if len(sub) == 0:
        return df.tail(min(len(df), 5000)).copy()
    return sub


def pearson(a: Iterable[float], b: Iterable[float]) -> float:
    x = np.asarray(a, dtype=float)
    y = np.asarray(b, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3:
        return np.nan
    if np.nanstd(x[m]) == 0 or np.nanstd(y[m]) == 0:
        return np.nan
    return float(np.corrcoef(x[m], y[m])[0, 1])


def detect_filter_columns(df: pd.DataFrame) -> list[str]:
    cols = []
    for name in FILTER_LAMBDA_A:
        if name in df.columns:
            vals = finite_values(df, name)
            if len(vals) > 10 and np.nanstd(vals) > 0:
                cols.append(name)
    cols.sort(key=lambda c: FILTER_LAMBDA_A[c])
    return cols


def analyze_seds(sed_dir: Path, results: list[CheckResult], max_files_per_filter: int = 200) -> dict:
    info: dict = {"exists": sed_dir.exists(), "filters": {}, "sample_table": None}
    if not sed_dir.exists():
        add(results, "SED directory", "WARN", f"{sed_dir} does not exist; cannot verify per-model SED output.")
        return info

    files = sorted(sed_dir.glob("*_SED_*.csv"))
    if not files:
        # Maybe sed_per_model was false.
        files = sorted(sed_dir.glob("*_SED.csv"))
        if files:
            add(results, "SED per model", "WARN", f"Found {len(files)} overwritten SED files but no per-model files (*_SED_########.csv).")
        else:
            add(results, "SED files", "WARN", f"No SED CSV files found in {sed_dir}.")
        return info

    filt_to_files: dict[str, list[Path]] = {}
    pat = re.compile(r"(.+)_SED_(\d+)\.csv$")
    for p in files:
        m = pat.match(p.name)
        if not m:
            continue
        filt_to_files.setdefault(m.group(1), []).append(p)
    for k in filt_to_files:
        filt_to_files[k].sort(key=lambda p: int(p.stem.split("_")[-1]))

    info["filters"] = {k: len(v) for k, v in filt_to_files.items()}
    if not filt_to_files:
        add(results, "SED files", "FAIL", f"Found {len(files)} SED-like files but could not parse filter/model names.")
        return info

    counts = np.array([len(v) for v in filt_to_files.values()])
    add(results, "SED files", "PASS", f"Found {len(files)} per-model SED files across {len(filt_to_files)} filters; median {int(np.median(counts))} files/filter.")

    # Inspect representative files for finite positive flux and convolved flux.
    sample_rows = []
    problems = []
    for filt, flist in sorted(filt_to_files.items(), key=lambda kv: FILTER_LAMBDA_A.get(kv[0], 1e9)):
        sample = flist[:: max(1, len(flist) // min(max_files_per_filter, len(flist)))]
        n_ok = 0
        fints = []
        synth_fluxes = []
        for p in sample:
            try:
                df = pd.read_csv(p)
                df.columns = df.columns.str.strip()
                if not {"wavelengths", "fluxes"}.issubset(df.columns):
                    problems.append(f"{p.name}: missing wavelengths/fluxes")
                    continue
                w = pd.to_numeric(df["wavelengths"], errors="coerce").to_numpy(dtype=float)
                f = pd.to_numeric(df["fluxes"], errors="coerce").to_numpy(dtype=float)
                m = np.isfinite(w) & np.isfinite(f) & (w > 0)
                if m.sum() < 10:
                    problems.append(f"{p.name}: too few finite SED points")
                    continue
                if np.nanmin(f[m]) < -1e-30:
                    problems.append(f"{p.name}: negative SED flux")
                    continue
                fint = float(np.trapezoid(f[m], w[m]))
                if not np.isfinite(fint) or fint <= 0:
                    problems.append(f"{p.name}: non-positive integrated SED flux")
                    continue
                fints.append(fint)
                # Optional pseudo synthetic flux using convolved_flux if present.
                if "convolved_flux" in df.columns:
                    cf = pd.to_numeric(df["convolved_flux"], errors="coerce").to_numpy(dtype=float)
                    mm = m & np.isfinite(cf)
                    if mm.sum() > 10:
                        # filter_on_grid = convolved_flux / fluxes where possible
                        denom_grid = np.zeros_like(cf)
                        good = mm & (np.abs(f) > 0)
                        denom_grid[good] = cf[good] / f[good]
                        num = float(np.trapezoid(cf[mm] * w[mm], w[mm]))
                        den = float(np.trapezoid(denom_grid[mm] * w[mm], w[mm]))
                        if np.isfinite(num) and np.isfinite(den) and den > 0 and num > 0:
                            synth_fluxes.append(num / den)
                n_ok += 1
            except Exception as e:  # noqa: BLE001
                problems.append(f"{p.name}: {e}")
        if n_ok:
            sample_rows.append({
                "filter": filt,
                "n_files": len(flist),
                "n_sampled": n_ok,
                "Fint_min": float(np.nanmin(fints)) if fints else np.nan,
                "Fint_max": float(np.nanmax(fints)) if fints else np.nan,
                "pseudo_flux_min": float(np.nanmin(synth_fluxes)) if synth_fluxes else np.nan,
                "pseudo_flux_max": float(np.nanmax(synth_fluxes)) if synth_fluxes else np.nan,
            })

    if problems:
        add(results, "SED integrity", "WARN", f"SED parsing found {len(problems)} issues; first: {problems[0]}")
    else:
        add(results, "SED integrity", "PASS", "Sampled SED files have finite positive wavelength grids and positive integrated fluxes.")

    if sample_rows:
        info["sample_table"] = pd.DataFrame(sample_rows)
    return info


def make_plots(df: pd.DataFrame, sub: pd.DataFrame, outdir: Path, filters: list[str], sed_info: dict) -> list[Path]:
    if plt is None:
        return []
    outdir.mkdir(parents=True, exist_ok=True)
    paths = []
    teff = get_teff(sub)
    phase = clean_numeric_series(sub["rsp_phase"]) if "rsp_phase" in sub.columns else None
    if phase is None or phase.notna().sum() < 10:
        if "star_age_day" in sub.columns and "rsp_period_in_days" in sub.columns:
            age = clean_numeric_series(sub["star_age_day"])
            p = np.nanmedian(clean_numeric_series(sub["rsp_period_in_days"]))
            phase = ((age - np.nanmin(age)) / p) % 1.0
        else:
            phase = pd.Series(np.linspace(0, 1, len(sub)), index=sub.index)

    # Teff and R over phase.
    if teff is not None:
        fig, ax = plt.subplots(figsize=(8, 4.5))
        order = np.argsort(phase.to_numpy(dtype=float))
        ax.scatter(phase.iloc[order], teff.iloc[order], s=8)
        ax.set_xlabel("Pulsation phase")
        ax.set_ylabel(r"$T_{\rm eff}$ (K)")
        ax.set_title("RR Lyrae Teff cycle sanity check")
        ax.grid(alpha=0.25)
        p = outdir / "verify_teff_phase.png"
        fig.tight_layout()
        fig.savefig(p, dpi=180)
        plt.close(fig)
        paths.append(p)

    # Filter light curves if present.
    if filters:
        fig, ax = plt.subplots(figsize=(9, 5))
        order = np.argsort(phase.to_numpy(dtype=float))
        for fcol in filters[:7]:
            y = clean_numeric_series(sub[fcol])
            if y.notna().sum() > 5:
                ax.plot(phase.iloc[order], y.iloc[order], ".", ms=2, label=fcol)
        ax.invert_yaxis()
        ax.set_xlabel("Pulsation phase")
        ax.set_ylabel("Magnitude")
        ax.set_title("Colors light curves in history")
        ax.legend(ncol=2, fontsize=8)
        ax.grid(alpha=0.25)
        p = outdir / "verify_lightcurves_history.png"
        fig.tight_layout()
        fig.savefig(p, dpi=180)
        plt.close(fig)
        paths.append(p)

        # Color-temperature if possible.
        if teff is not None and len(filters) >= 2:
            blue, red = filters[0], filters[-1]
            color = clean_numeric_series(sub[blue]) - clean_numeric_series(sub[red])
            fig, ax = plt.subplots(figsize=(6, 5))
            sc = ax.scatter(color, teff, c=phase, s=8, cmap="viridis")
            ax.set_xlabel(f"{blue} - {red} (mag)")
            ax.set_ylabel(r"$T_{\rm eff}$ (K)")
            ax.set_title("Color-temperature loop")
            ax.grid(alpha=0.25)
            fig.colorbar(sc, ax=ax, label="Pulsation phase")
            p = outdir / "verify_color_teff.png"
            fig.tight_layout()
            fig.savefig(p, dpi=180)
            plt.close(fig)
            paths.append(p)

    return paths


def main() -> int:
    ap = argparse.ArgumentParser(description="Verify that a MESA RSP RR Lyrae + Colors run is physically coherent.")
    ap.add_argument("--history", default=None, help="Path to history.data. Default: newest common LOGS*/history file.")
    ap.add_argument("--sed-dir", default="SED", help="Directory containing *_SED_########.csv files.")
    ap.add_argument("--last-periods", type=int, default=10, help="Number of final periods/cycles to analyze.")
    ap.add_argument("--grekm-threshold", type=float, default=1e-3, help="Settled threshold for abs(rsp_GREKM).")
    ap.add_argument("--outdir", default="verify_rrlyrae_physical", help="Directory for report tables/plots.")
    ap.add_argument("--no-plots", action="store_true", help="Disable PNG diagnostic plots.")
    args = ap.parse_args()

    results: list[CheckResult] = []
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    try:
        history_path = find_history(args.history)
        h = read_mesa_history(history_path)
    except Exception as e:  # noqa: BLE001
        print(f"FAIL: could not read history file: {e}", file=sys.stderr)
        return 2

    add(results, "History file", "INFO", f"Using {history_path} with {len(h)} rows and {len(h.columns)} columns.")
    if len(h) < 50:
        add(results, "History length", "FAIL", "Too few history rows for a pulsation-cycle sanity check.")
    else:
        add(results, "History length", "PASS", f"History has {len(h)} rows.")

    required = ["model_number", "log_L", "rsp_num_periods", "rsp_period_in_days", "rsp_GREKM"]
    missing = [c for c in required if c not in h.columns]
    if missing:
        add(results, "Required RSP columns", "FAIL", "Missing: " + ", ".join(missing))
    else:
        add(results, "Required RSP columns", "PASS", "Found model, luminosity, period, cycle, and GREKM diagnostics.")

    teff = get_teff(h)
    if teff is None or teff.notna().sum() < 10:
        add(results, "Teff column", "FAIL", "Could not find Teff/effective_T/log_Teff.")
    else:
        vals = teff.to_numpy(dtype=float)
        add(results, "Teff column", "PASS", f"Teff range {np.nanmin(vals):.0f}--{np.nanmax(vals):.0f} K.")

    sub = select_last_periods(h, args.last_periods)
    add(results, "Analysis window", "INFO", f"Analyzing {len(sub)} rows from the last {args.last_periods} period(s), where available.")

    # Period and RR Lyrae parameter sanity.
    if "rsp_period_in_days" in h.columns:
        p = clean_numeric_series(sub["rsp_period_in_days"])
        pvals = p.to_numpy(dtype=float)
        pvals = pvals[np.isfinite(pvals) & (pvals > 0)]
        if len(pvals):
            pmed = float(np.nanmedian(pvals))
            prel = float(np.nanstd(pvals) / pmed) if pmed > 0 else np.nan
            rr_ok = 0.2 <= pmed <= 1.2
            add(results, "RR Lyrae period range", "PASS" if rr_ok else "WARN", f"Median period {pmed:.5f} d; relative scatter in window {prel:.3e}.")
            add(results, "Period stability", "PASS" if prel < 2e-3 else "WARN", f"Relative period scatter {prel:.3e} in analysis window.")

    if "rsp_num_periods" in h.columns:
        nums = clean_numeric_series(h["rsp_num_periods"]).to_numpy(dtype=float)
        nums = nums[np.isfinite(nums)]
        if len(nums):
            n_cycles = int(np.nanmax(nums) - np.nanmin(nums) + 1)
            add(results, "Completed periods", "PASS" if np.nanmax(nums) >= args.last_periods else "WARN", f"rsp_num_periods range {np.nanmin(nums):.0f}--{np.nanmax(nums):.0f} ({n_cycles} integer cycles spanned).")

    if "rsp_GREKM" in h.columns:
        g = np.abs(clean_numeric_series(h["rsp_GREKM"]).to_numpy(dtype=float))
        g = g[np.isfinite(g)]
        if len(g):
            last = float(g[-1])
            loglast = math.log10(last) if last > 0 else -np.inf
            add(results, "RSP growth/stability", "PASS" if last <= args.grekm_threshold else "WARN", f"Final abs(rsp_GREKM)={last:.4e}; log10={loglast:.3f}; target <= {args.grekm_threshold:.1e}.")

    if teff is not None:
        tsub = teff.loc[sub.index]
        trange = float(np.nanmax(tsub) - np.nanmin(tsub)) if tsub.notna().sum() else np.nan
        add(results, "Teff pulsation amplitude", "PASS" if trange > 300 else "WARN", f"Teff peak-to-peak range in window: {trange:.0f} K.")
        tmed = float(np.nanmedian(tsub)) if tsub.notna().sum() else np.nan
        add(results, "RR Lyrae Teff regime", "PASS" if 5500 <= tmed <= 8500 else "WARN", f"Median Teff in window: {tmed:.0f} K.")

    if "log_L" in h.columns:
        log_l = clean_numeric_series(h["log_L"])
        lsub = log_l.loc[sub.index]
        lmed = float(np.nanmedian(lsub)) if lsub.notna().sum() else np.nan
        add(results, "RR Lyrae luminosity regime", "PASS" if 1.2 <= lmed <= 2.2 else "WARN", f"Median log_L/Lsun in window: {lmed:.3f}.")
    else:
        log_l = pd.Series(np.nan, index=h.index)

    if teff is not None and "log_L" in h.columns:
        r, r_name, sb_mad = infer_radius_rsun(h, teff, log_l)
        if r is not None and sb_mad is not None:
            rsub = r.loc[sub.index]
            rmed = float(np.nanmedian(rsub))
            ramp = float(np.nanmax(rsub) - np.nanmin(rsub))
            add(results, "Radius inference", "INFO", f"Using {r_name}; median R={rmed:.3f} Rsun; peak-to-peak range={ramp:.3f} Rsun.")
            add(results, "Stefan-Boltzmann consistency", "PASS" if sb_mad < 0.05 else "WARN", f"Median |log_L - (2logR+4logTeff)| = {sb_mad:.4f} dex.")
            add(results, "RR Lyrae radius regime", "PASS" if 3.0 <= rmed <= 9.0 else "WARN", f"Median radius {rmed:.3f} Rsun.")
        else:
            add(results, "Radius/Stellar consistency", "WARN", "Could not infer a radius column for Stefan-Boltzmann check.")

    # Dynamics: check shocks/supersonic velocities if available.
    for col in ["max_abs_v_div_cs", "v_div_cs"]:
        if col in h.columns:
            v = np.abs(clean_numeric_series(sub[col]).to_numpy(dtype=float))
            v = v[np.isfinite(v)]
            if len(v):
                vmax = float(np.nanmax(v))
                status = "PASS" if vmax < 1.0 else "WARN"
                add(results, "Mach/shock diagnostic", status, f"max |v/cs| in window from {col}: {vmax:.3f}. Supersonic RSP phases may be shocky but not automatically invalid.")
                break

    # Color checks from history, if present.
    filters = detect_filter_columns(h)
    if filters:
        add(results, "Filter magnitudes in history", "INFO", "Found: " + ", ".join(filters))
        amps = []
        for fcol in filters:
            y = clean_numeric_series(sub[fcol]).to_numpy(dtype=float)
            y = y[np.isfinite(y)]
            if len(y):
                amps.append((fcol, FILTER_LAMBDA_A[fcol], float(np.nanmax(y) - np.nanmin(y))))
        if amps:
            amp_msg = ", ".join(f"{name}:{amp:.3f}" for name, _, amp in amps)
            add(results, "Filter amplitudes", "INFO", amp_msg)
            if len(amps) >= 2:
                lam = np.array([a[1] for a in amps])
                ampv = np.array([a[2] for a in amps])
                rho = pearson(lam, ampv)
                add(results, "Amplitude-wavelength trend", "PASS" if rho < -0.4 else "WARN", f"Pearson r(amplitude, wavelength)={rho:.3f}; shorter wavelengths should generally have larger amplitudes.")
        if teff is not None and len(filters) >= 2:
            blue, red = filters[0], filters[-1]
            color = clean_numeric_series(sub[blue]) - clean_numeric_series(sub[red])
            rho = pearson(color, teff.loc[sub.index])
            add(results, "Color-temperature relation", "PASS" if rho < -0.8 else "WARN", f"Pearson r({blue}-{red}, Teff)={rho:.3f}; hotter phases should be bluer.")
    else:
        add(results, "Filter magnitudes in history", "INFO", "No recognized filter magnitude columns found in history; this is okay if relying on per-model SED files.")

    sed_info = analyze_seds(Path(args.sed_dir), results)

    # Save a machine-readable summary.
    report_df = pd.DataFrame([r.__dict__ for r in results])
    report_path = outdir / "verification_report.csv"
    report_df.to_csv(report_path, index=False)

    sed_table = sed_info.get("sample_table")
    if isinstance(sed_table, pd.DataFrame):
        sed_table.to_csv(outdir / "sed_sample_summary.csv", index=False)

    plot_paths: list[Path] = []
    if not args.no_plots:
        try:
            plot_paths = make_plots(h, sub, outdir, filters, sed_info)
        except Exception as e:  # noqa: BLE001
            add(results, "Plot generation", "WARN", f"Could not generate diagnostic plots: {e}")

    # Print report.
    print("\nRR LYRAE PHYSICAL VERIFICATION REPORT")
    print("=" * 72)
    print(f"workdir: {Path.cwd()}")
    print(f"history: {history_path}")
    print(f"analysis rows: {len(sub)}")
    print("-" * 72)
    for r in results:
        print(f"{status_symbol(r.status)} {r.status:4s} | {r.name}: {r.message}")
    print("-" * 72)
    print(f"wrote: {report_path}")
    if isinstance(sed_table, pd.DataFrame):
        print(f"wrote: {outdir / 'sed_sample_summary.csv'}")
    for p in plot_paths:
        print(f"wrote: {p}")

    n_fail = sum(r.status == "FAIL" for r in results)
    n_warn = sum(r.status == "WARN" for r in results)
    print("-" * 72)
    if n_fail:
        print(f"VERDICT: FAIL ({n_fail} failed checks, {n_warn} warnings). This run is not yet safe to present.")
        return 2
    if n_warn:
        print(f"VERDICT: PLAUSIBLE WITH CAVEATS ({n_warn} warnings). Read warnings before using in the paper.")
        return 1
    print("VERDICT: PASS. This looks physically coherent for a MESA RSP+Colors demonstration.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
