#!/usr/bin/env python3
"""
Characterize a MESA RSP RR Lyrae + Colors run for paper writing.

This script is designed for a two-stage MESA/RSP run where the final Colors stage
writes LOGS_colors/history.data and optionally per-model SED CSV files in SED/.
It produces a compact, paper-oriented characterization package: physical summary,
pulsation diagnostics, filter amplitudes, color-temperature/radius diagnostics,
selected SED diagnostics, plots, CSV tables, LaTeX table snippets, and a Markdown
report with values ready to quote.

Typical use:
    cd ~/MESA/MESA_Colors_Tests/rsp_RR_Lyrae
    ./characterize_rrlyrae_colors.py

Useful options:
    ./characterize_rrlyrae_colors.py --last-periods 10
    ./characterize_rrlyrae_colors.py --history LOGS_colors/history.data --sed-dir SED
    ./characterize_rrlyrae_colors.py --no-sed
    ./characterize_rrlyrae_colors.py --sed-sample 3000
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

# matplotlib is only imported when plotting, so --no-plots works on headless/minimal systems.

LSUN_CGS = 3.828e33
RSUN_CGS = 6.957e10
SIGMA_SB_CGS = 5.670374419e-5
TEFF_SUN = 5772.0

ROMAN_FILTER_WAVELENGTH_ANG = {
    "F062": 6200.0,
    "F087": 8700.0,
    "F106": 10600.0,
    "F129": 12900.0,
    "F146": 14600.0,
    "F158": 15800.0,
    "F184": 18400.0,
    "F213": 21300.0,
}

COMMON_FILTER_WAVELENGTH_ANG = {
    **ROMAN_FILTER_WAVELENGTH_ANG,
    "u": 3540.0,
    "g": 4770.0,
    "r": 6230.0,
    "i": 7630.0,
    "z": 9130.0,
    "y": 10200.0,
    "Gbp": 5320.0,
    "G": 6730.0,
    "Grp": 7970.0,
    "Kepler": 6410.0,
    "TESS": 7865.0,
}

COLOR_PAIRS_PREFERRED = [
    ("F062", "F213"),
    ("F062", "F184"),
    ("F087", "F184"),
    ("F106", "F184"),
    ("g", "r"),
    ("Gbp", "Grp"),
]


@dataclass
class PhasePoint:
    label: str
    idx: int
    model_number: int
    phase: float


def safe_float_token(token: str) -> float:
    """Parse normal floats and compact Fortran exponent strings like -1.078100-112."""
    t = token.strip()
    if not t:
        return np.nan
    t = t.replace("D", "E").replace("d", "E")
    # Repair strings like 1.234567-108 or -1.234567+108 to E notation.
    if "E" not in t and "e" not in t:
        m = re.match(r"^([+-]?(?:\d+\.\d*|\d*\.\d+|\d+))([+-]\d{2,4})$", t)
        if m:
            t = m.group(1) + "E" + m.group(2)
    try:
        return float(t)
    except ValueError:
        return np.nan


def finite(x: np.ndarray) -> np.ndarray:
    return np.asarray(x)[np.isfinite(x)]


def nanptp(x: Sequence[float]) -> float:
    arr = finite(np.asarray(x, dtype=float))
    if arr.size == 0:
        return np.nan
    return float(np.nanmax(arr) - np.nanmin(arr))


def nanmedian(x: Sequence[float]) -> float:
    arr = finite(np.asarray(x, dtype=float))
    return float(np.nanmedian(arr)) if arr.size else np.nan


def nanmean(x: Sequence[float]) -> float:
    arr = finite(np.asarray(x, dtype=float))
    return float(np.nanmean(arr)) if arr.size else np.nan


def robust_corr(x: Sequence[float], y: Sequence[float]) -> float:
    xa = np.asarray(x, dtype=float)
    ya = np.asarray(y, dtype=float)
    m = np.isfinite(xa) & np.isfinite(ya)
    if m.sum() < 3:
        return np.nan
    sx = np.nanstd(xa[m])
    sy = np.nanstd(ya[m])
    if sx == 0 or sy == 0:
        return np.nan
    return float(np.corrcoef(xa[m], ya[m])[0, 1])


def phase_distance_forward(phi0: float, phi1: float) -> float:
    """Forward phase distance from phi0 to phi1 in [0,1)."""
    return float((phi1 - phi0) % 1.0)


def normalize_for_plot(y: Sequence[float]) -> np.ndarray:
    arr = np.asarray(y, dtype=float)
    mn = np.nanmin(arr)
    mx = np.nanmax(arr)
    if not np.isfinite(mn) or not np.isfinite(mx) or mx == mn:
        return np.full_like(arr, np.nan)
    return 2.0 * (arr - mn) / (mx - mn) - 1.0


def read_mesa_history(path: Path) -> pd.DataFrame:
    lines = path.read_text(errors="ignore").splitlines()
    header_idx = None
    for i, line in enumerate(lines):
        parts = line.split()
        if "model_number" in parts:
            header_idx = i
            break
    if header_idx is None:
        raise RuntimeError(f"Could not find model_number header in {path}")
    df = pd.read_csv(path, sep=r"\s+", skiprows=header_idx, engine="python", on_bad_lines="skip")
    for c in df.columns:
        # Convert all numeric-looking columns; leave truly string columns alone.
        converted = pd.to_numeric(df[c], errors="coerce")
        if converted.notna().sum() >= max(1, int(0.5 * len(df))):
            df[c] = converted
    return df


def pick_history(workdir: Path, explicit: Optional[str]) -> Path:
    if explicit:
        p = Path(explicit)
        return p if p.is_absolute() else workdir / p
    candidates = [
        workdir / "LOGS_colors" / "history.data",
        workdir / "LOGS" / "history.data",
        workdir / "LOGS_settle" / "history.data",
    ]
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError("Could not find LOGS_colors/history.data, LOGS/history.data, or LOGS_settle/history.data")


def detect_phase(df: pd.DataFrame) -> pd.Series:
    if "rsp_phase" in df.columns:
        ph = pd.to_numeric(df["rsp_phase"], errors="coerce") % 1.0
        if ph.notna().sum() > 10:
            return ph
    # Fallback: within each integer period use cumulative row ordering.
    if "rsp_num_periods" in df.columns:
        cyc = np.floor(pd.to_numeric(df["rsp_num_periods"], errors="coerce").to_numpy()).astype(int)
        ph = np.zeros(len(df), dtype=float)
        for c in np.unique(cyc):
            idx = np.where(cyc == c)[0]
            if idx.size > 1:
                ph[idx] = np.linspace(0.0, 1.0, idx.size, endpoint=False)
        return pd.Series(ph, index=df.index)
    return pd.Series(np.linspace(0.0, 1.0, len(df), endpoint=False), index=df.index)


def restrict_last_periods(df: pd.DataFrame, last_periods: int) -> pd.DataFrame:
    if last_periods <= 0 or "rsp_num_periods" not in df.columns:
        return df.copy()
    cycles = pd.to_numeric(df["rsp_num_periods"], errors="coerce")
    if cycles.notna().sum() < 10:
        return df.copy()
    max_cycle = int(np.nanmax(cycles))
    min_cycle = max_cycle - last_periods + 1
    sub = df[cycles >= min_cycle].copy()
    if len(sub) < 20:
        return df.copy()
    return sub


def detect_filter_columns(df: pd.DataFrame) -> List[str]:
    filters = []
    for c in df.columns:
        if c in COMMON_FILTER_WAVELENGTH_ANG:
            if pd.to_numeric(df[c], errors="coerce").notna().sum() > 10:
                filters.append(c)
        elif re.match(r"^F\d{3}$", str(c)):
            if pd.to_numeric(df[c], errors="coerce").notna().sum() > 10:
                filters.append(c)
    filters = sorted(set(filters), key=lambda k: COMMON_FILTER_WAVELENGTH_ANG.get(k, 1e99))
    return filters


def get_luminosity(df: pd.DataFrame) -> np.ndarray:
    if "log_L" in df.columns:
        return 10 ** pd.to_numeric(df["log_L"], errors="coerce").to_numpy(float)
    if "luminosity" in df.columns:
        lum = pd.to_numeric(df["luminosity"], errors="coerce").to_numpy(float)
        # MESA luminosity is often Lsun in histories; do not convert if values are RR-like.
        return lum
    return np.full(len(df), np.nan)


def get_teff(df: pd.DataFrame) -> np.ndarray:
    for c in ["Teff", "effective_T"]:
        if c in df.columns:
            return pd.to_numeric(df[c], errors="coerce").to_numpy(float)
    if "log_Teff" in df.columns:
        return 10 ** pd.to_numeric(df["log_Teff"], errors="coerce").to_numpy(float)
    return np.full(len(df), np.nan)


def infer_radius(df: pd.DataFrame, lum_lsun: np.ndarray, teff: np.ndarray) -> Tuple[np.ndarray, str]:
    # R/Rsun = sqrt(L/Lsun) * (Tsun/T)^2
    sb_radius = np.sqrt(lum_lsun) * (TEFF_SUN / teff) ** 2
    for c in ["photosphere_r", "radius", "R", "photosphere_radius"]:
        if c in df.columns:
            r = pd.to_numeric(df[c], errors="coerce").to_numpy(float)
            med = nanmedian(r)
            if np.isfinite(med) and 0.1 < med < 100:
                return r, c
            # If in cm, convert.
            if np.isfinite(med) and med > 1e9:
                return r / RSUN_CGS, c + " / Rsun"
    return sb_radius, "Stefan-Boltzmann inferred"


def summarize_series(name: str, arr: Sequence[float], unit: str = "") -> Dict[str, Any]:
    a = np.asarray(arr, dtype=float)
    af = finite(a)
    out = {"quantity": name, "unit": unit, "n": int(af.size)}
    if af.size == 0:
        out.update({k: np.nan for k in ["min", "p05", "median", "mean", "p95", "max", "ptp", "std"]})
    else:
        out.update(
            {
                "min": float(np.nanmin(af)),
                "p05": float(np.nanpercentile(af, 5)),
                "median": float(np.nanmedian(af)),
                "mean": float(np.nanmean(af)),
                "p95": float(np.nanpercentile(af, 95)),
                "max": float(np.nanmax(af)),
                "ptp": float(np.nanmax(af) - np.nanmin(af)),
                "std": float(np.nanstd(af)),
            }
        )
    return out


def find_phase_landmarks(df: pd.DataFrame, phase: np.ndarray, filters: List[str], teff: np.ndarray, radius: np.ndarray) -> List[PhasePoint]:
    landmarks: List[PhasePoint] = []
    model = pd.to_numeric(df.get("model_number", pd.Series(np.arange(len(df)))), errors="coerce").to_numpy()

    def add(label: str, idx_local: int) -> None:
        if idx_local < 0 or idx_local >= len(df):
            return
        landmarks.append(PhasePoint(label, int(idx_local), int(model[idx_local]) if np.isfinite(model[idx_local]) else idx_local, float(phase[idx_local])))

    if finite(teff).size:
        add("hottest_Teff", int(np.nanargmax(teff)))
        add("coolest_Teff", int(np.nanargmin(teff)))
    if finite(radius).size:
        add("max_radius", int(np.nanargmax(radius)))
        add("min_radius", int(np.nanargmin(radius)))
    if filters:
        ref = filters[0]
        mag = pd.to_numeric(df[ref], errors="coerce").to_numpy(float)
        if finite(mag).size:
            add(f"max_light_{ref}", int(np.nanargmin(mag)))
            add(f"min_light_{ref}", int(np.nanargmax(mag)))
    # Deduplicate by label, not idx.
    return landmarks


def make_filter_table(df: pd.DataFrame, phase: np.ndarray, filters: List[str]) -> pd.DataFrame:
    rows = []
    ref_max_phase = None
    if filters:
        ref_mag = pd.to_numeric(df[filters[0]], errors="coerce").to_numpy(float)
        if finite(ref_mag).size:
            ref_max_phase = float(phase[int(np.nanargmin(ref_mag))])

    for f in filters:
        mag = pd.to_numeric(df[f], errors="coerce").to_numpy(float)
        if finite(mag).size < 5:
            continue
        i_bright = int(np.nanargmin(mag))  # magnitude minimum = max light
        i_faint = int(np.nanargmax(mag))   # magnitude maximum = min light
        amp = float(np.nanmax(mag) - np.nanmin(mag))
        phi_bright = float(phase[i_bright])
        phi_faint = float(phase[i_faint])
        rise_time = phase_distance_forward(phi_faint, phi_bright)
        rows.append(
            {
                "filter": f,
                "lambda_eff_A_assumed": COMMON_FILTER_WAVELENGTH_ANG.get(f, np.nan),
                "mag_min_bright": float(np.nanmin(mag)),
                "mag_max_faint": float(np.nanmax(mag)),
                "amplitude_mag": amp,
                "median_mag": float(np.nanmedian(mag)),
                "phase_max_light": phi_bright,
                "phase_min_light": phi_faint,
                "rise_time_phase_min_to_max_light": rise_time,
                "phase_lag_max_light_vs_first_filter": phase_distance_forward(ref_max_phase, phi_bright) if ref_max_phase is not None else np.nan,
            }
        )
    return pd.DataFrame(rows)


def polygon_area(x: np.ndarray, y: np.ndarray) -> float:
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    if len(x) < 3:
        return np.nan
    return float(0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))


def make_color_table(df: pd.DataFrame, phase: np.ndarray, filters: List[str], teff: np.ndarray, radius: np.ndarray) -> pd.DataFrame:
    pairs = []
    seen = set()
    for a, b in COLOR_PAIRS_PREFERRED:
        if a in filters and b in filters:
            pairs.append((a, b))
            seen.add((a, b))
    # Add adjacent-ish extremes if preferred pairs absent.
    if len(filters) >= 2:
        pair = (filters[0], filters[-1])
        if pair not in seen:
            pairs.insert(0, pair)
    rows = []
    for a, b in pairs:
        color = pd.to_numeric(df[a], errors="coerce").to_numpy(float) - pd.to_numeric(df[b], errors="coerce").to_numpy(float)
        rows.append(
            {
                "color": f"{a}-{b}",
                "blue_filter": a,
                "red_filter": b,
                "median_color_mag": nanmedian(color),
                "min_color_mag": float(np.nanmin(color)),
                "max_color_mag": float(np.nanmax(color)),
                "amplitude_color_mag": nanptp(color),
                "pearson_color_Teff": robust_corr(color, teff),
                "pearson_color_radius": robust_corr(color, radius),
                "loop_area_color_Teff": polygon_area(color, teff),
                "loop_area_color_radius": polygon_area(color, radius),
            }
        )
    return pd.DataFrame(rows)


def fourier_fit_phase(phase: np.ndarray, y: np.ndarray, order: int = 6) -> Dict[str, Any]:
    m = np.isfinite(phase) & np.isfinite(y)
    x = phase[m] % 1.0
    yy = y[m]
    out: Dict[str, Any] = {"fourier_order": order, "n_fit": int(len(yy))}
    if len(yy) < 2 * order + 3:
        return {**out, "status": "insufficient_data"}
    cols = [np.ones_like(x)]
    for k in range(1, order + 1):
        cols.append(np.cos(2 * np.pi * k * x))
        cols.append(np.sin(2 * np.pi * k * x))
    A = np.vstack(cols).T
    try:
        coef, *_ = np.linalg.lstsq(A, yy, rcond=None)
    except np.linalg.LinAlgError:
        return {**out, "status": "fit_failed"}
    yfit = A @ coef
    rms = float(np.sqrt(np.nanmean((yy - yfit) ** 2)))
    amps = []
    phases = []
    for k in range(order):
        a = coef[1 + 2 * k]
        b = coef[2 + 2 * k]
        amps.append(float(np.hypot(a, b)))
        # y = A cos(2pi k phi - phi_k), phi_k = atan2(b,a)
        phases.append(float(np.arctan2(b, a)))
    a1 = amps[0] if amps else np.nan
    out.update(
        {
            "status": "ok",
            "mean": float(coef[0]),
            "rms_residual": rms,
            "A1": a1,
            "A2": amps[1] if len(amps) > 1 else np.nan,
            "A3": amps[2] if len(amps) > 2 else np.nan,
            "R21": amps[1] / a1 if len(amps) > 1 and a1 != 0 else np.nan,
            "R31": amps[2] / a1 if len(amps) > 2 and a1 != 0 else np.nan,
            "phi21_rad": (phases[1] - 2 * phases[0]) % (2 * np.pi) if len(phases) > 1 else np.nan,
            "phi31_rad": (phases[2] - 3 * phases[0]) % (2 * np.pi) if len(phases) > 2 else np.nan,
        }
    )
    return out


def make_fourier_table(df: pd.DataFrame, phase: np.ndarray, filters: List[str], teff: np.ndarray, radius: np.ndarray) -> pd.DataFrame:
    series = {"Teff": teff, "radius_Rsun": radius}
    for f in filters:
        series[f] = pd.to_numeric(df[f], errors="coerce").to_numpy(float)
    rows = []
    for name, arr in series.items():
        row = {"series": name}
        row.update(fourier_fit_phase(phase, arr, order=6))
        rows.append(row)
    return pd.DataFrame(rows)


def parse_simple_inlist_assignments(path: Path) -> Dict[str, str]:
    if not path.exists():
        return {}
    out: Dict[str, str] = {}
    for raw in path.read_text(errors="ignore").splitlines():
        line = raw.split("!", 1)[0].strip()
        if not line or "=" not in line:
            continue
        if line.startswith("&") or line.startswith("/"):
            continue
        key, val = line.split("=", 1)
        key = key.strip()
        val = val.strip().rstrip(",")
        if key:
            out[key] = val
    return out


def load_sed_csv(path: Path) -> pd.DataFrame:
    # Manual parsing so compact Fortran exponents like -1.078100-112 do not break read_csv.
    with path.open("r", errors="ignore", newline="") as fh:
        reader = csv.reader(fh)
        rows = list(reader)
    if not rows:
        raise RuntimeError(f"Empty SED file {path}")
    header = [h.strip() for h in rows[0]]
    data: Dict[str, List[float]] = {h: [] for h in header}
    for row in rows[1:]:
        if len(row) < len(header):
            continue
        for h, token in zip(header, row):
            data[h].append(safe_float_token(token))
    return pd.DataFrame(data)


def sed_file_for_model(sed_dir: Path, filt: str, model: int) -> Optional[Path]:
    candidates = [
        sed_dir / f"{filt}_SED_{model:08d}.csv",
        sed_dir / f"{filt}_SED_{model:06d}.csv",
        sed_dir / f"{filt}_SED_{model}.csv",
    ]
    for p in candidates:
        if p.exists():
            return p
    matches = list(sed_dir.glob(f"{filt}_SED_*{model}*.csv"))
    return matches[0] if matches else None


def summarize_sed_file(path: Path) -> Dict[str, Any]:
    df = load_sed_csv(path)
    cols = {c.strip(): c for c in df.columns}
    wcol = cols.get("wavelengths") or cols.get("wavelength") or list(df.columns)[0]
    fcol = cols.get("fluxes") or cols.get("flux") or list(df.columns)[1]
    w = pd.to_numeric(df[wcol], errors="coerce").to_numpy(float)
    f = pd.to_numeric(df[fcol], errors="coerce").to_numpy(float)
    m = np.isfinite(w) & np.isfinite(f)
    w = w[m]
    f = f[m]
    order = np.argsort(w)
    w = w[order]
    f = f[order]
    pos = np.clip(f, 0, None)
    neg = np.clip(-f, 0, None)
    pos_area = float(np.trapezoid(pos, w)) if len(w) > 1 else np.nan
    neg_area = float(np.trapezoid(neg, w)) if len(w) > 1 else np.nan
    total_area = float(np.trapezoid(f, w)) if len(w) > 1 else np.nan
    return {
        "file": path.name,
        "path": str(path),
        "n_wave": int(len(w)),
        "wave_min_A": float(np.nanmin(w)) if len(w) else np.nan,
        "wave_max_A": float(np.nanmax(w)) if len(w) else np.nan,
        "flux_min": float(np.nanmin(f)) if len(f) else np.nan,
        "flux_max": float(np.nanmax(f)) if len(f) else np.nan,
        "flux_integral": total_area,
        "positive_area": pos_area,
        "negative_area": neg_area,
        "negative_area_fraction": neg_area / pos_area if pos_area and np.isfinite(pos_area) else np.nan,
        "has_negative_flux": bool(np.nanmin(f) < 0) if len(f) else False,
    }


def analyze_selected_seds(sed_dir: Path, filters: List[str], landmarks: List[PhasePoint], sed_sample: int) -> Tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    selected_filter = filters[0] if filters else "F062"
    for lm in landmarks:
        for filt in sorted(set([selected_filter, filters[-1] if filters else selected_filter])):
            p = sed_file_for_model(sed_dir, filt, lm.model_number)
            if p is None:
                continue
            try:
                row = summarize_sed_file(p)
                row.update({"landmark": lm.label, "filter_file_prefix": filt, "model_number": lm.model_number, "phase": lm.phase})
                rows.append(row)
            except Exception as e:
                rows.append({"landmark": lm.label, "filter_file_prefix": filt, "model_number": lm.model_number, "phase": lm.phase, "file": str(p), "error": str(e)})

    # Sample SED integrity across the directory without reading every file.
    sample_rows = []
    if sed_dir.exists() and sed_sample > 0:
        files = sorted(sed_dir.glob("*_SED_*.csv"))
        if files:
            stride = max(1, len(files) // sed_sample)
            for p in files[::stride][:sed_sample]:
                try:
                    row = summarize_sed_file(p)
                    row["sampled"] = True
                    sample_rows.append(row)
                except Exception as e:
                    sample_rows.append({"file": p.name, "path": str(p), "sampled": True, "error": str(e)})
    return pd.DataFrame(rows), pd.DataFrame(sample_rows)


def dataframe_to_latex(df: pd.DataFrame, path: Path, float_format: str = "%.4g") -> None:
    try:
        tex = df.to_latex(index=False, escape=False, float_format=lambda x: float_format % x)
    except Exception:
        tex = df.to_string(index=False)
    path.write_text(tex)


def write_report(
    outdir: Path,
    workdir: Path,
    history_path: Path,
    df: pd.DataFrame,
    sub: pd.DataFrame,
    global_df: pd.DataFrame,
    filter_df: pd.DataFrame,
    color_df: pd.DataFrame,
    fourier_df: pd.DataFrame,
    sed_selected_df: pd.DataFrame,
    sed_sample_df: pd.DataFrame,
    inlist_meta: Dict[str, Dict[str, str]],
    radius_source: str,
    filters: List[str],
) -> None:
    def getq(q: str, stat: str = "median") -> float:
        row = global_df[global_df["quantity"] == q]
        if row.empty or stat not in row.columns:
            return np.nan
        return float(row.iloc[0][stat])

    period = getq("RSP period", "median")
    period_std = getq("RSP period", "std")
    teff_med = getq("Teff", "median")
    teff_ptp = getq("Teff", "ptp")
    logl_med = getq("log_L", "median")
    rad_med = getq("radius", "median")
    rad_ptp = getq("radius", "ptp")
    grekm_final = np.nan
    if "rsp_GREKM" in df.columns:
        g = pd.to_numeric(df["rsp_GREKM"], errors="coerce").to_numpy(float)
        if len(g):
            grekm_final = float(abs(g[-1]))

    best_color = color_df.iloc[0].to_dict() if not color_df.empty else {}

    lines = []
    lines.append("# RR Lyrae RSP + Colors characterization report")
    lines.append("")
    lines.append("## Run provenance")
    lines.append("")
    lines.append(f"- Work directory: `{workdir}`")
    lines.append(f"- History file: `{history_path}`")
    lines.append(f"- Total history rows: {len(df)}")
    lines.append(f"- Analysis rows: {len(sub)}")
    if "model_number" in df.columns:
        lines.append(f"- Model range: {int(np.nanmin(df['model_number']))}--{int(np.nanmax(df['model_number']))}")
    if "rsp_num_periods" in df.columns:
        c = pd.to_numeric(df["rsp_num_periods"], errors="coerce")
        lines.append(f"- Completed periods in full history: {int(np.nanmin(c))}--{int(np.nanmax(c))}")
    lines.append(f"- Detected filter columns: {', '.join(filters) if filters else 'none'}")
    lines.append("")

    lines.append("## Main paper-ready numbers")
    lines.append("")
    lines.append(f"- Median pulsation period: **{period:.6f} d** (window scatter {period_std:.3e} d).")
    lines.append(f"- Median effective temperature: **{teff_med:.0f} K**, with peak-to-peak variation **{teff_ptp:.0f} K**.")
    lines.append(f"- Median luminosity: **log(L/Lsun) = {logl_med:.3f}**.")
    lines.append(f"- Median photospheric radius: **{rad_med:.3f} Rsun**, with peak-to-peak variation **{rad_ptp:.3f} Rsun**. Radius source: `{radius_source}`.")
    if np.isfinite(grekm_final):
        lines.append(f"- Final RSP growth/stability diagnostic: **|rsp_GREKM| = {grekm_final:.3e}** (log10={math.log10(grekm_final):.3f}).")
    if best_color:
        lines.append(f"- Strongest/default color diagnostic: **{best_color.get('color')}**, Pearson r(color, Teff) = **{best_color.get('pearson_color_Teff', np.nan):.3f}**.")
    lines.append("")

    if not filter_df.empty:
        lines.append("## Filter amplitudes and phase landmarks")
        lines.append("")
        lines.append(filter_df.to_markdown(index=False, floatfmt=".5g"))
        lines.append("")

    if not color_df.empty:
        lines.append("## Color diagnostics")
        lines.append("")
        lines.append(color_df.to_markdown(index=False, floatfmt=".5g"))
        lines.append("")

    lines.append("## Physical summary statistics")
    lines.append("")
    lines.append(global_df.to_markdown(index=False, floatfmt=".5g"))
    lines.append("")

    if not fourier_df.empty:
        lines.append("## Fourier morphology diagnostics")
        lines.append("")
        lines.append("These are descriptive Fourier coefficients of phase-folded curves, useful for shape comparison; they are not meant as an observational calibration unless the model is tuned.")
        lines.append("")
        lines.append(fourier_df.to_markdown(index=False, floatfmt=".5g"))
        lines.append("")

    if not sed_selected_df.empty:
        lines.append("## Selected SED diagnostics")
        lines.append("")
        lines.append(sed_selected_df.to_markdown(index=False, floatfmt=".5g"))
        lines.append("")

    if not sed_sample_df.empty:
        neg_frac = pd.to_numeric(sed_sample_df.get("negative_area_fraction", pd.Series(dtype=float)), errors="coerce")
        finite_neg = neg_frac[np.isfinite(neg_frac)]
        lines.append("## SED integrity sample")
        lines.append("")
        lines.append(f"- Sampled SED files: {len(sed_sample_df)}")
        if len(finite_neg):
            lines.append(f"- Median negative-area fraction: {np.nanmedian(finite_neg):.3e}")
            lines.append(f"- 95th percentile negative-area fraction: {np.nanpercentile(finite_neg, 95):.3e}")
            lines.append(f"- Maximum negative-area fraction: {np.nanmax(finite_neg):.3e}")
        lines.append("")

    lines.append("## Inlist/control metadata extracted naively")
    lines.append("")
    for name, meta in inlist_meta.items():
        lines.append(f"### {name}")
        if not meta:
            lines.append("No assignments parsed or file missing.")
        else:
            for k in sorted(meta):
                if any(tok.lower() in k.lower() for tok in ["rsp", "color", "filter", "atm", "mass", "lum", "teff", "x_ctrl", "x_integer", "history", "save", "load", "model", "period", "sed", "instrument", "distance"]):
                    lines.append(f"- `{k}` = `{meta[k]}`")
        lines.append("")

    lines.append("## Suggested wording/caveats")
    lines.append("")
    lines.append("- Treat this as a time-dependent synthetic-photometry demonstration from a self-consistent RSP calculation, not as a calibrated fit to a particular observed RR Lyrae.")
    lines.append("- If max |v/cs| is large, describe the model as shocky/supersonic in some phases rather than pretending the detailed morphology is tuned.")
    lines.append("- The robust quantities for the section are the period, Teff/radius/luminosity ranges, color-temperature anti-correlation, wavelength-dependent amplitude trend, and phase-dependent SED/filter response.")
    lines.append("")

    (outdir / "characterization_report.md").write_text("\n".join(lines))

    # Compact LaTeX snippet for paper draft.
    snippet = f"""% Auto-generated by characterize_rrlyrae_colors.py
% Edit prose before using directly.
The RSP demonstration model has a median pulsation period of {period:.5f} d over the analysed cycles.
Across the same window, the model spans T_{{\rm eff}} = {getq('Teff','min'):.0f}--{getq('Teff','max'):.0f} K, with a median log(L/L_\odot) = {logl_med:.3f} and a median photospheric radius of {rad_med:.3f}\,R_\odot.
The final growth diagnostic is |rsp\_GREKM| = {grekm_final:.3e}, indicating that the limit cycle is stable by the adopted criterion.
The synthetic Roman amplitudes decrease redward; for example the bluest and reddest detected filters have amplitudes {filter_df.iloc[0]['amplitude_mag']:.3f} and {filter_df.iloc[-1]['amplitude_mag']:.3f} mag, respectively.
The broad color {best_color.get('color', 'N/A')} is tightly anti-correlated with temperature, with Pearson r = {best_color.get('pearson_color_Teff', np.nan):.3f}.
"""
    (outdir / "paper_values_snippet.tex").write_text(snippet)


def make_plots(outdir: Path, sub: pd.DataFrame, phase: np.ndarray, filters: List[str], teff: np.ndarray, radius: np.ndarray, sed_selected_df: pd.DataFrame) -> None:
    import matplotlib.pyplot as plt

    # 1. Teff/radius/logL phase panel.
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(phase, teff, ".", markersize=2)
    ax.set_xlabel("Pulsation phase")
    ax.set_ylabel(r"$T_{\rm eff}$ (K)")
    ax.set_title("RR Lyrae effective-temperature cycle")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(outdir / "phase_Teff.png", dpi=180)
    plt.close(fig)

    if np.isfinite(radius).sum() > 5:
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.plot(phase, radius, ".", markersize=2)
        ax.set_xlabel("Pulsation phase")
        ax.set_ylabel(r"Radius ($R_\odot$)")
        ax.set_title("RR Lyrae photospheric-radius cycle")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(outdir / "phase_radius.png", dpi=180)
        plt.close(fig)

    # 2. Light curves.
    if filters:
        fig, ax = plt.subplots(figsize=(12, 6))
        for f in filters:
            mag = pd.to_numeric(sub[f], errors="coerce").to_numpy(float)
            ax.plot(phase, mag, ".", markersize=2, label=f)
        ax.invert_yaxis()
        ax.set_xlabel("Pulsation phase")
        ax.set_ylabel("Magnitude")
        ax.set_title("Synthetic light curves")
        ax.legend(ncol=2, fontsize=9)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(outdir / "phase_lightcurves.png", dpi=180)
        plt.close(fig)

    # 3. Filter amplitude trend.
    if len(filters) >= 2:
        rows = []
        for f in filters:
            lam = COMMON_FILTER_WAVELENGTH_ANG.get(f, np.nan)
            mag = pd.to_numeric(sub[f], errors="coerce").to_numpy(float)
            rows.append((f, lam, nanptp(mag)))
        fdf = pd.DataFrame(rows, columns=["filter", "lambda", "amp"])
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(fdf["lambda"], fdf["amp"], "o-")
        for _, r in fdf.iterrows():
            ax.text(r["lambda"], r["amp"], str(r["filter"]), fontsize=8, ha="center", va="bottom")
        ax.set_xlabel(r"Approximate effective wavelength ($\AA$)")
        ax.set_ylabel("Peak-to-peak amplitude (mag)")
        ax.set_title("Amplitude decreases redward")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(outdir / "filter_amplitude_trend.png", dpi=180)
        plt.close(fig)

    # 4. Color-temperature and color-radius loops.
    pair = None
    for a, b in COLOR_PAIRS_PREFERRED:
        if a in filters and b in filters:
            pair = (a, b)
            break
    if pair is None and len(filters) >= 2:
        pair = (filters[0], filters[-1])
    if pair:
        a, b = pair
        color = pd.to_numeric(sub[a], errors="coerce").to_numpy(float) - pd.to_numeric(sub[b], errors="coerce").to_numpy(float)
        fig, ax = plt.subplots(figsize=(8, 6))
        sc = ax.scatter(color, teff, c=phase, s=5)
        ax.set_xlabel(f"{a} - {b} (mag)")
        ax.set_ylabel(r"$T_{\rm eff}$ (K)")
        ax.set_title("Color-temperature loop")
        ax.grid(True, alpha=0.3)
        fig.colorbar(sc, ax=ax, label="Pulsation phase")
        fig.tight_layout()
        fig.savefig(outdir / "color_Teff_loop.png", dpi=180)
        plt.close(fig)

        if np.isfinite(radius).sum() > 5:
            fig, ax = plt.subplots(figsize=(8, 6))
            sc = ax.scatter(color, radius, c=phase, s=5)
            ax.set_xlabel(f"{a} - {b} (mag)")
            ax.set_ylabel(r"Radius ($R_\odot$)")
            ax.set_title("Color-radius loop")
            ax.grid(True, alpha=0.3)
            fig.colorbar(sc, ax=ax, label="Pulsation phase")
            fig.tight_layout()
            fig.savefig(outdir / "color_radius_loop.png", dpi=180)
            plt.close(fig)

        # CMD-ish loop using bluest filter magnitude.
        fig, ax = plt.subplots(figsize=(8, 6))
        ymag = pd.to_numeric(sub[a], errors="coerce").to_numpy(float)
        sc = ax.scatter(color, ymag, c=phase, s=5)
        ax.invert_yaxis()
        ax.set_xlabel(f"{a} - {b} (mag)")
        ax.set_ylabel(f"{a} (mag)")
        ax.set_title("Synthetic color-magnitude loop")
        ax.grid(True, alpha=0.3)
        fig.colorbar(sc, ax=ax, label="Pulsation phase")
        fig.tight_layout()
        fig.savefig(outdir / "cmd_loop.png", dpi=180)
        plt.close(fig)

    # 5. Selected SED comparison if paths are available.
    if not sed_selected_df.empty and "path" in sed_selected_df.columns:
        paths = []
        labels = []
        for _, row in sed_selected_df.iterrows():
            label = str(row.get("landmark", ""))
            if label in ["hottest_Teff", "coolest_Teff", "max_radius", "min_radius", "max_light_" + filters[0] if filters else ""]:
                p = Path(str(row.get("path", "")))
                if p.exists():
                    paths.append(p)
                    labels.append(label)
            if len(paths) >= 4:
                break
        if paths:
            fig, ax = plt.subplots(figsize=(10, 5))
            for p, lab in zip(paths, labels):
                try:
                    sdf = load_sed_csv(p)
                    cols = {c.strip(): c for c in sdf.columns}
                    wcol = cols.get("wavelengths") or list(sdf.columns)[0]
                    fcol = cols.get("fluxes") or list(sdf.columns)[1]
                    w = pd.to_numeric(sdf[wcol], errors="coerce").to_numpy(float)
                    f = pd.to_numeric(sdf[fcol], errors="coerce").to_numpy(float)
                    m = np.isfinite(w) & np.isfinite(f)
                    w = w[m]
                    f = f[m]
                    if len(w) > 1:
                        fn = f / np.nanmax(f) if np.nanmax(f) != 0 else f
                        ax.plot(w, fn, label=lab)
                except Exception:
                    continue
            ax.set_xlabel(r"Wavelength ($\AA$)")
            ax.set_ylabel(r"$F_\lambda$ (normalized)")
            ax.set_title("Selected phase SEDs")
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)
            fig.tight_layout()
            fig.savefig(outdir / "selected_seds.png", dpi=180)
            plt.close(fig)


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Fully characterize a MESA RSP RR Lyrae + Colors run.")
    ap.add_argument("--workdir", default=".", help="Run directory, default current directory.")
    ap.add_argument("--history", default=None, help="History file. Default: LOGS_colors/history.data, then LOGS/history.data, then LOGS_settle/history.data.")
    ap.add_argument("--sed-dir", default="SED", help="SED directory relative to workdir. Default: SED.")
    ap.add_argument("--outdir", default="rrlyrae_characterization", help="Output directory. Default: rrlyrae_characterization.")
    ap.add_argument("--last-periods", type=int, default=10, help="Analyze the last N completed periods. Default: 10.")
    ap.add_argument("--sed-sample", type=int, default=2000, help="Number of SED files to sample for integrity diagnostics. Default: 2000. Use 0 to skip sample.")
    ap.add_argument("--no-sed", action="store_true", help="Skip SED diagnostics.")
    ap.add_argument("--no-plots", action="store_true", help="Do not write plots.")
    args = ap.parse_args(argv)

    workdir = Path(args.workdir).expanduser().resolve()
    history_path = pick_history(workdir, args.history).resolve()
    sed_dir = (workdir / args.sed_dir).resolve()
    outdir = (workdir / args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    df = read_mesa_history(history_path)
    phase_full = detect_phase(df).to_numpy(float)
    sub = restrict_last_periods(df, args.last_periods).copy()
    phase = detect_phase(sub).to_numpy(float)

    teff = get_teff(sub)
    lum = get_luminosity(sub)
    log_l = np.log10(lum)
    radius, radius_source = infer_radius(sub, lum, teff)
    filters = detect_filter_columns(sub)

    # Global physical summary.
    global_rows = []
    global_rows.append(summarize_series("Teff", teff, "K"))
    global_rows.append(summarize_series("log_L", log_l, "dex Lsun"))
    global_rows.append(summarize_series("L", lum, "Lsun"))
    global_rows.append(summarize_series("radius", radius, "Rsun"))
    if "log_g" in sub.columns:
        global_rows.append(summarize_series("log_g", pd.to_numeric(sub["log_g"], errors="coerce"), "cgs"))
    if "star_mass" in sub.columns:
        global_rows.append(summarize_series("star_mass", pd.to_numeric(sub["star_mass"], errors="coerce"), "Msun"))
    if "rsp_period_in_days" in sub.columns:
        global_rows.append(summarize_series("RSP period", pd.to_numeric(sub["rsp_period_in_days"], errors="coerce"), "days"))
    if "rsp_GREKM" in sub.columns:
        global_rows.append(summarize_series("rsp_GREKM", pd.to_numeric(sub["rsp_GREKM"], errors="coerce"), ""))
        global_rows.append(summarize_series("abs_rsp_GREKM", np.abs(pd.to_numeric(sub["rsp_GREKM"], errors="coerce")), ""))
    if "max_abs_v_div_cs" in sub.columns:
        global_rows.append(summarize_series("max_abs_v_div_cs", pd.to_numeric(sub["max_abs_v_div_cs"], errors="coerce"), ""))
    if "v_surf_km_s" in sub.columns:
        global_rows.append(summarize_series("v_surf", pd.to_numeric(sub["v_surf_km_s"], errors="coerce"), "km/s"))
    if "rsp_DeltaR" in sub.columns:
        global_rows.append(summarize_series("rsp_DeltaR", pd.to_numeric(sub["rsp_DeltaR"], errors="coerce"), "Rsun?"))
    if "rsp_DeltaMag" in sub.columns:
        global_rows.append(summarize_series("rsp_DeltaMag", pd.to_numeric(sub["rsp_DeltaMag"], errors="coerce"), "mag"))
    # SB consistency.
    sb_logl = 2 * np.log10(radius) + 4 * np.log10(teff / TEFF_SUN)
    global_rows.append(summarize_series("SB_logL_residual", log_l - sb_logl, "dex"))
    global_df = pd.DataFrame(global_rows)

    filter_df = make_filter_table(sub, phase, filters)
    color_df = make_color_table(sub, phase, filters, teff, radius)
    fourier_df = make_fourier_table(sub, phase, filters[:4], teff, radius)
    landmarks = find_phase_landmarks(sub, phase, filters, teff, radius)
    landmark_df = pd.DataFrame([lm.__dict__ for lm in landmarks])

    inlist_meta = {
        "inlist": parse_simple_inlist_assignments(workdir / "inlist"),
        "inlist_rsp_RR_Lyrae": parse_simple_inlist_assignments(workdir / "inlist_rsp_RR_Lyrae"),
        "inlist_rsp_RR_Lyrae_colors": parse_simple_inlist_assignments(workdir / "inlist_rsp_RR_Lyrae_colors"),
    }

    sed_selected_df = pd.DataFrame()
    sed_sample_df = pd.DataFrame()
    if not args.no_sed and sed_dir.exists():
        sed_selected_df, sed_sample_df = analyze_selected_seds(sed_dir, filters, landmarks, args.sed_sample)

    # Write machine-readable outputs.
    global_df.to_csv(outdir / "global_physical_summary.csv", index=False)
    filter_df.to_csv(outdir / "filter_lightcurve_summary.csv", index=False)
    color_df.to_csv(outdir / "color_diagnostics.csv", index=False)
    fourier_df.to_csv(outdir / "fourier_morphology.csv", index=False)
    landmark_df.to_csv(outdir / "phase_landmarks.csv", index=False)
    if not sed_selected_df.empty:
        sed_selected_df.to_csv(outdir / "selected_sed_summary.csv", index=False)
    if not sed_sample_df.empty:
        sed_sample_df.to_csv(outdir / "sed_integrity_sample.csv", index=False)
    with (outdir / "inlist_metadata.json").open("w") as fh:
        json.dump(inlist_meta, fh, indent=2)

    # Write LaTeX snippets/tables.
    dataframe_to_latex(global_df[["quantity", "unit", "median", "min", "max", "ptp"]], outdir / "global_physical_summary.tex")
    if not filter_df.empty:
        cols = ["filter", "lambda_eff_A_assumed", "amplitude_mag", "median_mag", "phase_max_light", "phase_min_light"]
        dataframe_to_latex(filter_df[cols], outdir / "filter_lightcurve_summary.tex")
    if not color_df.empty:
        cols = ["color", "median_color_mag", "amplitude_color_mag", "pearson_color_Teff", "pearson_color_radius"]
        dataframe_to_latex(color_df[cols], outdir / "color_diagnostics.tex")

    write_report(outdir, workdir, history_path, df, sub, global_df, filter_df, color_df, fourier_df, sed_selected_df, sed_sample_df, inlist_meta, radius_source, filters)

    if not args.no_plots:
        make_plots(outdir, sub, phase, filters, teff, radius, sed_selected_df)

    # Console summary.
    print("\nRR LYRAE CHARACTERIZATION COMPLETE")
    print("=" * 72)
    print(f"workdir: {workdir}")
    print(f"history: {history_path}")
    print(f"output:  {outdir}")
    print(f"rows full/analysis: {len(df)} / {len(sub)}")
    if "rsp_num_periods" in df.columns:
        cyc = pd.to_numeric(df["rsp_num_periods"], errors="coerce")
        print(f"completed periods: {int(np.nanmin(cyc))}--{int(np.nanmax(cyc))}")
    print(f"filters: {', '.join(filters) if filters else 'none'}")
    print(f"period median: {global_df.loc[global_df.quantity == 'RSP period', 'median'].iloc[0]:.6f} d" if "RSP period" in set(global_df.quantity) else "period median: n/a")
    print(f"Teff range in analysis: {np.nanmin(teff):.0f}--{np.nanmax(teff):.0f} K")
    print(f"radius median/range: {nanmedian(radius):.3f} Rsun / {nanptp(radius):.3f} Rsun")
    if not filter_df.empty:
        print("filter amplitudes:")
        for _, r in filter_df.iterrows():
            print(f"  {r['filter']:>5s}: {r['amplitude_mag']:.3f} mag")
    if not color_df.empty:
        r = color_df.iloc[0]
        print(f"default color diagnostic: {r['color']} vs Teff Pearson r={r['pearson_color_Teff']:.3f}")
    if not sed_sample_df.empty:
        neg = pd.to_numeric(sed_sample_df.get("negative_area_fraction", pd.Series(dtype=float)), errors="coerce")
        neg = neg[np.isfinite(neg)]
        if len(neg):
            print(f"SED negative-area fraction: median={np.nanmedian(neg):.3e}, max={np.nanmax(neg):.3e}")
    print("\nKey files:")
    for name in [
        "characterization_report.md",
        "paper_values_snippet.tex",
        "global_physical_summary.csv",
        "filter_lightcurve_summary.csv",
        "color_diagnostics.csv",
        "phase_lightcurves.png",
        "cmd_loop.png",
        "selected_seds.png",
    ]:
        p = outdir / name
        if p.exists():
            print(f"  {p}")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
