#!/usr/bin/env python3
"""
Run a small FN Lyr RSP+Colors parameter grid.

Run 4: eddy-viscosity (RSP_alfam) sweep.

Stellar parameters (Teff, L) and kick are FIXED at the best converged
model from run 3 (T6700/L47.7/kick40, A ~= 0.95, period err ~1.2%).
The only axis varied here is RSP_alfam, the RSP eddy viscosity.

Rationale: T/L is exhausted as an amplitude lever -- going hotter buys
amplitude but pushes the period match the wrong way. Lowering alfam
raises the saturation amplitude at FIXED period, and is the same knob
governing the shock/shoulder damping. alfam=0.25 is the standard RSP
default and serves as the control: it must reproduce the run-3 amplitude.

This script:
  - edits inlist_rsp_RR_Lyrae
  - runs ./rn
  - saves LOGS/history.data, final.mod, inlist, and stdout log
  - measures late-cycle Kepler amplitude, period, and convergence
  - writes a summary CSV

Use from:
  ~/MESA/MESA_Colors_Tests/rsp_RR_Lyrae/FN_Lyr
"""

from __future__ import annotations

import csv
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


OBS_PERIOD = 0.527398471
OBS_AMP = 1.0464

# Amplitude floor for the convective-parameter sweeps. The goal is to clean
# the waveform (merge the shock-split peak) WITHOUT giving up the amplitude
# that alfam=0.15 bought. Any run whose CONVERGED amplitude (last-10 median)
# falls below this is flagged amp_above_floor=False and should be rejected:
# it cleaned the shape only by losing amplitude, which is just alfam in
# disguise. gammar is a genuine damping term so this floor is load-bearing
# here -- it is the test of whether clean-waveform AND 1.05-amplitude are
# simultaneously reachable at this M/L/Teff.
AMP_FLOOR = 1.03

INLIST = Path("inlist_rsp_RR_Lyrae")
RUN_CMD = ["./rn"]

RESULTS_DIR = Path("fn_lyr_grid_gammar")
RESULTS_DIR.mkdir(exist_ok=True)


@dataclass
class Trial:
    name: str
    teff: float
    lum: float
    kick: float
    alfam: float
    alfat: float
    gammar: float
    periods: int
    max_model_number: int = 280000


# Anchor: the alfam=0.15 model that matched FN Lyr's amplitude (1.05 mag) at
# the correct period, but shows a shock-split maximum (vmax/cs ~7.9) that the
# data lacks. The alfat sweep ruled out turbulent flux as the shock lever
# (vmax/cs and waveform were flat across alfat). The shock is a surface
# VELOCITY shock at maximum compression, so the lever is radiative damping:
# RSP_gammar (default 0.0 in the inherited radiative set; ~1.0 in the full
# convective RSP set). alfat is set back to 0 -- it demonstrably did nothing,
# so keep one knob moving at a time.
#
# gammar = 0.0 is the CONTROL: it must reproduce the shocky alfam=0.15 fit
# (amp ~1.05, double-peaked, r~0.982). Unlike alfat, gammar damps the driving,
# so amplitude is expected to drop as gammar rises -- watch amp_above_floor.
ANCHOR_TEFF = 6700.0
ANCHOR_L = 47.7
ANCHOR_KICK = 40.0
ANCHOR_ALFAM = 0.15
ANCHOR_ALFAT = 0.0

TRIALS = [
    Trial("T6700_L47p7_kick40_alfam0p15_gammar0p0", ANCHOR_TEFF, ANCHOR_L, ANCHOR_KICK, ANCHOR_ALFAM, ANCHOR_ALFAT, 0.0, 180),
    Trial("T6700_L47p7_kick40_alfam0p15_gammar0p5", ANCHOR_TEFF, ANCHOR_L, ANCHOR_KICK, ANCHOR_ALFAM, ANCHOR_ALFAT, 0.5, 180),
    Trial("T6700_L47p7_kick40_alfam0p15_gammar1p0", ANCHOR_TEFF, ANCHOR_L, ANCHOR_KICK, ANCHOR_ALFAM, ANCHOR_ALFAT, 1.0, 180),
    Trial("T6700_L47p7_kick40_alfam0p15_gammar2p0", ANCHOR_TEFF, ANCHOR_L, ANCHOR_KICK, ANCHOR_ALFAM, ANCHOR_ALFAT, 2.0, 180),
]


def replace_or_insert(text: str, pattern: str, replacement: str, anchor: str | None = None) -> str:
    """
    Replace a Fortran namelist assignment. If not found, insert after anchor.
    """
    new, n = re.subn(pattern, replacement, text, flags=re.MULTILINE)
    if n > 0:
        return new

    if anchor is not None and anchor in text:
        return text.replace(anchor, anchor + "\n   " + replacement)

    return text + "\n   " + replacement + "\n"


def set_inlist(trial: Trial) -> None:
    s = INLIST.read_text()

    # Force fresh model creation.
    s = replace_or_insert(
        s,
        r"^\s*create_RSP_model\s*=.*$",
        "   create_RSP_model = .true.",
    )
    s = replace_or_insert(
        s,
        r"^\s*load_saved_model\s*=.*$",
        "   load_saved_model = .false.",
    )

    # Remove saved-model filename from star_job if present.
    s = re.sub(r"^\s*load_model_filename\s*=.*\n", "", s, flags=re.MULTILINE)
    s = re.sub(r"^\s*saved_model_name\s*=.*\n", "", s, flags=re.MULTILINE)

    # Stellar parameters (held fixed at the anchor, but still written
    # explicitly so the inlist is self-documenting per run).
    s = replace_or_insert(s, r"^\s*RSP_Teff\s*=.*$", f"   RSP_Teff = {trial.teff:.0f}d0")
    s = replace_or_insert(s, r"^\s*RSP_L\s*=.*$", f"   RSP_L = {trial.lum:.4g}d0")

    # Kick (held fixed).
    s = replace_or_insert(
        s,
        r"^\s*RSP_kick_vsurf_km_per_sec\s*=.*$",
        f"   RSP_kick_vsurf_km_per_sec = {trial.kick:.4g}d0",
        anchor="RSP_Z",
    )

    # Eddy viscosity -- the swept axis for run 4.
    # If RSP_alfam already exists, overwrite it. If not (the usual case,
    # since the inlist inherits the 0.25 default rather than setting it),
    # insert it. We anchor on the FULL RSP_Z assignment line, not on a bare
    # keyword: replace_or_insert's anchor branch does a plain substring
    # replace, so anchoring on a bare token like RSP_kick_vsurf_km_per_sec
    # splits that token mid-line and produces a namelist read error
    # ("Equal sign must follow namelist object name"). Anchoring on the
    # whole "RSP_Z = 0.0004d0" line is safe because that exact string is not
    # a prefix of any other keyword, and it keeps the new line inside
    # &controls (before the closing '/').
    z_anchor = "RSP_Z = 0.0004d0"
    if "RSP_alfam" in s:
        s = replace_or_insert(
            s,
            r"^\s*RSP_alfam\s*=.*$",
            f"   RSP_alfam = {trial.alfam:.4g}d0",
        )
    else:
        s = replace_or_insert(
            s,
            r"^\s*RSP_alfam\s*=.*$",
            f"   RSP_alfam = {trial.alfam:.4g}d0",
            anchor=z_anchor,
        )

    # Turbulent flux -- the swept axis for this pass. Default in the inlist's
    # inherited RSP set is 0.0 (purely radiative). Turning it on smooths the
    # convective flux through the ionization zones and damps the shock-split
    # peak / pre-rise bump without touching the amplitude-setting alfam.
    # Same insert discipline as alfam: overwrite if present, else insert
    # after a full-line anchor. RSP_alfam was just written above and is
    # therefore always present now, so anchor on its full assignment line
    # (which contains '= ...d0' and cannot be split into a bare keyword).
    alfam_anchor = f"RSP_alfam = {trial.alfam:.4g}d0"
    if "RSP_alfat" in s:
        s = replace_or_insert(
            s,
            r"^\s*RSP_alfat\s*=.*$",
            f"   RSP_alfat = {trial.alfat:.4g}d0",
        )
    else:
        s = replace_or_insert(
            s,
            r"^\s*RSP_alfat\s*=.*$",
            f"   RSP_alfat = {trial.alfat:.4g}d0",
            anchor=alfam_anchor,
        )

    # Radiative damping -- the swept axis for this pass. Default in the
    # inherited radiative RSP set is 0.0. This is the lever for the surface
    # velocity shock that splits the maximum: raising it damps the shock at
    # maximum compression. Unlike alfat it acts on the driving, so amplitude
    # is expected to fall as gammar rises (the AMP_FLOOR gate catches that).
    # Same insert discipline: overwrite if present, else insert after a full
    # assignment-line anchor. RSP_alfat was just written above and is always
    # present now, so anchor on its full assignment line.
    alfat_anchor = f"RSP_alfat = {trial.alfat:.4g}d0"
    if "RSP_gammar" in s:
        s = replace_or_insert(
            s,
            r"^\s*RSP_gammar\s*=.*$",
            f"   RSP_gammar = {trial.gammar:.4g}d0",
        )
    else:
        s = replace_or_insert(
            s,
            r"^\s*RSP_gammar\s*=.*$",
            f"   RSP_gammar = {trial.gammar:.4g}d0",
            anchor=alfat_anchor,
        )

    # Target period and number of periods.
    s = replace_or_insert(
        s,
        r"^\s*x_ctrl\(1\)\s*=.*$",
        f"   x_ctrl(1) = {OBS_PERIOD:.9f}d0    ! observed FN Lyr period",
    )
    s = replace_or_insert(
        s,
        r"^\s*x_integer_ctrl\(1\)\s*=.*$",
        f"   x_integer_ctrl(1) = {trial.periods}   ! which period to check",
    )

    # Avoid stopping before the requested periods.
    s = replace_or_insert(
        s,
        r"^\s*max_model_number\s*=.*$",
        f"   max_model_number = {trial.max_model_number}",
    )

    # Keep output manageable.
    s = replace_or_insert(s, r"^\s*make_csv\s*=.*$", "   make_csv = .false.")
    s = replace_or_insert(s, r"^\s*sed_per_model\s*=.*$", "   sed_per_model = .false.")

    # Save a generic final.mod because wrappers expect it.
    s = replace_or_insert(
        s,
        r"^\s*save_model_filename\s*=.*$",
        "      save_model_filename = 'final.mod'",
    )

    INLIST.write_text(s)


def read_history(path: Path) -> pd.DataFrame:
    lines = path.read_text().splitlines()
    header_idx = next(i for i, line in enumerate(lines) if "model_number" in line.split())
    return pd.read_csv(path, sep=r"\s+", skiprows=header_idx)


def analyze_history(path: Path) -> dict:
    h = read_history(path)
    mag_col = "Kepler" if "Kepler" in h.columns else "K"

    if "rsp_num_periods" not in h.columns:
        raise RuntimeError("No rsp_num_periods in history.data")

    nmax = int(np.nanmax(h["rsp_num_periods"]))
    rows = []

    for n in range(max(0, nmax - 25), nmax + 1):
        m = h[h["rsp_num_periods"].astype(int) == n].copy()
        if len(m) < 200:
            continue

        phase_span = np.nanmax(m["rsp_phase"]) - np.nanmin(m["rsp_phase"]) if "rsp_phase" in m else np.nan
        amp = np.nanmax(m[mag_col]) - np.nanmin(m[mag_col])
        p_est = np.nanmax(m["star_age_day"]) - np.nanmin(m["star_age_day"])
        dteff = np.nanmax(m["effective_T"]) - np.nanmin(m["effective_T"])
        dlogr = np.nanmax(m["log_R"]) - np.nanmin(m["log_R"])
        vmax = np.nanmax(np.abs(m["max_abs_v_div_cs"])) if "max_abs_v_div_cs" in m else np.nan

        rows.append(
            {
                "cycle": n,
                "npts": len(m),
                "phase_span": phase_span,
                "period": p_est,
                "amp": amp,
                "dteff": dteff,
                "dlogr": dlogr,
                "vmax_cs": vmax,
                "deltaR_med": np.nanmedian(m["rsp_DeltaR"]) if "rsp_DeltaR" in m else np.nan,
            }
        )

    if not rows:
        return {
            "nmax": nmax,
            "mag_col": mag_col,
            "best_cycle": np.nan,
            "best_period": np.nan,
            "best_amp": np.nan,
            "best_score": np.inf,
            "median_amp_last": np.nan,
            "median_period_last": np.nan,
            "max_amp_last": np.nan,
            "amp_last10_median": np.nan,
            "amp_last10_range": np.nan,
            "amp_slope_last10": np.nan,
            "vmax_cs_last": np.nan,
            "amp_above_floor": False,
        }

    df = pd.DataFrame(rows)

    # Favor amplitude near FN Lyr, period near FN Lyr, complete cycles,
    # and avoid very shocky cycles where possible.
    df["amp_err"] = np.abs(df["amp"] - OBS_AMP) / OBS_AMP
    df["period_err"] = np.abs(df["period"] - OBS_PERIOD) / OBS_PERIOD
    df["shock_penalty"] = np.maximum(df["vmax_cs"] - 1.5, 0.0)
    df["score"] = 2.0 * df["amp_err"] + df["period_err"] + 0.25 * df["shock_penalty"]

    best = df.sort_values("score").iloc[0]

    # Convergence metrics over the last 10 available cycles.
    # These are the real acceptance gate: at a limit cycle the amplitude
    # range and slope over the last cycles should be ~0.
    last10 = df.sort_values("cycle").tail(10)
    amp_last10_median = float(np.nanmedian(last10["amp"]))
    amp_last10_range = float(np.nanmax(last10["amp"]) - np.nanmin(last10["amp"]))
    if len(last10) >= 2:
        amp_slope_last10 = float(
            np.polyfit(last10["cycle"].to_numpy(dtype=float),
                       last10["amp"].to_numpy(dtype=float), 1)[0]
        )
    else:
        amp_slope_last10 = np.nan
    vmax_cs_last = float(np.nanmax(last10["vmax_cs"])) if "vmax_cs" in last10 else np.nan

    # Amplitude-floor gate. Judge the CONVERGED amplitude (last-10 median),
    # not best_amp, so a single high cycle can't sneak past the floor.
    amp_above_floor = bool(amp_last10_median >= AMP_FLOOR)

    return {
        "nmax": nmax,
        "mag_col": mag_col,
        "best_cycle": int(best["cycle"]),
        "best_period": float(best["period"]),
        "best_amp": float(best["amp"]),
        "best_score": float(best["score"]),
        "median_amp_last": float(np.nanmedian(df["amp"])),
        "median_period_last": float(np.nanmedian(df["period"])),
        "max_amp_last": float(np.nanmax(df["amp"])),
        "amp_last10_median": amp_last10_median,
        "amp_last10_range": amp_last10_range,
        "amp_slope_last10": amp_slope_last10,
        "vmax_cs_last": vmax_cs_last,
        "amp_above_floor": amp_above_floor,
    }


def clean_run_products() -> None:
    for name in ["LOGS", "photos", "png", "SED"]:
        p = Path(name)
        if p.exists():
            shutil.rmtree(p)
    Path("SED").mkdir(exist_ok=True)
    Path("png").mkdir(exist_ok=True)
    if Path("final.mod").exists():
        Path("final.mod").unlink()


def save_trial_outputs(trial: Trial, stdout_text: str, summary: dict | None) -> None:
    out = RESULTS_DIR / trial.name
    out.mkdir(exist_ok=True)

    (out / "stdout.log").write_text(stdout_text)

    if INLIST.exists():
        shutil.copy2(INLIST, out / "inlist_rsp_RR_Lyrae")

    if Path("final.mod").exists():
        shutil.copy2("final.mod", out / "final.mod")

    if Path("LOGS/history.data").exists():
        shutil.copy2("LOGS/history.data", out / "history.data")

    if summary is not None:
        with (out / "summary.txt").open("w") as f:
            for k, v in summary.items():
                f.write(f"{k}: {v}\n")


def run_trial(trial: Trial) -> dict:
    print(f"\n=== Running {trial.name} ===")
    print(f"Teff={trial.teff} L={trial.lum} kick={trial.kick} "
          f"alfam={trial.alfam} alfat={trial.alfat} gammar={trial.gammar} "
          f"periods={trial.periods}")

    set_inlist(trial)
    clean_run_products()

    proc = subprocess.run(
        RUN_CMD,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )

    stdout_text = proc.stdout

    summary = None
    status = "failed"

    try:
        if Path("LOGS/history.data").exists():
            summary = analyze_history(Path("LOGS/history.data"))
            status = "ok" if proc.returncode == 0 else "nonzero_return"
        else:
            status = "no_history"
    except Exception as exc:
        summary = {"analysis_error": repr(exc)}
        status = "analysis_failed"

    save_trial_outputs(trial, stdout_text, summary)

    row = {
        "name": trial.name,
        "teff": trial.teff,
        "lum": trial.lum,
        "kick": trial.kick,
        "alfam": trial.alfam,
        "alfat": trial.alfat,
        "gammar": trial.gammar,
        "periods_requested": trial.periods,
        "returncode": proc.returncode,
        "status": status,
    }

    if summary:
        row.update(summary)

    print("status:", status)
    for key in ["best_cycle", "best_period", "best_amp",
                "median_amp_last", "max_amp_last",
                "amp_last10_median", "amp_last10_range", "amp_slope_last10",
                "vmax_cs_last", "amp_above_floor", "best_score"]:
        if key in row:
            print(f"{key}: {row[key]}")

    return row


def main() -> None:
    rows = []

    original = INLIST.read_text()
    (RESULTS_DIR / "inlist_original_before_grid").write_text(original)

    try:
        for trial in TRIALS:
            row = run_trial(trial)
            rows.append(row)

            with (RESULTS_DIR / "grid_summary.csv").open("w", newline="") as f:
                fieldnames = sorted({k for r in rows for k in r.keys()})
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(rows)

    finally:
        # Leave the last run in place, but keep original safely stored.
        pass

    print("\nWrote:", RESULTS_DIR / "grid_summary.csv")
    print("Sort with:")
    print("  column -s, -t fn_lyr_grid_gammar/grid_summary.csv | less -S")


if __name__ == "__main__":
    main()
