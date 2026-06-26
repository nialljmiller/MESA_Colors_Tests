#!/usr/bin/env python3
"""
Compare FN Lyr / KIC 6936115 Kepler photometry with a MESA RSP+Colors model.

Inputs expected in the current directory:
  FN_Lyr_KIC6936115_kepler_binned.csv
  LOGS/history.data

Outputs:
  FN_Lyr_MESA_vs_Kepler_normalized.png
  FN_Lyr_MESA_vs_Kepler_residuals.png
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


OBS_FILE = "FN_Lyr_KIC6936115_kepler_binned.csv"
MESA_FILE = "LOGS/history.data"

# Observed FN Lyr period and best MESA period from your successful run.
P_OBS = 0.527398471
P_MODEL = 0.52690624947998876

# Your MESA filter column is currently named "K".
# This likely corresponds to the Kepler filter file name.
FORCE_MAG_COL = "K"


def read_mesa_history(path):
    """
    Read a MESA history.data file into a pandas DataFrame.

    MESA history files have metadata at the top. The actual column-name line
    is the one containing 'model_number'.
    """
    path = Path(path)
    lines = path.read_text().splitlines()

    header_idx = None
    for i, line in enumerate(lines):
        toks = line.split()
        if "model_number" in toks:
            header_idx = i
            break

    if header_idx is None:
        raise RuntimeError("Could not find MESA history header line containing model_number")

    df = pd.read_csv(path, sep=r"\s+", skiprows=header_idx)
    return df


def periodic_interp(x, y, xnew):
    """
    Interpolate a periodic phase curve y(x) onto xnew.

    x, xnew are phases on [0, 1).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    xnew = np.asarray(xnew, dtype=float)

    good = np.isfinite(x) & np.isfinite(y)
    x = x[good]
    y = y[good]

    order = np.argsort(x)
    x = x[order]
    y = y[order]

    # Wrap curve to avoid edge artifacts.
    xx = np.concatenate([x - 1.0, x, x + 1.0])
    yy = np.concatenate([y, y, y])

    return np.interp(xnew, xx, yy)


def normalize_mag(m):
    """
    Normalize a magnitude curve by subtracting the mean and dividing by
    robust amplitude.

    Magnitudes are left as magnitudes; the plotting routine inverts the y-axis.
    """
    m = np.asarray(m, dtype=float)
    good = np.isfinite(m)
    m = m[good]

    amp = np.nanpercentile(m, 99) - np.nanpercentile(m, 1)
    norm = (m - np.nanmean(m)) / amp

    return norm, amp


def choose_time_column(mesa):
    """
    Select the best time column from MESA history.data.
    """
    candidates = [
        "star_age_day",
        "age_days",
        "star_age",
        "age",
    ]

    for col in candidates:
        if col in mesa.columns:
            return col

    raise RuntimeError(
        "No usable time column found. "
        f"Tried: {candidates}\n"
        f"Available columns: {list(mesa.columns)}"
    )


def choose_mag_column(mesa):
    """
    Select the photometric magnitude column.

    For your current run, this is expected to be 'K'.
    """
    if FORCE_MAG_COL is not None:
        if FORCE_MAG_COL not in mesa.columns:
            raise RuntimeError(
                f"Forced magnitude column '{FORCE_MAG_COL}' not found.\n"
                f"Available columns: {list(mesa.columns)}"
            )
        return FORCE_MAG_COL

    candidates = [
        c for c in mesa.columns
        if "kepler" in c.lower()
        or c.lower() in {"kp", "kplr", "k"}
    ]

    if not candidates:
        raise RuntimeError(
            "No obvious Kepler/Kp/K column found in MESA history.data.\n"
            f"Available columns: {list(mesa.columns)}"
        )

    return candidates[0]


def extract_final_cycle(mesa, time_col, mag_col, period):
    """
    Extract the final full pulsation cycle from MESA history.
    """
    t = mesa[time_col].to_numpy(dtype=float)
    m = mesa[mag_col].to_numpy(dtype=float)

    good = np.isfinite(t) & np.isfinite(m)
    t = t[good]
    m = m[good]

    if len(t) < 10:
        raise RuntimeError("Too few valid MESA points after filtering.")

    t_end = np.nanmax(t)
    mask = t >= (t_end - period)

    t_cycle = t[mask]
    m_cycle = m[mask]

    if len(t_cycle) < 20:
        raise RuntimeError(
            "Final-cycle extraction produced too few points. "
            "Check time units and period."
        )

    phase = ((t_cycle - t_cycle.min()) / period) % 1.0

    order = np.argsort(phase)
    return phase[order], m_cycle[order]


def best_phase_align(phase_obs, obs_norm, phase_model, model_norm):
    """
    Find phase shift that maximizes R^2 between normalized observed and model curves.
    """
    shifts = np.linspace(0.0, 1.0, 3000, endpoint=False)

    best = None

    for shift in shifts:
        shifted_phase = (phase_obs - shift) % 1.0
        model_at_obs = periodic_interp(phase_model, model_norm, shifted_phase)

        resid = obs_norm - model_at_obs
        ss_res = np.sum(resid**2)
        ss_tot = np.sum((obs_norm - np.mean(obs_norm))**2)

        r2 = 1.0 - ss_res / ss_tot
        corr = np.corrcoef(obs_norm, model_at_obs)[0, 1]
        rms = np.sqrt(np.mean(resid**2))

        if best is None or r2 > best["r2"]:
            best = {
                "shift": shift,
                "r2": r2,
                "corr": corr,
                "rms": rms,
                "model_at_obs": model_at_obs,
                "resid": resid,
            }

    return best


def main():
    obs_path = Path(OBS_FILE)
    mesa_path = Path(MESA_FILE)

    if not obs_path.exists():
        raise FileNotFoundError(f"Missing observed binned file: {obs_path}")

    if not mesa_path.exists():
        raise FileNotFoundError(f"Missing MESA history file: {mesa_path}")

    obs = pd.read_csv(obs_path)
    mesa = read_mesa_history(mesa_path)

    if "phase" not in obs.columns or "mag" not in obs.columns:
        raise RuntimeError(
            "Observed file must contain columns named 'phase' and 'mag'. "
            f"Found: {list(obs.columns)}"
        )

    time_col = choose_time_column(mesa)
    mag_col = choose_mag_column(mesa)

    print("Using:")
    print(f"  Observed file: {OBS_FILE}")
    print(f"  MESA file:     {MESA_FILE}")
    print(f"  Time column:   {time_col}")
    print(f"  Mag column:    {mag_col}")
    print()

    # Observed binned Kepler curve.
    phase_obs = obs["phase"].to_numpy(dtype=float)
    mag_obs = obs["mag"].to_numpy(dtype=float)

    good_obs = np.isfinite(phase_obs) & np.isfinite(mag_obs)
    phase_obs = phase_obs[good_obs]
    mag_obs = mag_obs[good_obs]

    order_obs = np.argsort(phase_obs)
    phase_obs = phase_obs[order_obs]
    mag_obs = mag_obs[order_obs]

    # MESA final cycle.
    phase_model, mag_model = extract_final_cycle(
        mesa=mesa,
        time_col=time_col,
        mag_col=mag_col,
        period=P_MODEL,
    )

    # Normalize both curves.
    obs_norm, obs_amp = normalize_mag(mag_obs)
    model_norm, model_amp = normalize_mag(mag_model)

    # Because normalize_mag removes bad points internally, keep obs arrays clean here.
    good_obs_norm = np.isfinite(mag_obs)
    phase_obs_clean = phase_obs[good_obs_norm]

    # Recompute obs_norm using same mask length logic.
    obs_norm = (mag_obs - np.nanmean(mag_obs)) / (
        np.nanpercentile(mag_obs, 99) - np.nanpercentile(mag_obs, 1)
    )

    # Align by phase.
    best = best_phase_align(
        phase_obs=phase_obs_clean,
        obs_norm=obs_norm,
        phase_model=phase_model,
        model_norm=model_norm,
    )

    # Useful direct amplitude diagnostics.
    obs_amp_maxmin = np.nanmax(mag_obs) - np.nanmin(mag_obs)
    model_amp_maxmin = np.nanmax(mag_model) - np.nanmin(mag_model)

    print("FN Lyr / KIC 6936115 comparison")
    print("--------------------------------")
    print(f"Observed period:          {P_OBS:.9f} d")
    print(f"Model period:             {P_MODEL:.9f} d")
    print(f"Period difference:        {P_MODEL - P_OBS:+.9f} d")
    print(f"Fractional difference:    {(P_MODEL - P_OBS) / P_OBS:+.4%}")
    print()
    print(f"Observed amp P99-P01:     {obs_amp:.4f} mag")
    print(f"Observed amp max-min:     {obs_amp_maxmin:.4f} mag")
    print(f"Model amp P99-P01:        {model_amp:.4f} mag")
    print(f"Model amp max-min:        {model_amp_maxmin:.4f} mag")
    print()
    print(f"Best phase shift:         {best['shift']:.5f}")
    print(f"Pearson r, shape:         {best['corr']:.5f}")
    print(f"R2, shape:                {best['r2']:.5f}")
    print(f"RMS normalized residual:  {best['rms']:.5f}")

    # Plot normalized overlay.
    plt.figure(figsize=(7.2, 4.4))
    plt.plot(
        phase_obs_clean,
        obs_norm,
        ".",
        ms=4,
        label="FN Lyr Kepler, binned",
    )
    plt.plot(
        phase_obs_clean,
        best["model_at_obs"],
        "-",
        lw=2,
        label="MESA RSP + Colors",
    )

    plt.gca().invert_yaxis()
    plt.xlabel("Pulsation phase")
    plt.ylabel("Normalized relative magnitude")
    plt.title("FN Lyr / KIC 6936115 observational validation")
    plt.legend()
    plt.tight_layout()
    plt.savefig("FN_Lyr_MESA_vs_Kepler_normalized.png", dpi=250)

    # Plot residuals.
    plt.figure(figsize=(7.2, 3.2))
    plt.axhline(0.0, lw=1)
    plt.plot(
        phase_obs_clean,
        best["resid"],
        ".",
        ms=4,
    )
    plt.xlabel("Pulsation phase")
    plt.ylabel("Observed - model")
    plt.title("FN Lyr normalized residuals")
    plt.tight_layout()
    plt.savefig("FN_Lyr_MESA_vs_Kepler_residuals.png", dpi=250)

    # Also save comparison table for later paper plotting if desired.
    out = pd.DataFrame(
        {
            "phase": phase_obs_clean,
            "obs_norm_mag": obs_norm,
            "model_norm_mag": best["model_at_obs"],
            "residual_obs_minus_model": best["resid"],
        }
    )
    out.to_csv("FN_Lyr_MESA_vs_Kepler_comparison.csv", index=False)

    print()
    print("Wrote:")
    print("  FN_Lyr_MESA_vs_Kepler_normalized.png")
    print("  FN_Lyr_MESA_vs_Kepler_residuals.png")
    print("  FN_Lyr_MESA_vs_Kepler_comparison.csv")


if __name__ == "__main__":
    main()
