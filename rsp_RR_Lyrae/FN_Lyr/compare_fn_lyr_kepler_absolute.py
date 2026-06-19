#!/usr/bin/env python3
"""
Absolute-amplitude FN Lyr comparison.

This intentionally does NOT normalize amplitudes.

It fits only:
  - phase shift
  - vertical magnitude offset

It plots observed FN Lyr over two phases and overlays the best model cycle.
Optionally overlays the previous and next model cycles too.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


OBS_FILE = "FN_Lyr_KIC6936115_kepler_binned.csv"
MESA_FILE = "LOGS/history.data"

P_OBS = 0.527398471
OBS_AMP = 1.0464

# Search only late cycles by default. Increase if needed.
SEARCH_LAST_N_CYCLES = 80

# Plot adjacent cycles around the selected cycle.
PLOT_ADJACENT = False

# Avoid very incomplete cycles.
MIN_POINTS_PER_CYCLE = 250
MIN_PHASE_SPAN = 0.80


def read_mesa_history(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    lines = path.read_text().splitlines()
    header_idx = next(i for i, line in enumerate(lines) if "model_number" in line.split())
    return pd.read_csv(path, sep=r"\s+", skiprows=header_idx)


def periodic_interp(phase: np.ndarray, mag: np.ndarray, phase_new: np.ndarray) -> np.ndarray:
    phase = np.asarray(phase, dtype=float) % 1.0
    mag = np.asarray(mag, dtype=float)
    phase_new = np.asarray(phase_new, dtype=float) % 1.0

    good = np.isfinite(phase) & np.isfinite(mag)
    phase = phase[good]
    mag = mag[good]

    order = np.argsort(phase)
    phase = phase[order]
    mag = mag[order]

    # Remove duplicate phase points enough for np.interp.
    # This is simple and robust for dense MESA output.
    _, unique_idx = np.unique(np.round(phase, 8), return_index=True)
    unique_idx = np.sort(unique_idx)
    phase = phase[unique_idx]
    mag = mag[unique_idx]

    x = np.concatenate([phase - 1.0, phase, phase + 1.0])
    y = np.concatenate([mag, mag, mag])

    return np.interp(phase_new, x, y)


def get_mag_col(h: pd.DataFrame) -> str:
    if "Kepler" in h.columns:
        return "Kepler"
    if "K" in h.columns:
        return "K"
    candidates = [c for c in h.columns if "kepler" in c.lower() or c.lower() in {"kp", "kplr"}]
    if candidates:
        return candidates[0]
    raise RuntimeError(f"No Kepler magnitude column found. Columns: {list(h.columns)}")


def get_obs() -> pd.DataFrame:
    obs = pd.read_csv(OBS_FILE)
    if not {"phase", "mag"}.issubset(obs.columns):
        raise RuntimeError(f"{OBS_FILE} needs phase and mag columns")
    obs = obs[np.isfinite(obs["phase"]) & np.isfinite(obs["mag"])].copy()
    obs = obs.sort_values("phase")
    return obs


def cycle_table(h: pd.DataFrame, mag_col: str) -> pd.DataFrame:
    rows = []
    nmax = int(np.nanmax(h["rsp_num_periods"]))

    for n in range(max(0, nmax - SEARCH_LAST_N_CYCLES), nmax + 1):
        cyc = h[h["rsp_num_periods"].astype(int) == n].copy()
        if len(cyc) < MIN_POINTS_PER_CYCLE:
            continue

        phase = np.asarray(cyc["rsp_phase"], dtype=float) % 1.0
        mag = np.asarray(cyc[mag_col], dtype=float)

        good = np.isfinite(phase) & np.isfinite(mag)
        phase = phase[good]
        mag = mag[good]

        if len(phase) < MIN_POINTS_PER_CYCLE:
            continue

        phase_span = np.nanmax(phase) - np.nanmin(phase)
        if phase_span < MIN_PHASE_SPAN:
            continue

        amp = np.nanmax(mag) - np.nanmin(mag)
        p_est = np.nanmax(cyc["star_age_day"]) - np.nanmin(cyc["star_age_day"])
        vmax = np.nanmax(np.abs(cyc["max_abs_v_div_cs"])) if "max_abs_v_div_cs" in cyc else np.nan

        rows.append(
            {
                "cycle": n,
                "npts": len(phase),
                "phase_span": phase_span,
                "amp": amp,
                "period_est": p_est,
                "vmax_cs": vmax,
            }
        )

    if not rows:
        raise RuntimeError("No complete cycles found")

    return pd.DataFrame(rows)


def get_cycle(h: pd.DataFrame, mag_col: str, cycle_num: int) -> tuple[np.ndarray, np.ndarray]:
    cyc = h[h["rsp_num_periods"].astype(int) == cycle_num].copy()
    phase = np.asarray(cyc["rsp_phase"], dtype=float) % 1.0
    mag = np.asarray(cyc[mag_col], dtype=float)

    good = np.isfinite(phase) & np.isfinite(mag)
    phase = phase[good]
    mag = mag[good]

    order = np.argsort(phase)
    return phase[order], mag[order]


def score_cycle(
    obs_phase: np.ndarray,
    obs_mag: np.ndarray,
    model_phase: np.ndarray,
    model_mag: np.ndarray,
) -> dict:
    """
    Search phase shift. Fit vertical offset analytically for each shift.

    We do not rescale amplitude.
    """
    shifts = np.linspace(0, 1, 2500, endpoint=False)

    best = None

    for shift in shifts:
        # model evaluated at observed phase after phase shift
        m_eval = periodic_interp(model_phase, model_mag, (obs_phase - shift) % 1.0)

        # best vertical offset: obs ~= model + offset
        offset = np.nanmedian(obs_mag - m_eval)
        pred = m_eval + offset
        resid = obs_mag - pred

        rms = float(np.sqrt(np.nanmean(resid**2)))
        mae = float(np.nanmedian(np.abs(resid)))

        # Shape metric after mean subtraction, but no amplitude rescaling.
        obs_dm = obs_mag - np.nanmean(obs_mag)
        pred_dm = pred - np.nanmean(pred)
        corr = float(np.corrcoef(obs_dm, pred_dm)[0, 1])

        amp_model = float(np.nanmax(model_mag) - np.nanmin(model_mag))
        amp_penalty = abs(amp_model - OBS_AMP) / OBS_AMP

        # Main score: residuals plus amplitude mismatch.
        score = rms + 0.25 * amp_penalty

        if best is None or score < best["score"]:
            best = {
                "score": score,
                "shift": float(shift),
                "offset": float(offset),
                "rms": rms,
                "mae": mae,
                "corr": corr,
                "amp_model": amp_model,
                "pred": pred,
                "resid": resid,
            }

    return best


def main() -> None:
    obs = get_obs()
    h = read_mesa_history(MESA_FILE)
    mag_col = get_mag_col(h)

    obs_phase = obs["phase"].to_numpy(dtype=float)
    obs_mag = obs["mag"].to_numpy(dtype=float)

    obs_amp = float(np.nanmax(obs_mag) - np.nanmin(obs_mag))

    ct = cycle_table(h, mag_col)

    results = []
    for _, row in ct.iterrows():
        cycle = int(row["cycle"])
        model_phase, model_mag = get_cycle(h, mag_col, cycle)
        best = score_cycle(obs_phase, obs_mag, model_phase, model_mag)

        results.append(
            {
                "cycle": cycle,
                "npts": int(row["npts"]),
                "period_est": float(row["period_est"]),
                "phase_span": float(row["phase_span"]),
                "vmax_cs": float(row["vmax_cs"]),
                **best,
            }
        )

    res = pd.DataFrame(results)
    res = res.sort_values("score")
    res.to_csv("FN_Lyr_cycle_scores_absolute.csv", index=False)

    best = res.iloc[0].to_dict()
    best_cycle = int(best["cycle"])

    print("Best absolute-amplitude cycle")
    print("-----------------------------")
    print(f"MESA mag column:      {mag_col}")
    print(f"Observed amp:         {obs_amp:.4f} mag")
    print(f"Best cycle:           {best_cycle}")
    print(f"Model amp:            {best['amp_model']:.4f} mag")
    print(f"Period estimate:      {best['period_est']:.5f} d")
    print(f"Best phase shift:     {best['shift']:.5f}")
    print(f"Best vertical offset: {best['offset']:.5f} mag")
    print(f"RMS residual:         {best['rms']:.5f} mag")
    print(f"Median abs residual:  {best['mae']:.5f} mag")
    print(f"Pearson r:            {best['corr']:.5f}")
    print()
    print("Top 10 cycles:")
    print(
        res[
            ["cycle", "period_est", "amp_model", "rms", "mae", "corr", "vmax_cs", "score"]
        ].head(10).to_string(index=False)
    )

    # Prepare two-phase observed curve.
    obs_phase_2 = np.concatenate([obs_phase, obs_phase + 1.0])
    obs_mag_2 = np.concatenate([obs_mag, obs_mag])

    plt.figure(figsize=(7.2, 4.2))

    plt.plot(
        obs_phase_2,
        obs_mag_2,
        ".",
        ms=3.5,
        label="FN Lyr Kepler, binned",
    )

    # Plot best cycle, and optionally previous/next cycles.
    cycles_to_plot = [best_cycle]
    if PLOT_ADJACENT:
        cycles_to_plot = [best_cycle - 1, best_cycle, best_cycle + 1]

    for cycnum in cycles_to_plot:
        if cycnum not in set(ct["cycle"].astype(int)):
            continue

        model_phase, model_mag = get_cycle(h, mag_col, cycnum)

        # Refit offset for this specific cycle using best cycle's phase shift for visual consistency.
        model_at_obs = periodic_interp(model_phase, model_mag, (obs_phase - best["shift"]) % 1.0)
        offset = float(np.nanmedian(obs_mag - model_at_obs))

        model_phase_shifted = (model_phase + best["shift"]) % 1.0
        order = np.argsort(model_phase_shifted)

        phase_one = model_phase_shifted[order]
        mag_one = model_mag[order] + offset

        # Plot over two phases.
        phase_two = np.concatenate([phase_one, phase_one + 1.0])
        mag_two = np.concatenate([mag_one, mag_one])

        if cycnum == best_cycle:
            lw = 2.2
            alpha = 1.0
            label = f"MESA cycle {cycnum}"
        else:
            lw = 1.2
            alpha = 0.45
            label = f"MESA cycle {cycnum}"

        plt.plot(phase_two, mag_two, "-", lw=lw, alpha=alpha, label=label)

    plt.xlim(0, 2)
    plt.gca().invert_yaxis()
    plt.xlabel("Pulsation phase")
    plt.ylabel("Relative Kepler magnitude")
    plt.title("FN Lyr / KIC 6936115 absolute-amplitude comparison")
    plt.legend(loc="best", fontsize=9)
    plt.tight_layout()
    plt.savefig("FN_Lyr_MESA_vs_Kepler_absolute_2phase.png", dpi=250)

    # Residuals over one phase for best cycle.
    model_phase, model_mag = get_cycle(h, mag_col, best_cycle)
    model_at_obs = periodic_interp(model_phase, model_mag, (obs_phase - best["shift"]) % 1.0)
    pred = model_at_obs + best["offset"]
    resid = obs_mag - pred

    plt.figure(figsize=(7.2, 2.8))
    plt.axhline(0, lw=1)
    plt.plot(obs_phase, resid, ".", ms=3.5)
    plt.xlabel("Pulsation phase")
    plt.ylabel("Obs. - model (mag)")
    plt.title(f"Residuals for MESA cycle {best_cycle}")
    plt.tight_layout()
    plt.savefig("FN_Lyr_MESA_vs_Kepler_absolute_residuals.png", dpi=250)

    out = pd.DataFrame(
        {
            "phase": obs_phase,
            "obs_mag": obs_mag,
            "model_mag_offset": pred,
            "residual_obs_minus_model": resid,
        }
    )
    out.to_csv("FN_Lyr_MESA_vs_Kepler_absolute_comparison.csv", index=False)

    print()
    print("Wrote:")
    print("  FN_Lyr_cycle_scores_absolute.csv")
    print("  FN_Lyr_MESA_vs_Kepler_absolute_2phase.png")
    print("  FN_Lyr_MESA_vs_Kepler_absolute_residuals.png")
    print("  FN_Lyr_MESA_vs_Kepler_absolute_comparison.csv")


if __name__ == "__main__":
    main()
