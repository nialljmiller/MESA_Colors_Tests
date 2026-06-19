#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd

OBS_FILE = Path("FN_Lyr_KIC6936115_kepler_binned.csv")
OBS_PERIOD = 0.527398471
SEARCH_LAST_N_CYCLES = 80
MIN_POINTS = 250
MIN_PHASE_SPAN = 0.80

HISTORY_GLOBS = [
    "fn_lyr_grid/T*/history.data",
    "fn_lyr_grid_second_pass/T*/history.data",
    "fn_lyr_grid_first_pass/T*/history.data",
]


def read_history(path):
    path = Path(path)
    lines = path.read_text(errors="ignore").splitlines()
    header_idx = next(i for i, l in enumerate(lines) if "model_number" in l.split())
    return pd.read_csv(path, sep=r"\s+", skiprows=header_idx, on_bad_lines="skip")


def periodic_interp(phase, mag, phase_new):
    phase = np.asarray(phase, dtype=float) % 1.0
    mag = np.asarray(mag, dtype=float)
    phase_new = np.asarray(phase_new, dtype=float) % 1.0

    good = np.isfinite(phase) & np.isfinite(mag)
    phase = phase[good]
    mag = mag[good]

    order = np.argsort(phase)
    phase = phase[order]
    mag = mag[order]

    # de-duplicate dense phases
    _, idx = np.unique(np.round(phase, 8), return_index=True)
    idx = np.sort(idx)
    phase = phase[idx]
    mag = mag[idx]

    x = np.concatenate([phase - 1.0, phase, phase + 1.0])
    y = np.concatenate([mag, mag, mag])

    return np.interp(phase_new, x, y)


def score_cycle(obs_phase, obs_mag, model_phase, model_mag):
    shifts = np.linspace(0.0, 1.0, 2000, endpoint=False)

    best = None

    for shift in shifts:
        m_eval = periodic_interp(model_phase, model_mag, (obs_phase - shift) % 1.0)

        # Fit only vertical offset. No amplitude scaling.
        offset = np.nanmedian(obs_mag - m_eval)
        pred = m_eval + offset
        resid = obs_mag - pred

        rms = float(np.sqrt(np.nanmean(resid**2)))
        mad = float(np.nanmedian(np.abs(resid)))

        obs_dm = obs_mag - np.nanmean(obs_mag)
        pred_dm = pred - np.nanmean(pred)
        pearson = float(np.corrcoef(obs_dm, pred_dm)[0, 1])

        if best is None or rms < best["rms"]:
            best = {
                "phase_shift": float(shift),
                "mag_offset": float(offset),
                "rms": rms,
                "mad": mad,
                "pearson": pearson,
            }

    return best


def cycle_metrics(h, mag_col, cycle):
    m = h[h["rsp_num_periods"].astype(int) == cycle].copy()

    phase = np.asarray(m["rsp_phase"], dtype=float) % 1.0
    mag = np.asarray(m[mag_col], dtype=float)

    good = np.isfinite(phase) & np.isfinite(mag)
    phase = phase[good]
    mag = mag[good]

    if len(phase) < MIN_POINTS:
        return None

    phase_span = np.nanmax(phase) - np.nanmin(phase)
    if phase_span < MIN_PHASE_SPAN:
        return None

    order = np.argsort(phase)
    phase = phase[order]
    mag = mag[order]

    amp = float(np.nanmax(mag) - np.nanmin(mag))
    period = float(np.nanmax(m["star_age_day"]) - np.nanmin(m["star_age_day"]))

    vmax = np.nan
    if "max_abs_v_div_cs" in m.columns:
        vmax = float(np.nanmax(np.abs(m["max_abs_v_div_cs"])))

    dteff = np.nan
    if "effective_T" in m.columns:
        dteff = float(np.nanmax(m["effective_T"]) - np.nanmin(m["effective_T"]))

    # crude morphology roughness metric
    total_variation = float(np.nansum(np.abs(np.diff(mag))))
    wiggle = total_variation / amp if amp > 0 else np.inf

    return {
        "phase": phase,
        "mag": mag,
        "amp": amp,
        "period": period,
        "vmax_cs": vmax,
        "dteff": dteff,
        "wiggle": wiggle,
        "npts": len(phase),
        "phase_span": phase_span,
    }


def convergence_metrics(h, mag_col):
    if "rsp_num_periods" not in h.columns:
        return {}

    nmax = int(np.nanmax(h["rsp_num_periods"]))
    rows = []

    for n in range(max(0, nmax - 25), nmax + 1):
        cm = cycle_metrics(h, mag_col, n)
        if cm is None:
            continue
        rows.append((n, cm["period"], cm["amp"], cm["vmax_cs"]))

    if len(rows) < 5:
        return {}

    arr = np.array(rows, dtype=float)
    cycles = arr[:, 0]
    periods = arr[:, 1]
    amps = arr[:, 2]
    vmax = arr[:, 3]

    last = min(10, len(amps))
    x = np.arange(last, dtype=float)
    amp_slope = float(np.polyfit(x, amps[-last:], 1)[0])

    out = {
        "conv_ncycles": int(len(rows)),
        "amp_last10_median": float(np.nanmedian(amps[-last:])),
        "amp_last10_std": float(np.nanstd(amps[-last:])),
        "amp_last10_range": float(np.nanmax(amps[-last:]) - np.nanmin(amps[-last:])),
        "amp_slope_last10": amp_slope,
        "period_last10_median": float(np.nanmedian(periods[-last:])),
        "period_last10_std": float(np.nanstd(periods[-last:])),
        "vmax_last10_median": float(np.nanmedian(vmax[-last:])),
    }

    grekm_cols = [c for c in h.columns if "grekm" in c.lower()]
    if grekm_cols:
        c = grekm_cols[0]
        out["grekm_col"] = c
        out["grekm_last_median"] = float(np.nanmedian(h[c].tail(500)))
        out["grekm_last_max"] = float(np.nanmax(h[c].tail(500)))
    else:
        out["grekm_col"] = ""

    return out


def main():
    obs = pd.read_csv(OBS_FILE)
    obs = obs[np.isfinite(obs["phase"]) & np.isfinite(obs["mag"])].copy()
    obs = obs.sort_values("phase")

    obs_phase = obs["phase"].to_numpy(dtype=float)
    obs_mag = obs["mag"].to_numpy(dtype=float)
    obs_amp = float(np.nanmax(obs_mag) - np.nanmin(obs_mag))

    histories = []
    for pattern in HISTORY_GLOBS:
        histories.extend(Path(".").glob(pattern))

    histories = sorted(set(histories))

    print(f"Found {len(histories)} histories")

    rows = []

    for hist in histories:
        run = hist.parent.name
        grid = hist.parent.parent.name

        print(f"Scoring {grid}/{run}")

        h = read_history(hist)
        mag_col = "Kepler" if "Kepler" in h.columns else "K"

        nmax = int(np.nanmax(h["rsp_num_periods"]))
        best = None
        cycle_rows = []

        for n in range(max(0, nmax - SEARCH_LAST_N_CYCLES), nmax + 1):
            cm = cycle_metrics(h, mag_col, n)
            if cm is None:
                continue

            fit = score_cycle(obs_phase, obs_mag, cm["phase"], cm["mag"])

            amp_err = abs(cm["amp"] - obs_amp)
            period_err_pct = 100.0 * abs(cm["period"] - OBS_PERIOD) / OBS_PERIOD

            # strict paper-facing score:
            # amplitude and RMS dominate; period and shockiness are secondary.
            shock_penalty = max(cm["vmax_cs"] - 6.0, 0.0) * 0.015 if np.isfinite(cm["vmax_cs"]) else 0.0
            score = fit["rms"] + 0.75 * amp_err + 0.015 * period_err_pct + shock_penalty

            row = {
                "grid": grid,
                "run": run,
                "history": str(hist),
                "cycle": int(n),
                "period": cm["period"],
                "period_err_pct": period_err_pct,
                "amp": cm["amp"],
                "obs_amp": obs_amp,
                "amp_frac": cm["amp"] / obs_amp,
                "amp_err_mag": amp_err,
                "rms": fit["rms"],
                "mad": fit["mad"],
                "pearson": fit["pearson"],
                "phase_shift": fit["phase_shift"],
                "mag_offset": fit["mag_offset"],
                "vmax_cs": cm["vmax_cs"],
                "dteff": cm["dteff"],
                "wiggle": cm["wiggle"],
                "score": score,
            }

            cycle_rows.append(row)

            if best is None or score < best["score"]:
                best = row

        if best is None:
            continue

        best.update(convergence_metrics(h, mag_col))
        rows.append(best)

        # Save all cycle scores for this run
        if cycle_rows:
            outdir = Path("model_score_tables")
            outdir.mkdir(exist_ok=True)
            pd.DataFrame(cycle_rows).sort_values("score").to_csv(
                outdir / f"{grid}__{run}__cycle_scores.csv",
                index=False,
            )

    df = pd.DataFrame(rows).sort_values("score")
    df.to_csv("FN_Lyr_model_scores.csv", index=False)

    cols = [
        "grid", "run", "cycle",
        "period", "period_err_pct",
        "amp", "amp_frac",
        "rms", "mad", "pearson",
        "vmax_cs",
        "amp_last10_median", "amp_last10_range", "amp_slope_last10",
        "score",
    ]

    print()
    print(df[cols].head(25).to_string(index=False))
    print()
    print("Wrote FN_Lyr_model_scores.csv")
    print("Wrote model_score_tables/*__cycle_scores.csv")


if __name__ == "__main__":
    main()
