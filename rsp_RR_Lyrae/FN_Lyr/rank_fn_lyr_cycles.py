#!/usr/bin/env python3
# rank_fn_lyr_cycles.py

from pathlib import Path
import sys
import numpy as np
import pandas as pd

OBS_AMP = 1.0508
OBS_PERIOD = 0.527398471

def read_history(path):
    path = Path(path)
    lines = path.read_text(errors="ignore").splitlines()
    header_idx = next(i for i, l in enumerate(lines) if "model_number" in l.split())
    return pd.read_csv(path, sep=r"\s+", skiprows=header_idx, on_bad_lines="skip")

def main():
    if len(sys.argv) < 2:
        raise SystemExit("usage: python rank_fn_lyr_cycles.py path/to/history.data")

    h = read_history(sys.argv[1])
    mag_col = "Kepler" if "Kepler" in h.columns else "K"

    rows = []
    nmax = int(np.nanmax(h["rsp_num_periods"]))

    for n in range(max(0, nmax - 80), nmax + 1):
        m = h[h["rsp_num_periods"].astype(int) == n]
        if len(m) < 300:
            continue

        phase_span = np.nanmax(m["rsp_phase"]) - np.nanmin(m["rsp_phase"])
        if phase_span < 0.8:
            continue

        amp = np.nanmax(m[mag_col]) - np.nanmin(m[mag_col])
        p_est = np.nanmax(m["star_age_day"]) - np.nanmin(m["star_age_day"])
        dteff = np.nanmax(m["effective_T"]) - np.nanmin(m["effective_T"])
        vmax = np.nanmax(np.abs(m["max_abs_v_div_cs"])) if "max_abs_v_div_cs" in m else np.nan

        # crude wiggle metric: total variation relative to amplitude
        phase = np.asarray(m["rsp_phase"], dtype=float)
        mag = np.asarray(m[mag_col], dtype=float)
        order = np.argsort(phase)
        mag_s = mag[order]
        tv = np.nansum(np.abs(np.diff(mag_s)))
        wiggle = tv / amp if amp > 0 else np.inf

        amp_err = abs(amp - OBS_AMP) / OBS_AMP
        per_err = abs(p_est - OBS_PERIOD) / OBS_PERIOD

        # This prefers high amplitude but penalizes obviously violent/shocky cycles.
        shock_penalty = max(vmax - 4.0, 0.0) / 4.0
        wiggle_penalty = max(wiggle - 2.5, 0.0)

        score = 2.0 * amp_err + 1.0 * per_err + 0.5 * shock_penalty + 0.25 * wiggle_penalty

        rows.append({
            "cycle": n,
            "period": p_est,
            "amp": amp,
            "amp_frac": amp / OBS_AMP,
            "dteff": dteff,
            "vmax_cs": vmax,
            "wiggle": wiggle,
            "score": score,
        })

    out = pd.DataFrame(rows).sort_values("score")
    print(out.head(20).to_string(index=False))
    out.to_csv("ranked_cycles.csv", index=False)

if __name__ == "__main__":
    main()
