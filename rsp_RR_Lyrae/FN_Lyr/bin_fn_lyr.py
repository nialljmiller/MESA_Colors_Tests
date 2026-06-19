# bin_fn_lyr.py

import numpy as np
import pandas as pd

df = pd.read_csv("FN_Lyr_KIC6936115_kepler_lc.csv")

nbin = 200
bins = np.linspace(0, 1, nbin + 1)
which = np.digitize(df["phase"], bins) - 1

rows = []
for k in range(nbin):
    x = df.loc[which == k, "rel_mag"].to_numpy()
    if len(x) < 3:
        continue
    rows.append({
        "phase": 0.5 * (bins[k] + bins[k+1]),
        "mag": np.nanmedian(x),
        "mag_std": np.nanstd(x),
        "n": len(x),
    })

out = pd.DataFrame(rows)
out.to_csv("FN_Lyr_KIC6936115_kepler_binned.csv", index=False)
print(out.head())
