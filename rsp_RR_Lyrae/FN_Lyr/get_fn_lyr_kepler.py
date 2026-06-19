# get_fn_lyr_kepler.py

import numpy as np
import pandas as pd
import lightkurve as lk

target = "KIC 6936115"   # FN Lyr
period = 0.527398471     # days, Nemec et al. 2011

# Start with long cadence because it is lighter and enough for a validation plot.
search = lk.search_lightcurve(target, mission="Kepler", cadence="long")
print(search)

lcs = search.download_all()
lc = lcs.stitch().remove_nans().remove_outliers(sigma=8)

t = np.asarray(lc.time.value, dtype=float)
flux = np.asarray(lc.flux.value, dtype=float)

good = np.isfinite(t) & np.isfinite(flux) & (flux > 0)
t = t[good]
flux = flux[good]

# Relative Kepler magnitude. This is enough for shape/amplitude comparison.
mag = -2.5 * np.log10(flux / np.nanmedian(flux))

# Fold. Absolute epoch does not matter because you will align phase later.
phase = ((t - t[0]) / period) % 1.0

df = pd.DataFrame({
    "time": t,
    "phase": phase,
    "rel_mag": mag,
    "flux": flux,
})

df = df.sort_values("phase")
df.to_csv("FN_Lyr_KIC6936115_kepler_lc.csv", index=False)

print(df.head())
print("N =", len(df))
print("Observed relative amplitude =", df["rel_mag"].max() - df["rel_mag"].min())
