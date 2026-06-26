#!/usr/bin/env python3
"""
Regenerate per-model SEDs for the settled overnight RR Lyrae run by
post-processing its history through SED_Model's forward model.

This uses the SAME flux cube and interpolation as the MESA Colors module
(sed_model.grid reads the identical flux_cube.bin; sed_model.forward calls
the same cc_api Fortran interpolation), so the SEDs produced here match what
Colors would have written at runtime had sed_per_model been enabled.

For each history row it calls run_forward(Teff, logg, [M/H], R, d) and writes
the interpolated SED to a per-model CSV in the format the plot scripts expect.

VERIFY BEFORE TRUSTING (see CONFIG block):
  - grid / filter / vega paths match the inlist
  - mag_system and interp_method match the Colors run
  - the SED output filename + columns match what plot_rsp_roman expects
  - the consistency check at the end passes (post-processed F-mags should
    match the F062..F213 already in the history within ~1e-3 mag)
"""

from __future__ import annotations
import numpy as np, pandas as pd, subprocess, io as _io
from pathlib import Path

from sed_model import load_grid, load_filters_from_instrument_dir, run_forward

# ===========================================================================
# CONFIG -- MATCH THESE TO YOUR inlist_rsp_RR_Lyrae COLORS BLOCK
# ===========================================================================
HISTORY   = "LOGS_overnight/history_last100.data"
# These three are RELATIVE in your inlist (data/colors_data/...). Use the same
# absolute paths the MESA run resolved them to. Adjust if your data dir differs.
DATA_ROOT = Path("/home/njm/MESA/mesa/data/colors_data/")
GRID_DIR  = DATA_ROOT / "stellar_models/Kurucz2003all"
FILTER_DIR= DATA_ROOT / "filters/Roman/WFI"
VEGA_SED  = DATA_ROOT / "stellar_models/vega_flam.csv"

MAG_SYSTEM    = "AB"        # inlist: mag_system = 'AB'
INTERP_METHOD = "hermite"   # <-- confirm vs colors.defaults; 'linear' if Colors used linear
# Grid axis is [M/H] in DEX (-2.5..+0.5), NOT Z mass fraction. The model is
# metal-poor: Z=0.0004 with a09 (Zsun~0.0134) -> [M/H] = log10(0.0004/0.0134) ~ -1.53.
# (Confirmed indirectly: runtime Colors Interp_rad ~0.06 is small, i.e. it sat near
# a grid node in the metal-poor region, not at solar.)
META          = -1.53
DISTANCE_CM   = 3.0857e19   # inlist: distance = 3.0857d19
RSUN_CM       = 6.957e10
# Only the Roman WFI science filters get SED files (skip Grism/Prism unless needed)
WANT_FILTERS  = ["F062","F087","F106","F129","F146","F158","F184","F213"]

OUTDIR = Path("SED")        # where per-model SEDs go (matches colors_results_directory)
OUTDIR.mkdir(exist_ok=True)

# ===========================================================================
# Read the (already-trimmed) settled history
# ===========================================================================
def read_history(path):
    with open(path) as f:
        for i, l in enumerate(f):
            if "model_number" in l.split():
                hdr = i; break
    return pd.read_csv(path, sep=r"\s+", skiprows=hdr)

h = read_history(HISTORY)
print(f"rows: {len(h)}  cycles: {int(h['rsp_num_periods'].min())}-{int(h['rsp_num_periods'].max())}")

# Trim to the last N COMPLETE settled cycles (skip the degenerate final cycle).
# The SED figures only need ~1 settled cycle (hot/cold/maxR/minR); a few gives margin.
N_CYCLES = 5
nmax = int(h['rsp_num_periods'].max())
# skip the final cycle if it's degenerate (very few points)
if (h['rsp_num_periods'].astype(int) == nmax).sum() < 200:
    nmax -= 1
h = h[h['rsp_num_periods'].astype(int).between(nmax - N_CYCLES + 1, nmax)].copy()
print(f"trimmed to cycles {nmax-N_CYCLES+1}-{nmax}: {len(h)} rows")

# radius in cm: prefer photosphere_r (Rsun) if present, else 10**log_R
if "photosphere_r" in h.columns:
    R_cm = h["photosphere_r"].to_numpy(float) * RSUN_CM
elif "log_R" in h.columns:
    R_cm = (10.0 ** h["log_R"].to_numpy(float)) * RSUN_CM
else:
    raise RuntimeError("no radius column (photosphere_r or log_R)")

teff = h["effective_T"].to_numpy(float) if "effective_T" in h.columns \
       else 10.0**h["log_Teff"].to_numpy(float)
logg = h["log_g"].to_numpy(float)
mnum = h["model_number"].to_numpy(int)

# ===========================================================================
# Load grid + filters ONCE (this is the expensive part)
# ===========================================================================
print("loading grid + filters...")
grid    = load_grid(str(GRID_DIR))
filters = load_filters_from_instrument_dir(str(FILTER_DIR),
            vega_sed_path=str(VEGA_SED) if Path(VEGA_SED).exists() else None)
print(f"  {grid}")
print(f"  filters: {[f.name for f in filters]}")

# ===========================================================================
# Forward-model every row, write per-model SED, collect mags for the check
# ===========================================================================
check_rows = []
for k in range(len(h)):
    res = run_forward(teff=float(teff[k]), logg=float(logg[k]), meta=float(META),
                      R=float(R_cm[k]), d=DISTANCE_CM,
                      grid=grid, filters=filters,
                      mag_system=MAG_SYSTEM, interp_method=INTERP_METHOD)
    # ---- write ONE file PER FILTER PER MODEL, matching the real Colors format:
    #   columns: wavelengths, fluxes, convolved_flux, filter_wavelengths, filter_trans
    #   filename: {FILTER}_SED_{model:08d}.csv   (e.g. F062_SED_00002441.csv)
    # convolved_flux = observed_flux * (filter transmission interpolated onto SED grid)
    sed_wl  = res.wavelengths
    obs_flx = res.observed_flux
    for filt in filters:
        if filt.name not in WANT_FILTERS:
            continue
        trans_on_sed = np.interp(sed_wl, filt.wavelengths, filt.transmission,
                                 left=0.0, right=0.0)
        convolved = obs_flx * trans_on_sed
        # Match MESA Colors synthetic.f90 CSV convention:
        # max(size(wavelengths), size(filter_wavelengths)) rows, zero-padded,
        # with the fluxes column containing the diluted observer-frame flux.
        nrow = max(len(sed_wl), len(filt.wavelengths))
        wl = np.zeros(nrow); fl = np.zeros(nrow); cf = np.zeros(nrow)
        fw = np.zeros(nrow); ft = np.zeros(nrow)
        wl[:len(sed_wl)] = sed_wl
        fl[:len(obs_flx)] = obs_flx
        cf[:len(convolved)] = convolved
        fw[:len(filt.wavelengths)] = filt.wavelengths
        ft[:len(filt.transmission)] = filt.transmission
        out = OUTDIR / f"{filt.name}_SED_{mnum[k]:08d}.csv"
        pd.DataFrame({
            "wavelengths":        wl,
            "fluxes":             fl,
            "convolved_flux":     cf,
            "filter_wavelengths": fw,
            "filter_trans":       ft,
        }).to_csv(out, index=False)
    # collect synthetic mags to compare to history's F062..F213
    row = {"model_number": mnum[k]}
    row.update({fn: m for fn, m in res.magnitudes.items()})
    check_rows.append(row)
    if k % 500 == 0:
        print(f"  {k}/{len(h)}  model={mnum[k]}  Teff={teff[k]:.0f}")

chk = pd.DataFrame(check_rows).set_index("model_number")
print(f"wrote {len(h)} SED files to {OUTDIR}/")

# ===========================================================================
# CONSISTENCY CHECK -- the whole point.
# The history ALREADY has F062..F213 from the runtime Colors module.
# Post-processed mags must match within ~1e-3 mag if the machinery is identical.
# ===========================================================================
print("\n=== consistency check: post-processed vs runtime Colors mags ===")
hi = h.set_index("model_number")
for fn in [f.name for f in filters]:
    if fn in hi.columns:
        d = (chk[fn] - hi.loc[chk.index, fn]).abs()
        print(f"  {fn}: median|Δmag|={d.median():.5f}  max={d.max():.5f}")
    else:
        print(f"  {fn}: not in history columns (skip)")
print("If these Δmag are ~<1e-3, SED_Model reproduces Colors. If large, the "
      "config (META/[M/H], interp_method, grid) does not match the MESA run.")
