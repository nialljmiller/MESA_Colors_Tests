from pathlib import Path
import numpy as np
import pandas as pd

MESA = Path("/home/njm/MESA/mesa")
FILTER_DIR = MESA / "data/colors_data/filters/GAIA/GAIA"
VEGA = MESA / "data/colors_data/stellar_models/vega_flam.csv"


def clean_num_array(s):
    s = s.astype(str).str.strip()
    s = s.str.replace(
        r"^([+-]?(?:\d+\.\d*|\d+|\.\d+))([+-]\d{2,4})$",
        r"\1e\2",
        regex=True,
    )
    return pd.to_numeric(s, errors="coerce").to_numpy()


def trapz(x, y):
    good = np.isfinite(x) & np.isfinite(y)
    x = np.asarray(x[good], dtype=float)
    y = np.asarray(y[good], dtype=float)
    order = np.argsort(x)
    return np.trapezoid(y[order], x[order])


def read_history(path):
    lines = Path(path).read_text().splitlines()

    header_idx = None
    for i, line in enumerate(lines):
        if "model_number" in line.split():
            header_idx = i

    if header_idx is None:
        raise RuntimeError(f"Could not find model_number header in {path}")

    cols = lines[header_idx].split()
    rows = []

    for line in lines[header_idx + 1:]:
        parts = line.split()
        if len(parts) != len(cols):
            continue
        try:
            rows.append([float(x) for x in parts])
        except ValueError:
            pass

    return pd.DataFrame(rows, columns=cols)


def read_filter(name):
    df = pd.read_csv(FILTER_DIR / f"{name}.dat")
    wave = clean_num_array(df.iloc[:, 0])
    trans = clean_num_array(df.iloc[:, 1])
    good = np.isfinite(wave) & np.isfinite(trans) & (wave > 0)
    return wave[good], trans[good]


def read_vega():
    df = pd.read_csv(VEGA)
    return clean_num_array(df["WAVELENGTH"]), clean_num_array(df["FLUX"])


def read_sed(name, model):
    path = Path(f"SED/{name}_SED_{model:08d}.csv")
    df = pd.read_csv(path)
    wave = clean_num_array(df["wavelengths"])
    flux = clean_num_array(df["fluxes"])
    return wave, flux


def vega_mag(name, model):
    filt_wave, filt_trans = read_filter(name)
    sed_wave, sed_flux = read_sed(name, model)
    vega_wave, vega_flux = read_vega()

    # Important: hard-zero filter outside passband.
    filt_on_sed = np.interp(
        sed_wave,
        filt_wave,
        filt_trans,
        left=0.0,
        right=0.0,
    )

    vega_on_filter = np.interp(
        filt_wave,
        vega_wave,
        vega_flux,
        left=0.0,
        right=0.0,
    )

    syn = trapz(sed_wave, sed_flux * filt_on_sed * sed_wave) / trapz(
        sed_wave, filt_on_sed * sed_wave
    )

    zp = trapz(filt_wave, vega_on_filter * filt_trans * filt_wave) / trapz(
        filt_wave, filt_trans * filt_wave
    )

    return -2.5 * np.log10(syn / zp)


hist = read_history("LOGS_wd_cool/history.data")

rows = []
for _, row in hist.iterrows():
    model = int(row["model_number"])

    try:
        G = vega_mag("G", model)
        Gbp = vega_mag("Gbp", model)
        Grp = vega_mag("Grp", model)
    except FileNotFoundError:
        continue

    rows.append(
        {
            "model_number": model,
            "star_age": row.get("star_age", np.nan),
            "Teff": row.get("Teff", np.nan),
            "log_Teff": row.get("log_Teff", np.nan),
            "log_g": row.get("log_g", np.nan),
            "log_L": row.get("log_L", np.nan),
            "G": G,
            "Gbp": Gbp,
            "Grp": Grp,
            "bp_rp": Gbp - Grp,
            "M_G": G,
        }
    )

out = pd.DataFrame(rows)

if len(out) == 0:
    raise RuntimeError("No per-model SED files found. Check SED/G_SED_00000001.csv exists.")

out.to_parquet("mesa_gaia_track_corrected.parquet", index=False)
out.to_csv("mesa_gaia_track_corrected.csv", index=False)

print(out[["model_number", "Teff", "G", "Gbp", "Grp", "bp_rp"]].head())
print(out[["model_number", "Teff", "G", "Gbp", "Grp", "bp_rp"]].tail())
print("wrote mesa_gaia_track_corrected.parquet")
print("wrote mesa_gaia_track_corrected.csv")
