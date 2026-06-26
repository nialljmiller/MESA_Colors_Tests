"""
plot_yrec.py
============
Plotting script for the YREC SPOTS section of the MESA Colors paper.

Three figure functions:

  make_modulation_figure   — top/middle/bottom: disc cartoon, SEDs, light curves
  make_amplitude_grid      — amplitude vs f_spot, one filter, multiple masses
  make_amplitude_multifilter — amplitude vs f_spot, one mass, multiple filters

Paths — edit the block at the bottom before running.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Publication style
# ---------------------------------------------------------------------------

APJ_SINGLE = 3.5
APJ_DOUBLE = 7.0

def setup_style():
    plt.rcParams.update({
        'font.size': 9,
        'font.family': 'serif',
        'font.serif': ['Times', 'Times New Roman', 'DejaVu Serif'],
        'axes.labelsize': 10,
        'axes.titlesize': 9,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'legend.fontsize': 7,
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'text.usetex': False,
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'xtick.minor.width': 0.5,
        'ytick.minor.width': 0.5,
        'lines.linewidth': 1.0,
    })

# ---------------------------------------------------------------------------
# Filter metadata
# ---------------------------------------------------------------------------

JOHNSON = {
    'U': {'wl': 3650,  'color': '#7B2D8B'},
    'B': {'wl': 4400,  'color': '#3366CC'},
    'V': {'wl': 5500,  'color': '#339933'},
    'R': {'wl': 6400,  'color': '#CC3300'},
    'I': {'wl': 8000,  'color': '#990000'},
    'J': {'wl': 12500, 'color': '#663300'},
    'M': {'wl': 47000, 'color': '#333333'},
}

PLOT_FILTERS = ['U', 'B', 'V', 'R', 'I', 'J']

# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def read_colors_track(track_dir: Path, stem: str) -> pd.DataFrame:
    path = track_dir / stem / f"{stem}_colors.track"
    with open(path) as fh:
        for line in fh:
            if line.startswith("# ") and "StepNum" in line:
                cols = line.lstrip("# ").split()
                break
    df = pd.read_csv(path, sep=r"\s+", comment="#", header=None, names=cols)
    return df


def read_sed(sed_dir: Path, step: int, kind: str) -> tuple[np.ndarray, np.ndarray]:
    path = sed_dir / f"step_{step:05d}_{kind}.csv"
    data = np.genfromtxt(path, delimiter=",", comments="#")
    return data[:, 0], data[:, 2]   # wavelength_AA, observed_flux


def find_ms_step(df: pd.DataFrame, x_cen_frac: float = 0.5) -> int:
    x0     = df["X_cen"].iloc[0]
    target = x0 * x_cen_frac
    idx    = (df["X_cen"] - target).abs().idxmin()
    return int(df.loc[idx, "StepNum"])


def interpolate_to_age(df: pd.DataFrame, age_gyr: float) -> pd.Series | None:
    if age_gyr < df["Age_Gyr"].min() or age_gyr > df["Age_Gyr"].max():
        return None
    idx = (df["Age_Gyr"] - age_gyr).abs().idxmin()
    return df.loc[idx]

# ---------------------------------------------------------------------------
# Physics helpers
# ---------------------------------------------------------------------------

def f_visible(phase: np.ndarray, fspot: float) -> np.ndarray:
    """
    Projected spot covering fraction vs rotational phase.
    φ=0: spot at disc centre; φ=0.25: spot at limb; φ≥0.25: spot hidden.
    """
    return fspot * np.maximum(0.0, np.cos(2.0 * np.pi * phase))


def compute_delta_mag(m_hot: float, m_cool: float, fspot: float,
                      phase: np.ndarray) -> np.ndarray:
    """
    Magnitude offset relative to unspotted star at each rotation phase.
    Positive Δm = star is fainter.
    """
    dm   = m_cool - m_hot
    fvis = f_visible(phase, fspot)
    return -2.5 * np.log10((1.0 - fvis) + fvis * 10.0 ** (-dm / 2.5))


def net_sed(wl_hot, flux_hot, wl_cool, flux_cool, fvis: float) -> np.ndarray:
    if not np.allclose(wl_hot, wl_cool):
        flux_cool = np.interp(wl_hot, wl_cool, flux_cool)
    return (1.0 - fvis) * flux_hot + fvis * flux_cool

# ---------------------------------------------------------------------------
# Stellar disc cartoon
# ---------------------------------------------------------------------------

STAR_COLOR = '#F5A623'
SPOT_COLOR = '#8B4513'

def draw_star_disc(ax, fspot: float, phase: float):
    ax.set_aspect('equal')
    ax.set_xlim(-1.35, 1.35)
    ax.set_ylim(-1.35, 1.35)
    ax.axis('off')

    star_patch = plt.Circle((0, 0.5), 1.0, color=STAR_COLOR, zorder=1)
    ax.add_patch(star_patch)

    # Subtle limb darkening
    for r in np.linspace(0.95, 0.1, 20):
        ax.add_patch(plt.Circle((0, 0.5), r, color='white',
                                alpha=0.1 * (1 - r), zorder=2))

    # Spot geometry:
    #   φ=0   → spot at disc centre (facing observer), full circular projection
    #   φ=0.25 → spot at limb, foreshortened to zero width
    #   φ>0.25 → spot on far side, not drawn
    # x-position on disc: sin(2π φ)  [0 at centre, 1 at limb]
    # x-width (foreshortening): spot_r * cos(2π φ)
    theta     = 2.0 * np.pi * phase
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    spot_r    = np.sqrt(fspot) * 0.8

    if cos_theta > 0:                        # spot on visible hemisphere
        spot_x  = sin_theta                  # centre moves from 0 → 1 as spot goes to limb
        spot_rx = spot_r * cos_theta         # foreshortened width
        if spot_rx > 0.01:
            spot_patch = mpatches.Ellipse(
                (spot_x, 0.5), 2 * spot_rx, 2 * spot_r,
                color=SPOT_COLOR, zorder=3
            )
            ax.add_patch(spot_patch)
            spot_patch.set_clip_path(star_patch)   # keep spot inside disc

    ax.text(-1.25, 1.2, rf'$\phi={phase:.2f}$',
            ha='center', va='bottom', fontsize=7)

# ---------------------------------------------------------------------------
# Figure 1: modulation figure  (top / middle / bottom)
# ---------------------------------------------------------------------------

REP_PHASES = [0.0, 1/6, 0.25]   # spot centred → halfway to limb → at limb (f_vis=0)

def make_modulation_figure(
    stem: str,
    track_dir: Path,
    output_dir: Path,
    step: int | None = None,
    n_phases: int = 200,
    filters_to_show: list = PLOT_FILTERS,
    tag: str = '',
):
    sed_dir = track_dir / stem / "SEDs"
    df      = read_colors_track(track_dir, stem)

    if step is None:
        step = find_ms_step(df, x_cen_frac=0.5)

    row    = df[df["StepNum"] == step].iloc[0]
    fspot  = float(row["Fspot"])
    mass   = float(row["Mass"])
    age    = float(row["Age_Gyr"])
    t_hot  = float(row["T_hot"])
    t_cool = float(row["T_cool"])

    wl_hot,  flux_hot  = read_sed(sed_dir, step, "hot")
    wl_cool, flux_cool = read_sed(sed_dir, step, "cool")

    phase  = np.linspace(0, 1, n_phases)
    n_rep  = len(REP_PHASES)

    # Layout
    fig     = plt.figure(figsize=(APJ_DOUBLE, APJ_DOUBLE * 1.1))
    gs_main = gridspec.GridSpec(3, 1, figure=fig,
                                height_ratios=[1.0, 1.4, 1.4], hspace=0.00)
    gs_top  = gridspec.GridSpecFromSubplotSpec(1, n_rep, subplot_spec=gs_main[0], wspace=0.05)
    gs_mid  = gridspec.GridSpecFromSubplotSpec(1, n_rep, subplot_spec=gs_main[1], wspace=0.00)
    ax_discs = [fig.add_subplot(gs_top[i]) for i in range(n_rep)]
    ax_seds  = [fig.add_subplot(gs_mid[i]) for i in range(n_rep)]
    ax_lc    = fig.add_subplot(gs_main[2])

    # Top: disc cartoons
    for i, phi in enumerate(REP_PHASES):
        draw_star_disc(ax_discs[i], fspot, phi)

    # Middle: SEDs
    wl_um    = wl_hot / 1e4
    wl_range = (wl_um > 0.3) & (wl_um < 2.5)

    for i, phi in enumerate(REP_PHASES):
        ax   = ax_seds[i]
        fvis = float(f_visible(np.array([phi]), fspot)[0])
        f_net = net_sed(wl_hot, flux_hot, wl_cool, flux_cool, fvis)

        ax.semilogy(wl_um[wl_range], flux_hot[wl_range],
                    color='#E87722', lw=1.0, ls='-',  label=r'$T_{\rm hot}$', zorder=3)
        ax.semilogy(wl_um[wl_range], flux_cool[wl_range],
                    color=SPOT_COLOR, lw=1.0, ls='--', label=r'$T_{\rm cool}$', zorder=3)
        ax.semilogy(wl_um[wl_range], f_net[wl_range],
                    color='k', lw=1.2, ls='-', label='Net', zorder=4)

        for fname in filters_to_show:
            if fname not in JOHNSON:
                continue
            wl_f = JOHNSON[fname]['wl'] / 1e4
            if 0.3 < wl_f < 2.5:
                ax.axvline(wl_f, color=JOHNSON[fname]['color'], alpha=0.5, lw=0.7, ls=':')

            ax.set_xlim(0.3, 2.5)

            ax.set_xlabel(r'$\lambda\;(\mu{\rm m})$', fontsize=8)
            ax.xaxis.set_label_position('top')
            ax.xaxis.tick_top()
            ax.tick_params(axis='x', which='both',
                           top=True, labeltop=True,
                           bottom=False, labelbottom=False)

        if i == 0:
            ax.set_ylabel(r'$F_\lambda$ (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)', fontsize=8)
            ax.legend(fontsize=6, loc='lower right')
        else:
            ax.set_yticklabels([])

        ax.text(0.97, 0.97, rf'$\phi={phi:.3f}$',
                transform=ax.transAxes, ha='right', va='top', fontsize=7)

    # Bottom: light curves
    for fname in filters_to_show:
        hot_col  = f"{fname}_hot"
        cool_col = f"{fname}_cool"
        if hot_col not in row.index or cool_col not in row.index:
            continue
        m_h = float(row[hot_col])
        m_c = float(row[cool_col])
        if not (np.isfinite(m_h) and np.isfinite(m_c)):
            continue
        ax_lc.plot(phase, compute_delta_mag(m_h, m_c, fspot, phase),
                   color=JOHNSON[fname]['color'], lw=1.1, label=fname)

    for phi in REP_PHASES:
        ax_lc.axvline(phi, color='grey', ls='--', lw=0.7, alpha=0.7)

    ax_lc.set_xlabel('Rotational phase', fontsize=10)
    ax_lc.set_ylabel(r'$\Delta m$ (mag, spotted $-$ unspotted)', fontsize=10)
    ax_lc.invert_yaxis()
    ax_lc.set_xlim(0, 1)
    ax_lc.legend(ncol=len(filters_to_show), fontsize=7,
                 loc='lower center', framealpha=0.8)
    ax_lc.xaxis.set_minor_locator(plt.MultipleLocator(0.05))

    fname_out = f"fig_yrec_modulation_{stem}_step{step:05d}{tag}.pdf"
    fig.savefig(output_dir / fname_out)
    plt.close(fig)
    print(f"  Saved {fname_out}")


# ---------------------------------------------------------------------------
# Figure 2: two-panel amplitude figure
#   Left  — amplitude vs f_spot, one filter, multiple masses (viridis colourmap)
#   Right — amplitude vs f_spot, one mass,   multiple filters (per-filter colours)
# ---------------------------------------------------------------------------
 
def compute_amplitude(df: pd.DataFrame, age_gyr: float,
                      fname: str) -> tuple[float, float] | None:
    row = interpolate_to_age(df, age_gyr)
    if row is None:
        return None
    hot_col  = f"{fname}_hot"
    cool_col = f"{fname}_cool"
    if hot_col not in row.index or cool_col not in row.index:
        return None
    m_hot  = float(row[hot_col])
    m_cool = float(row[cool_col])
    fspot  = float(row["Fspot"])
    if not (np.isfinite(m_hot) and np.isfinite(m_cool)) or fspot == 0:
        return None
    phase = np.linspace(0, 1, 500)
    dm    = compute_delta_mag(m_hot, m_cool, fspot, phase)
    return fspot, dm.max() - dm.min()
 
 
def make_amplitude_figure(
    track_dir:       Path,
    output_dir:      Path,
    age_gyr:         float,
    filter_name:     str          = 'V',    # filter for the left panel
    mass_multifilter: float       = 1.10,   # mass for the right panel
    masses:          list[float] | None = None,
    filters_to_show: list         = PLOT_FILTERS,
    tag:             str          = '',
):
    """
    Two-panel figure (full double-column width):
      Left  — peak-to-peak amplitude in `filter_name` vs f_spot for multiple masses.
      Right — peak-to-peak amplitude vs f_spot for multiple filters at `mass_multifilter`.
    Shared y-axis so panels align. No titles.
    """
    all_stems  = sorted(p.name for p in track_dir.iterdir() if p.is_dir())
    fspot_vals = sorted({_parse_stem(s)[1] for s in all_stems
                         if _parse_stem_safe(s) is not None and _parse_stem(s)[1] > 0})
 
    if masses is None:
        masses = sorted({_parse_stem(s)[0] for s in all_stems
                         if _parse_stem_safe(s) is not None})
 
    from matplotlib.gridspec import GridSpec
 
    fig = plt.figure(figsize=(APJ_SINGLE * 1.5, APJ_DOUBLE * 0.9))
    gs  = GridSpec(2, 2, figure=fig,
                   width_ratios=[1, 0.05],
                   hspace=0.0, wspace=0.0)
 
    ax_top = fig.add_subplot(gs[0, 0])
    ax_bot = fig.add_subplot(gs[1, 0], sharex=ax_top)
    cax    = fig.add_subplot(gs[0, 1])   # colorbar spans both rows
 
    plt.setp(ax_top.get_xticklabels(), visible=False)
 
    # ---- Top panel: one line per mass, colour = mass ----
    cmap = plt.cm.viridis
    norm = plt.Normalize(vmin=min(masses), vmax=max(masses))
 
    for mass in masses:
        pts = []
        for fspot in fspot_vals:
            stem = _mass_fspot_to_stem(mass, fspot)
            if not (track_dir / stem).exists():
                continue
            try:
                result = compute_amplitude(
                    read_colors_track(track_dir, stem), age_gyr, filter_name)
            except Exception:
                continue
            if result is not None:
                pts.append(result)
        if len(pts) < 2:
            continue
        pts.sort()
        fs_arr, amp_arr = zip(*pts)
        ax_top.plot(fs_arr, amp_arr, 'o-', color=cmap(norm(mass)),
                    lw=1.0, ms=3)
 
    ax_top.set_ylabel(r'Peak-to-peak amplitude (mag)', fontsize=10)
    ax_top.text(0.04, 0.96, rf'${filter_name}$ band',
                transform=ax_top.transAxes, ha='left', va='top', fontsize=8)
 
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, cax=cax, pad=0.0).set_label(r'$M\;(M_\odot)$', fontsize=8)
 
    # ---- Bottom panel: one line per filter, mass fixed ----
    for fname in filters_to_show:
        pts = []
        for fspot in fspot_vals:
            stem = _mass_fspot_to_stem(mass_multifilter, fspot)
            if not (track_dir / stem).exists():
                continue
            try:
                result = compute_amplitude(
                    read_colors_track(track_dir, stem), age_gyr, fname)
            except Exception:
                continue
            if result is not None:
                pts.append(result)
        if len(pts) < 2:
            continue
        pts.sort()
        fs_arr, amp_arr = zip(*pts)
        ax_bot.plot(fs_arr, amp_arr, 'o-',
                    color=JOHNSON[fname]['color'], lw=1.0, ms=3, label=fname)
 
    ax_bot.set_xlabel(r'$f_{\rm spot}$', fontsize=10)
    ax_bot.set_ylabel(r'Peak-to-peak amplitude (mag)', fontsize=10)
    ax_bot.legend(fontsize=7, title='Filter', loc='center left')
    ax_bot.text(0.04, 0.96, rf'$M = {mass_multifilter:.2f}\,M_\odot$',
                transform=ax_bot.transAxes, ha='left', va='top', fontsize=8)
 
    fname_out = (f"fig_yrec_amplitude_{filter_name}_"
                 f"m{int(mass_multifilter*100):03d}_age{age_gyr:.1f}gyr{tag}.pdf")
    fig.savefig(output_dir / fname_out, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {fname_out}")
 
# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
 
def _parse_stem(stem: str) -> tuple[float, float]:
    import re
    m = re.match(r"m(\d+)_f(\d+)", stem)
    if not m:
        raise ValueError(stem)
    return int(m.group(1)) / 100.0, int(m.group(2)) / 1000.0
 
def _parse_stem_safe(stem: str):
    try:
        return _parse_stem(stem)
    except ValueError:
        return None
 
def _mass_fspot_to_stem(mass: float, fspot: float) -> str:
    return f"m{int(round(mass * 100)):03d}_f{int(round(fspot * 1000)):03d}"
 
# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
 
if __name__ == "__main__":
 
    TRACK_DIR  = Path('/home/njm/MESA/MESA_Colors_Tests/YREC/')
    OUTPUT_DIR = Path('/home/njm/MESA/MESA_Colors_Tests/YREC/figures/')
    OUTPUT_DIR.mkdir(exist_ok=True)
 
    setup_style()
 
    # Figure 1 variants
    print("=== Modulation figures ===")
    mod_stems = ['m110_f034']
    for stem in mod_stems:
        if not (TRACK_DIR / stem).exists():
            continue
        make_modulation_figure(stem, TRACK_DIR, OUTPUT_DIR, tag='_midMS')
        df_tmp = read_colors_track(TRACK_DIR, stem)
        for frac, label in [(0.9, '_earlyMS'), (0.1, '_lateMS')]:
            s = find_ms_step(df_tmp, x_cen_frac=frac)
            make_modulation_figure(stem, TRACK_DIR, OUTPUT_DIR, step=s, tag=label)
 
    # Figure 2 variants: combined two-panel amplitude figure
    print("=== Amplitude figures ===")
    selected_masses = [0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20]
    for age in [3.0]:
        for filt in ['V', 'B', 'R']:
            for mass_mf in [1.10]:
                make_amplitude_figure(
                    TRACK_DIR, OUTPUT_DIR,
                    age_gyr=age,
                    filter_name=filt,
                    mass_multifilter=mass_mf,
                    masses=selected_masses,
                )
 
    print("\nAll figures written to", OUTPUT_DIR)
 
