#!/usr/bin/env python3
"""
plot_tpagb.py

Analysis and plotting for MESA Colors TP-AGB (Thermal Pulse AGB) demonstration.
Generates publication-quality figures showing photometric evolution during 
thermal pulses, with SED snapshots at key pulse phases.

Figures produced (each saved as .pdf AND .png, full range AND TP-phase crop):
- fig_tpagb_lightcurves     : Multi-band light curves with pulse markers
- fig_tpagb_cmd             : CMD with SED insets at pulse phases
- fig_tpagb_diagnostics     : He luminosity, core mass, Teff evolution
- fig_tpagb_summary         : 4-panel summary figure
- fig_tpagb_sed_movie.mp4   : Animated SED movie (SED + light curve + CMD per timestep)

Usage:
    python plot_tpagb.py [--logs_dir LOGS] [--sed_dir SED] [--output_dir figures]

Author: Miller, Joyce, Mocz et al.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from pathlib import Path
from collections import defaultdict
import argparse
import csv
import glob

from scipy.signal import find_peaks

from plot_utils import (
    setup_apj_style, read_mesa_history, find_magnitude_columns,
    get_model_number_at_index, add_em_spectrum_regions,
    JOHNSON_FILTERS, APJ_SINGLE_COL, APJ_DOUBLE_COL
)


# =============================================================================
# HELPERS
# =============================================================================

def save_fig(fig, output_dir, basename, dpi=300):
    """Save figure as both PDF and PNG."""
    for ext in ['pdf', 'png']:
        path = output_dir / f"{basename}.{ext}"
        fig.savefig(path, dpi=dpi)
        print(f"  Saved: {path}")


def find_tp_xlim(history, pre_pulse_margin_kyr=500):
    """
    Return (xmin, xmax) in kyr that covers just the thermal pulsing phase.
    xmin = first pulse time minus margin.
    xmax = end of run.
    """
    pulse_idx, pulse_times = identify_thermal_pulses(history)
    age = history['star_age']
    time_kyr = age / 1e3

    if len(pulse_times) == 0:
        return None  # Can't determine; caller will skip cropped version

    xmin = max(0, pulse_times[0] / 1e3 - pre_pulse_margin_kyr)
    xmax = time_kyr[-1]
    return (xmin, xmax)


def load_sed_csv(filepath):
    """
    Load SED data from CSV file.

    Reads wavelength/flux pairs row-by-row so arrays stay aligned.
    Rows where wavelength is zero (padding) are skipped.
    """
    try:
        wavelengths = []
        fluxes = []
        convolved = []

        with open(filepath, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    wl = float(row['wavelengths'])
                    fl = float(row['fluxes'])
                except (ValueError, KeyError):
                    continue
                if wl == 0:
                    continue  # skip zero-padded rows (keep fl=0 rows; those are real)
                wavelengths.append(wl)
                fluxes.append(fl)
                try:
                    convolved.append(float(row['convolved_flux']))
                except (ValueError, KeyError):
                    convolved.append(0.0)

        return {
            'wavelengths': np.array(wavelengths),
            'fluxes': np.array(fluxes),
            'convolved_flux': np.array(convolved),
        }
    except Exception as e:
        print(f"Warning: Could not read {filepath}: {e}")
        return {}


def identify_thermal_pulses(history, min_prominence=0.5):
    """
    Identify thermal pulse events from He-burning luminosity peaks.
    """
    log_LHe = history.get('log_LHe')
    age = history.get('star_age')

    if log_LHe is None or age is None:
        return np.array([]), np.array([])

    valid = np.isfinite(log_LHe)
    if np.sum(valid) < 10:
        return np.array([]), np.array([])

    peaks, _ = find_peaks(log_LHe[valid],
                          prominence=min_prominence,
                          distance=50)

    valid_indices = np.where(valid)[0]
    pulse_indices = valid_indices[peaks]
    pulse_times = age[pulse_indices]

    return pulse_indices, pulse_times


def find_poi_tpagb(history):
    """
    Find points of interest for TP-AGB evolution.
    Returns indices for pre_pulse, pulse_peak, interpulse, late_agb.
    """
    pulse_idx, pulse_times = identify_thermal_pulses(history)

    if len(pulse_idx) < 2:
        n = len(history['model_number'])
        return {'start': 0, 'middle': n // 2, 'end': n - 1}

    poi = {}

    if pulse_idx[0] > 10:
        poi['pre_pulse'] = pulse_idx[0] - 10
    else:
        poi['pre_pulse'] = 0

    mid_pulse = len(pulse_idx) // 2
    poi['pulse_peak'] = pulse_idx[mid_pulse]

    if mid_pulse + 1 < len(pulse_idx):
        interpulse_idx = (pulse_idx[mid_pulse] + pulse_idx[mid_pulse + 1]) // 2
        poi['interpulse'] = interpulse_idx

    poi['late_agb'] = len(history['model_number']) - 1

    return poi


# =============================================================================
# PLOT FUNCTIONS  (each accepts optional xlim for the TP-phase cropped version)
# =============================================================================

def plot_lightcurves(history, output_dir, basename='fig_tpagb_lightcurves', xlim=None):
    """
    Plot multi-band light curves during thermal pulses.

    Parameters
    ----------
    xlim : tuple or None
        (xmin, xmax) in kyr. If None, full range is used.
    """
    fig, axes = plt.subplots(2, 1, figsize=(APJ_DOUBLE_COL, 4.5), sharex=True)

    age = history['star_age']
    time_kyr = age / 1e3

    pulse_idx, pulse_times = identify_thermal_pulses(history)

    ax1 = axes[0]
    mag_cols = find_magnitude_columns(history, 'johnson')

    for band in ['B', 'V', 'R', 'I']:
        if band in mag_cols:
            mag = history[mag_cols[band]]
            valid = np.isfinite(mag)
            if np.sum(valid) > 10:
                ax1.plot(time_kyr[valid], mag[valid], '-',
                        color=JOHNSON_FILTERS[band]['color'],
                        label=band, lw=0.8, alpha=0.9)

    for pt in pulse_times:
        ax1.axvline(pt/1e3, color='gray', ls=':', alpha=0.4, lw=0.5)

    ax1.set_ylabel('Magnitude')
    ax1.invert_yaxis()
    ax1.legend(loc='upper right', ncol=4, fontsize=7, framealpha=0.9)
    ax1.set_title(f'TP-AGB Evolution ({len(pulse_idx)} thermal pulses)', fontsize=10)
    if xlim:
        ax1.set_xlim(xlim)

    ax2 = axes[1]
    if 'Mag_bol' in history:
        mag_bol = history['Mag_bol']
        valid = np.isfinite(mag_bol)
        ax2.plot(time_kyr[valid], mag_bol[valid], 'k-', lw=0.8, label='Bolometric')

    if 'log_L' in history:
        log_L = history['log_L']
        mbol_approx = 4.74 - 2.5 * log_L
        ax2.plot(time_kyr, mbol_approx, '--', color='gray', lw=0.6,
                label=r'From $\log L$', alpha=0.6)

    for pt in pulse_times:
        ax2.axvline(pt/1e3, color='gray', ls=':', alpha=0.4, lw=0.5)

    ax2.set_xlabel('Time (kyr)')
    ax2.set_ylabel('Bolometric Mag')
    ax2.invert_yaxis()
    ax2.legend(loc='upper right', fontsize=7, framealpha=0.9)
    if xlim:
        ax2.set_xlim(xlim)

    plt.tight_layout()
    save_fig(fig, output_dir, basename)
    plt.close()


def plot_cmd_with_seds(history, sed_dir, output_dir,
                       basename='fig_tpagb_cmd', xlim=None):
    """
    Plot CMD with SED insets at thermal pulse phases.

    Parameters
    ----------
    xlim : tuple or None
        (xmin, xmax) in kyr. Used to restrict which data points appear in CMD.
    """
    fig = plt.figure(figsize=(APJ_DOUBLE_COL, 4.5))

    gs = GridSpec(2, 3, width_ratios=[1.3, 1, 1], wspace=0.3, hspace=0.3,
                  left=0.08, right=0.97, top=0.92, bottom=0.12)

    ax_cmd = fig.add_subplot(gs[:, 0])

    mag_cols = find_magnitude_columns(history, 'johnson')

    if 'B' not in mag_cols or 'V' not in mag_cols:
        print("Warning: B and V not found, trying LSST filters")
        mag_cols = find_magnitude_columns(history, 'lsst')
        if 'g' in mag_cols and 'r' in mag_cols:
            color = history[mag_cols['g']] - history[mag_cols['r']]
            mag = history[mag_cols['g']]
            color_label = r'$g - r$'
            mag_label = r'$g$'
        else:
            print("Error: No suitable filter pairs found")
            plt.close()
            return
    else:
        B = history[mag_cols['B']]
        V = history[mag_cols['V']]
        color = B - V
        mag = V
        color_label = r'$B - V$'
        mag_label = r'$V$'

    age = history['star_age']
    time_kyr = age / 1e3

    valid = np.isfinite(color) & np.isfinite(mag)
    if xlim is not None:
        time_mask = (time_kyr >= xlim[0]) & (time_kyr <= xlim[1])
        valid = valid & time_mask

    scatter = ax_cmd.scatter(color[valid], mag[valid], c=time_kyr[valid],
                            cmap='viridis', s=3, alpha=0.7, rasterized=True)

    cbar = fig.colorbar(scatter, ax=ax_cmd, pad=0.02, aspect=25)
    cbar.set_label('Time (kyr)', fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    poi = find_poi_tpagb(history)

    poi_style = {
        'pre_pulse': {'marker': 'o', 'color': '#2ca02c', 'label': 'Pre-pulse'},
        'pulse_peak': {'marker': '*', 'color': '#d62728', 'label': 'Pulse peak'},
        'interpulse': {'marker': 's', 'color': '#1f77b4', 'label': 'Interpulse'},
        'late_agb': {'marker': 'D', 'color': '#9467bd', 'label': 'Late AGB'},
    }

    for name, idx in poi.items():
        if name in poi_style:
            style = poi_style[name]
            ax_cmd.scatter(color[idx], mag[idx], c=style['color'], s=80,
                          marker=style['marker'], edgecolors='black',
                          linewidths=0.8, zorder=5, label=style['label'])

    ax_cmd.set_xlabel(f'{color_label} (mag)')
    ax_cmd.set_ylabel(f'{mag_label} (mag)')
    ax_cmd.invert_yaxis()
    ax_cmd.legend(loc='lower left', fontsize=6, framealpha=0.9)

    # SED panels
    sed_pois = ['pre_pulse', 'pulse_peak', 'interpulse', 'late_agb']
    sed_titles = ['Pre-pulse', 'Pulse Peak', 'Interpulse', 'Late AGB']
    sed_positions = [(0, 1), (0, 2), (1, 1), (1, 2)]

    for (row, col), poi_name, title in zip(sed_positions, sed_pois, sed_titles):
        ax_sed = fig.add_subplot(gs[row, col])

        if poi_name in poi:
            idx = poi[poi_name]
            model_num = get_model_number_at_index(history, idx)

            sed_file = find_sed_file(sed_dir, model_num)

            if sed_file:
                sed_data = load_sed_csv(sed_file)

                if sed_data and 'wavelengths' in sed_data:
                    wl = sed_data['wavelengths']
                    fl = sed_data['fluxes']

                    ax_sed.plot(wl, fl, 'k-', lw=0.7)
                    ax_sed.set_xscale('log')

                    for filt, props in JOHNSON_FILTERS.items():
                        if filt in ['B', 'V', 'R', 'I']:
                            ax_sed.axvline(props['wavelength'],
                                          color=props['color'],
                                          alpha=0.4, lw=1.2, zorder=0)

                    ax_sed.set_xlim(3500, 25000)

                    mask_vis = (wl > 3500) & (wl < 25000)
                    if mask_vis.any():
                        fl_vis = fl[mask_vis]
                        if len(fl_vis) > 0 and fl_vis.max() > 0:
                            ax_sed.set_ylim(0, fl_vis.max() * 1.1)

                    add_em_spectrum_regions(ax_sed, alpha=0.03)

            teff = 10**history['log_Teff'][idx]
            ax_sed.text(0.97, 0.95, f'{teff:.0f} K',
                       transform=ax_sed.transAxes, fontsize=6,
                       ha='right', va='top',
                       bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                                alpha=0.7, edgecolor='none'))

        ax_sed.set_title(title, fontsize=8)
        ax_sed.tick_params(labelsize=6)

        if row == 1:
            ax_sed.set_xlabel(r'$\lambda$ ($\mathrm{\AA}$)', fontsize=7)
        else:
            ax_sed.set_xticklabels([])

    save_fig(fig, output_dir, basename)
    plt.close()


def plot_diagnostics(history, output_dir, basename='fig_tpagb_diagnostics', xlim=None):
    """
    Plot TP-AGB diagnostic panels.
    """
    fig, axes = plt.subplots(4, 1, figsize=(APJ_SINGLE_COL, 6), sharex=True)

    age = history['star_age']
    time_kyr = age / 1e3

    pulse_idx, pulse_times = identify_thermal_pulses(history)

    ax1 = axes[0]
    if 'log_LHe' in history:
        log_LHe = history['log_LHe']
        valid = np.isfinite(log_LHe)
        ax1.plot(time_kyr[valid], log_LHe[valid], 'r-', lw=0.7)

        for pi in pulse_idx:
            ax1.plot(time_kyr[pi], log_LHe[pi], 'ko', ms=3)

    ax1.set_ylabel(r'$\log(L_{\rm He}/L_\odot)$', fontsize=8)
    ax1.set_title(f'{len(pulse_idx)} thermal pulses', fontsize=9)
    if xlim:
        ax1.set_xlim(xlim)

    ax2 = axes[1]
    if 'log_Teff' in history:
        log_Teff = history['log_Teff']
        ax2.plot(time_kyr, 10**log_Teff, 'k-', lw=0.7)
    ax2.set_ylabel(r'$T_{\rm eff}$ (K)', fontsize=8)

    ax3 = axes[2]
    if 'he_core_mass' in history:
        ax3.plot(time_kyr, history['he_core_mass'], 'k-', lw=0.8, label='He core')
    if 'c_core_mass' in history:
        ax3.plot(time_kyr, history['c_core_mass'], 'b--', lw=0.6, label='C/O core')
    ax3.set_ylabel(r'Core mass ($M_\odot$)', fontsize=8)
    ax3.legend(loc='upper left', fontsize=6, framealpha=0.9)

    ax4 = axes[3]
    mag_cols = find_magnitude_columns(history, 'johnson')
    if 'V' in mag_cols:
        V = history[mag_cols['V']]
        valid = np.isfinite(V)
        ax4.plot(time_kyr[valid], V[valid], 'g-', lw=0.7)
    ax4.set_ylabel('$V$ (mag)', fontsize=8)
    ax4.invert_yaxis()
    ax4.set_xlabel('Time (kyr)')
    if xlim:
        ax4.set_xlim(xlim)

    for ax in axes:
        for pt in pulse_times:
            ax.axvline(pt/1e3, color='gray', ls=':', alpha=0.3, lw=0.5)

    plt.tight_layout()
    save_fig(fig, output_dir, basename)
    plt.close()


def plot_summary(history, sed_dir, output_dir, basename='fig_tpagb_summary', xlim=None):
    """
    Create 4-panel summary figure with SED evolution.
    """
    fig = plt.figure(figsize=(APJ_DOUBLE_COL, 5))
    gs = GridSpec(2, 2, hspace=0.3, wspace=0.3,
                  left=0.08, right=0.95, top=0.93, bottom=0.10)

    age = history['star_age']
    time_kyr = age / 1e3

    pulse_idx, pulse_times = identify_thermal_pulses(history)

    # (a) Multi-band light curves
    ax1 = fig.add_subplot(gs[0, 0])

    mag_cols = find_magnitude_columns(history, 'johnson')
    for band in ['B', 'V', 'R', 'I']:
        if band in mag_cols:
            mag = history[mag_cols[band]]
            valid = np.isfinite(mag)
            if np.sum(valid) > 10:
                ax1.plot(time_kyr[valid], mag[valid], '-',
                        color=JOHNSON_FILTERS[band]['color'],
                        label=band, lw=0.7, alpha=0.9)

    ax1.set_xlabel('Time (kyr)')
    ax1.set_ylabel('Magnitude')
    ax1.invert_yaxis()
    ax1.legend(loc='upper right', ncol=2, fontsize=6)
    ax1.set_title('(a) Multi-band Light Curves', fontsize=9)
    if xlim:
        ax1.set_xlim(xlim)

    # (b) CMD
    ax2 = fig.add_subplot(gs[0, 1])

    if 'B' in mag_cols and 'V' in mag_cols:
        B = history[mag_cols['B']]
        V = history[mag_cols['V']]
        BmV = B - V
        valid = np.isfinite(BmV) & np.isfinite(V)
        if xlim is not None:
            time_mask = (time_kyr >= xlim[0]) & (time_kyr <= xlim[1])
            valid = valid & time_mask

        scatter = ax2.scatter(BmV[valid], V[valid], c=time_kyr[valid],
                             cmap='viridis', s=2, alpha=0.6, rasterized=True)
        if valid.any():
            first_valid = np.where(valid)[0][0]
            last_valid = np.where(valid)[0][-1]
            ax2.plot(BmV[first_valid], V[first_valid], 'go', ms=6, zorder=10)
            ax2.plot(BmV[last_valid], V[last_valid], 'rs', ms=6, zorder=10)

    ax2.set_xlabel(r'$B - V$')
    ax2.set_ylabel('$V$')
    ax2.invert_yaxis()
    ax2.set_title('(b) Color-Magnitude Diagram', fontsize=9)

    # (c) HR diagram
    ax3 = fig.add_subplot(gs[1, 0])

    if 'log_Teff' in history and 'log_L' in history:
        log_Teff = history['log_Teff']
        log_L = history['log_L']

        if xlim is not None:
            time_mask = (time_kyr >= xlim[0]) & (time_kyr <= xlim[1])
            scatter = ax3.scatter(log_Teff[time_mask], log_L[time_mask],
                                 c=time_kyr[time_mask], cmap='viridis',
                                 s=2, alpha=0.6, rasterized=True)
            idx_arr = np.where(time_mask)[0]
            ax3.plot(log_Teff[idx_arr[0]], log_L[idx_arr[0]], 'go', ms=6, zorder=10)
            ax3.plot(log_Teff[idx_arr[-1]], log_L[idx_arr[-1]], 'rs', ms=6, zorder=10)
        else:
            scatter = ax3.scatter(log_Teff, log_L, c=time_kyr, cmap='viridis',
                                 s=2, alpha=0.6, rasterized=True)
            ax3.plot(log_Teff[0], log_L[0], 'go', ms=6, zorder=10)
            ax3.plot(log_Teff[-1], log_L[-1], 'rs', ms=6, zorder=10)

    ax3.set_xlabel(r'$\log T_{\rm eff}$ (K)')
    ax3.set_ylabel(r'$\log(L/L_\odot)$')
    ax3.invert_xaxis()
    ax3.set_title('(c) HR Diagram', fontsize=9)

    # (d) SED evolution
    ax4 = fig.add_subplot(gs[1, 1])

    poi = find_poi_tpagb(history)
    cmap_sed = plt.cm.viridis

    for i, (poi_name, label) in enumerate([
        ('pre_pulse', 'Pre-pulse'),
        ('pulse_peak', 'Pulse peak'),
        ('interpulse', 'Interpulse'),
        ('late_agb', 'Late AGB')
    ]):
        if poi_name not in poi:
            continue

        idx = poi[poi_name]
        model_num = get_model_number_at_index(history, idx)

        sed_file = find_sed_file(sed_dir, model_num)
        if sed_file:
            sed_data = load_sed_csv(sed_file)
            if sed_data and 'wavelengths' in sed_data:
                wl = sed_data['wavelengths']
                fl = sed_data['fluxes']

                color_sed = cmap_sed(i / 3)
                ax4.plot(wl, fl, color=color_sed, lw=0.8, label=label, alpha=0.9)

    ax4.set_xscale('log')
    ax4.set_xlabel(r'Wavelength ($\mathrm{\AA}$)')
    ax4.set_ylabel(r'$F_\lambda$')
    ax4.set_xlim(3000, 30000)
    ax4.legend(loc='upper right', fontsize=6, framealpha=0.9)
    ax4.set_title('(d) SED Evolution', fontsize=9)

    add_em_spectrum_regions(ax4, alpha=0.03)

    fig.suptitle('MESA Colors: Thermal Pulse AGB Evolution', fontsize=11)

    save_fig(fig, output_dir, basename)
    plt.close()


# =============================================================================
# SED MOVIE
# =============================================================================

def find_sed_file(sed_dir, model_num, preferred_prefix='V'):
    """
    Find SED file for a given model number, handling zero-padded filenames.
    Tries preferred_prefix first, then any available filter.
    """
    sed_dir = Path(sed_dir)
    # Try zero-padded 8-digit format first (standard MESA Colors output)
    for prefix in ([preferred_prefix] + ['B', 'R', 'I', 'U']):
        p = sed_dir / f'{prefix}_SED_{model_num:08d}.csv'
        if p.exists():
            return p
    # Fallback: glob with zero-padded
    matches = list(sed_dir.glob(f'*_SED_{model_num:08d}.csv'))
    if matches:
        return matches[0]
    # Last resort: unpadded
    matches = list(sed_dir.glob(f'*_SED_{model_num}.csv'))
    return matches[0] if matches else None
    """
    Load one SED CSV file per model number from sed_dir.

    MESA Colors writes one file per filter per timestep (e.g. U_SED_00000319.csv,
    B_SED_00000319.csv, …).  The stellar SED (wavelengths, fluxes) is identical
    across filter files for the same model, so we only need one.

    Parameters
    ----------
    sed_dir : Path
    preferred_prefix : str or None
        Prefer files starting with this prefix (e.g. 'V') when multiple filters
        exist for the same model.  If None, any file is accepted.

    Returns
    -------
    list of (model_number, sed_data) tuples, sorted by model_number
    """
    sed_dir = Path(sed_dir)
    all_files = sorted(sed_dir.glob('*_SED_*.csv'))

    # Group files by model number
    by_model = defaultdict(list)
    for f in all_files:
        stem = f.stem  # e.g. 'V_SED_00000319'
        parts = stem.split('_')
        try:
            model_num = int(parts[-1])
        except ValueError:
            continue
        by_model[model_num].append(f)

    seds = []
    for model_num in sorted(by_model.keys()):
        candidates = by_model[model_num]
        chosen = None
        if preferred_prefix:
            for c in candidates:
                if c.name.startswith(preferred_prefix + '_'):
                    chosen = c
                    break
        if chosen is None:
            chosen = candidates[0]

        data = load_sed_csv(chosen)
        if data and 'wavelengths' in data and len(data['wavelengths']) > 0:
            seds.append((model_num, data))

    return seds


def load_all_seds(sed_dir, preferred_prefix=None):
    """
    Load one SED CSV file per model number from sed_dir.

    MESA Colors writes one file per filter per timestep (e.g. U_SED_00000319.csv,
    B_SED_00000319.csv, …).  The stellar SED (wavelengths, fluxes) is identical
    across filter files for the same model, so we only need one.

    Parameters
    ----------
    sed_dir : Path
    preferred_prefix : str or None
        Prefer files starting with this prefix (e.g. 'V') when multiple filters
        exist for the same model.  If None, any file is accepted.

    Returns
    -------
    list of (model_number, sed_data) tuples, sorted by model_number
    """
    sed_dir = Path(sed_dir)
    all_files = sorted(sed_dir.glob('*_SED_*.csv'))

    # Group files by model number
    by_model = defaultdict(list)
    for f in all_files:
        stem = f.stem  # e.g. 'V_SED_00000319'
        parts = stem.split('_')
        try:
            model_num = int(parts[-1])
        except ValueError:
            continue
        by_model[model_num].append(f)

    seds = []
    for model_num in sorted(by_model.keys()):
        candidates = by_model[model_num]
        chosen = None
        if preferred_prefix:
            for c in candidates:
                if c.name.startswith(preferred_prefix + '_'):
                    chosen = c
                    break
        if chosen is None:
            chosen = candidates[0]

        data = load_sed_csv(chosen)
        if data and 'wavelengths' in data and len(data['wavelengths']) > 0:
            seds.append((model_num, data))

    return seds


def make_sed_movie(history, sed_dir, output_dir, filename='fig_tpagb_sed_movie.mp4',
                   fps=15, tp_only=False):
    """
    Create an animation showing SED evolution frame by frame, alongside
    the multi-band light curve (with time cursor) and CMD (with current point).

    Parameters
    ----------
    history : dict
        MESA history data
    sed_dir : Path
        Directory containing SED CSV files
    output_dir : Path
        Output directory
    filename : str
        Output movie filename
    fps : int
        Frames per second
    tp_only : bool
        If True, only animate the thermal pulsing phase
    """
    print("Loading all SED files for movie...")
    seds = load_all_seds(sed_dir, preferred_prefix='V')

    if not seds:
        print("Warning: No SED files found. Skipping movie.")
        return

    print(f"  Found {len(seds)} SED frames")

    age = history['star_age']
    time_kyr = age / 1e3
    model_nums = history['model_number'].astype(int)

    mag_cols = find_magnitude_columns(history, 'johnson')
    has_bv = ('B' in mag_cols and 'V' in mag_cols)

    # Map model number -> history index
    model_to_idx = {int(m): i for i, m in enumerate(model_nums)}

    # Filter SEDs to TP phase if requested
    if tp_only:
        tp_xlim = find_tp_xlim(history)
        if tp_xlim is not None:
            tp_t_min_yr = tp_xlim[0] * 1e3
            tp_t_max_yr = tp_xlim[1] * 1e3
            seds = [(m, d) for m, d in seds
                    if m in model_to_idx
                    and tp_t_min_yr <= age[model_to_idx[m]] <= tp_t_max_yr]
            print(f"  Filtered to TP phase: {len(seds)} frames")

    if not seds:
        print("Warning: No SED frames in TP phase. Skipping movie.")
        return

    # Determine SED wavelength range and flux scale
    wl_min, wl_max = 3000, 30000
    peak_fluxes = []
    for _, d in seds:
        wl = d['wavelengths']
        fl = d['fluxes']
        mask = (wl >= wl_min) & (wl <= wl_max) & np.isfinite(fl)
        if mask.any():
            peak_fluxes.append(fl[mask].max())
    fl_max = (np.percentile(peak_fluxes, 99) * 1.15) if peak_fluxes else 1.0
    fl_max = fl_max if fl_max > 0 else 1.0

    # Precompute history arrays for the movie (only within TP phase if needed)
    if tp_only and tp_xlim is not None:
        time_mask = (time_kyr >= tp_xlim[0]) & (time_kyr <= tp_xlim[1])
    else:
        time_mask = np.ones(len(time_kyr), dtype=bool)

    lc_time = time_kyr[time_mask]
    lc_bands = {}
    for band in ['B', 'V', 'R', 'I']:
        if band in mag_cols:
            m = history[mag_cols[band]][time_mask]
            valid = np.isfinite(m)
            if valid.sum() > 5:
                lc_bands[band] = (lc_time[valid], m[valid])

    if has_bv:
        B_all = history[mag_cols['B']][time_mask]
        V_all = history[mag_cols['V']][time_mask]
        BmV_all = B_all - V_all
        V_plot = V_all
        cmd_valid = np.isfinite(BmV_all) & np.isfinite(V_plot)
    else:
        cmd_valid = np.zeros(time_mask.sum(), dtype=bool)

    # Figure layout: SED top-left, light curve top-right, CMD bottom-left,
    #                He-core mass vs B-R bottom-right
    fig = plt.figure(figsize=(APJ_DOUBLE_COL, 6))
    gs = GridSpec(2, 2, hspace=0.35, wspace=0.32,
                  left=0.09, right=0.97, top=0.93, bottom=0.10)

    ax_sed  = fig.add_subplot(gs[0, 0])
    ax_lc   = fig.add_subplot(gs[0, 1])
    ax_cmd  = fig.add_subplot(gs[1, 0])
    ax_hec  = fig.add_subplot(gs[1, 1])

    # Static light curve background
    for band, (bt, bm) in lc_bands.items():
        ax_lc.plot(bt, bm, '-', color=JOHNSON_FILTERS[band]['color'],
                   label=band, lw=0.6, alpha=0.7)
    ax_lc.invert_yaxis()
    ax_lc.set_xlabel('Time (kyr)', fontsize=8)
    ax_lc.set_ylabel('Magnitude', fontsize=8)
    ax_lc.legend(loc='upper right', ncol=4, fontsize=6)
    ax_lc.set_title('Light Curves', fontsize=9)
    # Rescale x-axis to match the actual data range (important for tp_only mode)
    if lc_bands:
        all_times = np.concatenate([bt for bt, _ in lc_bands.values()])
        ax_lc.set_xlim(all_times.min(), all_times.max())

    # Static CMD background
    if cmd_valid.any():
        ax_cmd.scatter(BmV_all[cmd_valid], V_plot[cmd_valid],
                       c='lightgray', s=2, alpha=0.5, rasterized=True, zorder=1)
    ax_cmd.invert_yaxis()
    ax_cmd.set_xlabel(r'$B - V$ (mag)', fontsize=8)
    ax_cmd.set_ylabel('$V$ (mag)', fontsize=8)
    ax_cmd.set_title('CMD', fontsize=9)

    # Static He-core mass vs B-R background
    has_hec = 'he_core_mass' in history
    has_br  = has_bv and 'R' in mag_cols
    if has_hec and has_br:
        hec_all  = history['he_core_mass'][time_mask]
        R_all    = history[mag_cols['R']][time_mask]
        B_all_hec = history[mag_cols['B']][time_mask]
        BmR_all  = B_all_hec - R_all
        hec_br_valid = np.isfinite(hec_all) & np.isfinite(BmR_all)
        if hec_br_valid.any():
            ax_hec.scatter(BmR_all[hec_br_valid], hec_all[hec_br_valid],
                           c='lightgray', s=2, alpha=0.5, rasterized=True, zorder=1)
    ax_hec.set_xlabel(r'$B - R$ (mag)', fontsize=8)
    ax_hec.set_ylabel(r'$M_{\rm He\,core}\,(M_\odot)$', fontsize=8)
    ax_hec.set_title(r'He Core Mass vs $B-R$', fontsize=9)

    # Dynamic elements
    sed_line, = ax_sed.plot([], [], 'k-', lw=0.9)
    ax_sed.set_xscale('log')
    ax_sed.set_xlim(wl_min, wl_max)
    ax_sed.set_ylim(0, fl_max)
    ax_sed.set_xlabel(r'Wavelength ($\mathrm{\AA}$)', fontsize=8)
    ax_sed.set_ylabel(r'$F_\lambda$', fontsize=8)
    ax_sed.set_title('SED', fontsize=9)
    for filt, props in JOHNSON_FILTERS.items():
        if filt in ['B', 'V', 'R', 'I']:
            ax_sed.axvline(props['wavelength'], color=props['color'],
                          alpha=0.35, lw=1.0, zorder=0)
    add_em_spectrum_regions(ax_sed, alpha=0.03)

    vline    = ax_lc.axvline(0, color='red', lw=0.8, ls='--', alpha=0.8)
    cmd_dot, = ax_cmd.plot([], [], 'ro', ms=5, zorder=5)
    hec_dot, = ax_hec.plot([], [], 'ro', ms=5, zorder=5)

    fig.suptitle('MESA Colors: TP-AGB SED Evolution', fontsize=10)

    # Pulse markers
    pulse_idx, pulse_times = identify_thermal_pulses(history)
    for pt in pulse_times:
        t_kyr = pt / 1e3
        if not tp_only or (tp_xlim and tp_xlim[0] <= t_kyr <= tp_xlim[1]):
            ax_lc.axvline(t_kyr, color='gray', ls=':', alpha=0.3, lw=0.5)

    def init():
        sed_line.set_data([], [])
        vline.set_xdata([0])
        cmd_dot.set_data([], [])
        hec_dot.set_data([], [])
        return sed_line, vline, cmd_dot, hec_dot

    def update(frame_idx):
        model_num, sed_data = seds[frame_idx]
        wl = sed_data['wavelengths']
        fl = sed_data['fluxes']

        # Mask to plot range
        mask = (wl >= wl_min) & (wl <= wl_max)
        sed_line.set_data(wl[mask], fl[mask])

        # History index for this model
        hidx = model_to_idx.get(model_num)
        if hidx is not None:
            t_kyr = time_kyr[hidx]
            vline.set_xdata([t_kyr])

            if has_bv and cmd_valid.any():
                b_val = history[mag_cols['B']][hidx]
                v_val = history[mag_cols['V']][hidx]
                if np.isfinite(b_val) and np.isfinite(v_val):
                    cmd_dot.set_data([b_val - v_val], [v_val])

            if has_hec and has_br:
                b_val = history[mag_cols['B']][hidx]
                r_val = history[mag_cols['R']][hidx]
                hec_val = history['he_core_mass'][hidx]
                if np.isfinite(b_val) and np.isfinite(r_val) and np.isfinite(hec_val):
                    hec_dot.set_data([b_val - r_val], [hec_val])

        return sed_line, vline, cmd_dot, hec_dot

    n_frames = len(seds)
    print(f"  Rendering {n_frames} frames at {fps} fps...")

    ani = animation.FuncAnimation(fig, update, frames=n_frames,
                                  init_func=init, blit=False, interval=1000/fps)

    out_path = output_dir / filename
    writer = animation.FFMpegWriter(fps=fps, bitrate=1800,
                                    metadata={'title': 'TP-AGB SED Evolution'})
    try:
        ani.save(str(out_path), writer=writer)
        print(f"  Saved: {out_path}")
    except Exception as e:
        print(f"  Warning: Could not save movie with FFMpegWriter: {e}")
        print("  Trying PillowWriter (gif)...")
        gif_path = output_dir / filename.replace('.mp4', '.gif')
        writer_gif = animation.PillowWriter(fps=fps)
        ani.save(str(gif_path), writer=writer_gif)
        print(f"  Saved: {gif_path}")

    plt.close()


# =============================================================================
# SINGLE-PULSE ANALYSIS
# =============================================================================

def select_single_pulse(history, skip_first=5):
    """
    Select a well-formed thermal pulse for detailed examination.

    Skips the first `skip_first` pulses (which are still growing in amplitude)
    and returns the next fully-developed pulse together with four key phase
    indices:

        lead_in   – deep quiescent phase: 20 % of the way from the *previous*
                    pulse peak toward this one (well after the prior pulse has
                    decayed, well before the next rise begins)
        knee      – onset of rapid rise: first index where log_LHe crosses
                    25 % of the way from baseline to peak
        peak      – maximum log_LHe of the selected pulse
        return_q  – midpoint of the interpulse period *after* this pulse

    Window is padded by 15 % on each side so the pulse anatomy is visible
    with generous context.

    Returns dict or None if there are not enough pulses.
    """
    pulse_idx, pulse_times = identify_thermal_pulses(history)

    if len(pulse_idx) < skip_first + 2:
        if len(pulse_idx) < 2:
            return None
        chosen = len(pulse_idx) - 2
    else:
        chosen = skip_first

    log_LHe = history['log_LHe']
    age     = history['star_age']

    prev_peak_idx = pulse_idx[chosen - 1] if chosen > 0 else 0
    this_peak_idx = pulse_idx[chosen]
    next_peak_idx = (pulse_idx[chosen + 1]
                     if chosen + 1 < len(pulse_idx) else len(age) - 1)

    # --- lead_in: 20 % from previous peak toward this peak -----------------
    # This puts us well into the quiescent interpulse, far from the knee.
    lead_in_idx = prev_peak_idx + int(0.20 * (this_peak_idx - prev_peak_idx))

    # --- knee: first point on rising slope above 25 % of peak excursion ----
    baseline  = log_LHe[lead_in_idx]
    peak_val  = log_LHe[this_peak_idx]
    threshold = baseline + 0.25 * (peak_val - baseline)
    knee_idx  = lead_in_idx
    for k in range(lead_in_idx, this_peak_idx):
        if np.isfinite(log_LHe[k]) and log_LHe[k] >= threshold:
            knee_idx = k
            break

    # --- return_q: midpoint of interpulse after this peak -------------------
    return_q_idx = (this_peak_idx + next_peak_idx) // 2

    # --- window: lead_in → return_q with 15 % padding on each side ----------
    t_start = age[lead_in_idx]
    t_end   = age[return_q_idx]
    pad     = 0.15 * (t_end - t_start)
    t_window = (t_start - pad, t_end + pad)

    return {
        'pulse_number'   : chosen,
        'pulse_hist_idx' : this_peak_idx,
        'lead_in'        : lead_in_idx,
        'knee'           : knee_idx,
        'peak'           : this_peak_idx,
        'return_q'       : return_q_idx,
        't_window'       : t_window,
    }


# ---------------------------------------------------------------------------
# STATIC FIGURE  –  single pulse + four SED panels
# ---------------------------------------------------------------------------

def plot_single_pulse_with_seds(history, sed_dir, output_dir,
                                basename='fig_tpagb_single_pulse',
                                skip_first=5):
    """
    Left column  : Two stacked panels sharing the x-axis —
                     top:    log_LHe vs time  (He-burning driver)
                     bottom: V-band magnitude vs time  (observed response)
                   Both show colour-coded vertical markers at the four phases.
    Right panels : 2 × 2 grid of SEDs at those four phases, y-axis unified.
    """
    info = select_single_pulse(history, skip_first=skip_first)
    if info is None:
        print("Warning: Not enough thermal pulses to isolate a single pulse. Skipping.")
        return

    log_LHe  = history['log_LHe']
    age      = history['star_age']
    time_kyr = age / 1e3

    t_win_kyr = (info['t_window'][0] / 1e3, info['t_window'][1] / 1e3)
    win_mask  = (time_kyr >= t_win_kyr[0]) & (time_kyr <= t_win_kyr[1])

    # V-band magnitudes
    mag_cols = find_magnitude_columns(history, 'johnson')
    has_V    = 'V' in mag_cols
    if has_V:
        V_mag = history[mag_cols['V']]

    # ---- Phase definitions --------------------------------------------------
    phases = [
        ('lead_in',  'Lead-in',         '#4878d0'),
        ('knee',     'Onset (knee)',     '#ee854a'),
        ('peak',     'Pulse peak',       '#d62728'),
        ('return_q', 'Return to quiet',  '#6acc65'),
    ]

    # ---- Consistent SED y-axis across panels --------------------------------
    wl_min, wl_max = 3000, 30000
    fl_peaks = []
    for phase_key, _, _ in phases:
        idx = info[phase_key]
        mn  = get_model_number_at_index(history, idx)
        sf  = find_sed_file(sed_dir, mn)
        if sf:
            d = load_sed_csv(sf)
            if d and 'wavelengths' in d:
                m = (d['wavelengths'] >= wl_min) & (d['wavelengths'] <= wl_max)
                if m.any() and d['fluxes'][m].max() > 0:
                    fl_peaks.append(d['fluxes'][m].max())
    fl_max_global = (np.percentile(fl_peaks, 99) * 1.15) if fl_peaks else 1.0

    # ---- Figure layout ------------------------------------------------------
    # Left column: two stacked panels — V-band (larger, observable) on top,
    #              log_LHe (smaller, driver) on bottom
    # Right column: 2×2 SED grid, tightly packed
    fig = plt.figure(figsize=(APJ_DOUBLE_COL, 5.0))

    gs_outer = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.6],
                                 wspace=0.28,
                                 left=0.08, right=0.98,
                                 top=0.91, bottom=0.09)

    # V-band on top (larger) to emphasise the observable; log_LHe below
    gs_left = gs_outer[0].subgridspec(2, 1, hspace=0.06,
                                       height_ratios=[1.4, 1])
    ax_V   = fig.add_subplot(gs_left[0])
    ax_lhe = fig.add_subplot(gs_left[1], sharex=ax_V)

    gs_sed = gs_outer[1].subgridspec(2, 2, hspace=0.12, wspace=0.12)

    # ---- helper: draw phase markers on an axis ------------------------------
    def draw_phase_markers(ax, draw_legend=False):
        handles = []
        for phase_key, label, color in phases:
            idx  = info[phase_key]
            t_ph = time_kyr[idx]
            ax.axvline(t_ph, color=color, lw=1.2, ls='--', alpha=0.85, zorder=3)
            handles.append(Line2D([0], [0], color=color, lw=1.5, ls='--',
                                   marker='o', markersize=5,
                                   markeredgecolor='k', markeredgewidth=0.4,
                                   label=label))
        if draw_legend:
            ax.legend(handles=handles, fontsize=6.5,
                      loc='upper left', framealpha=0.9)
        return handles

    # ---- Top left: V-band photometry (observable — the focus) ---------------
    if has_V:
        valid_V = win_mask & np.isfinite(V_mag)
        ax_V.plot(time_kyr[valid_V], V_mag[valid_V],
                  color=JOHNSON_FILTERS['V']['color'], lw=1.3, zorder=2)

        for phase_key, label, color in phases:
            idx  = info[phase_key]
            t_ph = time_kyr[idx]
            v_pt = V_mag[idx]
            if np.isfinite(v_pt):
                ax_V.plot(t_ph, v_pt, 'o', color=color, ms=6, zorder=5,
                          markeredgecolor='k', markeredgewidth=0.5)
        draw_phase_markers(ax_V, draw_legend=True)

        ax_V.invert_yaxis()
        ax_V.set_ylabel(r'$V$ (mag)', fontsize=8)
    else:
        ax_V.text(0.5, 0.5, 'V-band not available',
                  ha='center', va='center', transform=ax_V.transAxes, fontsize=8)
        draw_phase_markers(ax_V, draw_legend=True)

    ax_V.set_xlim(t_win_kyr)
    ax_V.set_title(f'Pulse #{info["pulse_number"] + 1}', fontsize=9)
    plt.setp(ax_V.get_xticklabels(), visible=False)

    # ---- Bottom left: log_LHe (physical driver) -----------------------------
    valid_lhe = win_mask & np.isfinite(log_LHe)
    ax_lhe.plot(time_kyr[valid_lhe], log_LHe[valid_lhe],
                'k-', lw=1.1, zorder=2)

    for phase_key, label, color in phases:
        idx  = info[phase_key]
        t_ph = time_kyr[idx]
        lhe  = log_LHe[idx]
        ax_lhe.plot(t_ph, lhe, 'o', color=color, ms=6, zorder=5,
                    markeredgecolor='k', markeredgewidth=0.5)
    draw_phase_markers(ax_lhe)

    ax_lhe.set_xlim(t_win_kyr)
    ax_lhe.set_ylabel(r'$\log(L_{\rm He}/L_\odot)$', fontsize=8)
    ax_lhe.set_xlabel('Time (kyr)', fontsize=8)

    # ---- Right: four SED panels — share x and y axes -----------------------
    sed_axes = []
    for i, (phase_key, label, color) in enumerate(phases):
        row, col = divmod(i, 2)
        # Share x and y with first panel for uniform axes, no wasted space
        if i == 0:
            ax_sed = fig.add_subplot(gs_sed[row, col])
        elif col == 0:
            ax_sed = fig.add_subplot(gs_sed[row, col], sharey=sed_axes[0])
        else:
            ax_sed = fig.add_subplot(gs_sed[row, col],
                                     sharex=sed_axes[row * 2],
                                     sharey=sed_axes[0])
        sed_axes.append(ax_sed)

        idx = info[phase_key]
        mn  = get_model_number_at_index(history, idx)
        sf  = find_sed_file(sed_dir, mn)

        if sf:
            d = load_sed_csv(sf)
            if d and 'wavelengths' in d:
                wl = d['wavelengths']
                fl = d['fluxes']
                m  = (wl >= wl_min) & (wl <= wl_max)
                ax_sed.plot(wl[m], fl[m], '-', color=color, lw=0.85)

                for filt, props in JOHNSON_FILTERS.items():
                    if filt in ['B', 'V', 'R', 'I']:
                        ax_sed.axvline(props['wavelength'],
                                       color=props['color'],
                                       alpha=0.35, lw=0.8, zorder=0)
                add_em_spectrum_regions(ax_sed, alpha=0.03)

        ax_sed.set_xscale('log')
        ax_sed.set_xlim(wl_min, wl_max)
        ax_sed.set_ylim(0, fl_max_global)

        teff = 10**history['log_Teff'][idx] if 'log_Teff' in history else float('nan')
        t_ph = time_kyr[idx]
        ax_sed.text(0.97, 0.95, f'{teff:.0f} K\nt={t_ph:.1f} kyr',
                    transform=ax_sed.transAxes, fontsize=5.5,
                    ha='right', va='top', family='monospace',
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                              alpha=0.75, edgecolor='none'))

        for spine in ax_sed.spines.values():
            spine.set_linewidth(0.6)
        ax_sed.spines['top'].set_color(color)
        ax_sed.spines['top'].set_linewidth(2.0)

        ax_sed.set_title(label, fontsize=7.5, color=color, fontweight='bold', pad=2)
        ax_sed.tick_params(labelsize=5.5)

        # Only show axis labels on outer edges
        if row == 1:
            ax_sed.set_xlabel(r'$\lambda$ ($\mathrm{\AA}$)', fontsize=6.5)
        else:
            plt.setp(ax_sed.get_xticklabels(), visible=False)
        if col == 0:
            ax_sed.set_ylabel(r'$F_\lambda$', fontsize=6.5)
        else:
            plt.setp(ax_sed.get_yticklabels(), visible=False)

    fig.suptitle('TP-AGB: Single Thermal Pulse Anatomy', fontsize=10)

    save_fig(fig, output_dir, basename)
    plt.close()
    print(f"  Pulse #{info['pulse_number']+1} selected "
          f"(history index {info['pulse_hist_idx']}), "
          f"window {t_win_kyr[0]:.1f}–{t_win_kyr[1]:.1f} kyr")


# ---------------------------------------------------------------------------
# ANIMATION  –  single pulse traversal  +  live SED
# ---------------------------------------------------------------------------

def make_single_pulse_movie(history, sed_dir, output_dir,
                             filename='fig_tpagb_single_pulse_movie.mp4',
                             fps=20, skip_first=5):
    """
    Animate the traversal of one thermal pulse.

    Left column : Two stacked panels (log_LHe top, V-band bottom) with a
                  shared red cursor tracking the current model.
    Right panel : SED for the current model, updating each frame.
    """
    info = select_single_pulse(history, skip_first=skip_first)
    if info is None:
        print("Warning: Not enough thermal pulses for single-pulse movie. Skipping.")
        return

    log_LHe  = history['log_LHe']
    age      = history['star_age']
    time_kyr = age / 1e3

    t_win_yr  = info['t_window']
    t_win_kyr = (t_win_yr[0] / 1e3, t_win_yr[1] / 1e3)

    # V-band
    mag_cols = find_magnitude_columns(history, 'johnson')
    has_V    = 'V' in mag_cols
    if has_V:
        V_mag = history[mag_cols['V']]

    # ---- Load SED frames within the time window ----------------------------
    print("  Loading SED files for single-pulse movie...")
    all_seds = load_all_seds(sed_dir, preferred_prefix='V')
    if not all_seds:
        print("  Warning: No SED files found. Skipping single-pulse movie.")
        return

    model_nums_hist = history['model_number'].astype(int)
    model_to_idx    = {int(m): i for i, m in enumerate(model_nums_hist)}

    pulse_seds = [
        (mn, d) for mn, d in all_seds
        if mn in model_to_idx
        and t_win_yr[0] <= age[model_to_idx[mn]] <= t_win_yr[1]
    ]

    if not pulse_seds:
        print("  Warning: No SED frames found inside the pulse window. Skipping.")
        return

    print(f"  {len(pulse_seds)} frames in pulse window")

    # ---- SED display range -------------------------------------------------
    wl_min, wl_max = 3000, 30000
    fl_peaks = []
    for _, d in pulse_seds:
        wl = d['wavelengths']
        fl = d['fluxes']
        m  = (wl >= wl_min) & (wl <= wl_max) & np.isfinite(fl)
        if m.any() and fl[m].max() > 0:
            fl_peaks.append(fl[m].max())
    fl_max = (np.percentile(fl_peaks, 99) * 1.15) if fl_peaks else 1.0

    # ---- Phase reference lines (shared across both left panels) ------------
    phases = [
        ('lead_in',  'Lead-in',    '#4878d0'),
        ('knee',     'Onset',      '#ee854a'),
        ('peak',     'Peak',       '#d62728'),
        ('return_q', 'Quiescence', '#6acc65'),
    ]

    win_mask  = (time_kyr >= t_win_kyr[0]) & (time_kyr <= t_win_kyr[1])
    valid_lhe = win_mask & np.isfinite(log_LHe)

    # ---- Build figure -------------------------------------------------------
    fig = plt.figure(figsize=(APJ_DOUBLE_COL, 4.5))
    gs  = fig.add_gridspec(1, 2, width_ratios=[1, 1.5], wspace=0.32,
                            left=0.09, right=0.97, top=0.88, bottom=0.12)

    gs_left = gs[0].subgridspec(2, 1, hspace=0.08, height_ratios=[1, 1])
    ax_lhe  = fig.add_subplot(gs_left[0])
    ax_V    = fig.add_subplot(gs_left[1], sharex=ax_lhe)
    ax_sed  = fig.add_subplot(gs[1])

    # ---- Static: log_LHe background -----------------------------------------
    ax_lhe.plot(time_kyr[valid_lhe], log_LHe[valid_lhe], 'k-', lw=1.0, zorder=2)
    for phase_key, label, color in phases:
        idx = info[phase_key]
        ax_lhe.axvline(time_kyr[idx], color=color, lw=0.9, ls='--',
                       alpha=0.75, zorder=3)
        ax_lhe.plot(time_kyr[idx], log_LHe[idx], 'o', color=color,
                    ms=4, zorder=4, markeredgecolor='k', markeredgewidth=0.4)
    ax_lhe.set_xlim(t_win_kyr)
    ax_lhe.set_ylabel(r'$\log(L_{\rm He}/L_\odot)$', fontsize=7)
    ax_lhe.set_title(f'Pulse #{info["pulse_number"]+1}', fontsize=9)
    plt.setp(ax_lhe.get_xticklabels(), visible=False)

    legend_handles = [
        Line2D([0], [0], color=c, lw=1.1, ls='--', marker='o',
               markersize=4, markeredgecolor='k', markeredgewidth=0.3, label=lbl)
        for _, lbl, c in phases
    ]
    ax_lhe.legend(handles=legend_handles, fontsize=5.5,
                  loc='upper left', framealpha=0.9)

    # ---- Static: V-band background ------------------------------------------
    if has_V:
        valid_V = win_mask & np.isfinite(V_mag)
        ax_V.plot(time_kyr[valid_V], V_mag[valid_V],
                  color=JOHNSON_FILTERS['V']['color'], lw=1.0, zorder=2)
        for phase_key, label, color in phases:
            idx = info[phase_key]
            ax_V.axvline(time_kyr[idx], color=color, lw=0.9, ls='--',
                         alpha=0.75, zorder=3)
            v_pt = V_mag[idx]
            if np.isfinite(v_pt):
                ax_V.plot(time_kyr[idx], v_pt, 'o', color=color,
                          ms=4, zorder=4, markeredgecolor='k', markeredgewidth=0.4)
        ax_V.invert_yaxis()
        ax_V.set_ylabel(r'$V$ (mag)', fontsize=7)
    else:
        ax_V.text(0.5, 0.5, 'V-band not available',
                  ha='center', va='center', transform=ax_V.transAxes, fontsize=7)
        for phase_key, label, color in phases:
            ax_V.axvline(time_kyr[info[phase_key]], color=color,
                         lw=0.9, ls='--', alpha=0.75, zorder=3)

    ax_V.set_xlim(t_win_kyr)
    ax_V.set_xlabel('Time (kyr)', fontsize=7)

    # ---- Static SED axes decoration -----------------------------------------
    ax_sed.set_xscale('log')
    ax_sed.set_xlim(wl_min, wl_max)
    ax_sed.set_ylim(0, fl_max)
    ax_sed.set_xlabel(r'Wavelength ($\mathrm{\AA}$)', fontsize=8)
    ax_sed.set_ylabel(r'$F_\lambda$', fontsize=8)
    ax_sed.set_title('SED', fontsize=9)
    for filt, props in JOHNSON_FILTERS.items():
        if filt in ['B', 'V', 'R', 'I']:
            ax_sed.axvline(props['wavelength'], color=props['color'],
                           alpha=0.35, lw=0.9, zorder=0)
    add_em_spectrum_regions(ax_sed, alpha=0.03)

    # ---- Dynamic elements ---------------------------------------------------
    # Cursors on both left panels
    cursor_lhe, = ax_lhe.plot([], [], 'r|', ms=14, mew=2.0, zorder=6)
    cursor_V,   = ax_V.plot([], [], 'r|', ms=14, mew=2.0, zorder=6)
    sed_line,   = ax_sed.plot([], [], 'k-', lw=0.9)
    info_txt    = ax_sed.text(
        0.97, 0.95, '', transform=ax_sed.transAxes,
        fontsize=7, ha='right', va='top', family='monospace',
        bbox=dict(boxstyle='round,pad=0.25', facecolor='white',
                  alpha=0.8, edgecolor='none'))

    fig.suptitle('TP-AGB: Thermal Pulse Anatomy', fontsize=10)

    def init():
        cursor_lhe.set_data([], [])
        cursor_V.set_data([], [])
        sed_line.set_data([], [])
        info_txt.set_text('')
        return cursor_lhe, cursor_V, sed_line, info_txt

    def update(frame_idx):
        mn, d = pulse_seds[frame_idx]
        hidx  = model_to_idx.get(mn)

        if hidx is not None:
            t_now  = time_kyr[hidx]
            lhe_now = (log_LHe[hidx]
                       if np.isfinite(log_LHe[hidx])
                       else np.nanmean(log_LHe[valid_lhe]))
            cursor_lhe.set_data([t_now], [lhe_now])

            if has_V:
                v_now = V_mag[hidx]
                if np.isfinite(v_now):
                    cursor_V.set_data([t_now], [v_now])
                else:
                    cursor_V.set_data([], [])
            else:
                cursor_V.set_data([], [])

        wl = d['wavelengths']
        fl = d['fluxes']
        m  = (wl >= wl_min) & (wl <= wl_max)
        sed_line.set_data(wl[m], fl[m])

        if hidx is not None:
            teff  = (10**history['log_Teff'][hidx]
                     if 'log_Teff' in history else float('nan'))
            log_L = (history['log_L'][hidx]
                     if 'log_L'    in history else float('nan'))
            t_now = time_kyr[hidx]
            info_txt.set_text(
                f"t = {t_now:.2f} kyr\n"
                f"Teff = {teff:.0f} K\n"
                f"log L = {log_L:.2f}\n"
                f"Model {mn}"
            )

        return cursor_lhe, cursor_V, sed_line, info_txt

    n_frames = len(pulse_seds)
    print(f"  Rendering {n_frames} frames at {fps} fps...")

    ani = animation.FuncAnimation(fig, update, frames=n_frames,
                                  init_func=init, blit=True,
                                  interval=1000 / fps)

    out_path = output_dir / filename
    writer   = animation.FFMpegWriter(fps=fps, bitrate=2200,
                                       metadata={'title': 'TP-AGB Single Pulse'})
    try:
        ani.save(str(out_path), writer=writer)
        print(f"  Saved: {out_path}")
    except Exception as e:
        print(f"  Warning: FFMpegWriter failed ({e}). Trying PillowWriter...")
        gif_path = output_dir / filename.replace('.mp4', '.gif')
        ani.save(str(gif_path), writer=animation.PillowWriter(fps=fps))
        print(f"  Saved: {gif_path}")

    plt.close()


# =============================================================================
# SUMMARY PRINT
# =============================================================================

def print_summary(history):
    """Print TP-AGB evolution summary."""
    pulse_idx, pulse_times = identify_thermal_pulses(history)

    print("\n" + "=" * 60)
    print("TP-AGB EVOLUTION SUMMARY")
    print("=" * 60)

    print(f"Total models: {len(history['model_number'])}")
    print(f"Thermal pulses identified: {len(pulse_idx)}")

    if 'log_Teff' in history:
        teff = 10**history['log_Teff']
        print(f"Teff range: {teff.min():.0f} - {teff.max():.0f} K")

    if 'log_L' in history:
        log_L = history['log_L']
        print(f"log(L/Lsun) range: {log_L.min():.2f} - {log_L.max():.2f}")

    mag_cols = find_magnitude_columns(history, 'johnson')
    if 'V' in mag_cols:
        V = history[mag_cols['V']]
        valid = np.isfinite(V)
        if valid.any():
            print(f"V magnitude range: {V[valid].min():.2f} - {V[valid].max():.2f}")

    if len(pulse_times) > 1:
        interpulse = np.diff(pulse_times) / 1e3
        print(f"Mean interpulse period: {np.mean(interpulse):.1f} kyr")

    print("=" * 60 + "\n")


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Generate MESA Colors TP-AGB figures'
    )
    parser.add_argument('--logs_dir', type=str, default='../LOGS',
                       help='Path to MESA LOGS directory')
    parser.add_argument('--sed_dir', type=str, default='../SED',
                       help='Path to SED output directory')
    parser.add_argument('--output_dir', type=str, default='.',
                       help='Output directory for figures')
    parser.add_argument('--movie_fps', type=int, default=15,
                       help='Movie frames per second')
    args = parser.parse_args()

    setup_apj_style()

    logs_dir = Path(args.logs_dir)
    sed_dir = Path(args.sed_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)

    # Find history file
    history_file = logs_dir / 'history.data'
    if not history_file.exists():
        for alt in ['LOGS_TP/history.data', '../LOGS/history.data']:
            if Path(alt).exists():
                history_file = Path(alt)
                break

    if not history_file.exists():
        print(f"ERROR: {history_file} not found")
        return

    print(f"Reading: {history_file}")
    history = read_mesa_history(history_file)

    print_summary(history)

    # Determine TP-phase x-limits
    tp_xlim = find_tp_xlim(history)
    if tp_xlim:
        print(f"TP phase range: {tp_xlim[0]:.0f} - {tp_xlim[1]:.0f} kyr")
    else:
        print("Warning: Could not determine TP phase range; cropped versions will be skipped.")

    print("\nGenerating figures (full range)...")

    # --- Full range ---
    plot_lightcurves(history, output_dir,
                     basename='fig_tpagb_lightcurves')
    plot_diagnostics(history, output_dir,
                     basename='fig_tpagb_diagnostics')

    if sed_dir.exists():
        plot_cmd_with_seds(history, sed_dir, output_dir,
                           basename='fig_tpagb_cmd')
        plot_summary(history, sed_dir, output_dir,
                     basename='fig_tpagb_summary')
    else:
        print(f"Warning: SED directory {sed_dir} not found – skipping CMD and summary.")

    # --- TP-phase cropped ---
    if tp_xlim:
        print("\nGenerating figures (TP phase only)...")

        plot_lightcurves(history, output_dir,
                         basename='fig_tpagb_lightcurves_tp', xlim=tp_xlim)
        plot_diagnostics(history, output_dir,
                         basename='fig_tpagb_diagnostics_tp', xlim=tp_xlim)

        if sed_dir.exists():
            plot_cmd_with_seds(history, sed_dir, output_dir,
                               basename='fig_tpagb_cmd_tp', xlim=tp_xlim)
            plot_summary(history, sed_dir, output_dir,
                         basename='fig_tpagb_summary_tp', xlim=tp_xlim)

    # --- Single-pulse anatomy figure + movie ---
    if sed_dir.exists():
        print("\nGenerating single-pulse anatomy figure...")
        plot_single_pulse_with_seds(history, sed_dir, output_dir,
                                    basename='fig_tpagb_single_pulse')

        print("\nGenerating single-pulse anatomy movie...")
        make_single_pulse_movie(history, sed_dir, output_dir,
                                filename='fig_tpagb_single_pulse_movie.mp4',
                                fps=args.movie_fps)
    else:
        print(f"Warning: SED directory {sed_dir} not found - skipping single-pulse outputs.")

    # --- SED Movie ---
    if sed_dir.exists():
        print("\nGenerating SED movie (full range)...")
        make_sed_movie(history, sed_dir, output_dir,
                       filename='fig_tpagb_sed_movie.mp4',
                       fps=args.movie_fps, tp_only=False)

        if tp_xlim:
            print("\nGenerating SED movie (TP phase only)...")
            make_sed_movie(history, sed_dir, output_dir,
                           filename='fig_tpagb_sed_movie_tp.mp4',
                           fps=args.movie_fps, tp_only=True)
    else:
        print(f"Warning: SED directory {sed_dir} not found – skipping movie.")

    print(f"\nDone! All outputs saved to: {output_dir}")


if __name__ == '__main__':
    main()