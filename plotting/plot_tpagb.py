#!/usr/bin/env python3
"""
plot_tpagb.py

Analysis and plotting for MESA Colors TP-AGB (Thermal Pulse AGB) demonstration.
Generates publication-quality figures showing photometric evolution during 
thermal pulses, with SED snapshots at key pulse phases.

Figures produced:
- fig_tpagb_lightcurves.pdf : Multi-band light curves with pulse markers
- fig_tpagb_cmd.pdf         : CMD with SED insets at pulse phases
- fig_tpagb_diagnostics.pdf : He luminosity, core mass, Teff evolution
- fig_tpagb_summary.pdf     : 4-panel summary figure

Usage:
    python plot_tpagb.py [--logs_dir LOGS] [--sed_dir SED] [--output_dir figures]

Author: Miller, Joyce, Mocz et al.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path
import argparse
import csv

from scipy.signal import find_peaks

from plot_utils import (
    setup_apj_style, read_mesa_history, find_magnitude_columns,
    get_model_number_at_index, get_sed_counter_at_index,
    add_em_spectrum_regions,
    JOHNSON_FILTERS, APJ_SINGLE_COL, APJ_DOUBLE_COL
)


def load_sed_csv(filepath):
    """Load SED data from CSV file."""
    try:
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f)
            data = {'wavelengths': [], 'fluxes': [], 'convolved_flux': []}
            
            for row in reader:
                try:
                    wl = float(row['wavelengths'])
                    fl = float(row['fluxes'])
                    # Only keep rows where both are valid and non-zero
                    if wl > 0 and fl > 0:
                        data['wavelengths'].append(wl)
                        data['fluxes'].append(fl)
                        try:
                            cf = float(row.get('convolved_flux', 0))
                            data['convolved_flux'].append(cf)
                        except (ValueError, TypeError):
                            data['convolved_flux'].append(0)
                except (ValueError, KeyError):
                    continue
            
            for key in data:
                data[key] = np.array(data[key])
            
            return data
    except Exception as e:
        print(f"Warning: Could not read {filepath}: {e}")
        return {}


def identify_thermal_pulses(history, min_prominence=0.5):
    """
    Identify thermal pulse events from He-burning luminosity peaks.
    
    Parameters
    ----------
    history : dict
        MESA history data
    min_prominence : float
        Minimum peak prominence in log_LHe
        
    Returns
    -------
    tuple
        (pulse_indices, pulse_times) arrays
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
    
    Returns indices for:
    - pre_pulse: Just before first major pulse
    - pulse_peak: At maximum He luminosity
    - post_pulse: During interpulse phase
    - late_agb: Near end of AGB
    
    Parameters
    ----------
    history : dict
        MESA history data
        
    Returns
    -------
    dict
        Named indices for points of interest
    """
    pulse_idx, pulse_times = identify_thermal_pulses(history)
    
    if len(pulse_idx) < 2:
        # Not enough pulses, return basic points
        n = len(history['model_number'])
        return {
            'start': 0,
            'middle': n // 2,
            'end': n - 1
        }
    
    # Select specific pulses for SED display
    poi = {}
    
    # Just before first pulse
    if pulse_idx[0] > 10:
        poi['pre_pulse'] = pulse_idx[0] - 10
    else:
        poi['pre_pulse'] = 0
    
    # At pulse peak (use middle pulse for cleaner display)
    mid_pulse = len(pulse_idx) // 2
    poi['pulse_peak'] = pulse_idx[mid_pulse]
    
    # Post-pulse (interpulse phase)
    if mid_pulse + 1 < len(pulse_idx):
        interpulse_idx = (pulse_idx[mid_pulse] + pulse_idx[mid_pulse + 1]) // 2
        poi['interpulse'] = interpulse_idx
    
    # Late AGB
    poi['late_agb'] = len(history['model_number']) - 1
    
    return poi


def plot_lightcurves(history, output_dir, filename='fig_tpagb_lightcurves.pdf'):
    """
    Plot multi-band light curves during thermal pulses.
    
    Parameters
    ----------
    history : dict
        MESA history data
    output_dir : Path
        Output directory
    filename : str
        Output filename
    """
    fig, axes = plt.subplots(2, 1, figsize=(APJ_DOUBLE_COL, 4.5), sharex=True)
    
    age = history['star_age']
    time_kyr = age / 1e3
    
    pulse_idx, pulse_times = identify_thermal_pulses(history)
    
    # Top panel: Filter magnitudes
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
    
    # Mark thermal pulses
    for pt in pulse_times:
        ax1.axvline(pt/1e3, color='gray', ls=':', alpha=0.4, lw=0.5)
    
    ax1.set_ylabel('Magnitude')
    ax1.invert_yaxis()
    ax1.legend(loc='upper right', ncol=4, fontsize=7, framealpha=0.9)
    ax1.set_title(f'TP-AGB Evolution ({len(pulse_idx)} thermal pulses)', fontsize=10)
    
    # Bottom panel: Bolometric
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
    
    plt.tight_layout()
    plt.savefig(output_dir / filename, dpi=300)
    plt.close()
    print(f"Saved: {output_dir / filename}")


def plot_cmd_with_seds(history, sed_dir, output_dir, 
                       filename='fig_tpagb_cmd.pdf'):
    """
    Plot CMD with SED insets at thermal pulse phases.
    
    Parameters
    ----------
    history : dict
        MESA history data
    sed_dir : Path
        SED directory
    output_dir : Path
        Output directory
    filename : str
        Output filename
    """
    fig = plt.figure(figsize=(APJ_DOUBLE_COL, 4.5))
    
    gs = GridSpec(2, 3, width_ratios=[1.3, 1, 1], wspace=0.3, hspace=0.3,
                  left=0.08, right=0.97, top=0.92, bottom=0.12)
    
    ax_cmd = fig.add_subplot(gs[:, 0])
    
    # Get magnitudes
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
    
    # Plot CMD
    scatter = ax_cmd.scatter(color[valid], mag[valid], c=time_kyr[valid],
                            cmap='viridis', s=3, alpha=0.7, rasterized=True)
    
    cbar = fig.colorbar(scatter, ax=ax_cmd, pad=0.02, aspect=25)
    cbar.set_label('Time (kyr)', fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    
    # POI markers
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
            sed_counter = get_sed_counter_at_index(history, idx)
            
            sed_files = list(sed_dir.glob(f'*_SED_{sed_counter}.csv'))
            
            if sed_files:
                sed_data = load_sed_csv(sed_files[0])
                
                if sed_data and 'wavelengths' in sed_data:
                    wl = sed_data['wavelengths']
                    fl = sed_data['fluxes']
                    
                    ax_sed.plot(wl, fl, 'k-', lw=0.7)
                    ax_sed.set_xscale('log')
                    
                    # Filter markers
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
            
            # Stellar params
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
    
    plt.savefig(output_dir / filename, dpi=300)
    plt.close()
    print(f"Saved: {output_dir / filename}")


def plot_diagnostics(history, output_dir, filename='fig_tpagb_diagnostics.pdf'):
    """
    Plot TP-AGB diagnostic panels.
    
    Parameters
    ----------
    history : dict
        MESA history data
    output_dir : Path
        Output directory
    filename : str
        Output filename
    """
    fig, axes = plt.subplots(4, 1, figsize=(APJ_SINGLE_COL, 6), sharex=True)
    
    age = history['star_age']
    time_kyr = age / 1e3
    
    pulse_idx, pulse_times = identify_thermal_pulses(history)
    
    # Panel 1: He luminosity
    ax1 = axes[0]
    if 'log_LHe' in history:
        log_LHe = history['log_LHe']
        valid = np.isfinite(log_LHe)
        ax1.plot(time_kyr[valid], log_LHe[valid], 'r-', lw=0.7)
        
        for pi in pulse_idx:
            ax1.plot(time_kyr[pi], log_LHe[pi], 'ko', ms=3)
    
    ax1.set_ylabel(r'$\log(L_{\rm He}/L_\odot)$', fontsize=8)
    ax1.set_title(f'{len(pulse_idx)} thermal pulses', fontsize=9)
    
    # Panel 2: Teff
    ax2 = axes[1]
    if 'log_Teff' in history:
        log_Teff = history['log_Teff']
        ax2.plot(time_kyr, 10**log_Teff, 'k-', lw=0.7)
    ax2.set_ylabel(r'$T_{\rm eff}$ (K)', fontsize=8)
    
    # Panel 3: Core masses
    ax3 = axes[2]
    if 'he_core_mass' in history:
        ax3.plot(time_kyr, history['he_core_mass'], 'k-', lw=0.8, label='He core')
    if 'c_core_mass' in history:
        ax3.plot(time_kyr, history['c_core_mass'], 'b--', lw=0.6, label='C/O core')
    ax3.set_ylabel(r'Core mass ($M_\odot$)', fontsize=8)
    ax3.legend(loc='upper left', fontsize=6, framealpha=0.9)
    
    # Panel 4: V magnitude
    ax4 = axes[3]
    mag_cols = find_magnitude_columns(history, 'johnson')
    if 'V' in mag_cols:
        V = history[mag_cols['V']]
        valid = np.isfinite(V)
        ax4.plot(time_kyr[valid], V[valid], 'g-', lw=0.7)
    ax4.set_ylabel('$V$ (mag)', fontsize=8)
    ax4.invert_yaxis()
    ax4.set_xlabel('Time (kyr)')
    
    # Mark pulses on all panels
    for ax in axes:
        for pt in pulse_times:
            ax.axvline(pt/1e3, color='gray', ls=':', alpha=0.3, lw=0.5)
    
    plt.tight_layout()
    plt.savefig(output_dir / filename, dpi=300)
    plt.close()
    print(f"Saved: {output_dir / filename}")


def plot_summary(history, sed_dir, output_dir, filename='fig_tpagb_summary.pdf'):
    """
    Create 4-panel summary figure with SED evolution.
    
    Parameters
    ----------
    history : dict
        MESA history data
    sed_dir : Path
        SED directory
    output_dir : Path
        Output directory
    filename : str
        Output filename
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
    
    # (b) CMD
    ax2 = fig.add_subplot(gs[0, 1])
    
    if 'B' in mag_cols and 'V' in mag_cols:
        B = history[mag_cols['B']]
        V = history[mag_cols['V']]
        BmV = B - V
        valid = np.isfinite(BmV) & np.isfinite(V)
        
        scatter = ax2.scatter(BmV[valid], V[valid], c=time_kyr[valid],
                             cmap='viridis', s=2, alpha=0.6, rasterized=True)
        ax2.plot(BmV[valid][0], V[valid][0], 'go', ms=6, zorder=10)
        ax2.plot(BmV[valid][-1], V[valid][-1], 'rs', ms=6, zorder=10)
    
    ax2.set_xlabel(r'$B - V$')
    ax2.set_ylabel('$V$')
    ax2.invert_yaxis()
    ax2.set_title('(b) Color-Magnitude Diagram', fontsize=9)
    
    # (c) HR diagram
    ax3 = fig.add_subplot(gs[1, 0])
    
    if 'log_Teff' in history and 'log_L' in history:
        log_Teff = history['log_Teff']
        log_L = history['log_L']
        
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
    cmap = plt.cm.viridis
    
    for i, (poi_name, label) in enumerate([
        ('pre_pulse', 'Pre-pulse'),
        ('pulse_peak', 'Pulse peak'),
        ('interpulse', 'Interpulse'),
        ('late_agb', 'Late AGB')
    ]):
        if poi_name not in poi:
            continue
        
        idx = poi[poi_name]
        sed_counter = get_sed_counter_at_index(history, idx)
        
        sed_files = list(sed_dir.glob(f'*_SED_{sed_counter}.csv'))
        if sed_files:
            sed_data = load_sed_csv(sed_files[0])
            if sed_data and 'wavelengths' in sed_data:
                wl = sed_data['wavelengths']
                fl = sed_data['fluxes']
                
                color = cmap(i / 3)
                ax4.plot(wl, fl, color=color, lw=0.8, label=label, alpha=0.9)
    
    ax4.set_xscale('log')
    ax4.set_xlabel(r'Wavelength ($\mathrm{\AA}$)')
    ax4.set_ylabel(r'$F_\lambda$')
    ax4.set_xlim(3000, 30000)
    ax4.legend(loc='upper right', fontsize=6, framealpha=0.9)
    ax4.set_title('(d) SED Evolution', fontsize=9)
    
    add_em_spectrum_regions(ax4, alpha=0.03)
    
    fig.suptitle('MESA Colors: Thermal Pulse AGB Evolution', fontsize=11)
    
    plt.savefig(output_dir / filename, dpi=300)
    plt.close()
    print(f"Saved: {output_dir / filename}")


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
        interpulse = np.diff(pulse_times) / 1e3  # kyr
        print(f"Mean interpulse period: {np.mean(interpulse):.1f} kyr")
    
    print("=" * 60 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description='Generate MESA Colors TP-AGB figures'
    )
    parser.add_argument('--logs_dir', type=str, default='LOGS',
                       help='Path to MESA LOGS directory')
    parser.add_argument('--sed_dir', type=str, default='SED',
                       help='Path to SED output directory')
    parser.add_argument('--output_dir', type=str, default='figures',
                       help='Output directory for figures')
    args = parser.parse_args()
    
    setup_apj_style()
    
    logs_dir = Path(args.logs_dir)
    sed_dir = Path(args.sed_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Find history file
    history_file = logs_dir / 'history.data'
    if not history_file.exists():
        # Try alternative locations
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
    
    print("Generating figures...")
    
    plot_lightcurves(history, output_dir)
    
    if sed_dir.exists():
        plot_cmd_with_seds(history, sed_dir, output_dir)
        plot_summary(history, sed_dir, output_dir)
    else:
        print(f"Warning: SED directory {sed_dir} not found")
    
    plot_diagnostics(history, output_dir)
    
    print(f"\nDone! Figures saved to: {output_dir}")


if __name__ == '__main__':
    main()
