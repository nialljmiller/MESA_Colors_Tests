#!/usr/bin/env python3
"""
plot_rsp.py

Generate figures for MESA Colors RSP (Radial Stellar Pulsation) demonstration.
Shows phase-resolved photometry for RR Lyrae / Cepheid pulsations with
SED snapshots at key pulsation phases.

Figures produced:
- fig_rsp_lightcurve.pdf : Multi-band light curves with SED insets
- fig_rsp_colorloop.pdf  : Color-magnitude loop (B-V vs V)
- fig_rsp_sed_phases.pdf : SED comparison at different phases

Usage:
    python plot_rsp.py LOGS/history.data [--sed_dir SED] [--output_dir figures]

Author: Miller, Joyce, Mocz et al.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from pathlib import Path
import argparse
import csv

from plot_utils import (
    setup_apj_style, read_mesa_history, find_magnitude_columns,
    find_poi_rsp, get_model_number_at_index, add_em_spectrum_regions,
    JOHNSON_FILTERS, APJ_SINGLE_COL, APJ_DOUBLE_COL
)


def plot_lightcurve_with_seds(history, sed_dir, output_dir,
                               filename='fig_rsp_lightcurve.pdf'):
    """
    Plot multi-band light curves with SED insets at key phases.
    
    Parameters
    ----------
    history : dict
        MESA history data
    sed_dir : Path
        Directory containing SED files
    output_dir : Path
        Output directory
    filename : str
        Output filename
    """
    fig = plt.figure(figsize=(APJ_DOUBLE_COL, 4.5))
    
    gs = GridSpec(2, 3, width_ratios=[1.5, 1, 1], height_ratios=[1, 1],
                  wspace=0.35, hspace=0.25,
                  left=0.08, right=0.97, top=0.92, bottom=0.12)
    
    # Main light curve panel (spans both rows on left)
    ax_lc = fig.add_subplot(gs[:, 0])
    
    # Get phase
    if 'rsp_phase' in history:
        phase = history['rsp_phase']
    else:
        time = history.get('star_age_day', history['star_age'])
        period = 0.71  # days, approximate
        phase = ((time - time[0]) / period) % 1.0
    
    # Find magnitude columns
    mag_cols = find_magnitude_columns(history, 'johnson')
    
    # Plotting parameters
    band_info = {
        'B': {'color': '#0066CC', 'offset': 0.0},
        'V': {'color': '#33AA33', 'offset': 1.0},
        'R': {'color': '#CC6600', 'offset': 2.0},
        'I': {'color': '#990033', 'offset': 3.0},
    }
    
    # Select data after settling (skip first few cycles)
    if 'rsp_num_periods' in history:
        mask = (history['rsp_num_periods'] >= 4) & (history['rsp_num_periods'] < 6)
    else:
        n = len(phase)
        mask = np.arange(n) > n // 3
    
    if mask.sum() < 50:
        mask = np.ones(len(phase), dtype=bool)
    
    # Plot each band
    for band in ['B', 'V', 'R', 'I']:
        if band in mag_cols:
            col = mag_cols[band]
            mag = history[col][mask]
            ph = phase[mask]
            info = band_info[band]
            
            sort_idx = np.argsort(ph)
            ax_lc.plot(ph[sort_idx], mag[sort_idx] + info['offset'],
                      '.', ms=1.5, alpha=0.7, color=info['color'], 
                      label=f"${band}$", rasterized=True)
    
    ax_lc.set_xlabel('Pulsation Phase')
    ax_lc.set_ylabel('Magnitude + offset')
    ax_lc.set_xlim(0, 1)
    ax_lc.invert_yaxis()
    ax_lc.legend(loc='upper right', markerscale=4, fontsize=8, framealpha=0.9)
    
    # Annotation about offsets
    ax_lc.text(0.02, 0.02, 'Offsets: +1.0 mag/band',
              transform=ax_lc.transAxes, fontsize=7, alpha=0.6)
    
    # Find POIs for SED panels
    poi = find_poi_rsp(history)
    
    # SED panels (2x2 grid on right)
    sed_positions = [(0, 1), (0, 2), (1, 1), (1, 2)]
    sed_pois = ['maximum_light', 'rising', 'minimum_light', 'falling']
    sed_titles = ['Maximum Light', 'Rising', 'Minimum Light', 'Falling']
    phase_labels = [r'$\phi \approx 0$', r'$\phi \approx 0.25$', 
                    r'$\phi \approx 0.5$', r'$\phi \approx 0.75$']
    
    for (row, col), poi_name, title, ph_label in zip(sed_positions, sed_pois, 
                                                       sed_titles, phase_labels):
        ax_sed = fig.add_subplot(gs[row, col])
        
        if poi_name in poi:
            idx = poi[poi_name]
            model_num = get_model_number_at_index(history, idx)
            
            # Load SED
            sed_files = list(sed_dir.glob(f'*_SED_{model_num}.csv'))
            
            if sed_files:
                sed_data = load_sed_csv(sed_files[0])
                
                if sed_data and len(sed_data.get('wavelengths', [])) > 0:
                    wl = sed_data['wavelengths']
                    fl = sed_data['fluxes']
                    
                    ax_sed.plot(wl, fl, 'k-', lw=0.7)
                    ax_sed.set_xscale('log')
                    
                    # Filter markers
                    for filt, props in JOHNSON_FILTERS.items():
                        if filt in ['B', 'V', 'R', 'I']:
                            ax_sed.axvline(props['wavelength'], 
                                          color=props['color'],
                                          alpha=0.5, lw=1.2, zorder=0)
                    
                    ax_sed.set_xlim(3500, 10000)
                    
                    mask_vis = (wl > 3500) & (wl < 10000)
                    if mask_vis.any():
                        fl_vis = fl[mask_vis]
                        ax_sed.set_ylim(0, fl_vis.max() * 1.15)
                    
                    add_em_spectrum_regions(ax_sed, alpha=0.03)
            
            # Mark this phase on light curve
            ax_lc.axvline(phase[idx] % 1, color='gray', ls=':', alpha=0.4, lw=0.7)
        
        ax_sed.set_title(f'{title}\n{ph_label}', fontsize=8)
        ax_sed.tick_params(labelsize=6)
        
        if row == 1:
            ax_sed.set_xlabel(r'$\lambda$ ($\mathrm{\AA}$)', fontsize=7)
        else:
            ax_sed.set_xticklabels([])
        
        if col == 1:
            ax_sed.set_ylabel(r'$F_\lambda$', fontsize=7)
    
    plt.savefig(output_dir / filename, dpi=300)
    plt.close()
    print(f"Saved: {output_dir / filename}")


def plot_colorloop(history, sed_dir, output_dir, 
                   filename='fig_rsp_colorloop.pdf'):
    """
    Plot color-magnitude loop (B-V vs V) with SED comparison inset.
    
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
    fig = plt.figure(figsize=(APJ_SINGLE_COL * 1.3, 4))
    
    gs = GridSpec(2, 1, height_ratios=[2, 1], hspace=0.3,
                  left=0.15, right=0.95, top=0.95, bottom=0.12)
    
    ax_loop = fig.add_subplot(gs[0])
    ax_sed = fig.add_subplot(gs[1])
    
    # Get magnitudes
    mag_cols = find_magnitude_columns(history, 'johnson')
    
    if 'B' not in mag_cols or 'V' not in mag_cols:
        print("Error: Need B and V for color loop")
        return
    
    B = history[mag_cols['B']]
    V = history[mag_cols['V']]
    color = B - V
    
    # Get phase
    if 'rsp_phase' in history:
        phase = history['rsp_phase']
    else:
        time = history.get('star_age_day', history['star_age'])
        period = 0.71
        phase = ((time - time[0]) / period) % 1.0
    
    # Select one cycle
    if 'rsp_num_periods' in history:
        mask = history['rsp_num_periods'] >= 4
    else:
        mask = np.ones(len(phase), dtype=bool)
    
    # Plot color loop
    sc = ax_loop.scatter(color[mask], V[mask], c=phase[mask], 
                        cmap='twilight_shifted', s=6, alpha=0.8,
                        edgecolors='none', rasterized=True)
    
    cbar = fig.colorbar(sc, ax=ax_loop, pad=0.02, aspect=20)
    cbar.set_label('Phase', fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    
    # Mark max/min light
    poi = find_poi_rsp(history)
    for name, marker in [('maximum_light', '*'), ('minimum_light', 'o')]:
        if name in poi:
            idx = poi[name]
            ax_loop.plot(color[idx], V[idx], marker, ms=10, mfc='white',
                        mec='black', mew=1.5, zorder=10)
    
    # Arrow for loop direction
    c_arr = color[mask]
    v_arr = V[mask]
    p_arr = phase[mask]
    
    idx1 = np.argmin(np.abs(p_arr - 0.2))
    idx2 = np.argmin(np.abs(p_arr - 0.3))
    
    if idx1 != idx2:
        ax_loop.annotate('', xy=(c_arr[idx2], v_arr[idx2]),
                        xytext=(c_arr[idx1], v_arr[idx1]),
                        arrowprops=dict(arrowstyle='->', color='black', lw=1))
    
    ax_loop.set_xlabel(r'$B - V$ (mag)')
    ax_loop.set_ylabel(r'$V$ (mag)')
    ax_loop.invert_yaxis()
    
    # SED comparison panel: max vs min light
    if 'maximum_light' in poi and 'minimum_light' in poi:
        idx_max = poi['maximum_light']
        idx_min = poi['minimum_light']
        
        model_max = get_model_number_at_index(history, idx_max)
        model_min = get_model_number_at_index(history, idx_min)
        
        for model_num, label, color_line, ls in [
            (model_max, 'Max light', '#1f77b4', '-'),
            (model_min, 'Min light', '#d62728', '--')
        ]:
            sed_files = list(sed_dir.glob(f'*_SED_{model_num}.csv'))
            if sed_files:
                sed_data = load_sed_csv(sed_files[0])
                if sed_data and 'wavelengths' in sed_data:
                    wl = sed_data['wavelengths']
                    fl = sed_data['fluxes']
                    
                    # Normalize for comparison
                    fl_norm = fl / fl.max()
                    
                    ax_sed.plot(wl, fl_norm, color=color_line, ls=ls, 
                               lw=0.9, label=label, alpha=0.9)
    
    ax_sed.set_xscale('log')
    ax_sed.set_xlim(3500, 10000)
    ax_sed.set_ylim(0, 1.1)
    ax_sed.set_xlabel(r'Wavelength ($\mathrm{\AA}$)')
    ax_sed.set_ylabel(r'Normalized $F_\lambda$')
    ax_sed.legend(loc='upper right', fontsize=7, framealpha=0.9)
    
    add_em_spectrum_regions(ax_sed, alpha=0.03)
    
    plt.savefig(output_dir / filename, dpi=300)
    plt.close()
    print(f"Saved: {output_dir / filename}")


def plot_sed_phase_comparison(history, sed_dir, output_dir,
                               filename='fig_rsp_sed_phases.pdf'):
    """
    Direct SED comparison across pulsation phases.
    
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
    fig, ax = plt.subplots(figsize=(APJ_SINGLE_COL, 3))
    
    poi = find_poi_rsp(history)
    
    # Color map for phases
    cmap = plt.cm.twilight_shifted
    phase_colors = {
        'maximum_light': cmap(0.0),
        'rising': cmap(0.25),
        'minimum_light': cmap(0.5),
        'falling': cmap(0.75),
    }
    phase_labels = {
        'maximum_light': r'$\phi = 0$',
        'rising': r'$\phi = 0.25$',
        'minimum_light': r'$\phi = 0.5$',
        'falling': r'$\phi = 0.75$',
    }
    
    for poi_name in ['maximum_light', 'rising', 'minimum_light', 'falling']:
        if poi_name not in poi:
            continue
        
        idx = poi[poi_name]
        model_num = get_model_number_at_index(history, idx)
        
        sed_files = list(sed_dir.glob(f'*_SED_{model_num}.csv'))
        if sed_files:
            sed_data = load_sed_csv(sed_files[0])
            if sed_data and 'wavelengths' in sed_data:
                wl = sed_data['wavelengths']
                fl = sed_data['fluxes']
                
                ax.plot(wl, fl, color=phase_colors[poi_name], 
                       lw=0.9, label=phase_labels[poi_name], alpha=0.85)
    
    ax.set_xscale('log')
    ax.set_xlabel(r'Wavelength ($\mathrm{\AA}$)')
    ax.set_ylabel(r'$F_\lambda$ (erg s$^{-1}$ cm$^{-2}$ $\mathrm{\AA}^{-1}$)')
    ax.set_xlim(3500, 10000)
    ax.legend(loc='upper right', fontsize=7, framealpha=0.9)
    
    add_em_spectrum_regions(ax, alpha=0.03)
    
    plt.tight_layout()
    plt.savefig(output_dir / filename, dpi=300)
    plt.close()
    print(f"Saved: {output_dir / filename}")


def load_sed_csv(filepath):
    """Load SED data from CSV file."""
    try:
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f)
            data = {'wavelengths': [], 'fluxes': [], 'convolved_flux': []}
            
            for row in reader:
                for key in data:
                    try:
                        val = float(row[key])
                        if val != 0:
                            data[key].append(val)
                    except (ValueError, KeyError):
                        pass
            
            for key in data:
                data[key] = np.array(data[key])
            
            return data
    except Exception as e:
        print(f"Warning: Could not read {filepath}: {e}")
        return {}


def print_statistics(history):
    """Print RSP run statistics."""
    mag_cols = find_magnitude_columns(history, 'johnson')
    
    print("\n" + "=" * 60)
    print("RSP + COLORS RUN STATISTICS")
    print("=" * 60)
    
    print(f"Total models: {len(history['model_number'])}")
    
    if 'rsp_num_periods' in history:
        n_periods = int(np.max(history['rsp_num_periods']))
        print(f"Pulsation cycles: {n_periods}")
    
    if 'rsp_period_in_days' in history:
        periods = history['rsp_period_in_days']
        valid = periods > 0
        if valid.any():
            print(f"Mean period: {np.mean(periods[valid]):.5f} days")
    
    print("\nMagnitude ranges:")
    for band in ['B', 'V', 'R', 'I']:
        if band in mag_cols:
            mag = history[mag_cols[band]]
            amp = np.max(mag) - np.min(mag)
            print(f"  {band}: {np.min(mag):.3f} to {np.max(mag):.3f} "
                  f"(Δ = {amp:.3f} mag)")
    
    if 'B' in mag_cols and 'V' in mag_cols:
        color = history[mag_cols['B']] - history[mag_cols['V']]
        print(f"\n(B-V) range: {np.min(color):.3f} to {np.max(color):.3f}")
    
    print("=" * 60 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description='Generate MESA Colors RSP figures'
    )
    parser.add_argument('history_file', type=str, nargs='?',
                       default='../rsp_RR_Lyrae/LOGS/history.data',
                       help='Path to history.data file')
    parser.add_argument('--sed_dir', type=str, default='../rsp_RR_Lyrae/SED',
                       help='Path to SED output directory')
    parser.add_argument('--output_dir', type=str, default='../rsp_RR_Lyrae/figures',
                       help='Output directory for figures')
    args = parser.parse_args()
    
    setup_apj_style()
    
    history_file = Path(args.history_file)
    sed_dir = Path(args.sed_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    if not history_file.exists():
        print(f"ERROR: {history_file} not found")
        return
    
    print(f"Reading: {history_file}")
    history = read_mesa_history(history_file)
    
    print_statistics(history)
    
    print("Generating figures...")
    
    if sed_dir.exists():
        plot_lightcurve_with_seds(history, sed_dir, output_dir)
        plot_colorloop(history, sed_dir, output_dir)
        plot_sed_phase_comparison(history, sed_dir, output_dir)
    else:
        print(f"Warning: SED directory {sed_dir} not found")
        print("Run with --sed_dir to specify SED location")
    
    print(f"\nDone! Figures saved to: {output_dir}")


if __name__ == '__main__':
    main()
