#!/usr/bin/env python3
"""
plot_blue_loop.py

Analysis and plotting for MESA Colors Blue Loop demonstration.
Generates publication-quality figures showing chromatic evolution during
blue loop phase of a 5 Msun intermediate-mass star, with SED insets at
key evolutionary points.

Figures produced:
- fig_blueloop_cmd.pdf      : CMD track with SED insets at key points
- fig_blueloop_hr.pdf       : HR diagram with color-coded track
- fig_blueloop_evolution.pdf: Multi-panel magnitude/color evolution

Usage:
    python plot_blue_loop.py [--logs_dir LOGS] [--sed_dir SED] [--output_dir figures]

Author: Miller, Joyce, Mocz et al.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path
import argparse

from plot_utils import (
    setup_apj_style, read_mesa_history, find_magnitude_columns,
    get_sed_for_model, find_poi_blue_loop, get_model_number_at_index,
    add_em_spectrum_regions, LSST_FILTERS,
    APJ_SINGLE_COL, APJ_DOUBLE_COL
)


def plot_cmd_with_seds(history, sed_dir, output_dir, 
                       filename='fig_blueloop_cmd.pdf'):
    """
    Plot CMD track with SED insets at key evolutionary points.
    
    Parameters
    ----------
    history : dict
        MESA history data
    sed_dir : Path
        Directory containing SED files
    output_dir : Path
        Output directory for figures
    filename : str
        Output filename
    """
    fig = plt.figure(figsize=(APJ_DOUBLE_COL, 5))
    
    # Main CMD panel on left, SED panels on right
    gs = GridSpec(3, 2, width_ratios=[1.2, 1], wspace=0.3, hspace=0.35,
                  left=0.08, right=0.95, top=0.92, bottom=0.10)
    
    ax_cmd = fig.add_subplot(gs[:, 0])
    
    # Get magnitudes
    g = history['g']
    r = history['r']
    color = g - r
    age = history['star_age'] / 1e6  # Myr
    
    # Plot CMD track
    scatter = ax_cmd.scatter(color, g, c=age, cmap='viridis', s=3, alpha=0.8,
                             rasterized=True)
    ax_cmd.plot(color, g, 'k-', lw=0.3, alpha=0.3, zorder=0)
    
    # Find points of interest
    poi = find_poi_blue_loop(history)
    
    # Define markers and colors for POIs
    poi_style = {
        'rgb_tip': {'marker': 'o', 'color': '#d62728', 'label': 'RGB Tip'},
        'loop_ascending': {'marker': '^', 'color': '#2ca02c', 'label': 'Ascending'},
        'loop_maximum': {'marker': 's', 'color': '#1f77b4', 'label': 'Blue Maximum'},
        'end': {'marker': 'D', 'color': '#9467bd', 'label': 'AGB'},
    }
    
    # Plot POI markers
    for name, idx in poi.items():
        style = poi_style[name]
        ax_cmd.scatter(color[idx], g[idx], c=style['color'], s=80, 
                      marker=style['marker'], edgecolors='black', 
                      linewidths=0.8, zorder=5, label=style['label'])
    
    # Invert y-axis (magnitudes)
    ax_cmd.invert_yaxis()
    ax_cmd.set_xlabel(r'$g - r$ (mag)')
    ax_cmd.set_ylabel(r'$g$ (mag)')
    
    # Add instability strip region
    ax_cmd.axvspan(0.25, 0.65, alpha=0.1, color='blue', zorder=0)
    ax_cmd.text(0.45, ax_cmd.get_ylim()[0] - 0.3, 'IS', fontsize=8, 
               ha='center', color='blue', alpha=0.7)
    
    # Colorbar
    cbar = fig.colorbar(scatter, ax=ax_cmd, pad=0.02, aspect=30)
    cbar.set_label('Age (Myr)', fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    
    ax_cmd.legend(loc='lower left', fontsize=7, framealpha=0.9)
    
    # SED panels on right
    sed_axes = [fig.add_subplot(gs[i, 1]) for i in range(3)]
    
    # Select 3 POIs for SED display
    sed_pois = ['rgb_tip', 'loop_maximum', 'end']
    sed_titles = ['RGB Tip', 'Blue Maximum', 'Early AGB']
    
    for ax_sed, poi_name, title in zip(sed_axes, sed_pois, sed_titles):
        idx = poi[poi_name]
        model_num = get_model_number_at_index(history, idx)
        
        # Try to load SED
        sed_files = list(sed_dir.glob(f'*_SED_{model_num}.csv'))
        
        if sed_files:
            # Read first filter's SED (they all have the same underlying SED)
            sed_data = {}
            with open(sed_files[0], 'r') as f:
                import csv
                reader = csv.DictReader(f)
                for key in ['wavelengths', 'fluxes', 'convolved_flux']:
                    sed_data[key] = []
                
                for row in reader:
                    for key in sed_data:
                        try:
                            val = float(row[key])
                            if val != 0:
                                sed_data[key].append(val)
                        except (ValueError, KeyError):
                            pass
                
                for key in sed_data:
                    sed_data[key] = np.array(sed_data[key])
            
            if len(sed_data.get('wavelengths', [])) > 0:
                wl = sed_data['wavelengths']
                fl = sed_data['fluxes']
                
                # Plot SED
                ax_sed.plot(wl, fl, 'k-', lw=0.8)
                ax_sed.set_xscale('log')
                
                # Add filter bands as colored regions
                for filt, props in LSST_FILTERS.items():
                    wl_center = props['wavelength']
                    ax_sed.axvline(wl_center, color=props['color'], 
                                  alpha=0.5, lw=1.5, zorder=0)
                
                # Set limits to optical/NIR range
                ax_sed.set_xlim(3000, 11000)
                
                # Find flux range in this wavelength range
                mask = (wl > 3000) & (wl < 11000)
                if mask.any():
                    fl_vis = fl[mask]
                    ax_sed.set_ylim(0, fl_vis.max() * 1.1)
                
                add_em_spectrum_regions(ax_sed, alpha=0.03)
        
        # Stellar parameters annotation
        teff = 10**history['log_Teff'][idx]
        logg = history['log_g'][idx]
        ax_sed.text(0.97, 0.95, f'$T_{{\\rm eff}}$ = {teff:.0f} K\n$\\log g$ = {logg:.2f}',
                   transform=ax_sed.transAxes, fontsize=7,
                   ha='right', va='top',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                            alpha=0.8, edgecolor='gray', linewidth=0.5))
        
        ax_sed.set_title(title, fontsize=9)
        ax_sed.tick_params(labelsize=7)
        
        if ax_sed == sed_axes[-1]:
            ax_sed.set_xlabel(r'Wavelength ($\mathrm{\AA}$)', fontsize=8)
        else:
            ax_sed.set_xticklabels([])
    
    # Common y-label for SED panels
    fig.text(0.52, 0.5, r'$F_\lambda$ (erg s$^{-1}$ cm$^{-2}$ $\mathrm{\AA}^{-1}$)',
             rotation=90, va='center', fontsize=8)
    
    plt.savefig(output_dir / filename, dpi=300)
    plt.close()
    print(f"Saved: {output_dir / filename}")


def plot_hr_diagram(history, output_dir, filename='fig_blueloop_hr.pdf'):
    """
    Plot HR diagram with color-coded evolutionary track.
    
    Parameters
    ----------
    history : dict
        MESA history data
    output_dir : Path
        Output directory
    filename : str
        Output filename
    """
    fig, ax = plt.subplots(figsize=(APJ_SINGLE_COL, 3.5))
    
    log_teff = history['log_Teff']
    log_l = history['log_L']
    g_r = history['g'] - history['r']
    
    # Plot track colored by g-r
    scatter = ax.scatter(log_teff, log_l, c=g_r, cmap='coolwarm', 
                        s=2, alpha=0.8, rasterized=True)
    ax.plot(log_teff, log_l, 'k-', lw=0.2, alpha=0.3, zorder=0)
    
    # POI markers
    poi = find_poi_blue_loop(history)
    markers = {'rgb_tip': 'o', 'loop_maximum': 's', 'end': 'D'}
    
    for name, idx in poi.items():
        if name in markers:
            ax.scatter(log_teff[idx], log_l[idx], c='black', s=50,
                      marker=markers[name], edgecolors='white', 
                      linewidths=0.5, zorder=5)
    
    ax.invert_xaxis()
    ax.set_xlabel(r'$\log(T_{\rm eff}/{\rm K})$')
    ax.set_ylabel(r'$\log(L/L_\odot)$')
    
    # Instability strip edges
    ax.axvline(np.log10(5500), color='red', ls='--', alpha=0.4, lw=0.8)
    ax.axvline(np.log10(6500), color='blue', ls='--', alpha=0.4, lw=0.8)
    
    cbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    cbar.set_label(r'$g - r$ (mag)', fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    
    plt.tight_layout()
    plt.savefig(output_dir / filename, dpi=300)
    plt.close()
    print(f"Saved: {output_dir / filename}")


def plot_evolution_panels(history, output_dir, 
                          filename='fig_blueloop_evolution.pdf'):
    """
    Plot multi-panel magnitude and color evolution.
    
    Parameters
    ----------
    history : dict
        MESA history data
    output_dir : Path
        Output directory
    filename : str
        Output filename
    """
    fig, axes = plt.subplots(3, 1, figsize=(APJ_SINGLE_COL, 5), sharex=True)
    
    age = history['star_age'] / 1e6
    
    # Panel 1: Multi-band magnitudes
    ax1 = axes[0]
    bands = ['u', 'g', 'r', 'i', 'z', 'y']
    for band in bands:
        if band in history:
            props = LSST_FILTERS[band]
            ax1.plot(age, history[band], color=props['color'], 
                    lw=0.8, label=band, alpha=0.9)
    
    ax1.invert_yaxis()
    ax1.set_ylabel('Magnitude (AB)')
    ax1.legend(loc='upper right', ncol=3, fontsize=7, framealpha=0.9)
    
    # Panel 2: g-r color
    ax2 = axes[1]
    g_r = history['g'] - history['r']
    ax2.plot(age, g_r, 'k-', lw=0.8)
    ax2.axhspan(0.25, 0.65, alpha=0.1, color='blue', zorder=0)
    ax2.set_ylabel(r'$g - r$ (mag)')
    
    # Panel 3: Teff
    ax3 = axes[2]
    teff = 10**history['log_Teff']
    ax3.plot(age, teff, 'k-', lw=0.8)
    ax3.axhspan(5500, 6500, alpha=0.1, color='blue', zorder=0)
    ax3.set_ylabel(r'$T_{\rm eff}$ (K)')
    ax3.set_xlabel('Age (Myr)')
    
    # Mark POIs
    poi = find_poi_blue_loop(history)
    for name, idx in poi.items():
        for ax in axes:
            ax.axvline(age[idx], color='gray', ls=':', alpha=0.5, lw=0.5)
    
    plt.tight_layout()
    plt.savefig(output_dir / filename, dpi=300)
    plt.close()
    print(f"Saved: {output_dir / filename}")


def print_summary(history):
    """Print summary statistics for the evolution."""
    print("\n" + "=" * 60)
    print("BLUE LOOP EVOLUTION SUMMARY")
    print("=" * 60)
    
    teff = 10**history['log_Teff']
    g_r = history['g'] - history['r']
    
    print(f"Teff range: {teff.min():.0f} - {teff.max():.0f} K")
    print(f"g-r range: {g_r.min():.3f} - {g_r.max():.3f} mag")
    print(f"Color amplitude: {g_r.max() - g_r.min():.3f} mag")
    
    age_span = (history['star_age'][-1] - history['star_age'][0]) / 1e6
    print(f"Time span: {age_span:.2f} Myr")
    print(f"Models: {len(history['model_number'])}")
    
    poi = find_poi_blue_loop(history)
    print("\nPoints of Interest:")
    for name, idx in poi.items():
        model = int(history['model_number'][idx])
        print(f"  {name}: model {model}, Teff = {teff[idx]:.0f} K")
    
    print("=" * 60 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description='Generate MESA Colors Blue Loop figures'
    )
    parser.add_argument('--logs_dir', type=str, default='../cepheid_blue_loop/LOGS',
                       help='Path to MESA LOGS directory')
    parser.add_argument('--sed_dir', type=str, default='../cepheid_blue_loop/SED',
                       help='Path to SED output directory')
    parser.add_argument('--output_dir', type=str, default='../cepheid_blue_loop/figures',
                       help='Output directory for figures')
    args = parser.parse_args()
    
    # Setup
    setup_apj_style()
    
    logs_dir = Path(args.logs_dir)
    sed_dir = Path(args.sed_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Read history
    history_file = logs_dir / 'history.data'
    if not history_file.exists():
        print(f"ERROR: {history_file} not found")
        return
    
    print(f"Reading: {history_file}")
    history = read_mesa_history(history_file)
    
    print_summary(history)
    
    # Generate figures
    print("Generating figures...")
    
    if sed_dir.exists():
        plot_cmd_with_seds(history, sed_dir, output_dir)
    else:
        print(f"Warning: SED directory {sed_dir} not found, skipping CMD+SED figure")
    
    plot_hr_diagram(history, output_dir)
    plot_evolution_panels(history, output_dir)
    
    print(f"\nDone! Figures saved to: {output_dir}")


if __name__ == '__main__':
    main()
