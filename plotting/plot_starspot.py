#!/usr/bin/env python3
"""
plot_starspot.py

Compare MESA Colors output between spotted and unspotted stellar models.
Demonstrates chromatic signature of starspots with SED comparison.

Figures produced:
- fig_starspot_chromatic.pdf : Magnitude difference vs wavelength with SEDs
- fig_starspot_sed.pdf       : Direct SED comparison spotted vs unspotted

Usage:
    python plot_starspot.py [--logs_spotted LOGS_spotted] [--logs_unspotted LOGS_unspotted]

Author: Miller, Joyce, Mocz et al.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path
import argparse
import csv

from plot_utils import (
    setup_apj_style, read_mesa_history, find_magnitude_columns,
    get_model_number_at_index, get_sed_counter_at_index,
    add_em_spectrum_regions,
    LSST_FILTERS, APJ_SINGLE_COL, APJ_DOUBLE_COL
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


def plot_chromatic_signature_with_seds(history_spotted, history_unspotted,
                                        sed_dir_spotted, sed_dir_unspotted,
                                        output_dir, filename='fig_starspot_chromatic.pdf'):
    """
    Plot chromatic signature with SED comparison panel.
    
    Parameters
    ----------
    history_spotted : dict
        Spotted star history data
    history_unspotted : dict
        Unspotted star history data
    sed_dir_spotted : Path
        SED directory for spotted model
    sed_dir_unspotted : Path  
        SED directory for unspotted model
    output_dir : Path
        Output directory
    filename : str
        Output filename
    """
    fig = plt.figure(figsize=(APJ_DOUBLE_COL, 3.5))
    
    gs = GridSpec(1, 2, width_ratios=[1, 1.3], wspace=0.25,
                  left=0.08, right=0.97, top=0.90, bottom=0.15)
    
    ax_chrom = fig.add_subplot(gs[0])
    ax_sed = fig.add_subplot(gs[1])
    
    # Find magnitude columns
    mag_cols_s = find_magnitude_columns(history_spotted, 'lsst')
    mag_cols_u = find_magnitude_columns(history_unspotted, 'lsst')
    
    # Extract final timestep magnitudes
    wavelengths = []
    delta_mags = []
    filters_used = []
    
    for filt in ['u', 'g', 'r', 'i', 'z', 'y']:
        if filt in mag_cols_s and filt in mag_cols_u:
            mag_s = history_spotted[mag_cols_s[filt]][-1]
            mag_u = history_unspotted[mag_cols_u[filt]][-1]
            
            wavelengths.append(LSST_FILTERS[filt]['wavelength'])
            delta_mags.append(mag_s - mag_u)
            filters_used.append(filt)
    
    wavelengths = np.array(wavelengths)
    delta_mags = np.array(delta_mags)
    
    # Plot chromatic signature
    for w, dm, f in zip(wavelengths, delta_mags, filters_used):
        color = LSST_FILTERS[f]['color']
        ax_chrom.scatter(w, dm, s=100, c=color, edgecolors='black',
                        linewidths=0.8, zorder=5)
        ax_chrom.annotate(f, (w, dm), xytext=(5, 5), textcoords='offset points',
                         fontsize=9, fontweight='bold')
    
    ax_chrom.plot(wavelengths, delta_mags, 'k--', alpha=0.4, lw=0.8, zorder=1)
    ax_chrom.axhline(0, color='gray', ls=':', alpha=0.5, lw=0.8)
    
    ax_chrom.set_xlabel('Wavelength (nm)')
    ax_chrom.set_ylabel(r'$\Delta m = m_{\rm spotted} - m_{\rm unspotted}$ (mag)')
    ax_chrom.set_title('Chromatic Signature', fontsize=10)
    
    # Annotation
    ax_chrom.text(0.97, 0.97, 
                 'Spots affect blue\nbands more than red',
                 transform=ax_chrom.transAxes, ha='right', va='top',
                 fontsize=8, style='italic',
                 bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', 
                          alpha=0.7, edgecolor='none'))
    
    # SED comparison panel
    idx_s = len(history_spotted['model_number']) - 1
    idx_u = len(history_unspotted['model_number']) - 1
    
    model_s = get_model_number_at_index(history_spotted, idx_s)
    model_u = get_model_number_at_index(history_unspotted, idx_u)
    
    sed_counter_s = get_sed_counter_at_index(history_spotted, idx_s)
    sed_counter_u = get_sed_counter_at_index(history_unspotted, idx_u)
    
    sed_plotted = False
    
    # Load and plot spotted SED
    sed_files_s = list(sed_dir_spotted.glob(f'*_SED_{sed_counter_s}.csv'))
    sed_files_u = list(sed_dir_unspotted.glob(f'*_SED_{sed_counter_u}.csv'))
    
    if sed_files_s:
        sed_data_s = load_sed_csv(sed_files_s[0])
        if sed_data_s and 'wavelengths' in sed_data_s:
            wl = sed_data_s['wavelengths']
            fl = sed_data_s['fluxes']
            ax_sed.plot(wl, fl, 'r-', lw=0.9, label='Spotted', alpha=0.9)
            sed_plotted = True
    
    if sed_files_u:
        sed_data_u = load_sed_csv(sed_files_u[0])
        if sed_data_u and 'wavelengths' in sed_data_u:
            wl = sed_data_u['wavelengths']
            fl = sed_data_u['fluxes']
            ax_sed.plot(wl, fl, 'b-', lw=0.9, label='Unspotted', alpha=0.9)
            sed_plotted = True
    
    if sed_plotted:
        ax_sed.set_xscale('log')
        ax_sed.set_xlabel(r'Wavelength ($\mathrm{\AA}$)')
        ax_sed.set_ylabel(r'$F_\lambda$ (erg s$^{-1}$ cm$^{-2}$ $\mathrm{\AA}^{-1}$)')
        ax_sed.set_title('SED Comparison', fontsize=10)
        ax_sed.set_xlim(3000, 12000)
        ax_sed.legend(loc='upper right', fontsize=8, framealpha=0.9)
        
        # Add filter markers
        for filt, props in LSST_FILTERS.items():
            ax_sed.axvline(props['wavelength'], color=props['color'],
                          alpha=0.3, lw=1.5, zorder=0)
        
        add_em_spectrum_regions(ax_sed, alpha=0.03)
    else:
        ax_sed.text(0.5, 0.5, 'SED files not found',
                   transform=ax_sed.transAxes, ha='center', va='center')
    
    plt.savefig(output_dir / filename, dpi=300)
    plt.close()
    print(f"Saved: {output_dir / filename}")


def plot_sed_ratio(history_spotted, history_unspotted,
                   sed_dir_spotted, sed_dir_unspotted,
                   output_dir, filename='fig_starspot_sed_ratio.pdf'):
    """
    Plot SED ratio (spotted/unspotted) to highlight wavelength dependence.
    
    Parameters
    ----------
    history_spotted : dict
        Spotted star history data
    history_unspotted : dict
        Unspotted star history data
    sed_dir_spotted : Path
        SED directory for spotted model
    sed_dir_unspotted : Path
        SED directory for unspotted model
    output_dir : Path
        Output directory
    filename : str
        Output filename
    """
    fig, axes = plt.subplots(2, 1, figsize=(APJ_SINGLE_COL, 4), 
                             sharex=True, height_ratios=[2, 1])
    
    ax_sed = axes[0]
    ax_ratio = axes[1]
    
    # Get model numbers
    idx_s = len(history_spotted['model_number']) - 1
    idx_u = len(history_unspotted['model_number']) - 1
    
    model_s = get_model_number_at_index(history_spotted, idx_s)
    model_u = get_model_number_at_index(history_unspotted, idx_u)
    
    sed_counter_s = get_sed_counter_at_index(history_spotted, idx_s)
    sed_counter_u = get_sed_counter_at_index(history_unspotted, idx_u)
    
    sed_files_s = list(sed_dir_spotted.glob(f'*_SED_{sed_counter_s}.csv'))
    sed_files_u = list(sed_dir_unspotted.glob(f'*_SED_{sed_counter_u}.csv'))
    
    if not sed_files_s or not sed_files_u:
        print("SED files not found for ratio plot")
        plt.close()
        return
    
    sed_data_s = load_sed_csv(sed_files_s[0])
    sed_data_u = load_sed_csv(sed_files_u[0])
    
    if not sed_data_s or not sed_data_u:
        plt.close()
        return
    
    wl_s = sed_data_s['wavelengths']
    fl_s = sed_data_s['fluxes']
    wl_u = sed_data_u['wavelengths']
    fl_u = sed_data_u['fluxes']
    
    # Plot SEDs
    ax_sed.plot(wl_u, fl_u, 'b-', lw=0.9, label='Unspotted', alpha=0.9)
    ax_sed.plot(wl_s, fl_s, 'r-', lw=0.9, label='Spotted', alpha=0.9)
    
    ax_sed.set_xscale('log')
    ax_sed.set_ylabel(r'$F_\lambda$')
    ax_sed.legend(loc='upper right', fontsize=7, framealpha=0.9)
    ax_sed.set_xlim(3000, 12000)
    
    # Interpolate to common grid and compute ratio
    from scipy.interpolate import interp1d
    
    wl_common = np.logspace(np.log10(3000), np.log10(12000), 500)
    
    interp_s = interp1d(wl_s, fl_s, bounds_error=False, fill_value=np.nan)
    interp_u = interp1d(wl_u, fl_u, bounds_error=False, fill_value=np.nan)
    
    fl_s_interp = interp_s(wl_common)
    fl_u_interp = interp_u(wl_common)
    
    ratio = fl_s_interp / fl_u_interp
    
    ax_ratio.plot(wl_common, ratio, 'k-', lw=0.9)
    ax_ratio.axhline(1.0, color='gray', ls=':', alpha=0.5)
    ax_ratio.set_xscale('log')
    ax_ratio.set_xlabel(r'Wavelength ($\mathrm{\AA}$)')
    ax_ratio.set_ylabel(r'$F_{\rm spotted}/F_{\rm unspotted}$')
    ax_ratio.set_ylim(0.85, 1.02)
    
    # Add filter markers to ratio panel
    for filt, props in LSST_FILTERS.items():
        ax_ratio.axvline(props['wavelength'], color=props['color'],
                        alpha=0.4, lw=1.2, zorder=0)
    
    add_em_spectrum_regions(ax_sed, alpha=0.03)
    
    plt.tight_layout()
    plt.savefig(output_dir / filename, dpi=300)
    plt.close()
    print(f"Saved: {output_dir / filename}")


def print_summary(history_spotted, history_unspotted):
    """Print starspot comparison summary."""
    print("\n" + "=" * 60)
    print("STARSPOT CHROMATIC SIGNATURE SUMMARY")
    print("=" * 60)
    
    mag_cols_s = find_magnitude_columns(history_spotted, 'lsst')
    mag_cols_u = find_magnitude_columns(history_unspotted, 'lsst')
    
    print(f"\n{'Filter':<8} {'λ (nm)':<10} {'Δm (mag)':<10}")
    print("-" * 35)
    
    for filt in ['u', 'g', 'r', 'i', 'z', 'y']:
        if filt in mag_cols_s and filt in mag_cols_u:
            mag_s = history_spotted[mag_cols_s[filt]][-1]
            mag_u = history_unspotted[mag_cols_u[filt]][-1]
            dm = mag_s - mag_u
            wl = LSST_FILTERS[filt]['wavelength']
            print(f"{filt:<8} {wl:<10.0f} {dm:<+10.4f}")
    
    # Stellar parameters
    print("\nStellar Parameters (unspotted):")
    if 'log_Teff' in history_unspotted:
        print(f"  log Teff = {history_unspotted['log_Teff'][-1]:.4f}")
    if 'log_L' in history_unspotted:
        print(f"  log L/Lsun = {history_unspotted['log_L'][-1]:.4f}")
    if 'log_g' in history_unspotted:
        print(f"  log g = {history_unspotted['log_g'][-1]:.4f}")
    
    print("=" * 60 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description='Generate MESA Colors starspot figures'
    )
    parser.add_argument('--logs_spotted', type=str, default='LOGS_spotted',
                       help='LOGS directory for spotted model')
    parser.add_argument('--logs_unspotted', type=str, default='LOGS_unspotted',
                       help='LOGS directory for unspotted model')
    parser.add_argument('--sed_spotted', type=str, default='SED_spotted',
                       help='SED directory for spotted model')
    parser.add_argument('--sed_unspotted', type=str, default='SED_unspotted',
                       help='SED directory for unspotted model')
    parser.add_argument('--output_dir', type=str, default='figures',
                       help='Output directory for figures')
    args = parser.parse_args()
    
    setup_apj_style()
    
    logs_spotted = Path(args.logs_spotted)
    logs_unspotted = Path(args.logs_unspotted)
    sed_spotted = Path(args.sed_spotted)
    sed_unspotted = Path(args.sed_unspotted)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Read history files
    history_file_s = logs_spotted / 'history.data'
    history_file_u = logs_unspotted / 'history.data'
    
    if not history_file_s.exists():
        print(f"ERROR: {history_file_s} not found")
        return
    if not history_file_u.exists():
        print(f"ERROR: {history_file_u} not found")
        return
    
    print(f"Reading: {history_file_s}")
    history_spotted = read_mesa_history(history_file_s)
    
    print(f"Reading: {history_file_u}")
    history_unspotted = read_mesa_history(history_file_u)
    
    print_summary(history_spotted, history_unspotted)
    
    print("Generating figures...")
    
    plot_chromatic_signature_with_seds(
        history_spotted, history_unspotted,
        sed_spotted, sed_unspotted,
        output_dir
    )
    
    if sed_spotted.exists() and sed_unspotted.exists():
        plot_sed_ratio(
            history_spotted, history_unspotted,
            sed_spotted, sed_unspotted,
            output_dir
        )
    
    print(f"\nDone! Figures saved to: {output_dir}")


if __name__ == '__main__':
    main()
