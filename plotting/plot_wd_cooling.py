#!/usr/bin/env python3
"""
plot_wd_cooling.py

Analysis and plotting for MESA Colors White Dwarf cooling sequence.
Generates publication-quality figures showing photometric evolution during
WD cooling, with SED snapshots at key cooling ages.

Figures produced:
- fig_wd_cmd.pdf         : CMD (g-r vs g) with SED insets at key ages
- fig_wd_lightcurves.pdf : Multi-band magnitude evolution vs cooling age
- fig_wd_colors.pdf      : Color evolution vs Teff
- fig_wd_sed_evolution.pdf : SED sequence showing cooling progression
- fig_wd_summary.pdf     : 4-panel summary figure

Usage:
    python plot_wd_cooling.py [--logs_dir LOGS] [--sed_dir SED] [--output_dir figures]

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


def find_poi_wd_cooling(history):
    """
    Find points of interest for WD cooling sequence.
    
    Returns indices for:
    - hot_wd: Initial hot white dwarf
    - mid_cooling: Middle of cooling sequence
    - cool_wd: Cool white dwarf
    - crystallization: Near crystallization onset (if applicable)
    
    Parameters
    ----------
    history : dict
        MESA history data
        
    Returns
    -------
    dict
        Named indices for points of interest
    """
    n = len(history['model_number'])
    
    if 'log_Teff' in history:
        log_Teff = history['log_Teff']
        teff = 10**log_Teff
        
        # Hot WD: highest Teff
        hot_idx = np.argmax(teff)
        
        # Cool WD: lowest Teff
        cool_idx = np.argmin(teff)
        
        # Mid cooling: ~15000 K if available, else midpoint
        mid_temp_target = 15000
        mid_idx = np.argmin(np.abs(teff - mid_temp_target))
        
        # Crystallization: ~6000 K for DA WDs
        cryst_temp_target = 6000
        cryst_idx = np.argmin(np.abs(teff - cryst_temp_target))
        
    else:
        # Fallback to evenly spaced
        hot_idx = 0
        mid_idx = n // 3
        cryst_idx = 2 * n // 3
        cool_idx = n - 1
    
    return {
        'hot_wd': hot_idx,
        'mid_cooling': mid_idx,
        'crystallization': cryst_idx,
        'cool_wd': cool_idx
    }


def plot_cmd_with_seds(history, sed_dir, output_dir, 
                       filename='fig_wd_cmd.pdf'):
    """
    Plot CMD with SED insets at key cooling ages.
    
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
    mag_cols = find_magnitude_columns(history, 'lsst')
    
    if 'g' not in mag_cols or 'r' not in mag_cols:
        print("Error: g and r bands not found")
        plt.close()
        return
    
    g = history[mag_cols['g']]
    r = history[mag_cols['r']]
    color = g - r
    
    # Color by age
    if 'star_age' in history:
        age_gyr = history['star_age'] / 1e9
        
        valid = np.isfinite(color) & np.isfinite(g)
        
        scatter = ax_cmd.scatter(color[valid], g[valid], c=age_gyr[valid],
                                cmap='plasma', s=4, alpha=0.8,
                                rasterized=True)
        
        cbar = fig.colorbar(scatter, ax=ax_cmd, pad=0.02, aspect=25)
        cbar.set_label('Age (Gyr)', fontsize=8)
        cbar.ax.tick_params(labelsize=7)
    else:
        ax_cmd.plot(color, g, 'b-', lw=0.8, alpha=0.7)
    
    # POI markers
    poi = find_poi_wd_cooling(history)
    
    poi_style = {
        'hot_wd': {'marker': 'o', 'color': '#d62728', 'label': 'Hot WD'},
        'mid_cooling': {'marker': 's', 'color': '#ff7f0e', 'label': 'Mid cooling'},
        'crystallization': {'marker': '^', 'color': '#2ca02c', 'label': 'Crystallization'},
        'cool_wd': {'marker': 'D', 'color': '#1f77b4', 'label': 'Cool WD'},
    }
    
    for name, idx in poi.items():
        if name in poi_style and np.isfinite(color[idx]) and np.isfinite(g[idx]):
            style = poi_style[name]
            ax_cmd.scatter(color[idx], g[idx], c=style['color'], s=80,
                          marker=style['marker'], edgecolors='black',
                          linewidths=0.8, zorder=5, label=style['label'])
    
    ax_cmd.set_xlabel(r'$g - r$ (mag)')
    ax_cmd.set_ylabel(r'$M_g$ (mag)')
    ax_cmd.invert_yaxis()
    ax_cmd.legend(loc='lower right', fontsize=6, framealpha=0.9)
    ax_cmd.set_title('WD Cooling Sequence', fontsize=10)
    
    # SED panels
    sed_pois = ['hot_wd', 'mid_cooling', 'crystallization', 'cool_wd']
    sed_titles = ['Hot WD', 'Mid Cooling', 'Crystallization', 'Cool WD']
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
                    for filt, props in LSST_FILTERS.items():
                        ax_sed.axvline(props['wavelength'],
                                      color=props['color'],
                                      alpha=0.4, lw=1.2, zorder=0)
                    
                    ax_sed.set_xlim(3000, 12000)
                    
                    mask_vis = (wl > 3000) & (wl < 12000)
                    if mask_vis.any():
                        fl_vis = fl[mask_vis]
                        if len(fl_vis) > 0 and fl_vis.max() > 0:
                            ax_sed.set_ylim(0, fl_vis.max() * 1.1)
                    
                    add_em_spectrum_regions(ax_sed, alpha=0.03)
            
            # Stellar params
            if 'log_Teff' in history:
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


def plot_lightcurves(history, output_dir, filename='fig_wd_lightcurves.pdf'):
    """
    Plot multi-band magnitude evolution vs cooling age.
    
    Parameters
    ----------
    history : dict
        MESA history data
    output_dir : Path
        Output directory
    filename : str
        Output filename
    """
    fig, ax = plt.subplots(figsize=(APJ_SINGLE_COL, 3))
    
    mag_cols = find_magnitude_columns(history, 'lsst')
    
    if not mag_cols:
        print("Warning: No LSST magnitudes found")
        plt.close()
        return
    
    # Time axis
    if 'star_age' in history:
        time = history['star_age'] / 1e9  # Gyr
        xlabel = 'Age (Gyr)'
    else:
        time = np.arange(len(history['model_number']))
        xlabel = 'Model Number'
    
    # Plot each filter
    for filt in ['u', 'g', 'r', 'i', 'z', 'y']:
        if filt in mag_cols:
            mag = history[mag_cols[filt]]
            valid = np.isfinite(mag)
            ax.plot(time[valid], mag[valid], 
                   color=LSST_FILTERS[filt]['color'],
                   lw=0.9, label=filt, alpha=0.9)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Absolute Magnitude (AB)')
    ax.invert_yaxis()
    ax.legend(loc='upper left', ncol=3, fontsize=7, framealpha=0.9)
    
    plt.tight_layout()
    plt.savefig(output_dir / filename, dpi=300)
    plt.close()
    print(f"Saved: {output_dir / filename}")


def plot_color_evolution(history, output_dir, filename='fig_wd_colors.pdf'):
    """
    Plot color evolution vs effective temperature.
    
    Parameters
    ----------
    history : dict
        MESA history data
    output_dir : Path
        Output directory
    filename : str
        Output filename
    """
    if 'log_Teff' not in history:
        print("Warning: log_Teff not found, skipping color evolution plot")
        return
    
    teff = 10**history['log_Teff']
    valid = (teff > 1000) & (teff < 200000)
    
    mag_cols = find_magnitude_columns(history, 'lsst')
    
    fig, axes = plt.subplots(2, 2, figsize=(APJ_DOUBLE_COL, 4), sharex=True)
    axes = axes.flatten()
    
    color_pairs = [
        ('u', 'g', r'$u - g$'),
        ('g', 'r', r'$g - r$'),
        ('r', 'i', r'$r - i$'),
        ('i', 'z', r'$i - z$'),
    ]
    
    for ax, (blue, red, label) in zip(axes, color_pairs):
        if blue in mag_cols and red in mag_cols:
            color = history[mag_cols[blue]] - history[mag_cols[red]]
            
            sc = ax.scatter(teff[valid], color[valid], c=teff[valid], 
                           cmap='coolwarm', s=2, alpha=0.7, rasterized=True)
            ax.set_ylabel(label)
            ax.set_xscale('log')
            ax.invert_xaxis()
        else:
            ax.text(0.5, 0.5, f'{label}\nnot available',
                   transform=ax.transAxes, ha='center', va='center', fontsize=8)
    
    for ax in axes[2:]:
        ax.set_xlabel(r'$T_{\rm eff}$ (K)')
    
    plt.tight_layout()
    plt.savefig(output_dir / filename, dpi=300)
    plt.close()
    print(f"Saved: {output_dir / filename}")


def plot_sed_evolution(history, sed_dir, output_dir, 
                       filename='fig_wd_sed_evolution.pdf'):
    """
    Plot SED sequence showing cooling progression.
    
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
    fig, ax = plt.subplots(figsize=(APJ_SINGLE_COL, 3.5))
    
    # Get all SED files
    sed_files = sorted(sed_dir.glob('*_SED_*.csv'))
    
    if not sed_files:
        print("Warning: No SED files found for evolution plot")
        plt.close()
        return
    
    # Extract model numbers and select subset
    model_nums = []
    for f in sed_files:
        try:
            # Extract model number from filename like "filter_SED_123.csv"
            parts = f.stem.split('_')
            model_num = int(parts[-1])
            if model_num not in model_nums:
                model_nums.append(model_num)
        except (ValueError, IndexError):
            continue
    
    model_nums = sorted(set(model_nums))
    
    if len(model_nums) == 0:
        print("Warning: Could not parse model numbers from SED files")
        plt.close()
        return
    
    # Select ~8 evenly spaced SEDs
    n_seds = min(8, len(model_nums))
    indices = np.linspace(0, len(model_nums) - 1, n_seds, dtype=int)
    selected_models = [model_nums[i] for i in indices]
    
    # Colormap
    cmap = plt.cm.plasma
    colors = cmap(np.linspace(0.1, 0.9, len(selected_models)))
    
    for model_num, color in zip(selected_models, colors):
        sed_files_model = list(sed_dir.glob(f'*_SED_{model_num}.csv'))
        
        if sed_files_model:
            sed_data = load_sed_csv(sed_files_model[0])
            
            if sed_data and 'wavelengths' in sed_data:
                wl = sed_data['wavelengths']
                fl = sed_data['fluxes']
                
                # Get Teff for this model if available
                if 'model_number' in history:
                    model_idx = np.where(history['model_number'] == model_num)[0]
                    if len(model_idx) > 0 and 'log_Teff' in history:
                        teff = 10**history['log_Teff'][model_idx[0]]
                        label = f'{teff:.0f} K'
                    else:
                        label = f'Model {model_num}'
                else:
                    label = f'Model {model_num}'
                
                ax.plot(wl, fl, color=color, lw=0.8, alpha=0.9, label=label)
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'Wavelength ($\mathrm{\AA}$)')
    ax.set_ylabel(r'$F_\lambda$ (erg s$^{-1}$ cm$^{-2}$ $\mathrm{\AA}^{-1}$)')
    ax.set_xlim(1000, 30000)
    
    # Colorbar for cooling progression
    sm = ScalarMappable(cmap=cmap, norm=Normalize(0, 1))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02, aspect=20)
    cbar.set_label('Cooling →', fontsize=8)
    cbar.set_ticks([])
    
    ax.legend(loc='upper right', fontsize=6, framealpha=0.9, ncol=2)
    
    add_em_spectrum_regions(ax, alpha=0.03)
    
    plt.tight_layout()
    plt.savefig(output_dir / filename, dpi=300)
    plt.close()
    print(f"Saved: {output_dir / filename}")


def plot_hr_diagram(history, output_dir, filename='fig_wd_hr.pdf'):
    """
    Plot HR diagram with cooling track.
    
    Parameters
    ----------
    history : dict
        MESA history data
    output_dir : Path
        Output directory
    filename : str
        Output filename
    """
    if 'log_Teff' not in history or 'log_L' not in history:
        print("Warning: Missing log_Teff or log_L, skipping HR diagram")
        return
    
    fig, ax = plt.subplots(figsize=(APJ_SINGLE_COL, 3.5))
    
    log_teff = history['log_Teff']
    log_L = history['log_L']
    
    if 'star_age' in history:
        age_gyr = history['star_age'] / 1e9
        sc = ax.scatter(log_teff, log_L, c=age_gyr, cmap='plasma', 
                       s=3, alpha=0.8, rasterized=True)
        cbar = fig.colorbar(sc, ax=ax, pad=0.02)
        cbar.set_label('Age (Gyr)', fontsize=8)
        cbar.ax.tick_params(labelsize=7)
    else:
        ax.plot(log_teff, log_L, 'b-', lw=0.8)
    
    # Mark start/end
    ax.plot(log_teff[0], log_L[0], 'ro', ms=6, zorder=10, label='Start')
    ax.plot(log_teff[-1], log_L[-1], 'bs', ms=6, zorder=10, label='End')
    
    ax.set_xlabel(r'$\log T_{\rm eff}$ (K)')
    ax.set_ylabel(r'$\log(L/L_\odot)$')
    ax.invert_xaxis()
    ax.legend(loc='lower left', fontsize=7, framealpha=0.9)
    
    plt.tight_layout()
    plt.savefig(output_dir / filename, dpi=300)
    plt.close()
    print(f"Saved: {output_dir / filename}")


def plot_summary(history, sed_dir, output_dir, filename='fig_wd_summary.pdf'):
    """
    Create 4-panel summary figure.
    
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
    gs = GridSpec(2, 2, hspace=0.35, wspace=0.3,
                  left=0.08, right=0.95, top=0.93, bottom=0.10)
    
    mag_cols = find_magnitude_columns(history, 'lsst')
    
    # (a) CMD
    ax1 = fig.add_subplot(gs[0, 0])
    
    if 'g' in mag_cols and 'r' in mag_cols:
        g = history[mag_cols['g']]
        r = history[mag_cols['r']]
        color = g - r
        
        if 'star_age' in history:
            age_gyr = history['star_age'] / 1e9
            valid = np.isfinite(color) & np.isfinite(g)
            
            sc = ax1.scatter(color[valid], g[valid], c=age_gyr[valid],
                            cmap='plasma', s=3, alpha=0.7, rasterized=True)
            ax1.plot(color[0], g[0], 'ro', ms=5, zorder=10)
            ax1.plot(color[-1], g[-1], 'bs', ms=5, zorder=10)
    
    ax1.set_xlabel(r'$g - r$')
    ax1.set_ylabel(r'$M_g$')
    ax1.invert_yaxis()
    ax1.set_title('(a) Color-Magnitude Diagram', fontsize=9)
    
    # (b) Light curves
    ax2 = fig.add_subplot(gs[0, 1])
    
    if 'star_age' in history:
        time = history['star_age'] / 1e9
    else:
        time = np.arange(len(history['model_number']))
    
    for filt in ['u', 'g', 'r', 'i', 'z', 'y']:
        if filt in mag_cols:
            mag = history[mag_cols[filt]]
            valid = np.isfinite(mag)
            ax2.plot(time[valid], mag[valid],
                    color=LSST_FILTERS[filt]['color'],
                    lw=0.7, label=filt, alpha=0.9)
    
    ax2.set_xlabel('Age (Gyr)')
    ax2.set_ylabel('Magnitude')
    ax2.invert_yaxis()
    ax2.legend(loc='upper left', ncol=3, fontsize=6)
    ax2.set_title('(b) Multi-band Light Curves', fontsize=9)
    
    # (c) HR diagram
    ax3 = fig.add_subplot(gs[1, 0])
    
    if 'log_Teff' in history and 'log_L' in history:
        log_teff = history['log_Teff']
        log_L = history['log_L']
        
        if 'star_age' in history:
            age_gyr = history['star_age'] / 1e9
            sc = ax3.scatter(log_teff, log_L, c=age_gyr, cmap='plasma',
                            s=2, alpha=0.7, rasterized=True)
        ax3.plot(log_teff[0], log_L[0], 'ro', ms=5, zorder=10)
        ax3.plot(log_teff[-1], log_L[-1], 'bs', ms=5, zorder=10)
    
    ax3.set_xlabel(r'$\log T_{\rm eff}$ (K)')
    ax3.set_ylabel(r'$\log(L/L_\odot)$')
    ax3.invert_xaxis()
    ax3.set_title('(c) HR Diagram', fontsize=9)
    
    # (d) SED evolution
    ax4 = fig.add_subplot(gs[1, 1])
    
    poi = find_poi_wd_cooling(history)
    cmap = plt.cm.plasma
    
    for i, (poi_name, label) in enumerate([
        ('hot_wd', 'Hot WD'),
        ('mid_cooling', 'Mid cooling'),
        ('crystallization', 'Crystallization'),
        ('cool_wd', 'Cool WD')
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
    ax4.set_yscale('log')
    ax4.set_xlabel(r'Wavelength ($\mathrm{\AA}$)')
    ax4.set_ylabel(r'$F_\lambda$')
    ax4.set_xlim(1000, 30000)
    ax4.legend(loc='upper right', fontsize=6, framealpha=0.9)
    ax4.set_title('(d) SED Evolution', fontsize=9)
    
    add_em_spectrum_regions(ax4, alpha=0.03)
    
    fig.suptitle('MESA Colors: White Dwarf Cooling Sequence', fontsize=11)
    
    plt.savefig(output_dir / filename, dpi=300)
    plt.close()
    print(f"Saved: {output_dir / filename}")


def print_summary(history):
    """Print WD cooling summary statistics."""
    print("\n" + "=" * 60)
    print("WHITE DWARF COOLING SEQUENCE SUMMARY")
    print("=" * 60)
    
    print(f"Total models: {len(history['model_number'])}")
    
    if 'log_Teff' in history:
        teff = 10**history['log_Teff']
        print(f"Teff range: {teff.min():.0f} - {teff.max():.0f} K")
    
    if 'log_L' in history:
        log_L = history['log_L']
        print(f"log(L/Lsun) range: {log_L.min():.2f} - {log_L.max():.2f}")
    
    if 'star_age' in history:
        age = history['star_age']
        print(f"Age range: {age.min()/1e9:.3f} - {age.max()/1e9:.3f} Gyr")
    
    mag_cols = find_magnitude_columns(history, 'lsst')
    print(f"\nAvailable LSST filters: {list(mag_cols.keys())}")
    
    if 'g' in mag_cols and 'r' in mag_cols:
        g = history[mag_cols['g']]
        r = history[mag_cols['r']]
        color = g - r
        valid = np.isfinite(color)
        if valid.any():
            print(f"g-r color range: {color[valid].min():.3f} - {color[valid].max():.3f}")
    
    print("=" * 60 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description='Generate MESA Colors WD cooling figures'
    )
    parser.add_argument('--logs_dir', type=str, default='../WD/LOGS',
                       help='Path to MESA LOGS directory')
    parser.add_argument('--sed_dir', type=str, default='../WD/SED',
                       help='Path to SED output directory')
    parser.add_argument('--output_dir', type=str, default='../WD/figures',
                       help='Output directory for figures')
    args = parser.parse_args()
    
    setup_apj_style()
    
    logs_dir = Path(args.logs_dir)
    sed_dir = Path(args.sed_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    history_file = logs_dir / 'history.data'
    if not history_file.exists():
        print(f"ERROR: {history_file} not found")
        return
    
    print(f"Reading: {history_file}")
    history = read_mesa_history(history_file)
    
    print_summary(history)
    
    print("Generating figures...")
    
    if sed_dir.exists():
        plot_cmd_with_seds(history, sed_dir, output_dir)
        plot_sed_evolution(history, sed_dir, output_dir)
        plot_summary(history, sed_dir, output_dir)
    else:
        print(f"Warning: SED directory {sed_dir} not found")
    
    plot_lightcurves(history, output_dir)
    plot_color_evolution(history, output_dir)
    plot_hr_diagram(history, output_dir)
    
    print(f"\nDone! Figures saved to: {output_dir}")


if __name__ == '__main__':
    main()
