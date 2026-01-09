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
    get_sed_counter_at_index, add_em_spectrum_regions, LSST_FILTERS,
    APJ_SINGLE_COL, APJ_DOUBLE_COL
)
from matplotlib.ticker import NullLocator, NullFormatter, ScalarFormatter,FixedLocator


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
    # MODIFIED: hspace set to 0.0 so SEDs touch. 
    # Adjusted 'right' to 0.88 to make room for right-sided Y labels.
    gs = GridSpec(3, 2, width_ratios=[1.2, 1], wspace=0.1, hspace=0.0,
                  left=0.08, right=0.88, top=0.95, bottom=0.10)
    
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
    cbar = fig.colorbar(scatter, ax=ax_cmd, pad=0.0, aspect=30)
    cbar.set_label('Age (Myr)', fontsize=7)
    cbar.ax.tick_params(labelsize=5)
    
    ax_cmd.legend(loc='upper right', fontsize=7, framealpha=0.9)
    
    # SED panels on right
    sed_axes = [fig.add_subplot(gs[i, 1]) for i in range(3)]
    
    # Select 3 POIs for SED display
    sed_pois = ['rgb_tip', 'loop_maximum', 'end']
    sed_titles = ['RGB Tip', 'Blue Maximum', 'Early AGB']
    
    for i, (ax_sed, poi_name, title) in enumerate(zip(sed_axes, sed_pois, sed_titles)):
        idx = poi[poi_name]
        model_num = get_model_number_at_index(history, idx)
        sed_counter = get_sed_counter_at_index(history, idx)
        
        # Try to load SED (files use counter, not MESA model number)
        sed_files = list(sed_dir.glob(f'*_SED_{sed_counter}.csv'))
        
        print(f"  {poi_name}: idx={idx}, sed_counter={sed_counter}, found {len(sed_files)} files")
        
        # MODIFIED: Move Y-axis ticks to the right side
        ax_sed.yaxis.tick_right()
        
        if sed_files:
            # Read first filter's SED (they all have the same underlying SED)
            sed_data = {'wavelengths': [], 'fluxes': [], 'convolved_flux': []}
            with open(sed_files[0], 'r') as f:
                import csv
                reader = csv.DictReader(f)
                
                for row in reader:
                    try:
                        wl = float(row['wavelengths'])
                        fl = float(row['fluxes'])
                        # Only keep rows where both wavelength and flux are valid and non-zero
                        if wl > 0 and fl > 0:
                            sed_data['wavelengths'].append(wl)
                            sed_data['fluxes'].append(fl)
                            # Convolved flux is optional
                            try:
                                cf = float(row.get('convolved_flux', 0))
                                sed_data['convolved_flux'].append(cf)
                            except (ValueError, TypeError):
                                sed_data['convolved_flux'].append(0)
                    except (ValueError, KeyError):
                        continue
                
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
        ax_sed.text(0.97, 0.35, f'$T_{{\\rm eff}}$ = {teff:.0f} K\n$\\log g$ = {logg:.2f}',
                   transform=ax_sed.transAxes, fontsize=7,
                   ha='right', va='top',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                             alpha=0.6, edgecolor='none', linewidth=0.5))
        
        # MODIFIED: Retrieve color from poi_style and apply to text
        title_color = poi_style[poi_name]['color']
        
        ax_sed.text(0.97, 0.92, title, transform=ax_sed.transAxes,
                    fontsize=9, fontweight='bold', ha='right', va='top',
                    color=title_color,  # Apply the POI color here
                    bbox=dict(facecolor='white', alpha=0.6, pad=0, edgecolor='none'))

        ax_sed.text(0.98, 0.93, title, transform=ax_sed.transAxes,
                    fontsize=9, fontweight='bold', ha='right', va='top',
                    color='k',  # Apply the POI color here
                    bbox=dict(facecolor='white', alpha=0.6, pad=0, edgecolor='none'))
        
        
        ax_sed.tick_params(labelsize=7)
                
        # Shared X-axis logic
        ax_sed.xaxis.get_offset_text().set_visible(False)

        if ax_sed != sed_axes[-1]:
            ax_sed.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
        else:
            ax_sed.set_xlabel(r'Wavelength ($\mathrm{\AA}$)', fontsize=8)

    
    # MODIFIED: Common y-label moved to the far Right
    # x coordinate increased to ~0.98 to sit on the right side
    fig.text(0.93, 0.5, r'$F_\lambda$ (erg s$^{-1}$ cm$^{-2}$ $\mathrm{\AA}^{-1}$)',
             rotation=270, va='center', ha='center', fontsize=8)
    
    plt.savefig(output_dir / filename, dpi=300)
    plt.close()
    print(f"Saved: {output_dir / filename}")






def plot_evolution_panels(history, sed_dir, output_dir, 
                          filename='fig_blueloop_evolution.pdf'):

    sed_dir = Path(sed_dir)
    fig = plt.figure(figsize=(APJ_DOUBLE_COL, 6))

    # Tighter GridSpec
    gs = GridSpec(
        2, 3,
        height_ratios=[1.5, 1],
        hspace=0.0,
        wspace=0.0,
        left=0.10,
        right=0.97,
        top=0.95,
        bottom=0.08
    )


    global_fl_max = 0.0

    age = history['star_age'] / 1e6

    # -------------------------
    # Top panel: light curves
    # -------------------------
    ax_lc = fig.add_subplot(gs[0, :])

    bands = ['u', 'g', 'r', 'i', 'z', 'y']
    for band in bands:
        if band in history:
            props = LSST_FILTERS[band]
            ax_lc.plot(
                age, history[band],
                color=props['color'],
                lw=0.8,
                label=band,
                alpha=0.9
            )

    ax_lc.invert_yaxis()
    ax_lc.set_ylabel('Magnitude (AB)')

    # X-axis on top
    ax_lc.xaxis.set_ticks_position('top')
    ax_lc.xaxis.set_label_position('top')
    ax_lc.set_xlabel('Age (Myr)')
    ax_lc.tick_params(labelbottom=False, labeltop=True)

    ax_lc.legend(
        loc='upper right',
        ncol=3,
        fontsize=7,
        framealpha=0.9
    )

    poi = find_poi_blue_loop(history)
    for idx in poi.values():
        ax_lc.axvline(age[idx], color='gray', ls=':', alpha=0.4, lw=0.5)

    sed_pois   = ['loop_ascending', 'loop_maximum', 'end']
    sed_titles = ['Entering Loop', 'Blue Maximum', 'Returning']
    sed_colors = ['#2ca02c', '#1f77b4', '#9467bd']

    shared_ylim = None

    # -------------------------
    # Bottom row: SEDs
    # -------------------------
    for i, (poi_name, title, marker_color) in enumerate(
        zip(sed_pois, sed_titles, sed_colors)
    ):
        ax_sed = fig.add_subplot(gs[1, i])

        if poi_name in poi:
            idx = poi[poi_name]
            sed_counter = get_sed_counter_at_index(history, idx)

            ax_lc.axvline(age[idx], color=marker_color, lw=1.2, alpha=0.6)
            ax_lc.scatter(
                age[idx], history['g'][idx],
                c=marker_color, s=30,
                zorder=5, edgecolors='black', linewidths=0.5
            )

            sed_files = list(sed_dir.glob(f'*_SED_{sed_counter}.csv')) if sed_dir.exists() else []

            if sed_files:
                wl, fl = [], []
                with open(sed_files[0]) as f:
                    import csv
                    for row in csv.DictReader(f):
                        try:
                            w = float(row['wavelengths'])
                            f_ = float(row['fluxes'])
                            if w > 0 and f_ > 0:
                                wl.append(w)
                                fl.append(f_)
                        except Exception:
                            pass

                wl = np.asarray(wl)
                fl = np.asarray(fl)

                if len(wl):
                    ax_sed.plot(wl, fl, 'k-', lw=0.7)
                    ax_sed.set_xscale('log')
                    xticks = [3500, 5000, 7000, 11000]
                    ax_sed.xaxis.set_major_locator(FixedLocator(xticks))
                    ax_sed.xaxis.set_minor_locator(NullLocator())
                    ax_sed.set_xticklabels([str(x) for x in xticks])
                    ax_sed.xaxis.get_offset_text().set_visible(False)



                    ax_sed.set_xlim(3000, 15000)

                    for filt, props in LSST_FILTERS.items():
                        ax_sed.axvline(
                            props['wavelength'],
                            color=props['color'],
                            alpha=0.35,
                            lw=1.0
                        )

                mask = (wl > 3000) & (wl < 11000)
                if mask.any():
                    fl_max = fl[mask].max()
                    if fl_max > global_fl_max:
                        global_fl_max = fl_max


                    add_em_spectrum_regions(ax_sed, alpha=0.03)


        # Subtitle INSIDE subplot
        ax_sed.text(
            0.0254, 0.9455,
            title,
            transform=ax_sed.transAxes,
            fontsize=8,
            color=marker_color,
            ha='left',
            va='top',
            weight='bold'
        )

        ax_sed.text(
            0.03, 0.95,
            title,
            transform=ax_sed.transAxes,
            fontsize=8,
            color='k',
            ha='left',
            va='top',
            weight='bold'
        )


        if i != 0:
            ax_sed.tick_params(axis='y', which='both',
                               labelbottom=False)
        else:
            ax_sed.tick_params(axis='y', which='both',
                               labelsize=6)

        if global_fl_max > 0:
            for ax in fig.axes:
                if ax is not ax_lc:
                    ax.set_ylim(0, global_fl_max * 1.1)


        if i == 0:
            ax_sed.set_ylabel(r'$F_\lambda$', fontsize=8)
        else:
            ax_sed.set_yticklabels([])

        if i == 1:
            ax_sed.set_xlabel(r'Wavelength ($\mathrm{\AA}$)', fontsize=8)

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
    
    print(f"Reading: {history_file}")
    history = read_mesa_history(history_file)
    
    print_summary(history)
    
    # Generate figures
    print("Generating figures...")
    
    plot_cmd_with_seds(history, sed_dir, output_dir)
    
    plot_evolution_panels(history, sed_dir, output_dir)
    
    print(f"\nDone! Figures saved to: {output_dir}")


if __name__ == '__main__':
    main()