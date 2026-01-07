#!/usr/bin/env python3
"""
plot.py - Analysis script for MESA Colors nova burst demonstration

Generates figures demonstrating the Colors module's ability to capture
rapid photometric evolution during classical nova outburst.

Produces:
1. Multi-band light curves showing magnitude evolution
2. Color-magnitude diagram showing temperature evolution
3. Timestep analysis showing MESA's adaptive sampling

Author: Niall Miller (2025)
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Try to use mesa_reader if available, otherwise fall back to manual parsing
try:
    import mesa_reader as mr
    USE_MESA_READER = True
except ImportError:
    USE_MESA_READER = False
    print("mesa_reader not found, using manual history parsing")


def read_history_manual(filepath):
    """Read MESA history file without mesa_reader."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Find header line (starts with column names)
    header_idx = None
    for i, line in enumerate(lines):
        if line.strip().startswith('model_number'):
            header_idx = i
            break
    
    if header_idx is None:
        # Try alternative: look for line after blank line following header info
        for i, line in enumerate(lines):
            parts = line.strip().split()
            if len(parts) > 0 and parts[0].isdigit():
                header_idx = i - 1
                break
    
    if header_idx is None:
        raise ValueError("Could not find header in history file")
    
    # Parse column names
    names = lines[header_idx].strip().split()
    
    # Parse data
    data = {}
    for name in names:
        data[name] = []
    
    for line in lines[header_idx + 1:]:
        parts = line.strip().split()
        if len(parts) == len(names):
            for i, name in enumerate(names):
                try:
                    data[name].append(float(parts[i]))
                except ValueError:
                    data[name].append(np.nan)
    
    # Convert to numpy arrays
    for name in names:
        data[name] = np.array(data[name])
    
    return data


def load_history(logs_dir='LOGS'):
    """Load MESA history data."""
    history_file = Path(logs_dir) / 'history.data'
    
    if USE_MESA_READER:
        h = mr.MesaData(str(history_file))
        return h
    else:
        return read_history_manual(history_file)


def get_column(history, name):
    """Get column from history object (works with mesa_reader or dict)."""
    if USE_MESA_READER:
        return history.data(name)
    else:
        return history.get(name, np.array([]))


def has_column(history, name):
    """Check if column exists."""
    if USE_MESA_READER:
        return name in history.bulk_names
    else:
        return name in history


def plot_multiband_lightcurves(history, outdir='plots'):
    """
    Plot multi-band light curves during nova outburst.
    
    Demonstrates: MESA Colors captures rapid magnitude evolution
    at native timestep resolution.
    """
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # Get time axis
    age = get_column(history, 'star_age')  # years
    # Convert to days from start of nova phase
    time_days = (age - age[0]) * 365.25
    
    # Top panel: LSST bands (or whatever filters are available)
    ax1 = axes[0]
    
    # Check which filters are available
    lsst_bands = ['u', 'g', 'r', 'i', 'z', 'y']
    johnson_bands = ['U', 'B', 'V', 'R', 'I']
    
    colors_lsst = ['purple', 'green', 'orange', 'red', 'brown', 'black']
    colors_johnson = ['purple', 'blue', 'green', 'red', 'darkred']
    
    bands_found = []
    
    # Try LSST bands first
    for band, color in zip(lsst_bands, colors_lsst):
        if has_column(history, band):
            mag = get_column(history, band)
            valid = np.isfinite(mag)
            if np.sum(valid) > 10:
                ax1.plot(time_days[valid], mag[valid], '-', color=color, 
                        label=band, linewidth=1.5)
                bands_found.append(band)
    
    # Try Johnson bands if no LSST
    if not bands_found:
        for band, color in zip(johnson_bands, colors_johnson):
            if has_column(history, band):
                mag = get_column(history, band)
                valid = np.isfinite(mag)
                if np.sum(valid) > 10:
                    ax1.plot(time_days[valid], mag[valid], '-', color=color,
                            label=band, linewidth=1.5)
                    bands_found.append(band)
    
    ax1.set_ylabel('Absolute Magnitude', fontsize=12)
    ax1.invert_yaxis()  # Brighter = up
    ax1.legend(loc='upper right', ncol=2)
    ax1.set_title('Nova Outburst: Multi-band Light Curves', fontsize=14)
    ax1.grid(True, alpha=0.3)
    
    # Bottom panel: Bolometric
    ax2 = axes[1]
    
    if has_column(history, 'Mag_bol'):
        mag_bol = get_column(history, 'Mag_bol')
        valid = np.isfinite(mag_bol)
        ax2.plot(time_days[valid], mag_bol[valid], 'k-', linewidth=2, 
                label='Bolometric')
    
    # Also show log_L for comparison
    if has_column(history, 'log_L'):
        log_L = get_column(history, 'log_L')
        # Convert to approximate Mbol: Mbol_sun = 4.74, L_sun = 1
        # Mbol = 4.74 - 2.5 * log_L
        mbol_approx = 4.74 - 2.5 * log_L
        ax2.plot(time_days, mbol_approx, '--', color='gray', linewidth=1,
                label=r'From $\log L$', alpha=0.7)
    
    ax2.set_xlabel('Time (days)', fontsize=12)
    ax2.set_ylabel('Bolometric Magnitude', fontsize=12)
    ax2.invert_yaxis()
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    Path(outdir).mkdir(exist_ok=True)
    plt.savefig(f'{outdir}/nova_lightcurves.png', dpi=150, bbox_inches='tight')
    plt.savefig(f'{outdir}/nova_lightcurves.pdf', bbox_inches='tight')
    print(f"Saved: {outdir}/nova_lightcurves.png")
    plt.close()


def plot_color_magnitude(history, outdir='plots'):
    """
    Plot color-magnitude diagram showing temperature evolution.
    
    Demonstrates: Color evolution during rapid transient phase.
    """
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Try different color combinations based on available bands
    color_pairs = [
        ('g', 'r', 'g-r'),      # LSST
        ('r', 'i', 'r-i'),      # LSST
        ('B', 'V', 'B-V'),      # Johnson
        ('V', 'R', 'V-R'),      # Johnson
    ]
    
    mag_band = None
    color_name = None
    color_data = None
    
    for blue, red, cname in color_pairs:
        if has_column(history, blue) and has_column(history, red):
            blue_mag = get_column(history, blue)
            red_mag = get_column(history, red)
            
            valid = np.isfinite(blue_mag) & np.isfinite(red_mag)
            if np.sum(valid) > 10:
                color_data = blue_mag - red_mag
                mag_band = red
                color_name = cname
                break
    
    if color_data is None:
        print("Warning: No suitable filter pairs found for CMD")
        plt.close()
        return
    
    mag_data = get_column(history, mag_band)
    time = get_column(history, 'star_age')
    time_days = (time - time[0]) * 365.25
    
    valid = np.isfinite(color_data) & np.isfinite(mag_data)
    
    # Color by time
    scatter = ax.scatter(color_data[valid], mag_data[valid], 
                        c=time_days[valid], cmap='viridis', 
                        s=20, alpha=0.7, edgecolors='none')
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Time (days)', fontsize=12)
    
    # Mark start and end
    ax.plot(color_data[valid][0], mag_data[valid][0], 'go', markersize=15,
           label='Start', zorder=10)
    ax.plot(color_data[valid][-1], mag_data[valid][-1], 'rs', markersize=15,
           label='End', zorder=10)
    
    ax.set_xlabel(f'{color_name}', fontsize=14)
    ax.set_ylabel(f'{mag_band}', fontsize=14)
    ax.invert_yaxis()
    ax.legend(loc='upper left')
    ax.set_title('Nova Color-Magnitude Evolution', fontsize=14)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    Path(outdir).mkdir(exist_ok=True)
    plt.savefig(f'{outdir}/nova_cmd.png', dpi=150, bbox_inches='tight')
    plt.savefig(f'{outdir}/nova_cmd.pdf', bbox_inches='tight')
    print(f"Saved: {outdir}/nova_cmd.png")
    plt.close()


def plot_timestep_analysis(history, outdir='plots'):
    """
    Analyze timestep resolution during rapid evolution.
    
    Demonstrates: MESA Colors samples at native timestep resolution,
    capturing rapid phases that would be missed by coarse post-processing.
    """
    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    
    age = get_column(history, 'star_age')
    time_days = (age - age[0]) * 365.25
    
    # Panel 1: Timestep evolution
    ax1 = axes[0]
    if has_column(history, 'time_step_sec'):
        dt_sec = get_column(history, 'time_step_sec')
        ax1.semilogy(time_days, dt_sec / 3600, 'b-', linewidth=1)
        ax1.set_ylabel('Timestep (hours)', fontsize=12)
    elif has_column(history, 'log_dt'):
        log_dt = get_column(history, 'log_dt')
        dt_years = 10**log_dt
        dt_hours = dt_years * 365.25 * 24
        ax1.semilogy(time_days, dt_hours, 'b-', linewidth=1)
        ax1.set_ylabel('Timestep (hours)', fontsize=12)
    ax1.set_title('MESA Adaptive Timestep During Nova', fontsize=14)
    ax1.grid(True, alpha=0.3)
    
    # Mark typical survey cadence
    ax1.axhline(y=24, color='r', linestyle='--', alpha=0.5, label='1 day cadence')
    ax1.axhline(y=1, color='orange', linestyle='--', alpha=0.5, label='1 hour cadence')
    ax1.legend(loc='upper right')
    
    # Panel 2: Luminosity evolution
    ax2 = axes[1]
    if has_column(history, 'log_L'):
        log_L = get_column(history, 'log_L')
        ax2.plot(time_days, log_L, 'k-', linewidth=1.5)
        ax2.set_ylabel(r'$\log(L/L_\odot)$', fontsize=12)
    ax2.grid(True, alpha=0.3)
    
    # Panel 3: Temperature evolution
    ax3 = axes[2]
    if has_column(history, 'log_Teff'):
        log_Teff = get_column(history, 'log_Teff')
        ax3.plot(time_days, log_Teff, 'r-', linewidth=1.5)
        ax3.set_ylabel(r'$\log T_{\rm eff}$ (K)', fontsize=12)
    ax3.set_xlabel('Time (days)', fontsize=12)
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    Path(outdir).mkdir(exist_ok=True)
    plt.savefig(f'{outdir}/nova_timestep_analysis.png', dpi=150, bbox_inches='tight')
    plt.savefig(f'{outdir}/nova_timestep_analysis.pdf', bbox_inches='tight')
    print(f"Saved: {outdir}/nova_timestep_analysis.png")
    plt.close()


def plot_summary_grid(history, outdir='plots'):
    """
    Create a summary grid figure suitable for publication.
    """
    fig = plt.figure(figsize=(12, 10))
    
    gs = fig.add_gridspec(2, 2, hspace=0.25, wspace=0.25)
    
    age = get_column(history, 'star_age')
    time_days = (age - age[0]) * 365.25
    
    # Panel A: Multi-band light curves
    ax1 = fig.add_subplot(gs[0, 0])
    
    bands = ['u', 'g', 'r', 'i', 'z', 'y', 'U', 'B', 'V', 'R', 'I']
    colors = {'u': 'purple', 'g': 'green', 'r': 'orange', 'i': 'red', 
              'z': 'brown', 'y': 'black',
              'U': 'purple', 'B': 'blue', 'V': 'green', 'R': 'red', 'I': 'darkred'}
    
    for band in bands:
        if has_column(history, band):
            mag = get_column(history, band)
            valid = np.isfinite(mag)
            if np.sum(valid) > 10:
                ax1.plot(time_days[valid], mag[valid], '-', color=colors.get(band, 'gray'),
                        label=band, linewidth=1.5)
    
    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('Absolute Magnitude')
    ax1.invert_yaxis()
    ax1.legend(loc='upper right', ncol=2, fontsize=8)
    ax1.set_title('(a) Multi-band Light Curves')
    ax1.grid(True, alpha=0.3)
    
    # Panel B: HR diagram track
    ax2 = fig.add_subplot(gs[0, 1])
    
    if has_column(history, 'log_Teff') and has_column(history, 'log_L'):
        log_Teff = get_column(history, 'log_Teff')
        log_L = get_column(history, 'log_L')
        
        scatter = ax2.scatter(log_Teff, log_L, c=time_days, cmap='viridis',
                            s=10, alpha=0.7)
        cbar = plt.colorbar(scatter, ax=ax2)
        cbar.set_label('Time (days)', fontsize=10)
        
        ax2.plot(log_Teff[0], log_L[0], 'go', markersize=10, zorder=10)
        ax2.plot(log_Teff[-1], log_L[-1], 'rs', markersize=10, zorder=10)
    
    ax2.set_xlabel(r'$\log T_{\rm eff}$ (K)')
    ax2.set_ylabel(r'$\log(L/L_\odot)$')
    ax2.invert_xaxis()
    ax2.set_title('(b) HR Diagram Track')
    ax2.grid(True, alpha=0.3)
    
    # Panel C: Color evolution
    ax3 = fig.add_subplot(gs[1, 0])
    
    color_pairs = [('g', 'r'), ('B', 'V')]
    for blue, red in color_pairs:
        if has_column(history, blue) and has_column(history, red):
            blue_mag = get_column(history, blue)
            red_mag = get_column(history, red)
            color = blue_mag - red_mag
            valid = np.isfinite(color)
            if np.sum(valid) > 10:
                ax3.plot(time_days[valid], color[valid], '-', linewidth=1.5,
                        label=f'{blue}-{red}')
    
    ax3.set_xlabel('Time (days)')
    ax3.set_ylabel('Color (mag)')
    ax3.legend(loc='best')
    ax3.set_title('(c) Color Evolution')
    ax3.grid(True, alpha=0.3)
    
    # Panel D: Timestep resolution
    ax4 = fig.add_subplot(gs[1, 1])
    
    if has_column(history, 'time_step_sec'):
        dt_sec = get_column(history, 'time_step_sec')
        dt_hours = dt_sec / 3600
        ax4.semilogy(time_days, dt_hours, 'b-', linewidth=1)
    elif has_column(history, 'log_dt'):
        log_dt = get_column(history, 'log_dt')
        dt_hours = (10**log_dt) * 365.25 * 24
        ax4.semilogy(time_days, dt_hours, 'b-', linewidth=1)
    
    ax4.axhline(y=24, color='r', linestyle='--', alpha=0.5, label='Daily cadence')
    ax4.axhline(y=1, color='orange', linestyle='--', alpha=0.5, label='Hourly cadence')
    ax4.set_xlabel('Time (days)')
    ax4.set_ylabel('Timestep (hours)')
    ax4.legend(loc='upper right', fontsize=8)
    ax4.set_title('(d) Adaptive Timestep')
    ax4.grid(True, alpha=0.3)
    
    plt.suptitle('MESA Colors: Classical Nova Outburst', fontsize=14, y=1.02)
    
    Path(outdir).mkdir(exist_ok=True)
    plt.savefig(f'{outdir}/nova_summary.png', dpi=150, bbox_inches='tight')
    plt.savefig(f'{outdir}/nova_summary.pdf', bbox_inches='tight')
    print(f"Saved: {outdir}/nova_summary.png")
    plt.close()


def main():
    """Main analysis routine."""
    print("Loading MESA history data...")
    
    # Look for LOGS directory
    logs_dirs = ['LOGS', 'LOGS_nova', '../LOGS']
    history = None
    
    for logs_dir in logs_dirs:
        if Path(logs_dir).exists():
            try:
                history = load_history(logs_dir)
                print(f"Loaded history from {logs_dir}")
                break
            except Exception as e:
                print(f"Could not load from {logs_dir}: {e}")
    
    if history is None:
        print("ERROR: Could not find LOGS directory with history.data")
        print("Make sure you have run the MESA model first.")
        return
    
    # List available columns
    if USE_MESA_READER:
        cols = history.bulk_names
    else:
        cols = list(history.keys())
    
    print(f"\nAvailable columns ({len(cols)} total):")
    colors_cols = [c for c in cols if c in ['Mag_bol', 'Flux_bol', 'Interp_rad',
                                            'u', 'g', 'r', 'i', 'z', 'y',
                                            'U', 'B', 'V', 'R', 'I', 'J', 'H', 'K']]
    if colors_cols:
        print(f"  Colors columns: {', '.join(colors_cols)}")
    else:
        print("  WARNING: No Colors columns found - check that colors module ran")
    
    # Generate plots
    print("\nGenerating plots...")
    outdir = 'plots'
    
    plot_multiband_lightcurves(history, outdir)
    plot_color_magnitude(history, outdir)
    plot_timestep_analysis(history, outdir)
    plot_summary_grid(history, outdir)
    
    print("\nAnalysis complete!")


if __name__ == '__main__':
    main()
