#!/usr/bin/env python3
"""
plot.py - Analysis script for MESA Colors TP-AGB demonstration

Generates figures showing photometric evolution during thermal pulses
on the Asymptotic Giant Branch.

Produces:
1. Multi-band light curves with thermal pulse identification
2. Color-magnitude diagram showing pulse loops
3. HR diagram track with time coloring
4. Thermal pulse diagnostics (He luminosity, core mass)

Author: Niall Miller (2025)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path
from scipy.signal import find_peaks

# Try to use mesa_reader if available
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
    
    # Find header line
    header_idx = None
    for i, line in enumerate(lines):
        if 'model_number' in line and not line.strip().startswith('!'):
            header_idx = i
            break
    
    if header_idx is None:
        raise ValueError("Could not find header in history file")
    
    names = lines[header_idx].strip().split()
    
    data = {name: [] for name in names}
    
    for line in lines[header_idx + 1:]:
        parts = line.strip().split()
        if len(parts) == len(names):
            for i, name in enumerate(names):
                try:
                    data[name].append(float(parts[i]))
                except ValueError:
                    data[name].append(np.nan)
    
    for name in names:
        data[name] = np.array(data[name])
    
    return data


def load_history(logs_dir='LOGS'):
    """Load MESA history data."""
    history_file = Path(logs_dir) / 'history.data'
    
    if USE_MESA_READER:
        return mr.MesaData(str(history_file))
    else:
        return read_history_manual(history_file)


def get_column(history, name):
    """Get column from history object."""
    if USE_MESA_READER:
        try:
            return history.data(name)
        except:
            return None
    else:
        return history.get(name)


def has_column(history, name):
    """Check if column exists."""
    if USE_MESA_READER:
        return name in history.bulk_names
    else:
        return name in history


def identify_thermal_pulses(history, min_prominence=0.5):
    """
    Identify thermal pulse events from He-burning luminosity peaks.
    
    Returns indices and times of pulse maxima.
    """
    log_LHe = get_column(history, 'log_LHe')
    age = get_column(history, 'star_age')
    
    if log_LHe is None or age is None:
        return [], []
    
    # Find peaks in He luminosity
    # Thermal pulses show log_LHe spikes of 2+ orders of magnitude
    valid = np.isfinite(log_LHe)
    
    if np.sum(valid) < 10:
        return [], []
    
    peaks, properties = find_peaks(log_LHe[valid], 
                                   prominence=min_prominence,
                                   distance=50)  # Minimum separation
    
    # Map back to original indices
    valid_indices = np.where(valid)[0]
    pulse_indices = valid_indices[peaks]
    pulse_times = age[pulse_indices]
    
    return pulse_indices, pulse_times


def plot_lightcurves(history, outdir='plots'):
    """
    Plot multi-band light curves during thermal pulses.
    """
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    
    age = get_column(history, 'star_age')
    time_kyr = age / 1e3  # Convert to kyr
    
    # Identify thermal pulses
    pulse_idx, pulse_times = identify_thermal_pulses(history)
    
    # Top panel: Filter magnitudes
    ax1 = axes[0]
    
    bands = ['B', 'V', 'R', 'I', 'U', 'J', 'H', 'K']
    colors = {'U': 'violet', 'B': 'blue', 'V': 'green', 'R': 'orange', 
              'I': 'red', 'J': 'darkred', 'H': 'brown', 'K': 'black'}
    
    plotted = []
    for band in bands:
        if has_column(history, band):
            mag = get_column(history, band)
            valid = np.isfinite(mag)
            if np.sum(valid) > 10:
                ax1.plot(time_kyr[valid], mag[valid], '-', 
                        color=colors.get(band, 'gray'),
                        label=band, linewidth=1)
                plotted.append(band)
    
    # Mark thermal pulses
    for pt in pulse_times:
        ax1.axvline(pt/1e3, color='gray', linestyle=':', alpha=0.5, linewidth=0.5)
    
    ax1.set_ylabel('Magnitude', fontsize=12)
    ax1.invert_yaxis()
    ax1.legend(loc='upper right', ncol=len(plotted)//2 + 1, fontsize=9)
    ax1.set_title('TP-AGB: Multi-band Light Curves', fontsize=14)
    ax1.grid(True, alpha=0.3)
    
    # Bottom panel: Bolometric magnitude
    ax2 = axes[1]
    
    if has_column(history, 'Mag_bol'):
        mag_bol = get_column(history, 'Mag_bol')
        valid = np.isfinite(mag_bol)
        ax2.plot(time_kyr[valid], mag_bol[valid], 'k-', linewidth=1.5, 
                label='Bolometric')
    
    # Also show log_L for comparison
    if has_column(history, 'log_L'):
        log_L = get_column(history, 'log_L')
        mbol_approx = 4.74 - 2.5 * log_L
        ax2.plot(time_kyr, mbol_approx, '--', color='gray', linewidth=1,
                label=r'From $\log L$', alpha=0.7)
    
    for pt in pulse_times:
        ax2.axvline(pt/1e3, color='gray', linestyle=':', alpha=0.5, linewidth=0.5)
    
    ax2.set_xlabel('Time (kyr)', fontsize=12)
    ax2.set_ylabel('Bolometric Magnitude', fontsize=12)
    ax2.invert_yaxis()
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)
    
    # Add pulse count annotation
    if len(pulse_times) > 0:
        ax1.text(0.02, 0.98, f'{len(pulse_times)} thermal pulses identified',
                transform=ax1.transAxes, fontsize=10,
                verticalalignment='top', 
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    
    Path(outdir).mkdir(exist_ok=True)
    plt.savefig(f'{outdir}/tp_agb_lightcurves.png', dpi=150, bbox_inches='tight')
    plt.savefig(f'{outdir}/tp_agb_lightcurves.pdf', bbox_inches='tight')
    print(f"Saved: {outdir}/tp_agb_lightcurves.png")
    plt.close()


def plot_color_magnitude(history, outdir='plots'):
    """
    Plot color-magnitude diagram showing thermal pulse loops.
    """
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Try B-V vs V first, then V-R vs V
    color_pairs = [('B', 'V'), ('V', 'R'), ('V', 'I'), ('R', 'I')]
    
    color_name = None
    color_data = None
    mag_data = None
    mag_name = None
    
    for blue, red in color_pairs:
        if has_column(history, blue) and has_column(history, red):
            blue_mag = get_column(history, blue)
            red_mag = get_column(history, red)
            
            valid = np.isfinite(blue_mag) & np.isfinite(red_mag)
            if np.sum(valid) > 10:
                color_data = blue_mag - red_mag
                mag_data = red_mag
                color_name = f'{blue}-{red}'
                mag_name = red
                break
    
    if color_data is None:
        print("Warning: No suitable filter pairs found for CMD")
        plt.close()
        return
    
    age = get_column(history, 'star_age')
    time_kyr = age / 1e3
    
    valid = np.isfinite(color_data) & np.isfinite(mag_data)
    
    # Color points by time
    scatter = ax.scatter(color_data[valid], mag_data[valid],
                        c=time_kyr[valid], cmap='viridis',
                        s=10, alpha=0.6, edgecolors='none')
    
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Time (kyr)', fontsize=12)
    
    # Mark start and end
    ax.plot(color_data[valid][0], mag_data[valid][0], 'go', markersize=12,
           label='Start', zorder=10)
    ax.plot(color_data[valid][-1], mag_data[valid][-1], 'rs', markersize=12,
           label='End', zorder=10)
    
    ax.set_xlabel(f'{color_name}', fontsize=14)
    ax.set_ylabel(f'{mag_name}', fontsize=14)
    ax.invert_yaxis()
    ax.legend(loc='upper left')
    ax.set_title('TP-AGB Color-Magnitude Evolution', fontsize=14)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    Path(outdir).mkdir(exist_ok=True)
    plt.savefig(f'{outdir}/tp_agb_cmd.png', dpi=150, bbox_inches='tight')
    plt.savefig(f'{outdir}/tp_agb_cmd.pdf', bbox_inches='tight')
    print(f"Saved: {outdir}/tp_agb_cmd.png")
    plt.close()


def plot_pulse_diagnostics(history, outdir='plots'):
    """
    Plot thermal pulse diagnostics: He luminosity, core mass, etc.
    """
    fig, axes = plt.subplots(4, 1, figsize=(12, 12), sharex=True)
    
    age = get_column(history, 'star_age')
    time_kyr = age / 1e3
    
    pulse_idx, pulse_times = identify_thermal_pulses(history)
    
    # Panel 1: He-burning luminosity (the pulse indicator)
    ax1 = axes[0]
    if has_column(history, 'log_LHe'):
        log_LHe = get_column(history, 'log_LHe')
        valid = np.isfinite(log_LHe)
        ax1.plot(time_kyr[valid], log_LHe[valid], 'r-', linewidth=1)
        
        # Mark identified pulses
        for pi in pulse_idx:
            ax1.plot(time_kyr[pi], log_LHe[pi], 'ko', markersize=5)
    
    ax1.set_ylabel(r'$\log(L_{\rm He}/L_\odot)$', fontsize=12)
    ax1.set_title('Thermal Pulse Diagnostics', fontsize=14)
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Surface luminosity and Teff
    ax2 = axes[1]
    if has_column(history, 'log_L'):
        log_L = get_column(history, 'log_L')
        ax2.plot(time_kyr, log_L, 'b-', linewidth=1, label=r'$\log L$')
    ax2.set_ylabel(r'$\log(L/L_\odot)$', fontsize=12, color='blue')
    ax2.tick_params(axis='y', labelcolor='blue')
    
    ax2b = ax2.twinx()
    if has_column(history, 'log_Teff'):
        log_Teff = get_column(history, 'log_Teff')
        ax2b.plot(time_kyr, log_Teff, 'r-', linewidth=1, label=r'$\log T_{\rm eff}$')
    ax2b.set_ylabel(r'$\log T_{\rm eff}$', fontsize=12, color='red')
    ax2b.tick_params(axis='y', labelcolor='red')
    ax2.grid(True, alpha=0.3)
    
    # Panel 3: Core masses
    ax3 = axes[2]
    if has_column(history, 'he_core_mass'):
        he_core = get_column(history, 'he_core_mass')
        ax3.plot(time_kyr, he_core, 'k-', linewidth=1.5, label='He core')
    if has_column(history, 'c_core_mass'):
        c_core = get_column(history, 'c_core_mass')
        ax3.plot(time_kyr, c_core, 'b--', linewidth=1, label='C/O core')
    ax3.set_ylabel(r'Core mass ($M_\odot$)', fontsize=12)
    ax3.legend(loc='upper left')
    ax3.grid(True, alpha=0.3)
    
    # Panel 4: Visual magnitude (if available)
    ax4 = axes[3]
    if has_column(history, 'V'):
        V = get_column(history, 'V')
        valid = np.isfinite(V)
        ax4.plot(time_kyr[valid], V[valid], 'g-', linewidth=1)
    ax4.set_ylabel('V magnitude', fontsize=12)
    ax4.invert_yaxis()
    ax4.set_xlabel('Time (kyr)', fontsize=12)
    ax4.grid(True, alpha=0.3)
    
    # Mark pulses on all panels
    for ax in axes:
        for pt in pulse_times:
            ax.axvline(pt/1e3, color='gray', linestyle=':', alpha=0.3, linewidth=0.5)
    
    plt.tight_layout()
    
    Path(outdir).mkdir(exist_ok=True)
    plt.savefig(f'{outdir}/tp_agb_diagnostics.png', dpi=150, bbox_inches='tight')
    plt.savefig(f'{outdir}/tp_agb_diagnostics.pdf', bbox_inches='tight')
    print(f"Saved: {outdir}/tp_agb_diagnostics.png")
    plt.close()


def plot_hr_diagram(history, outdir='plots'):
    """
    Plot HR diagram track colored by time.
    """
    fig, ax = plt.subplots(figsize=(8, 8))
    
    log_Teff = get_column(history, 'log_Teff')
    log_L = get_column(history, 'log_L')
    age = get_column(history, 'star_age')
    
    if log_Teff is None or log_L is None:
        print("Warning: Cannot plot HR diagram - missing log_Teff or log_L")
        plt.close()
        return
    
    time_kyr = age / 1e3
    
    scatter = ax.scatter(log_Teff, log_L, c=time_kyr, cmap='viridis',
                        s=10, alpha=0.6, edgecolors='none')
    
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Time (kyr)', fontsize=12)
    
    # Mark start and end
    ax.plot(log_Teff[0], log_L[0], 'go', markersize=12, label='Start', zorder=10)
    ax.plot(log_Teff[-1], log_L[-1], 'rs', markersize=12, label='End', zorder=10)
    
    ax.set_xlabel(r'$\log T_{\rm eff}$ (K)', fontsize=14)
    ax.set_ylabel(r'$\log(L/L_\odot)$', fontsize=14)
    ax.invert_xaxis()
    ax.legend(loc='lower left')
    ax.set_title('TP-AGB HR Diagram Track', fontsize=14)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    Path(outdir).mkdir(exist_ok=True)
    plt.savefig(f'{outdir}/tp_agb_hr.png', dpi=150, bbox_inches='tight')
    plt.savefig(f'{outdir}/tp_agb_hr.pdf', bbox_inches='tight')
    print(f"Saved: {outdir}/tp_agb_hr.png")
    plt.close()


def plot_summary_grid(history, outdir='plots'):
    """
    Create publication-ready 4-panel summary figure.
    """
    fig = plt.figure(figsize=(12, 10))
    gs = GridSpec(2, 2, hspace=0.25, wspace=0.3)
    
    age = get_column(history, 'star_age')
    time_kyr = age / 1e3
    
    pulse_idx, pulse_times = identify_thermal_pulses(history)
    
    # Panel A: Multi-band light curves
    ax1 = fig.add_subplot(gs[0, 0])
    
    for band, color in [('B', 'blue'), ('V', 'green'), ('R', 'orange'), ('I', 'red')]:
        if has_column(history, band):
            mag = get_column(history, band)
            valid = np.isfinite(mag)
            if np.sum(valid) > 10:
                ax1.plot(time_kyr[valid], mag[valid], '-', color=color,
                        label=band, linewidth=1)
    
    ax1.set_xlabel('Time (kyr)')
    ax1.set_ylabel('Magnitude')
    ax1.invert_yaxis()
    ax1.legend(loc='upper right', ncol=2, fontsize=9)
    ax1.set_title('(a) Multi-band Light Curves')
    ax1.grid(True, alpha=0.3)
    
    # Panel B: Color-magnitude diagram
    ax2 = fig.add_subplot(gs[0, 1])
    
    if has_column(history, 'B') and has_column(history, 'V'):
        B = get_column(history, 'B')
        V = get_column(history, 'V')
        BmV = B - V
        valid = np.isfinite(BmV) & np.isfinite(V)
        
        scatter = ax2.scatter(BmV[valid], V[valid], c=time_kyr[valid], 
                             cmap='viridis', s=10, alpha=0.6)
        cbar = plt.colorbar(scatter, ax=ax2)
        cbar.set_label('Time (kyr)', fontsize=10)
        
        ax2.plot(BmV[valid][0], V[valid][0], 'go', markersize=10, zorder=10)
        ax2.plot(BmV[valid][-1], V[valid][-1], 'rs', markersize=10, zorder=10)
    
    ax2.set_xlabel('B-V')
    ax2.set_ylabel('V')
    ax2.invert_yaxis()
    ax2.set_title('(b) Color-Magnitude Diagram')
    ax2.grid(True, alpha=0.3)
    
    # Panel C: HR diagram
    ax3 = fig.add_subplot(gs[1, 0])
    
    log_Teff = get_column(history, 'log_Teff')
    log_L = get_column(history, 'log_L')
    
    if log_Teff is not None and log_L is not None:
        scatter = ax3.scatter(log_Teff, log_L, c=time_kyr, cmap='viridis',
                             s=10, alpha=0.6)
        ax3.plot(log_Teff[0], log_L[0], 'go', markersize=10, zorder=10)
        ax3.plot(log_Teff[-1], log_L[-1], 'rs', markersize=10, zorder=10)
    
    ax3.set_xlabel(r'$\log T_{\rm eff}$ (K)')
    ax3.set_ylabel(r'$\log(L/L_\odot)$')
    ax3.invert_xaxis()
    ax3.set_title('(c) HR Diagram Track')
    ax3.grid(True, alpha=0.3)
    
    # Panel D: He luminosity showing pulses
    ax4 = fig.add_subplot(gs[1, 1])
    
    if has_column(history, 'log_LHe'):
        log_LHe = get_column(history, 'log_LHe')
        valid = np.isfinite(log_LHe)
        ax4.plot(time_kyr[valid], log_LHe[valid], 'r-', linewidth=1)
        
        for pi in pulse_idx:
            ax4.plot(time_kyr[pi], log_LHe[pi], 'ko', markersize=4)
    
    ax4.set_xlabel('Time (kyr)')
    ax4.set_ylabel(r'$\log(L_{\rm He}/L_\odot)$')
    ax4.set_title(f'(d) He Luminosity ({len(pulse_idx)} pulses)')
    ax4.grid(True, alpha=0.3)
    
    plt.suptitle('MESA Colors: Thermal Pulse AGB Evolution', fontsize=14, y=1.02)
    
    Path(outdir).mkdir(exist_ok=True)
    plt.savefig(f'{outdir}/tp_agb_summary.png', dpi=150, bbox_inches='tight')
    plt.savefig(f'{outdir}/tp_agb_summary.pdf', bbox_inches='tight')
    print(f"Saved: {outdir}/tp_agb_summary.png")
    plt.close()


def main():
    """Main analysis routine."""
    print("=" * 60)
    print("MESA Colors TP-AGB Analysis")
    print("=" * 60)
    
    # Try to find LOGS directory
    logs_dirs = ['LOGS', 'LOGS_TP', '../LOGS']
    history = None
    
    for logs_dir in logs_dirs:
        if Path(logs_dir).exists() and (Path(logs_dir) / 'history.data').exists():
            try:
                history = load_history(logs_dir)
                print(f"Loaded history from {logs_dir}")
                break
            except Exception as e:
                print(f"Could not load from {logs_dir}: {e}")
    
    if history is None:
        print("\nERROR: Could not find LOGS directory with history.data")
        print("Make sure you have run the MESA model first (./rn)")
        return
    
    # Check for Colors columns
    if USE_MESA_READER:
        cols = history.bulk_names
    else:
        cols = list(history.keys())
    
    colors_cols = [c for c in cols if c in 
                   ['Mag_bol', 'Flux_bol', 'Interp_rad',
                    'U', 'B', 'V', 'R', 'I', 'J', 'H', 'K',
                    'u', 'g', 'r', 'i', 'z', 'y']]
    
    print(f"\nFound {len(cols)} total columns")
    if colors_cols:
        print(f"Colors columns: {', '.join(colors_cols)}")
    else:
        print("WARNING: No Colors columns found!")
        print("Check that the Colors module was enabled in Stage 2")
    
    # Identify thermal pulses
    pulse_idx, pulse_times = identify_thermal_pulses(history)
    print(f"\nIdentified {len(pulse_idx)} thermal pulses")
    
    # Generate plots
    print("\nGenerating plots...")
    outdir = 'plots'
    
    plot_lightcurves(history, outdir)
    plot_color_magnitude(history, outdir)
    plot_hr_diagram(history, outdir)
    plot_pulse_diagnostics(history, outdir)
    plot_summary_grid(history, outdir)
    
    print("\n" + "=" * 60)
    print("Analysis complete!")
    print("=" * 60)


if __name__ == '__main__':
    main()
