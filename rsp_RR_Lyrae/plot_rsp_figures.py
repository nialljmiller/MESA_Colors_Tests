#!/usr/bin/env python3
"""
Generate figures for MESA Colors RSP section (ApJ paper)

Usage:
    python plot_rsp_figures.py LOGS/history.data

Generates:
    1. rsp_lightcurve.pdf - Multi-band light curves vs phase
    2. rsp_colorloop.pdf  - Color-magnitude loop (B-V vs V)

Author: Generated for Miller, Joyce, Mocz et al. MESA Colors paper
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import sys
import os

# Publication-quality settings
plt.rcParams.update({
    'font.size': 11,
    'font.family': 'serif',
    'axes.labelsize': 12,
    'axes.titlesize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.figsize': (3.5, 3.5),  # Single column ApJ
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})


def load_mesa_history(filepath):
    """
    Load MESA history.data file
    Returns dict-like object with column names as attributes
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Find header line (contains column names)
    header_idx = None
    for i, line in enumerate(lines):
        if line.strip().startswith('model_number') or 'model_number' in line:
            header_idx = i
            break
    
    if header_idx is None:
        # Try standard MESA format: line 6 has headers
        header_idx = 5
    
    # Parse column names
    cols = lines[header_idx].split()
    
    # Load data (skip header rows)
    data = np.genfromtxt(filepath, skip_header=header_idx+1, names=cols)
    
    return data


def find_magnitude_columns(data):
    """
    Find magnitude column names in the history data
    MESA Colors may use different naming conventions
    """
    colnames = data.dtype.names
    
    # Possible naming patterns for Johnson bands
    patterns = {
        'B': ['B', 'abs_mag_B', 'mag_B', 'Johnson_B'],
        'V': ['V', 'abs_mag_V', 'mag_V', 'Johnson_V'],
        'R': ['R', 'abs_mag_R', 'mag_R', 'Johnson_R'],
        'I': ['I', 'abs_mag_I', 'mag_I', 'Johnson_I'],
        'U': ['U', 'abs_mag_U', 'mag_U', 'Johnson_U'],
    }
    
    found = {}
    for band, candidates in patterns.items():
        for cand in candidates:
            if cand in colnames:
                found[band] = cand
                break
    
    return found


def get_column(data, name, default=None):
    """Safely get a column from data"""
    if name in data.dtype.names:
        return data[name]
    return default


def plot_lightcurve(data, output_file='rsp_lightcurve.pdf'):
    """
    Plot multi-band light curves vs pulsation phase
    """
    fig, ax = plt.subplots(figsize=(7, 5))
    
    # Get time/phase data
    if 'rsp_phase' in data.dtype.names:
        phase = data['rsp_phase']
    else:
        # Compute phase from time
        time = get_column(data, 'star_age_day', get_column(data, 'star_age', None))
        if time is None:
            print("Error: No time column found")
            return
        period = 0.71  # days
        phase = ((time - time[0]) / period) % 1.0
    
    # Find magnitude columns
    mag_cols = find_magnitude_columns(data)
    
    if not mag_cols:
        print("Warning: No magnitude columns found. Available columns:")
        print(data.dtype.names)
        return
    
    print(f"Found magnitude columns: {mag_cols}")
    
    # Define plotting parameters
    band_info = {
        'B': {'color': '#0066CC', 'offset': 0.0, 'label': '$B$'},
        'V': {'color': '#33AA33', 'offset': 1.2, 'label': '$V$'},
        'R': {'color': '#CC6600', 'offset': 2.4, 'label': '$R$'},
        'I': {'color': '#990033', 'offset': 3.6, 'label': '$I$'},
    }
    
    # Select data for ~2-3 complete cycles after initial settling
    # Skip first ~3 cycles to avoid transients
    if 'rsp_num_periods' in data.dtype.names:
        mask = (data['rsp_num_periods'] >= 3) & (data['rsp_num_periods'] <= 5)
    else:
        # Use time-based selection
        time = get_column(data, 'star_age_day', data['star_age'])
        mask = (time > 2.0) & (time < 4.5)  # ~3 cycles
    
    if mask.sum() < 100:
        print("Warning: Few points in selected range, using all data")
        mask = np.ones(len(data), dtype=bool)
    
    # Plot each band
    for band in ['B', 'V', 'R', 'I']:
        if band in mag_cols:
            col = mag_cols[band]
            mag = data[col][mask]
            ph = phase[mask]
            info = band_info[band]
            
            # Sort by phase for cleaner plotting
            sort_idx = np.argsort(ph)
            
            ax.plot(ph[sort_idx], mag[sort_idx] + info['offset'], 
                   '.', ms=1.5, alpha=0.6, color=info['color'], 
                   label=info['label'])
    
    ax.set_xlabel('Pulsation Phase')
    ax.set_ylabel('Magnitude + offset')
    ax.set_xlim(0, 1)
    ax.invert_yaxis()
    ax.legend(loc='upper right', markerscale=4, framealpha=0.9)
    
    # Add annotation about offsets
    ax.text(0.02, 0.02, 'Offsets: +1.2 mag between bands',
            transform=ax.transAxes, fontsize=8, alpha=0.7)
    
    ax.grid(True, alpha=0.2)
    
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    print(f"Saved {output_file}")


def plot_colorloop(data, output_file='rsp_colorloop.pdf'):
    """
    Plot color-magnitude loop (B-V vs V)
    """
    fig, ax = plt.subplots(figsize=(5, 5))
    
    # Find magnitude columns
    mag_cols = find_magnitude_columns(data)
    
    if 'B' not in mag_cols or 'V' not in mag_cols:
        print("Error: Need B and V magnitudes for color loop")
        return
    
    B = data[mag_cols['B']]
    V = data[mag_cols['V']]
    color = B - V
    
    # Get phase for coloring
    if 'rsp_phase' in data.dtype.names:
        phase = data['rsp_phase']
    else:
        time = get_column(data, 'star_age_day', data['star_age'])
        period = 0.71
        phase = ((time - time[0]) / period) % 1.0
    
    # Select one complete cycle (after settling)
    if 'rsp_num_periods' in data.dtype.names:
        mask = (data['rsp_num_periods'] >= 4) & (data['rsp_num_periods'] < 5)
    else:
        time = get_column(data, 'star_age_day', data['star_age'])
        mask = (time > 3.0) & (time < 3.71)
    
    if mask.sum() < 50:
        # Fall back to using phase
        mask = np.ones(len(data), dtype=bool)
        # Take points from one cycle based on phase
        unique_phases = []
        unique_idx = []
        for i, p in enumerate(phase):
            if len(unique_phases) == 0 or abs(p - unique_phases[-1]) > 0.001:
                unique_phases.append(p)
                unique_idx.append(i)
            if len(unique_idx) > 1000:
                break
        mask = np.zeros(len(data), dtype=bool)
        mask[unique_idx[:1000]] = True
    
    # Create scatter plot colored by phase
    sc = ax.scatter(color[mask], V[mask], 
                   c=phase[mask], cmap='twilight_shifted', 
                   s=8, alpha=0.8, edgecolors='none')
    
    # Add colorbar
    cbar = plt.colorbar(sc, ax=ax, label='Pulsation Phase', shrink=0.8)
    
    # Add arrow to show loop direction
    c_arr = color[mask]
    v_arr = V[mask]
    p_arr = phase[mask]
    
    # Find point near phase 0.25 to place arrow
    idx_arrow = np.argmin(np.abs(p_arr - 0.25))
    idx_arrow2 = np.argmin(np.abs(p_arr - 0.30))
    
    if idx_arrow != idx_arrow2:
        ax.annotate('', 
                    xy=(c_arr[idx_arrow2], v_arr[idx_arrow2]),
                    xytext=(c_arr[idx_arrow], v_arr[idx_arrow]),
                    arrowprops=dict(arrowstyle='->', color='black', lw=1.5))
    
    ax.set_xlabel('$B - V$ (mag)')
    ax.set_ylabel('$V$ (mag)')
    ax.invert_yaxis()
    
    # Mark key phases
    for target_phase, label in [(0.0, 'max'), (0.5, 'min')]:
        idx = np.argmin(np.abs(p_arr - target_phase))
        ax.plot(c_arr[idx], v_arr[idx], 'k*', ms=10, zorder=10)
    
    ax.grid(True, alpha=0.2)
    
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    print(f"Saved {output_file}")


def plot_amplitude_spectrum(data, output_file='rsp_amplitudes.pdf'):
    """
    Optional: Plot pulsation amplitude vs wavelength
    """
    fig, ax = plt.subplots(figsize=(5, 3.5))
    
    mag_cols = find_magnitude_columns(data)
    
    # Effective wavelengths (Angstroms)
    wavelengths = {'U': 3650, 'B': 4400, 'V': 5500, 'R': 6400, 'I': 8000}
    colors_plot = {'U': '#6600CC', 'B': '#0066CC', 'V': '#33AA33', 
                   'R': '#CC6600', 'I': '#990033'}
    
    wl_plot = []
    amp_plot = []
    col_plot = []
    labels = []
    
    for band in ['U', 'B', 'V', 'R', 'I']:
        if band in mag_cols:
            mag = data[mag_cols[band]]
            amp = np.max(mag) - np.min(mag)
            wl_plot.append(wavelengths[band])
            amp_plot.append(amp)
            col_plot.append(colors_plot[band])
            labels.append(band)
    
    ax.scatter(wl_plot, amp_plot, c=col_plot, s=100, zorder=5, edgecolors='k')
    ax.plot(wl_plot, amp_plot, 'k--', alpha=0.4, zorder=1)
    
    for w, a, l in zip(wl_plot, amp_plot, labels):
        ax.annotate(l, (w, a), xytext=(5, 5), textcoords='offset points', fontsize=10)
    
    ax.set_xlabel('Effective Wavelength (Å)')
    ax.set_ylabel('Pulsation Amplitude (mag)')
    ax.grid(True, alpha=0.2)
    
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    print(f"Saved {output_file}")


def print_statistics(data):
    """Print useful statistics about the run"""
    mag_cols = find_magnitude_columns(data)
    
    print("\n" + "="*60)
    print("RSP + Colors Run Statistics")
    print("="*60)
    
    n_models = len(data)
    print(f"Total models: {n_models}")
    
    if 'rsp_num_periods' in data.dtype.names:
        n_periods = int(np.max(data['rsp_num_periods']))
        print(f"Pulsation cycles completed: {n_periods}")
        print(f"Points per cycle: ~{n_models // max(n_periods, 1)}")
    
    if 'rsp_period_in_days' in data.dtype.names:
        periods = data['rsp_period_in_days']
        valid = periods > 0
        if valid.any():
            print(f"Mean period: {np.mean(periods[valid]):.5f} days")
    
    print("\nMagnitude ranges:")
    for band in ['B', 'V', 'R', 'I']:
        if band in mag_cols:
            mag = data[mag_cols[band]]
            print(f"  {band}: {np.min(mag):.3f} to {np.max(mag):.3f} "
                  f"(Δ = {np.max(mag)-np.min(mag):.3f} mag)")
    
    if 'B' in mag_cols and 'V' in mag_cols:
        color = data[mag_cols['B']] - data[mag_cols['V']]
        print(f"\n(B-V) color range: {np.min(color):.3f} to {np.max(color):.3f}")
    
    print("="*60 + "\n")


def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_rsp_figures.py LOGS/history.data")
        print("\nThis script generates publication-quality figures for the")
        print("MESA Colors paper RSP demonstration section.")
        sys.exit(1)
    
    history_file = sys.argv[1]
    
    if not os.path.exists(history_file):
        print(f"Error: File not found: {history_file}")
        sys.exit(1)
    
    print(f"Loading {history_file}...")
    data = load_mesa_history(history_file)
    
    print(f"Loaded {len(data)} models")
    print(f"Columns: {', '.join(data.dtype.names[:10])}...")
    
    # Print statistics
    print_statistics(data)
    
    # Generate figures
    print("Generating figures...")
    plot_lightcurve(data, 'rsp_lightcurve.pdf')
    plot_colorloop(data, 'rsp_colorloop.pdf')
    
    # Optional: amplitude spectrum
    # plot_amplitude_spectrum(data, 'rsp_amplitudes.pdf')
    
    print("\nDone! Generated:")
    print("  - rsp_lightcurve.pdf")
    print("  - rsp_colorloop.pdf")


if __name__ == '__main__':
    main()
