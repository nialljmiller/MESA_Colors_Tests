#!/usr/bin/env python3
"""
White Dwarf Cooling Sequence Analysis Script
=============================================

Generates publication-quality figures from MESA WD cooling runs
with the colors module enabled.

Output figures:
    1. WD Cooling Track CMD (g vs g-r)
    2. Multi-band Light Curves (all LSST filters)
    3. Color Evolution vs Teff
    4. SED Evolution at selected ages (if SED files available)

Usage:
    python plot_wd_cooling.py

Requirements:
    - matplotlib
    - numpy
    - mesa_reader
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

# Try to import mesa_reader, fall back to manual parsing if not available
try:
    import mesa_reader as mr
    HAS_MESA_READER = True
except ImportError:
    HAS_MESA_READER = False
    print("Warning: mesa_reader not found. Using manual history parsing.")


# =============================================================================
# Configuration
# =============================================================================
HISTORY_FILE = 'LOGS/history.data'
SED_DIR = 'SED'
OUTPUT_DIR = 'plots'

# LSST filter names (as they appear in history.data)
LSST_FILTERS = ['u', 'g', 'r', 'i', 'z', 'y']

# Plot styling
plt.style.use('default')
mpl.rcParams['font.size'] = 12
mpl.rcParams['axes.labelsize'] = 14
mpl.rcParams['axes.titlesize'] = 14
mpl.rcParams['legend.fontsize'] = 10
mpl.rcParams['figure.figsize'] = (10, 8)
mpl.rcParams['figure.dpi'] = 150


# =============================================================================
# Data Loading
# =============================================================================
def load_history_manual(filepath):
    """
    Manually parse MESA history.data file.
    Returns dict with column names as keys.
    """
    data = {}
    with open(filepath, 'r') as f:
        # Skip header lines until we find column names
        for line in f:
            if line.strip().startswith('1'):
                # This is the column number line, next is names
                break
        
        # Column names
        col_line = f.readline()
        col_names = col_line.split()
        
        # Initialize arrays
        for col in col_names:
            data[col] = []
        
        # Read data
        for line in f:
            vals = line.split()
            if len(vals) == len(col_names):
                for i, col in enumerate(col_names):
                    try:
                        data[col].append(float(vals[i]))
                    except ValueError:
                        data[col].append(np.nan)
    
    # Convert to numpy arrays
    for col in data:
        data[col] = np.array(data[col])
    
    return data


def load_history(filepath=HISTORY_FILE):
    """
    Load MESA history file using mesa_reader or manual parsing.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"History file not found: {filepath}")
    
    if HAS_MESA_READER:
        h = mr.MesaData(filepath)
        # Convert to dict for consistent interface
        data = {}
        for col in h.bulk_names:
            data[col] = getattr(h, col)
        return data
    else:
        return load_history_manual(filepath)


def get_available_filters(data):
    """
    Detect which photometric filters are present in the data.
    """
    available = []
    for f in LSST_FILTERS:
        if f in data:
            available.append(f)
    return available


def compute_colors(data, filters):
    """
    Compute color indices from filter magnitudes.
    Returns dict of color arrays.
    """
    colors = {}
    
    # Standard colors
    color_pairs = [
        ('u', 'g'),
        ('g', 'r'),
        ('r', 'i'),
        ('i', 'z'),
        ('z', 'y'),
    ]
    
    for blue, red in color_pairs:
        if blue in data and red in data:
            color_name = f'{blue}_minus_{red}'
            colors[color_name] = data[blue] - data[red]
    
    return colors


# =============================================================================
# Figure 1: CMD - Cooling Track
# =============================================================================
def plot_cmd(data, colors, output_dir=OUTPUT_DIR):
    """
    Color-Magnitude Diagram showing WD cooling sequence.
    """
    fig, ax = plt.subplots(figsize=(8, 10))
    
    # Get data
    if 'g' not in data:
        print("Warning: g-band not found in data. Skipping CMD plot.")
        return
    
    g_mag = data['g']
    
    # Use g-r color if available
    if 'g_minus_r' in colors:
        color = colors['g_minus_r']
        xlabel = 'g - r (AB mag)'
    elif 'g' in data and 'r' in data:
        color = data['g'] - data['r']
        xlabel = 'g - r (AB mag)'
    else:
        print("Warning: Cannot compute g-r color. Skipping CMD plot.")
        return
    
    # Color by cooling time (or age)
    if 'star_age' in data:
        age = data['star_age']
        # Normalize to Gyr
        age_gyr = age / 1e9
        
        # Create colormap
        norm = Normalize(vmin=age_gyr.min(), vmax=age_gyr.max())
        cmap = plt.cm.plasma
        
        # Plot with color gradient
        sc = ax.scatter(color, g_mag, c=age_gyr, cmap=cmap, norm=norm,
                       s=5, alpha=0.8, edgecolors='none')
        
        # Colorbar
        cbar = plt.colorbar(sc, ax=ax, label='Age (Gyr)')
    else:
        ax.plot(color, g_mag, 'b-', lw=1, alpha=0.7)
    
    # Labels and formatting
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r'$M_g$ (AB mag)')
    ax.set_title('White Dwarf Cooling Sequence')
    
    # Invert y-axis (brighter = smaller magnitude)
    ax.invert_yaxis()
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    # Annotate key points
    if len(g_mag) > 10:
        # Mark start (hot WD)
        ax.annotate('Hot WD', xy=(color[0], g_mag[0]), 
                   xytext=(color[0]+0.1, g_mag[0]-0.5),
                   fontsize=10, ha='left',
                   arrowprops=dict(arrowstyle='->', color='gray'))
        
        # Mark end (cool WD)
        ax.annotate('Cool WD', xy=(color[-1], g_mag[-1]),
                   xytext=(color[-1]-0.1, g_mag[-1]+0.5),
                   fontsize=10, ha='right',
                   arrowprops=dict(arrowstyle='->', color='gray'))
    
    plt.tight_layout()
    
    # Save
    os.makedirs(output_dir, exist_ok=True)
    outfile = os.path.join(output_dir, 'wd_cmd.png')
    fig.savefig(outfile, dpi=200, bbox_inches='tight')
    print(f"Saved: {outfile}")
    plt.close(fig)


# =============================================================================
# Figure 2: Multi-band Light Curves
# =============================================================================
def plot_lightcurves(data, output_dir=OUTPUT_DIR):
    """
    Magnitude evolution in all available filters vs time.
    """
    available = get_available_filters(data)
    
    if not available:
        print("Warning: No filter magnitudes found. Skipping light curve plot.")
        return
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Get time axis
    if 'star_age' in data:
        time = data['star_age'] / 1e9  # Convert to Gyr
        xlabel = 'Age (Gyr)'
    else:
        time = np.arange(len(data[available[0]]))
        xlabel = 'Model Number'
    
    # Color palette for filters (blue to red)
    filter_colors = {
        'u': '#4B0082',  # Indigo
        'g': '#006400',  # Dark green
        'r': '#FF4500',  # Orange red
        'i': '#8B0000',  # Dark red
        'z': '#4A0000',  # Very dark red
        'y': '#2F0000',  # Near-black red
    }
    
    # Plot each filter
    for f in available:
        color = filter_colors.get(f, 'gray')
        ax.plot(time, data[f], label=f'{f}-band', color=color, lw=1.5, alpha=0.8)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Absolute Magnitude (AB)')
    ax.set_title('WD Cooling: Multi-band Light Curves')
    ax.invert_yaxis()
    ax.legend(loc='upper left')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    os.makedirs(output_dir, exist_ok=True)
    outfile = os.path.join(output_dir, 'wd_lightcurves.png')
    fig.savefig(outfile, dpi=200, bbox_inches='tight')
    print(f"Saved: {outfile}")
    plt.close(fig)


# =============================================================================
# Figure 3: Color Evolution vs Teff
# =============================================================================
def plot_color_evolution(data, colors, output_dir=OUTPUT_DIR):
    """
    Color indices as a function of effective temperature.
    """
    if 'log_Teff' not in data and 'Teff' not in data:
        print("Warning: Teff not found. Skipping color evolution plot.")
        return
    
    # Get Teff
    if 'log_Teff' in data:
        teff = 10**data['log_Teff']
    else:
        teff = data['Teff']
    
    # Filter out non-physical values
    valid = (teff > 1000) & (teff < 200000)
    teff = teff[valid]
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()
    
    color_pairs = [
        ('u_minus_g', 'u - g'),
        ('g_minus_r', 'g - r'),
        ('r_minus_i', 'r - i'),
        ('i_minus_z', 'i - z'),
    ]
    
    for ax, (color_key, label) in zip(axes, color_pairs):
        if color_key in colors:
            c = colors[color_key][valid]
            ax.scatter(teff, c, c=teff, cmap='coolwarm', s=3, alpha=0.7)
            ax.set_ylabel(label)
            ax.set_xlabel(r'$T_{\rm eff}$ (K)')
            ax.set_xscale('log')
            ax.invert_xaxis()  # Hot on left
            ax.grid(True, alpha=0.3)
        else:
            ax.text(0.5, 0.5, f'{label}\nnot available', 
                   transform=ax.transAxes, ha='center', va='center')
            ax.set_xlabel(r'$T_{\rm eff}$ (K)')
    
    fig.suptitle('WD Color Evolution', fontsize=14)
    plt.tight_layout()
    
    os.makedirs(output_dir, exist_ok=True)
    outfile = os.path.join(output_dir, 'wd_color_evolution.png')
    fig.savefig(outfile, dpi=200, bbox_inches='tight')
    print(f"Saved: {outfile}")
    plt.close(fig)


# =============================================================================
# Figure 4: SED Evolution
# =============================================================================
def plot_sed_evolution(sed_dir=SED_DIR, output_dir=OUTPUT_DIR):
    """
    Plot SEDs at different cooling ages if SED files are available.
    """
    # Find SED files
    sed_files = sorted(glob.glob(os.path.join(sed_dir, '*_SED.csv')))
    
    if not sed_files:
        print("Warning: No SED files found. Skipping SED evolution plot.")
        return
    
    # Select a subset of SEDs to plot (e.g., every 10th or specific ages)
    n_seds = len(sed_files)
    if n_seds > 10:
        # Select ~10 evenly spaced
        indices = np.linspace(0, n_seds-1, 10, dtype=int)
        selected = [sed_files[i] for i in indices]
    else:
        selected = sed_files
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Colormap for time evolution
    cmap = plt.cm.plasma
    colors = cmap(np.linspace(0.1, 0.9, len(selected)))
    
    for i, (sed_file, color) in enumerate(zip(selected, colors)):
        try:
            # Load SED (assumes wavelength, flux columns)
            sed_data = np.genfromtxt(sed_file, delimiter=',', skip_header=1)
            wavelength = sed_data[:, 0]  # Angstroms
            flux = sed_data[:, 1]        # erg/s/cm2/A
            
            # Convert to microns for plotting
            wave_um = wavelength / 1e4
            
            # Plot
            label = os.path.basename(sed_file).replace('_SED.csv', '')
            ax.plot(wave_um, flux, color=color, lw=1, alpha=0.8, label=label)
            
        except Exception as e:
            print(f"Warning: Could not load {sed_file}: {e}")
            continue
    
    ax.set_xlabel(r'Wavelength ($\mu$m)')
    ax.set_ylabel(r'Flux (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)')
    ax.set_title('WD SED Evolution')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)
    
    # Add colorbar to indicate time progression
    sm = ScalarMappable(cmap=cmap, norm=Normalize(0, 1))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Cooling Progression (early → late)')
    cbar.set_ticks([])
    
    plt.tight_layout()
    
    os.makedirs(output_dir, exist_ok=True)
    outfile = os.path.join(output_dir, 'wd_sed_evolution.png')
    fig.savefig(outfile, dpi=200, bbox_inches='tight')
    print(f"Saved: {outfile}")
    plt.close(fig)


# =============================================================================
# Figure 5: HR Diagram (bonus)
# =============================================================================
def plot_hr_diagram(data, output_dir=OUTPUT_DIR):
    """
    Standard HR diagram with cooling track.
    """
    if 'log_Teff' not in data or 'log_L' not in data:
        print("Warning: Missing log_Teff or log_L. Skipping HR diagram.")
        return
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    log_teff = data['log_Teff']
    log_L = data['log_L']
    
    # Color by age if available
    if 'star_age' in data:
        age = data['star_age'] / 1e9
        sc = ax.scatter(log_teff, log_L, c=age, cmap='plasma', s=5, alpha=0.8)
        cbar = plt.colorbar(sc, ax=ax, label='Age (Gyr)')
    else:
        ax.plot(log_teff, log_L, 'b-', lw=1)
    
    ax.set_xlabel(r'$\log T_{\rm eff}$ (K)')
    ax.set_ylabel(r'$\log L/L_\odot$')
    ax.set_title('WD Cooling Track - HR Diagram')
    ax.invert_xaxis()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    os.makedirs(output_dir, exist_ok=True)
    outfile = os.path.join(output_dir, 'wd_hr_diagram.png')
    fig.savefig(outfile, dpi=200, bbox_inches='tight')
    print(f"Saved: {outfile}")
    plt.close(fig)


# =============================================================================
# Main
# =============================================================================
def main():
    """
    Generate all analysis plots for WD cooling sequence.
    """
    print("=" * 60)
    print("White Dwarf Cooling Sequence Analysis")
    print("=" * 60)
    
    # Load data
    print(f"\nLoading history data from {HISTORY_FILE}...")
    try:
        data = load_history(HISTORY_FILE)
        print(f"  Loaded {len(data)} columns")
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        print("Run the MESA test suite first (./rn2)")
        return
    
    # Detect available filters
    available = get_available_filters(data)
    print(f"  Available filters: {available}")
    
    # Compute colors
    colors = compute_colors(data, available)
    print(f"  Computed colors: {list(colors.keys())}")
    
    # Generate plots
    print("\nGenerating plots...")
    
    print("  1. CMD (cooling track)...")
    plot_cmd(data, colors)
    
    print("  2. Multi-band light curves...")
    plot_lightcurves(data)
    
    print("  3. Color evolution vs Teff...")
    plot_color_evolution(data, colors)
    
    print("  4. SED evolution...")
    plot_sed_evolution()
    
    print("  5. HR diagram...")
    plot_hr_diagram(data)
    
    print("\n" + "=" * 60)
    print(f"Analysis complete! Figures saved to {OUTPUT_DIR}/")
    print("=" * 60)


if __name__ == '__main__':
    main()
