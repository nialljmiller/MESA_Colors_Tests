#!/usr/bin/env python3
"""
Plot color-magnitude diagram from starspots grid.
Recreates MESA VI Figure 15 as CMD using LSST g-r vs M_g.
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Style
plt.style.use('seaborn-v0_8-paper')
plt.rcParams['figure.figsize'] = (8, 10)
plt.rcParams['font.size'] = 11


def read_history(filepath):
    """Read MESA history.data file."""
    with open(filepath, 'r') as f:
        # Skip header lines
        for _ in range(5):
            f.readline()
        header = f.readline().split()
    
    data = np.genfromtxt(filepath, skip_header=6, names=header)
    return data


def main():
    # Find all LOGS directories
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(script_dir)
    
    log_dirs = sorted(glob.glob(os.path.join(parent_dir, "LOGS_M*_f*")))
    
    if not log_dirs:
        print("ERROR: No LOGS directories found.")
        print("Run the grid first: python run_grid.py")
        return
    
    print(f"Found {len(log_dirs)} model directories")
    
    # Color scheme for masses
    masses = [0.3, 0.6, 0.7, 0.9, 1.1]
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(masses)))
    mass_colors = {m: c for m, c in zip(masses, colors)}
    
    # Linestyles for fspot
    fspots = [0.2, 0.4, 0.6, 0.8]
    linestyles = ['-', '--', '-.', ':']
    fspot_styles = {f: ls for f, ls in zip(fspots, linestyles)}
    
    # Create figure
    fig, ax = plt.subplots()
    
    for logdir in log_dirs:
        # Parse mass and fspot from directory name
        dirname = os.path.basename(logdir)
        parts = dirname.replace("LOGS_M", "").replace("_f", " ").split()
        mass = float(parts[0])
        fspot = float(parts[1])
        
        history_file = os.path.join(logdir, "history.data")
        if not os.path.exists(history_file):
            print(f"  Warning: {history_file} not found, skipping")
            continue
        
        data = read_history(history_file)
        
        # Extract photometry columns
        try:
            g_mag = data['g']
            r_mag = data['r']
        except ValueError:
            print(f"  Warning: LSST magnitudes not found in {dirname}")
            continue
        
        # Color index
        g_r = g_mag - r_mag
        
        # Filter valid data (on main sequence: center_h1 > 0.01)
        if 'center_h1' in data.dtype.names:
            mask = data['center_h1'] > 0.01
            g_r = g_r[mask]
            g_mag = g_mag[mask]
        
        # Plot track
        ax.plot(g_r, g_mag, 
                color=mass_colors[mass],
                linestyle=fspot_styles[fspot],
                linewidth=1.5,
                alpha=0.8)
        
        print(f"  Plotted M={mass}, fspot={fspot}: {len(g_mag)} points")
    
    # Invert y-axis (brighter at top)
    ax.invert_yaxis()
    
    # Labels
    ax.set_xlabel(r'$g - r$ (mag)', fontsize=12)
    ax.set_ylabel(r'$M_g$ (mag)', fontsize=12)
    ax.set_title('Starspot Effects on CMD\n(MESA VI Figure 15 analog)', fontsize=13)
    
    # Legends
    # Mass legend
    mass_handles = [Line2D([0], [0], color=mass_colors[m], linewidth=2, 
                           label=f'{m} M$_\\odot$') for m in masses]
    leg1 = ax.legend(handles=mass_handles, loc='upper left', 
                     title='Mass', fontsize=9)
    ax.add_artist(leg1)
    
    # fspot legend
    fspot_handles = [Line2D([0], [0], color='gray', linestyle=fspot_styles[f],
                            linewidth=2, label=f'fspot={f}') for f in fspots]
    ax.legend(handles=fspot_handles, loc='lower right',
              title='Spot filling', fontsize=9)
    
    ax.grid(True, alpha=0.3)
    
    # Save
    outfile = os.path.join(script_dir, "starspots_cmd.pdf")
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    print(f"\nSaved: {outfile}")
    
    plt.savefig(outfile.replace('.pdf', '.png'), dpi=150, bbox_inches='tight')
    print(f"Saved: {outfile.replace('.pdf', '.png')}")


if __name__ == "__main__":
    main()
