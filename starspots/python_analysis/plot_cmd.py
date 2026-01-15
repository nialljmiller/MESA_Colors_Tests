#!/usr/bin/env python3
"""
Plot color-magnitude diagram AND HR diagram from starspots grid.
CMD: LSST g-r vs M_g
HR: log L vs log Teff
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Style
plt.style.use('seaborn-v0_8-paper')
plt.rcParams['font.size'] = 11


def read_history(filepath):
    """Read MESA history.data file."""
    with open(filepath, 'r') as f:
        for _ in range(5):
            f.readline()
        header = f.readline().split()
    data = np.genfromtxt(filepath, skip_header=6, names=header)
    return data


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(script_dir)

    log_dirs = sorted(glob.glob(os.path.join(parent_dir, "LOGS_M*_f*")))
    if not log_dirs:
        print("ERROR: No LOGS directories found.")
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

    # Figure: CMD + HR
    fig, (ax_cmd, ax_hr) = plt.subplots(1, 2, figsize=(14, 6))

    for logdir in log_dirs:
        dirname = os.path.basename(logdir)
        parts = dirname.replace("LOGS_M", "").replace("_f", " ").split()
        mass = float(parts[0])
        fspot = float(parts[1])

        history_file = os.path.join(logdir, "history.data")
        if not os.path.exists(history_file):
            continue

        data = read_history(history_file)

        try:
            g_mag = data['g']
            r_mag = data['r']
            log_teff = data['log_Teff']
            log_l = data['log_L']
        except ValueError:
            continue

        g_r = g_mag - r_mag

        if 'center_h1' in data.dtype.names:
            mask = data['center_h1'] > 0.01
            g_r = g_r[mask]
            g_mag = g_mag[mask]
            log_teff = 10 ** log_teff[mask]
            log_l = log_l[mask]

        # CMD
        ax_cmd.plot(
            g_r, g_mag,
            color=mass_colors[mass],
            linestyle=fspot_styles[fspot],
            linewidth=1.5,
            alpha=0.8
        )

        # HR diagram
        ax_hr.plot(
            log_teff, log_l,
            color=mass_colors[mass],
            linestyle=fspot_styles[fspot],
            linewidth=1.5,
            alpha=0.8
        )

        print(f"Plotted M={mass}, fspot={fspot}")

    # CMD formatting
    ax_cmd.invert_yaxis()
    ax_cmd.set_xlabel(r'$g - r$ (mag)')
    ax_cmd.set_ylabel(r'$M_g$ (mag)')
    ax_cmd.set_title('CMD')

    # HR formatting
    ax_hr.invert_xaxis()
    ax_hr.set_xlabel(r'$\log T_{\mathrm{eff}}$')
    ax_hr.set_ylabel(r'$\log L / L_\odot$')
    ax_hr.set_title('HR Diagram')

    # Legends (CMD only)
    mass_handles = [
        Line2D([0], [0], color=mass_colors[m], linewidth=2,
               label=f'{m} M$_\\odot$') for m in masses
    ]
    leg1 = ax_cmd.legend(handles=mass_handles, loc='upper left',
                         title='Mass', fontsize=9)
    ax_cmd.add_artist(leg1)

    fspot_handles = [
        Line2D([0], [0], color='gray', linestyle=fspot_styles[f],
               linewidth=2, label=f'fspot={f}') for f in fspots
    ]
    ax_cmd.legend(handles=fspot_handles, loc='lower right',
                  title='Spot filling', fontsize=9)

    ax_cmd.grid(True, alpha=0.3)
    ax_hr.grid(True, alpha=0.3)

    outfile = os.path.join(script_dir, "starspots_cmd_hr.pdf")
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.savefig(outfile.replace('.pdf', '.png'), dpi=150, bbox_inches='tight')
    print(f"Saved: {outfile}")


if __name__ == "__main__":
    main()
