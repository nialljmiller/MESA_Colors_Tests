#!/usr/bin/env python3
"""
plot_wd_analysis.py - Comprehensive WD Cooling Analysis with MESA Colors
=========================================================================

Publication-quality figures coupling photometric evolution with SEDs.
SEDs are matched to history via model number (1:1 correspondence).

Figures produced:
    1. CMD (g vs g-r) colored by Teff
    2. Multi-band light curves vs cooling age  
    3. Color-Teff relations (multiple colors)
    4. SED evolution with filter bandpasses overlaid
    5. SED + photometry validation panel
    6. HR diagram with cooling track
    7. Color-color diagrams
    8. Publication summary figure

Author: Niall Miller (2025)
"""

import os
import re
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize, LogNorm
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Polygon
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import warnings

# Optional imports
try:
    import mesa_reader as mr
    HAS_MESA_READER = True
except ImportError:
    HAS_MESA_READER = False

try:
    from scipy.interpolate import interp1d
    from scipy.integrate import trapezoid
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

# =============================================================================
# Configuration
# =============================================================================
HISTORY_FILE = 'LOGS/history.data'
SED_DIR = 'SED'
OUTPUT_DIR = 'plots'

# LSST filter configuration
LSST_FILTERS = ['u', 'g', 'r', 'i', 'z', 'y']
LSST_COLORS = {
    'u': '#56B4E9',  # Sky blue
    'g': '#009E73',  # Bluish green
    'r': '#F0E442',  # Yellow
    'i': '#E69F00',  # Orange
    'z': '#D55E00',  # Vermillion
    'y': '#CC79A7',  # Reddish purple
}

# Effective wavelengths (Angstroms) for LSST filters
LSST_LAMBDA_EFF = {
    'u': 3671,
    'g': 4827,
    'r': 6223,
    'i': 7546,
    'z': 8691,
    'y': 9712,
}

# Plot styling
plt.style.use('default')
mpl.rcParams.update({
    'font.size': 11,
    'font.family': 'serif',
    'axes.labelsize': 13,
    'axes.titlesize': 13,
    'legend.fontsize': 9,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'figure.figsize': (10, 8),
    'figure.dpi': 150,
    'savefig.dpi': 200,
    'savefig.bbox': 'tight',
    'axes.grid': True,
    'grid.alpha': 0.3,
})


# =============================================================================
# Data Loading
# =============================================================================
def load_history_manual(filepath: str) -> Dict[str, np.ndarray]:
    """Parse MESA history.data file manually."""
    data = {}
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Find header (line with column numbers, followed by names)
    header_idx = None
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped and stripped[0].isdigit():
            parts = stripped.split()
            if all(p.isdigit() for p in parts[:5]):
                header_idx = i
                break
    
    if header_idx is None:
        raise ValueError("Could not locate header in history file")
    
    col_names = lines[header_idx + 1].split()
    for col in col_names:
        data[col] = []
    
    for line in lines[header_idx + 2:]:
        vals = line.split()
        if len(vals) == len(col_names):
            for j, col in enumerate(col_names):
                try:
                    data[col].append(float(vals[j]))
                except ValueError:
                    data[col].append(np.nan)
    
    return {k: np.array(v) for k, v in data.items()}


def load_history(filepath: str = HISTORY_FILE) -> Dict[str, np.ndarray]:
    """Load MESA history using mesa_reader or manual parsing."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"History file not found: {filepath}")
    
    if HAS_MESA_READER:
        h = mr.MesaData(filepath)
        return {col: h.data(col) for col in h.bulk_names}
    return load_history_manual(filepath)


def get_available_filters(data: Dict) -> List[str]:
    """Detect available photometric filters."""
    return [f for f in LSST_FILTERS if f in data]


def compute_colors(data: Dict) -> Dict[str, np.ndarray]:
    """Compute color indices from magnitudes."""
    colors = {}
    pairs = [('u', 'g'), ('g', 'r'), ('r', 'i'), ('i', 'z'), ('z', 'y')]
    
    for blue, red in pairs:
        if blue in data and red in data:
            colors[f'{blue}_{red}'] = data[blue] - data[red]
    
    return colors


# =============================================================================
# SED Handling
# =============================================================================
def parse_sed_filename(filename: str) -> Tuple[str, int]:
    """
    Extract filter name and model number from SED filename.
    
    Patterns:
        {filter}_SED_{model}.csv -> filter, model
        {filter}_filter_SED_{model}.csv -> filter (integrated), model
    """
    basename = os.path.basename(filename)
    
    # Try pattern: X_SED_N.csv
    match = re.match(r'^([a-zA-Z]+)_SED_(\d+)\.csv$', basename)
    if match:
        return match.group(1), int(match.group(2))
    
    # Try pattern: X_filter_SED_N.csv  
    match = re.match(r'^([a-zA-Z]+)_filter_SED_(\d+)\.csv$', basename)
    if match:
        return f"{match.group(1)}_filter", int(match.group(2))
    
    return None, None


def discover_sed_files(sed_dir: str = SED_DIR) -> Dict[int, Dict[str, str]]:
    """
    Build mapping: model_number -> {filter: filepath}
    
    Returns dict where keys are model numbers and values are dicts
    mapping filter names to their SED file paths.
    """
    sed_map = {}
    
    if not os.path.isdir(sed_dir):
        return sed_map
    
    for fpath in glob.glob(os.path.join(sed_dir, '*.csv')):
        filt, model = parse_sed_filename(fpath)
        if model is not None:
            if model not in sed_map:
                sed_map[model] = {}
            sed_map[model][filt] = fpath
    
    return sed_map


def load_sed(filepath: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load SED CSV file.
    
    Returns (wavelength_angstrom, flux) arrays.
    Handles both comma and whitespace delimiters.
    """
    try:
        # Try comma delimiter first
        data = np.genfromtxt(filepath, delimiter=',', skip_header=1)
        if data.ndim == 1 or data.shape[1] < 2:
            # Fall back to whitespace
            data = np.genfromtxt(filepath, skip_header=1)
    except:
        data = np.genfromtxt(filepath, skip_header=1)
    
    return data[:, 0], data[:, 1]


def get_model_numbers_with_seds(sed_map: Dict, data: Dict) -> np.ndarray:
    """Get model numbers that exist in both history and SED files."""
    if 'model_number' not in data:
        return np.array([])
    
    history_models = set(data['model_number'].astype(int))
    sed_models = set(sed_map.keys())
    common = sorted(history_models & sed_models)
    
    return np.array(common)


# =============================================================================
# Figure 1: Color-Magnitude Diagram
# =============================================================================
def plot_cmd(data: Dict, colors: Dict, output_dir: str = OUTPUT_DIR):
    """CMD showing WD cooling sequence colored by Teff."""
    if 'g' not in data:
        print("Skipping CMD: g-band not found")
        return
    
    fig, ax = plt.subplots(figsize=(8, 10))
    
    g_mag = data['g']
    
    if 'g_r' in colors:
        color = colors['g_r']
        xlabel = '$g - r$ (AB mag)'
    elif 'r' in data:
        color = data['g'] - data['r']
        xlabel = '$g - r$ (AB mag)'
    else:
        print("Skipping CMD: cannot compute g-r")
        return
    
    # Color by Teff
    if 'log_Teff' in data:
        teff = 10**data['log_Teff']
        norm = LogNorm(vmin=max(3000, teff.min()), vmax=min(100000, teff.max()))
        
        sc = ax.scatter(color, g_mag, c=teff, cmap='plasma', norm=norm,
                       s=8, alpha=0.8, edgecolors='none')
        cbar = plt.colorbar(sc, ax=ax, pad=0.02)
        cbar.set_label(r'$T_{\rm eff}$ (K)')
    else:
        ax.scatter(color, g_mag, c='steelblue', s=5, alpha=0.7)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r'$M_g$ (AB mag)')
    ax.set_title('White Dwarf Cooling Sequence')
    ax.invert_yaxis()
    
    # Mark evolutionary direction
    if len(g_mag) > 20:
        ax.annotate('', xy=(color[-1], g_mag[-1]), xytext=(color[0], g_mag[0]),
                   arrowprops=dict(arrowstyle='->', color='gray', lw=1.5))
        ax.text(color[0], g_mag[0]-0.3, 'Hot', fontsize=9, ha='center', color='gray')
        ax.text(color[-1], g_mag[-1]+0.3, 'Cool', fontsize=9, ha='center', color='gray')
    
    plt.tight_layout()
    _save_figure(fig, output_dir, 'wd_cmd')


# =============================================================================
# Figure 2: Multi-band Light Curves
# =============================================================================
def plot_lightcurves(data: Dict, output_dir: str = OUTPUT_DIR):
    """Magnitude evolution in all filters vs cooling age."""
    available = get_available_filters(data)
    if not available:
        print("Skipping lightcurves: no filters found")
        return
    
    fig, axes = plt.subplots(2, 1, figsize=(12, 9), sharex=True,
                             gridspec_kw={'height_ratios': [2, 1]})
    
    # Time axis
    if 'star_age' in data:
        time = data['star_age'] / 1e9
        xlabel = 'Age (Gyr)'
    else:
        time = np.arange(len(data[available[0]]))
        xlabel = 'Model Number'
    
    # Top panel: magnitudes
    ax1 = axes[0]
    for f in available:
        ax1.plot(time, data[f], label=f, color=LSST_COLORS[f], lw=1.5, alpha=0.9)
    
    ax1.set_ylabel('Absolute Magnitude (AB)')
    ax1.invert_yaxis()
    ax1.legend(loc='upper left', ncol=len(available))
    ax1.set_title('WD Cooling: Multi-band Photometry')
    
    # Bottom panel: Teff and log_L
    ax2 = axes[1]
    if 'log_Teff' in data:
        ax2.plot(time, 10**data['log_Teff']/1000, 'k-', lw=1.5, label=r'$T_{\rm eff}$/1000 K')
    if 'log_L' in data:
        ax2b = ax2.twinx()
        ax2b.plot(time, data['log_L'], 'r--', lw=1, alpha=0.7, label=r'$\log L/L_\odot$')
        ax2b.set_ylabel(r'$\log L/L_\odot$', color='red')
        ax2b.tick_params(axis='y', labelcolor='red')
    
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel(r'$T_{\rm eff}$/1000 K')
    ax2.legend(loc='upper right')
    
    plt.tight_layout()
    _save_figure(fig, output_dir, 'wd_lightcurves')


# =============================================================================
# Figure 3: Color-Teff Relations
# =============================================================================
def plot_color_teff(data: Dict, colors: Dict, output_dir: str = OUTPUT_DIR):
    """Color indices as functions of effective temperature."""
    if 'log_Teff' not in data:
        print("Skipping color-Teff: Teff not found")
        return
    
    teff = 10**data['log_Teff']
    valid = (teff > 1000) & (teff < 200000)
    
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    axes = axes.flatten()
    
    color_configs = [
        ('u_g', '$u - g$'),
        ('g_r', '$g - r$'),
        ('r_i', '$r - i$'),
        ('i_z', '$i - z$'),
        ('z_y', '$z - y$'),
    ]
    
    # Use cooling age for color if available
    if 'star_age' in data:
        c_data = data['star_age'][valid] / 1e9
        cmap = 'viridis'
        clabel = 'Age (Gyr)'
    else:
        c_data = teff[valid]
        cmap = 'plasma'
        clabel = r'$T_{\rm eff}$ (K)'
    
    for ax, (ckey, clabel_ax) in zip(axes[:5], color_configs):
        if ckey in colors:
            c = colors[ckey][valid]
            sc = ax.scatter(teff[valid], c, c=c_data, cmap=cmap, s=5, alpha=0.7)
            ax.set_ylabel(clabel_ax)
        else:
            ax.text(0.5, 0.5, f'{clabel_ax}\nnot available',
                   transform=ax.transAxes, ha='center', va='center')
        
        ax.set_xlabel(r'$T_{\rm eff}$ (K)')
        ax.set_xscale('log')
        ax.invert_xaxis()
    
    # Use last panel for colorbar
    axes[5].axis('off')
    sm = ScalarMappable(cmap=cmap, norm=Normalize(c_data.min(), c_data.max()))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes[5], orientation='vertical', fraction=0.8)
    cbar.set_label(clabel)
    
    fig.suptitle('WD Color Evolution with Temperature', fontsize=14)
    plt.tight_layout()
    _save_figure(fig, output_dir, 'wd_color_teff')


# =============================================================================
# Figure 4: SED Evolution with Filter Bandpasses
# =============================================================================
def plot_sed_evolution(data: Dict, sed_map: Dict, output_dir: str = OUTPUT_DIR,
                       n_seds: int = 8):
    """
    Plot SED evolution at selected epochs with LSST filter bandpasses.
    Selects models spanning the cooling sequence.
    """
    common_models = get_model_numbers_with_seds(sed_map, data)
    
    if len(common_models) == 0:
        print("Skipping SED evolution: no matching models")
        return
    
    # Select evenly-spaced models
    if len(common_models) > n_seds:
        indices = np.linspace(0, len(common_models)-1, n_seds, dtype=int)
        selected_models = common_models[indices]
    else:
        selected_models = common_models
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Get Teff for coloring
    model_to_teff = {}
    if 'model_number' in data and 'log_Teff' in data:
        for m, t in zip(data['model_number'], data['log_Teff']):
            model_to_teff[int(m)] = 10**t
    
    teffs = [model_to_teff.get(m, None) for m in selected_models]
    valid_teffs = [t for t in teffs if t is not None]
    
    if valid_teffs:
        norm = LogNorm(vmin=min(valid_teffs), vmax=max(valid_teffs))
        cmap = plt.cm.plasma
    else:
        norm = Normalize(0, 1)
        cmap = plt.cm.viridis
    
    # Plot SEDs
    for i, model in enumerate(selected_models):
        # Try to load g-band SED (or any available)
        sed_files = sed_map.get(model, {})
        sed_file = sed_files.get('g') or sed_files.get(list(sed_files.keys())[0] if sed_files else None)
        
        if sed_file is None:
            continue
        
        try:
            wave, flux = load_sed(sed_file)
            wave_um = wave / 1e4  # Convert to microns
            
            teff = model_to_teff.get(model)
            if teff:
                color = cmap(norm(teff))
                label = f'Model {model} ({teff/1000:.1f} kK)'
            else:
                color = cmap(i / len(selected_models))
                label = f'Model {model}'
            
            ax.plot(wave_um, flux, color=color, lw=1, alpha=0.8, label=label)
        except Exception as e:
            print(f"Warning: Could not load SED for model {model}: {e}")
    
    # Add filter bandpasses (schematic)
    ymin, ymax = ax.get_ylim()
    band_height = (ymax - ymin) * 0.08
    
    for filt, wave_eff in LSST_LAMBDA_EFF.items():
        # Approximate filter width
        width = wave_eff * 0.15
        wave_um = wave_eff / 1e4
        width_um = width / 1e4
        
        rect = plt.Rectangle((wave_um - width_um/2, ymin), width_um, band_height,
                             facecolor=LSST_COLORS[filt], alpha=0.3, edgecolor='none')
        ax.add_patch(rect)
        ax.text(wave_um, ymin + band_height*1.2, filt, ha='center', va='bottom',
               fontsize=10, color=LSST_COLORS[filt], fontweight='bold')
    
    ax.set_xlabel(r'Wavelength ($\mu$m)')
    ax.set_ylabel(r'$F_\lambda$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)')
    ax.set_title('WD SED Evolution During Cooling')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(0.1, 2.0)
    ax.legend(loc='upper right', fontsize=8, ncol=2)
    
    # Colorbar
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label(r'$T_{\rm eff}$ (K)')
    
    plt.tight_layout()
    _save_figure(fig, output_dir, 'wd_sed_evolution')


# =============================================================================
# Figure 5: SED + Photometry Validation
# =============================================================================
def plot_sed_photometry_comparison(data: Dict, sed_map: Dict, 
                                   output_dir: str = OUTPUT_DIR,
                                   model_number: int = None):
    """
    Compare SED-integrated magnitudes with history values for validation.
    Shows a single model's SED with photometric points overlaid.
    """
    common_models = get_model_numbers_with_seds(sed_map, data)
    
    if len(common_models) == 0:
        print("Skipping SED-photometry comparison: no matching models")
        return
    
    # Select model (middle of sequence if not specified)
    if model_number is None:
        model_number = common_models[len(common_models)//2]
    elif model_number not in common_models:
        print(f"Model {model_number} not found, using {common_models[len(common_models)//2]}")
        model_number = common_models[len(common_models)//2]
    
    sed_files = sed_map.get(model_number, {})
    if not sed_files:
        print(f"No SED files for model {model_number}")
        return
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left panel: SED with photometric points
    # Load and plot SED
    sed_file = sed_files.get('g') or list(sed_files.values())[0]
    try:
        wave, flux = load_sed(sed_file)
        wave_um = wave / 1e4
        ax1.plot(wave_um, flux, 'k-', lw=1, alpha=0.7, label='SED')
    except:
        print("Could not load SED for comparison plot")
        plt.close()
        return
    
    # Get history row for this model
    if 'model_number' in data:
        model_idx = np.where(data['model_number'].astype(int) == model_number)[0]
        if len(model_idx) > 0:
            idx = model_idx[0]
            
            # Plot photometric points
            for filt in get_available_filters(data):
                mag = data[filt][idx]
                wave_eff = LSST_LAMBDA_EFF[filt] / 1e4  # microns
                
                # Convert magnitude to flux (approximate)
                # AB mag: m = -2.5 log10(f_nu) - 48.6
                # f_nu in erg/s/cm2/Hz
                f_nu = 10**(-(mag + 48.6)/2.5)
                # f_lambda = f_nu * c / lambda^2
                c = 3e10  # cm/s
                wave_cm = wave_eff * 1e-4
                f_lambda = f_nu * c / (wave_cm**2) / 1e8  # convert to per Angstrom
                
                ax1.scatter(wave_eff, f_lambda, s=100, marker='o',
                           color=LSST_COLORS[filt], edgecolors='black',
                           linewidths=1, zorder=10, label=f'{filt} = {mag:.2f}')
    
    ax1.set_xlabel(r'Wavelength ($\mu$m)')
    ax1.set_ylabel(r'$F_\lambda$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)')
    ax1.set_title(f'Model {model_number}: SED + Photometry')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(0.2, 1.5)
    ax1.legend(loc='upper right', fontsize=9)
    
    # Right panel: Filter-by-filter residuals across all models
    ax2.text(0.5, 0.5, 'Residual analysis\n(requires scipy)', 
            transform=ax2.transAxes, ha='center', va='center', fontsize=12)
    ax2.set_xlabel('Model Number')
    ax2.set_ylabel('Magnitude Residual (computed - history)')
    ax2.set_title('Photometry Validation')
    
    plt.tight_layout()
    _save_figure(fig, output_dir, 'wd_sed_photometry')


# =============================================================================
# Figure 6: HR Diagram
# =============================================================================
def plot_hr_diagram(data: Dict, output_dir: str = OUTPUT_DIR):
    """HR diagram with cooling track."""
    if 'log_Teff' not in data or 'log_L' not in data:
        print("Skipping HR diagram: missing Teff or L")
        return
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    log_teff = data['log_Teff']
    log_L = data['log_L']
    
    if 'star_age' in data:
        age_gyr = data['star_age'] / 1e9
        sc = ax.scatter(log_teff, log_L, c=age_gyr, cmap='plasma', s=8, alpha=0.8)
        cbar = plt.colorbar(sc, ax=ax)
        cbar.set_label('Age (Gyr)')
    else:
        ax.plot(log_teff, log_L, 'b-', lw=1.5)
    
    # Mark start/end
    ax.plot(log_teff[0], log_L[0], 'go', ms=12, mew=2, mfc='none', label='Start')
    ax.plot(log_teff[-1], log_L[-1], 'rs', ms=12, mew=2, mfc='none', label='End')
    
    ax.set_xlabel(r'$\log T_{\rm eff}$ (K)')
    ax.set_ylabel(r'$\log L/L_\odot$')
    ax.set_title('White Dwarf Cooling Track')
    ax.invert_xaxis()
    ax.legend(loc='lower left')
    
    plt.tight_layout()
    _save_figure(fig, output_dir, 'wd_hr_diagram')


# =============================================================================
# Figure 7: Color-Color Diagrams
# =============================================================================
def plot_color_color(data: Dict, colors: Dict, output_dir: str = OUTPUT_DIR):
    """Color-color diagrams for WD classification."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    pairs = [
        (('u_g', 'g_r'), ('$u - g$', '$g - r$')),
        (('g_r', 'r_i'), ('$g - r$', '$r - i$')),
        (('r_i', 'i_z'), ('$r - i$', '$i - z$')),
    ]
    
    if 'log_Teff' in data:
        teff = 10**data['log_Teff']
        c_data = teff
        cmap = 'plasma'
        clabel = r'$T_{\rm eff}$ (K)'
        norm = LogNorm(vmin=max(3000, teff.min()), vmax=min(100000, teff.max()))
    else:
        c_data = np.arange(len(data[list(data.keys())[0]]))
        cmap = 'viridis'
        clabel = 'Model index'
        norm = Normalize()
    
    for ax, ((cx, cy), (lx, ly)) in zip(axes, pairs):
        if cx in colors and cy in colors:
            sc = ax.scatter(colors[cx], colors[cy], c=c_data, cmap=cmap,
                           norm=norm, s=8, alpha=0.7)
            
            # Mark start/end
            ax.plot(colors[cx][0], colors[cy][0], 'go', ms=10, mew=2, mfc='none')
            ax.plot(colors[cx][-1], colors[cy][-1], 'rs', ms=10, mew=2, mfc='none')
            
            ax.set_xlabel(lx)
            ax.set_ylabel(ly)
        else:
            ax.text(0.5, 0.5, 'Colors not available',
                   transform=ax.transAxes, ha='center')
    
    # Add colorbar
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, orientation='vertical', fraction=0.02, pad=0.02)
    cbar.set_label(clabel)
    
    fig.suptitle('WD Color-Color Diagrams', fontsize=14)
    plt.tight_layout()
    _save_figure(fig, output_dir, 'wd_color_color')


# =============================================================================
# Figure 8: Publication Summary
# =============================================================================
def plot_summary(data: Dict, colors: Dict, sed_map: Dict, 
                 output_dir: str = OUTPUT_DIR):
    """Publication-ready 4-panel summary figure."""
    fig = plt.figure(figsize=(12, 11))
    gs = GridSpec(2, 2, hspace=0.28, wspace=0.28)
    
    # Panel A: CMD
    ax1 = fig.add_subplot(gs[0, 0])
    if 'g' in data and 'g_r' in colors:
        g_mag = data['g']
        color = colors['g_r']
        
        if 'log_Teff' in data:
            teff = 10**data['log_Teff']
            sc = ax1.scatter(color, g_mag, c=teff, cmap='plasma',
                            norm=LogNorm(vmin=3000, vmax=100000),
                            s=5, alpha=0.8)
            cbar = plt.colorbar(sc, ax=ax1, pad=0.02)
            cbar.set_label(r'$T_{\rm eff}$ (K)', fontsize=10)
        else:
            ax1.scatter(color, g_mag, c='steelblue', s=5, alpha=0.7)
        
        ax1.invert_yaxis()
        ax1.set_xlabel('$g - r$')
        ax1.set_ylabel('$M_g$')
    ax1.set_title('(a) Color-Magnitude Diagram')
    
    # Panel B: Multi-band lightcurves
    ax2 = fig.add_subplot(gs[0, 1])
    available = get_available_filters(data)
    
    if available and 'star_age' in data:
        time = data['star_age'] / 1e9
        for f in available:
            ax2.plot(time, data[f], color=LSST_COLORS[f], lw=1.2, label=f)
        ax2.invert_yaxis()
        ax2.set_xlabel('Age (Gyr)')
        ax2.set_ylabel('Magnitude (AB)')
        ax2.legend(loc='upper left', ncol=2, fontsize=8)
    ax2.set_title('(b) Cooling Light Curves')
    
    # Panel C: HR diagram
    ax3 = fig.add_subplot(gs[1, 0])
    if 'log_Teff' in data and 'log_L' in data:
        log_teff = data['log_Teff']
        log_L = data['log_L']
        
        if 'star_age' in data:
            sc = ax3.scatter(log_teff, log_L, c=data['star_age']/1e9,
                            cmap='viridis', s=5, alpha=0.8)
            cbar = plt.colorbar(sc, ax=ax3, pad=0.02)
            cbar.set_label('Age (Gyr)', fontsize=10)
        else:
            ax3.plot(log_teff, log_L, 'b-', lw=1)
        
        ax3.invert_xaxis()
        ax3.set_xlabel(r'$\log T_{\rm eff}$')
        ax3.set_ylabel(r'$\log L/L_\odot$')
    ax3.set_title('(c) HR Diagram')
    
    # Panel D: Color-Teff for g-r
    ax4 = fig.add_subplot(gs[1, 1])
    if 'g_r' in colors and 'log_Teff' in data:
        teff = 10**data['log_Teff']
        valid = (teff > 1000) & (teff < 200000)
        
        if 'star_age' in data:
            sc = ax4.scatter(teff[valid], colors['g_r'][valid],
                            c=data['star_age'][valid]/1e9,
                            cmap='viridis', s=5, alpha=0.8)
        else:
            ax4.scatter(teff[valid], colors['g_r'][valid], c='steelblue', s=5)
        
        ax4.set_xscale('log')
        ax4.invert_xaxis()
        ax4.set_xlabel(r'$T_{\rm eff}$ (K)')
        ax4.set_ylabel('$g - r$')
    ax4.set_title('(d) Color Evolution')
    
    fig.suptitle('MESA Colors: White Dwarf Cooling Sequence', fontsize=14, y=0.98)
    
    plt.tight_layout()
    _save_figure(fig, output_dir, 'wd_summary')


# =============================================================================
# Utilities
# =============================================================================
def _save_figure(fig, output_dir: str, name: str):
    """Save figure in PNG and PDF formats."""
    os.makedirs(output_dir, exist_ok=True)
    
    png_path = os.path.join(output_dir, f'{name}.png')
    pdf_path = os.path.join(output_dir, f'{name}.pdf')
    
    fig.savefig(png_path, dpi=200, bbox_inches='tight')
    fig.savefig(pdf_path, bbox_inches='tight')
    print(f"Saved: {png_path}")
    
    plt.close(fig)


def print_summary(data: Dict, sed_map: Dict, colors: Dict):
    """Print data summary."""
    print("\n" + "="*60)
    print("Data Summary")
    print("="*60)
    
    print(f"\nHistory columns: {len(data)}")
    
    available = get_available_filters(data)
    print(f"Photometric filters: {available}")
    print(f"Colors computed: {list(colors.keys())}")
    
    if 'model_number' in data:
        models = data['model_number']
        print(f"Model range: {int(models.min())} - {int(models.max())}")
    
    if 'star_age' in data:
        age = data['star_age']
        print(f"Age range: {age.min()/1e9:.3f} - {age.max()/1e9:.3f} Gyr")
    
    if 'log_Teff' in data:
        teff = 10**data['log_Teff']
        print(f"Teff range: {teff.min():.0f} - {teff.max():.0f} K")
    
    print(f"\nSED files: {len(sed_map)} models with SEDs")
    common = get_model_numbers_with_seds(sed_map, data)
    print(f"Models with both history + SED: {len(common)}")


# =============================================================================
# Main
# =============================================================================
def main():
    """Generate all analysis plots."""
    print("="*60)
    print("MESA Colors: White Dwarf Cooling Analysis")
    print("="*60)
    
    # Load history
    print(f"\nLoading history from {HISTORY_FILE}...")
    try:
        data = load_history(HISTORY_FILE)
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        return
    
    # Discover SEDs
    print(f"Scanning SED directory: {SED_DIR}...")
    sed_map = discover_sed_files(SED_DIR)
    
    # Compute colors
    available = get_available_filters(data)
    colors = compute_colors(data)
    
    # Print summary
    print_summary(data, sed_map, colors)
    
    # Generate figures
    print("\n" + "="*60)
    print("Generating Figures")
    print("="*60)
    
    print("\n1. Color-Magnitude Diagram...")
    plot_cmd(data, colors)
    
    print("\n2. Multi-band Light Curves...")
    plot_lightcurves(data)
    
    print("\n3. Color-Teff Relations...")
    plot_color_teff(data, colors)
    
    print("\n4. SED Evolution...")
    plot_sed_evolution(data, sed_map)
    
    print("\n5. SED-Photometry Comparison...")
    plot_sed_photometry_comparison(data, sed_map)
    
    print("\n6. HR Diagram...")
    plot_hr_diagram(data)
    
    print("\n7. Color-Color Diagrams...")
    plot_color_color(data, colors)
    
    print("\n8. Publication Summary...")
    plot_summary(data, colors, sed_map)
    
    print("\n" + "="*60)
    print(f"Analysis complete! Figures saved to {OUTPUT_DIR}/")
    print("="*60)


if __name__ == '__main__':
    main()
