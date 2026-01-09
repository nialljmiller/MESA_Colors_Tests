#!/usr/bin/env python3
"""
plot_utils.py

Shared utilities for MESA Colors ApJ paper figures.
Provides consistent styling, SED plotting, and history file handling.

Author: Miller, Joyce, Mocz et al.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import csv

# =============================================================================
# PUBLICATION STYLE SETTINGS
# =============================================================================

# ApJ single-column width: 3.5 inches, double-column: 7 inches
APJ_SINGLE_COL = 3.5
APJ_DOUBLE_COL = 7.0

def setup_apj_style():
    """Configure matplotlib for ApJ publication quality."""
    plt.rcParams.update({
        'font.size': 9,
        'font.family': 'serif',
        'font.serif': ['Times', 'Times New Roman', 'DejaVu Serif'],
        'axes.labelsize': 10,
        'axes.titlesize': 10,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'legend.fontsize': 8,
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'text.usetex': False,
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'xtick.minor.width': 0.5,
        'ytick.minor.width': 0.5,
        'lines.linewidth': 1.0,
    })


# =============================================================================
# FILTER DEFINITIONS
# =============================================================================

# LSST filter properties
LSST_FILTERS = {
    'u': {'wavelength': 3670, 'color': '#56106e'},
    'g': {'wavelength': 4830, 'color': '#3b528b'},
    'r': {'wavelength': 6220, 'color': '#21918c'},
    'i': {'wavelength': 7550, 'color': '#5ec962'},
    'z': {'wavelength': 8690, 'color': '#addc30'},
    'y': {'wavelength': 9710, 'color': '#fde725'},
}

# Johnson-Cousins filter properties
JOHNSON_FILTERS = {
    'U': {'wavelength': 3650, 'color': '#4B0082'},
    'B': {'wavelength': 4400, 'color': '#0000FF'},
    'V': {'wavelength': 5500, 'color': '#228B22'},
    'R': {'wavelength': 6400, 'color': '#FF4500'},
    'I': {'wavelength': 8000, 'color': '#8B0000'},
}

# Gaia filter properties
GAIA_FILTERS = {
    'G': {'wavelength': 6730, 'color': '#228B22'},
    'Gbp': {'wavelength': 5320, 'color': '#0000FF'},
    'Grp': {'wavelength': 7970, 'color': '#FF0000'},
}

# EM spectrum regions for SED plots
EM_REGIONS = [
    (100, 2000, 'FUV', '#4B0082'),
    (2000, 4000, 'NUV', '#9400D3'),
    (4000, 4500, '', '#0000FF'),
    (4500, 5500, '', '#00FF00'),
    (5500, 7000, '', '#FFFF00'),
    (7000, 10000, 'NIR', '#FF0000'),
    (10000, 25000, '', '#8B0000'),
]


# =============================================================================
# MESA HISTORY FILE HANDLING
# =============================================================================

def read_mesa_history(filepath):
    """
    Read MESA history.data file into a dictionary.
    
    Parameters
    ----------
    filepath : str or Path
        Path to history.data file
        
    Returns
    -------
    dict
        Column names as keys, numpy arrays as values
    """
    filepath = Path(filepath)
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Find header line (contains column names)
    header_idx = None
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith('model_number') or 'model_number' in stripped.split():
            header_idx = i
            break
    
    if header_idx is None:
        # Standard MESA format: line 6 has headers (0-indexed: 5)
        header_idx = 5
    
    # Parse column names
    cols = lines[header_idx].split()
    
    # Read data
    data = {col: [] for col in cols}
    
    for line in lines[header_idx + 1:]:
        if line.strip() and not line.strip().startswith('#'):
            vals = line.split()
            if len(vals) == len(cols):
                for col, val in zip(cols, vals):
                    try:
                        data[col].append(float(val))
                    except ValueError:
                        data[col].append(np.nan)
    
    # Convert to numpy arrays
    for col in data:
        data[col] = np.array(data[col])
    
    return data


def find_magnitude_columns(data, filter_set='lsst'):
    """
    Find magnitude column names in history data.
    
    Parameters
    ----------
    data : dict
        MESA history data
    filter_set : str
        'lsst', 'johnson', or 'gaia'
        
    Returns
    -------
    dict
        Filter name -> column name mapping
    """
    colnames = list(data.keys())
    
    if filter_set == 'lsst':
        filters = list(LSST_FILTERS.keys())
    elif filter_set == 'johnson':
        filters = list(JOHNSON_FILTERS.keys())
    elif filter_set == 'gaia':
        filters = list(GAIA_FILTERS.keys())
    else:
        filters = []
    
    found = {}
    for filt in filters:
        # Try exact match first
        if filt in colnames:
            found[filt] = filt
            continue
        # Try common patterns
        for col in colnames:
            col_lower = col.lower()
            filt_lower = filt.lower()
            if filt_lower == col_lower or f'mag_{filt_lower}' in col_lower:
                found[filt] = col
                break
    
    return found


# =============================================================================
# SED FILE HANDLING
# =============================================================================

def read_sed_file(filepath):
    """
    Read SED CSV file from MESA Colors output.
    
    Parameters
    ----------
    filepath : str or Path
        Path to SED CSV file
        
    Returns
    -------
    dict
        wavelengths, fluxes, convolved_flux arrays
    """
    filepath = Path(filepath)
    
    try:
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f)
            
            data = {'wavelengths': [], 'fluxes': [], 'convolved_flux': []}
            
            for row in reader:
                try:
                    wl = float(row['wavelengths'])
                    fl = float(row['fluxes'])
                    # Only keep rows where both wavelength and flux are valid and non-zero
                    if wl > 0 and fl > 0:
                        data['wavelengths'].append(wl)
                        data['fluxes'].append(fl)
                        # Convolved flux is optional
                        try:
                            cf = float(row.get('convolved_flux', 0))
                            data['convolved_flux'].append(cf)
                        except (ValueError, TypeError):
                            data['convolved_flux'].append(0)
                except (ValueError, KeyError):
                    continue
            
            # Convert to arrays
            result = {}
            for field in data:
                if data[field]:
                    result[field] = np.array(data[field])
            
            return result
    except Exception as e:
        print(f"Warning: Could not read {filepath}: {e}")
        return {}


def find_sed_files(sed_dir, model_number, filter_name=None):
    """
    Find SED files for a given model number.
    
    Parameters
    ----------
    sed_dir : str or Path
        Directory containing SED files
    model_number : int
        Model number to find
    filter_name : str, optional
        Specific filter to find, or None for all
        
    Returns
    -------
    list of Path
        Matching SED file paths
    """
    sed_dir = Path(sed_dir)
    pattern = f'*_SED_{model_number}.csv'
    
    if filter_name:
        pattern = f'*{filter_name}*_SED_{model_number}.csv'
    
    return list(sed_dir.glob(pattern))


def get_sed_for_model(sed_dir, model_number):
    """
    Load all SED data for a given model number.
    
    Parameters
    ----------
    sed_dir : str or Path
        Directory containing SED files
    model_number : int
        Model number to load
        
    Returns
    -------
    dict
        Filter name -> SED data dict mapping
    """
    sed_dir = Path(sed_dir)
    files = find_sed_files(sed_dir, model_number)
    
    result = {}
    for f in files:
        # Extract filter name from filename (e.g., "lsst_g_SED_123.csv" -> "g")
        name = f.stem.replace(f'_SED_{model_number}', '')
        # Clean up filter name
        for prefix in ['lsst_', 'johnson_', 'gaia_', 'LSST_']:
            name = name.replace(prefix, '')
        result[name] = read_sed_file(f)
    
    return result


# =============================================================================
# SED PLOTTING
# =============================================================================

def add_em_spectrum_regions(ax, alpha=0.05):
    """
    Add faint colored background regions for EM spectrum.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to add regions to
    alpha : float
        Transparency of regions
    """
    xlim = ax.get_xlim()
    
    for min_wave, max_wave, label, color in EM_REGIONS:
        if max_wave < xlim[0] or min_wave > xlim[1]:
            continue
        
        visible_min = max(min_wave, xlim[0])
        visible_max = min(max_wave, xlim[1])
        
        ax.axvspan(visible_min, visible_max, alpha=alpha, color=color, zorder=0)


def plot_sed(ax, wavelengths, fluxes, convolved_data=None, 
             label='SED', color='black', show_em_regions=True,
             filter_set='lsst'):
    """
    Plot an SED with optional convolved filter curves.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to plot on
    wavelengths : array
        Wavelength array in Angstroms
    fluxes : array
        Flux array
    convolved_data : dict, optional
        Filter name -> convolved flux array mapping
    label : str
        Label for the SED line
    color : str
        Color for the SED line
    show_em_regions : bool
        Whether to show EM spectrum background
    filter_set : str
        'lsst', 'johnson', or 'gaia' for filter colors
    """
    # Plot main SED
    ax.plot(wavelengths, fluxes, color=color, lw=1.2, label=label, zorder=2)
    
    # Get filter definitions
    if filter_set == 'lsst':
        filters = LSST_FILTERS
    elif filter_set == 'johnson':
        filters = JOHNSON_FILTERS
    elif filter_set == 'gaia':
        filters = GAIA_FILTERS
    else:
        filters = {}
    
    # Plot convolved curves if provided
    if convolved_data:
        for filt_name, conv_flux in convolved_data.items():
            if filt_name in filters:
                filt_color = filters[filt_name]['color']
            else:
                filt_color = 'gray'
            
            # Use same wavelength grid as SED
            if len(conv_flux) == len(wavelengths):
                ax.fill_between(wavelengths, 0, conv_flux, 
                               alpha=0.3, color=filt_color, 
                               label=filt_name, zorder=1)
    
    # Formatting
    ax.set_xlabel(r'Wavelength ($\mathrm{\AA}$)')
    ax.set_ylabel(r'$F_\lambda$ (erg s$^{-1}$ cm$^{-2}$ $\mathrm{\AA}^{-1}$)')
    ax.set_xscale('log')
    
    # Add EM regions after setting scale
    if show_em_regions:
        add_em_spectrum_regions(ax)


def create_sed_inset(ax_main, sed_data, position, size, 
                     title=None, filter_set='lsst'):
    """
    Create an inset SED plot on a main axis.
    
    Parameters
    ----------
    ax_main : matplotlib.axes.Axes
        Main axes to add inset to
    sed_data : dict
        SED data with 'wavelengths', 'fluxes' keys
    position : tuple
        (x, y, width, height) in axes coordinates
    size : tuple
        (width, height) of inset
    title : str, optional
        Title for inset
    filter_set : str
        Filter set for colors
        
    Returns
    -------
    matplotlib.axes.Axes
        The inset axes
    """
    ax_inset = ax_main.inset_axes(position)
    
    if 'wavelengths' in sed_data and 'fluxes' in sed_data:
        wavelengths = sed_data['wavelengths']
        fluxes = sed_data['fluxes']
        
        # Get filter info
        if filter_set == 'lsst':
            filters = LSST_FILTERS
        elif filter_set == 'johnson':
            filters = JOHNSON_FILTERS
        else:
            filters = GAIA_FILTERS
        
        # Plot SED
        ax_inset.plot(wavelengths, fluxes, 'k-', lw=0.8)
        
        # Plot convolved if available
        if 'convolved_flux' in sed_data:
            ax_inset.fill_between(wavelengths, 0, sed_data['convolved_flux'],
                                  alpha=0.3, color='steelblue')
        
        ax_inset.set_xscale('log')
        ax_inset.set_xlim(2000, 12000)
        
        # Minimal tick labels for inset
        ax_inset.tick_params(labelsize=6)
        ax_inset.set_xlabel('')
        ax_inset.set_ylabel('')
        
        if title:
            ax_inset.set_title(title, fontsize=7, pad=2)
    
    return ax_inset


# =============================================================================
# POINT OF INTEREST SELECTION
# =============================================================================

def find_extrema(data, column, n_points=1):
    """
    Find indices of local extrema in a data column.
    
    Parameters
    ----------
    data : dict
        MESA history data
    column : str
        Column name to search
    n_points : int
        Number of extrema to find
        
    Returns
    -------
    dict
        'maxima' and 'minima' index arrays
    """
    values = data[column]
    
    # Find local maxima
    maxima = []
    minima = []
    
    for i in range(1, len(values) - 1):
        if values[i] > values[i-1] and values[i] > values[i+1]:
            maxima.append(i)
        elif values[i] < values[i-1] and values[i] < values[i+1]:
            minima.append(i)
    
    return {'maxima': np.array(maxima), 'minima': np.array(minima)}


def find_poi_blue_loop(data):
    """
    Find points of interest for blue loop evolution.
    
    Returns indices for:
    - RGB tip (minimum Teff before loop)
    - Blue loop maximum (maximum Teff during loop)
    - AGB start (returning to cool temperatures)
    
    Parameters
    ----------
    data : dict
        MESA history data
        
    Returns
    -------
    dict
        Named indices for points of interest
    """
    teff = 10**data['log_Teff']
    
    # Find RGB tip (minimum Teff in first half)
    half = len(teff) // 2
    rgb_tip_idx = np.argmin(teff[:half])
    
    # Find blue loop maximum (maximum Teff after RGB tip)
    blue_max_idx = rgb_tip_idx + np.argmax(teff[rgb_tip_idx:])
    
    # Find a point during the loop (between RGB tip and max)
    loop_mid_idx = (rgb_tip_idx + blue_max_idx) // 2
    
    # Find end point
    end_idx = len(teff) - 1
    
    return {
        'rgb_tip': rgb_tip_idx,
        'loop_ascending': loop_mid_idx,
        'loop_maximum': blue_max_idx,
        'end': end_idx
    }


def find_poi_rsp(data):
    """
    Find points of interest for RSP pulsation.
    
    Returns indices for:
    - Maximum light (minimum magnitude)
    - Minimum light (maximum magnitude)
    - Rising branch
    - Falling branch
    
    Parameters
    ----------
    data : dict
        MESA history data
        
    Returns
    -------
    dict
        Named indices for points of interest
    """
    # Use V or g band magnitude
    mag_col = None
    for col in ['V', 'g', 'Johnson_V', 'lsst_g']:
        if col in data:
            mag_col = col
            break
    
    if mag_col is None:
        return {}
    
    mag = data[mag_col]
    
    # Get phase if available
    if 'rsp_phase' in data:
        phase = data['rsp_phase']
        # Find one complete cycle after settling
        if 'rsp_num_periods' in data:
            mask = (data['rsp_num_periods'] >= 4) & (data['rsp_num_periods'] < 5)
        else:
            mask = np.ones(len(mag), dtype=bool)
    else:
        mask = np.ones(len(mag), dtype=bool)
        phase = np.linspace(0, 1, len(mag))
    
    indices = np.where(mask)[0]
    if len(indices) < 10:
        indices = np.arange(len(mag))
    
    mag_subset = mag[indices]
    
    # Find extrema in the selected range
    max_light_local = np.argmin(mag_subset)  # Brightest = minimum mag
    min_light_local = np.argmax(mag_subset)  # Faintest = maximum mag
    
    max_light_idx = indices[max_light_local]
    min_light_idx = indices[min_light_local]
    
    # Find rising and falling points (quarter phases)
    rising_idx = indices[len(indices) // 4]
    falling_idx = indices[3 * len(indices) // 4]
    
    return {
        'maximum_light': max_light_idx,
        'minimum_light': min_light_idx,
        'rising': rising_idx,
        'falling': falling_idx
    }


def find_poi_starspot(data_spotted, data_unspotted):
    """
    Find points of interest for starspot comparison.
    
    Parameters
    ----------
    data_spotted : dict
        Spotted star history data
    data_unspotted : dict
        Unspotted star history data
        
    Returns
    -------
    dict
        Named indices (typically just final timestep)
    """
    return {
        'spotted_final': len(data_spotted['model_number']) - 1,
        'unspotted_final': len(data_unspotted['model_number']) - 1
    }


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def get_model_number_at_index(data, idx):
    """Get model number at a given array index."""
    return int(data['model_number'][idx])


def get_sed_counter_at_index(data, idx):
    """
    Get SED file counter for a given array index.
    
    SED files are numbered sequentially starting at 1,
    corresponding to history row indices 0, 1, 2, ...
    """
    return idx + 1


def annotate_point(ax, x, y, label, offset=(5, 5), fontsize=8, **kwargs):
    """Add annotation to a plot point."""
    ax.annotate(label, (x, y), 
                xytext=offset, textcoords='offset points',
                fontsize=fontsize, **kwargs)


def add_colorbar(fig, mappable, ax, label='', **kwargs):
    """Add a colorbar with consistent styling."""
    cbar = fig.colorbar(mappable, ax=ax, **kwargs)
    cbar.set_label(label, fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    return cbar
