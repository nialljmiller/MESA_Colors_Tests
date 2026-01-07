#!/usr/bin/env python3
"""
compare_photometry.py

Compare MESA Colors output between spotted and unspotted stellar models.
Demonstrates chromatic signature of starspots: blue bands affected more than red.

Usage:
    python compare_photometry.py

Outputs:
    plots/chromatic_signature.png - Magnitude difference vs wavelength
    plots/color_comparison.png    - Color-magnitude comparison
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# LSST filter properties (central wavelengths in nm)
LSST_FILTERS = {
    'u': 367,
    'g': 483,
    'r': 622,
    'i': 755,
    'z': 869,
    'y': 971
}

def read_mesa_history(filepath):
    """Read MESA history file, handling header format."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Find header line (column names)
    header_idx = None
    for i, line in enumerate(lines):
        if line.strip().startswith('model_number') or 'log_Teff' in line:
            header_idx = i
            break
    
    if header_idx is None:
        # Try standard MESA format: line 6 is header
        header_idx = 5
    
    # Parse header
    header = lines[header_idx].split()
    
    # Read data (skip header rows)
    data = {}
    for col in header:
        data[col] = []
    
    for line in lines[header_idx + 1:]:
        if line.strip() and not line.strip().startswith('#'):
            vals = line.split()
            if len(vals) == len(header):
                for col, val in zip(header, vals):
                    try:
                        data[col].append(float(val))
                    except ValueError:
                        data[col].append(np.nan)
    
    # Convert to numpy arrays
    for col in data:
        data[col] = np.array(data[col])
    
    return data


def find_magnitude_columns(data):
    """Identify magnitude columns in the data (filter-specific)."""
    mag_cols = {}
    for col in data.keys():
        col_lower = col.lower()
        for filt in LSST_FILTERS.keys():
            # Look for patterns like 'lsst_u', 'u_mag', 'abs_mag_u', etc.
            if filt in col_lower and ('mag' in col_lower or 'lsst' in col_lower):
                mag_cols[filt] = col
                break
    return mag_cols


def get_final_photometry(log_dir, history_name):
    """Extract final timestep photometry from a MESA run."""
    history_path = Path(log_dir) / history_name
    
    if not history_path.exists():
        print(f"Warning: {history_path} not found")
        return None
    
    data = read_mesa_history(history_path)
    
    # Get the final timestep values
    result = {}
    for col, vals in data.items():
        if len(vals) > 0:
            result[col] = vals[-1]
    
    return result


def main():
    # Paths
    base_dir = Path(__file__).parent.parent
    plots_dir = base_dir / 'plots'
    plots_dir.mkdir(exist_ok=True)
    
    # Read spotted and unspotted data
    spotted = get_final_photometry(
        base_dir / 'LOGS_spotted', 
        'history_spotted.data'
    )
    unspotted = get_final_photometry(
        base_dir / 'LOGS_unspotted',
        'history_unspotted.data'
    )
    
    if spotted is None or unspotted is None:
        print("Error: Could not read history files.")
        print("Make sure you have run both inlist_spotted and inlist_unspotted")
        sys.exit(1)
    
    # Find magnitude columns
    mag_cols_spotted = find_magnitude_columns(spotted)
    mag_cols_unspotted = find_magnitude_columns(unspotted)
    
    if not mag_cols_spotted:
        print("Warning: No magnitude columns found in history files.")
        print("Available columns:", list(spotted.keys())[:20], "...")
        print("\nCreating placeholder plot with expected behavior...")
        create_expected_plot(plots_dir)
        return
    
    # Extract magnitudes
    wavelengths = []
    delta_mags = []
    filters_used = []
    
    for filt in ['u', 'g', 'r', 'i', 'z', 'y']:
        if filt in mag_cols_spotted and filt in mag_cols_unspotted:
            col_s = mag_cols_spotted[filt]
            col_u = mag_cols_unspotted[filt]
            
            mag_spotted = spotted.get(col_s)
            mag_unspotted = unspotted.get(col_u)
            
            if mag_spotted is not None and mag_unspotted is not None:
                wavelengths.append(LSST_FILTERS[filt])
                delta_mags.append(mag_spotted - mag_unspotted)
                filters_used.append(filt)
    
    if len(wavelengths) == 0:
        print("No valid magnitude data found.")
        create_expected_plot(plots_dir)
        return
    
    wavelengths = np.array(wavelengths)
    delta_mags = np.array(delta_mags)
    
    # Create chromatic signature plot
    fig, ax = plt.subplots(figsize=(8, 6))
    
    ax.scatter(wavelengths, delta_mags, s=100, c='darkred', zorder=5)
    ax.plot(wavelengths, delta_mags, 'k--', alpha=0.5)
    
    for w, dm, f in zip(wavelengths, delta_mags, filters_used):
        ax.annotate(f, (w, dm), xytext=(5, 5), textcoords='offset points',
                   fontsize=12, fontweight='bold')
    
    ax.axhline(0, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('Wavelength (nm)', fontsize=12)
    ax.set_ylabel(r'$\Delta m$ = $m_{\rm spotted} - m_{\rm unspotted}$ (mag)', fontsize=12)
    ax.set_title('Chromatic Signature of Starspots\n'
                 r'($f_{\rm spot}=0.34$, $x_{\rm spot}=0.85$)', fontsize=14)
    
    # Add annotation explaining the trend
    ax.text(0.95, 0.95, 
            'Spots affect blue bands\nmore than red bands',
            transform=ax.transAxes, ha='right', va='top',
            fontsize=10, style='italic',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(plots_dir / 'chromatic_signature.png', dpi=150)
    print(f"Saved: {plots_dir / 'chromatic_signature.png'}")
    plt.close()
    
    # Print summary
    print("\n" + "="*50)
    print("STARSPOT CHROMATIC SIGNATURE SUMMARY")
    print("="*50)
    print(f"Spot coverage: 34%")
    print(f"Temperature contrast: 0.85 (spots 15% cooler)")
    print("-"*50)
    print(f"{'Filter':<8} {'λ (nm)':<10} {'Δm (mag)':<10}")
    print("-"*50)
    for f, w, dm in zip(filters_used, wavelengths, delta_mags):
        print(f"{f:<8} {w:<10.0f} {dm:<+10.4f}")
    print("="*50)
    
    # Stellar parameters
    print("\nStellar Parameters (unspotted):")
    if 'log_Teff' in unspotted:
        print(f"  log Teff = {unspotted['log_Teff']:.4f}")
    if 'log_L' in unspotted:
        print(f"  log L/Lsun = {unspotted['log_L']:.4f}")
    if 'log_g' in unspotted:
        print(f"  log g = {unspotted['log_g']:.4f}")


def create_expected_plot(plots_dir):
    """Create a plot showing expected behavior when real data unavailable."""
    
    # Expected behavior: spots make star fainter, more so in blue
    # Using approximate values based on blackbody + spot model
    wavelengths = np.array([367, 483, 622, 755, 869, 971])
    filters = ['u', 'g', 'r', 'i', 'z', 'y']
    
    # Approximate delta_mag for fspot=0.34, xspot=0.85
    # Using Stefan-Boltzmann and Wien approximation
    T_phot = 4500  # K, approximate for 0.7 Msun
    T_spot = 0.85 * T_phot
    fspot = 0.34
    
    # Simplified chromatic model
    delta_mags = []
    for lam in wavelengths:
        # Planck function ratio at this wavelength
        h, c, k = 6.626e-34, 3e8, 1.38e-23
        lam_m = lam * 1e-9
        
        def planck(T, lam):
            return 1.0 / (np.exp(h*c/(lam*k*T)) - 1) * (1/lam**5)
        
        B_phot = planck(T_phot, lam_m)
        B_spot = planck(T_spot, lam_m)
        
        # Flux from spotted star relative to unspotted
        F_spotted = (1 - fspot) * B_phot + fspot * B_spot
        F_unspotted = B_phot
        
        # Magnitude difference
        dm = -2.5 * np.log10(F_spotted / F_unspotted)
        delta_mags.append(dm)
    
    delta_mags = np.array(delta_mags)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    ax.scatter(wavelengths, delta_mags, s=100, c='darkred', zorder=5)
    ax.plot(wavelengths, delta_mags, 'k--', alpha=0.5)
    
    for w, dm, f in zip(wavelengths, delta_mags, filters):
        ax.annotate(f, (w, dm), xytext=(5, 5), textcoords='offset points',
                   fontsize=12, fontweight='bold')
    
    ax.axhline(0, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('Wavelength (nm)', fontsize=12)
    ax.set_ylabel(r'$\Delta m$ = $m_{\rm spotted} - m_{\rm unspotted}$ (mag)', fontsize=12)
    ax.set_title('Expected Chromatic Signature of Starspots\n'
                 r'($f_{\rm spot}=0.34$, $x_{\rm spot}=0.85$, $T_{\rm phot}=4500$ K)',
                 fontsize=14)
    
    ax.text(0.95, 0.95, 
            'EXPECTED BEHAVIOR\n(actual data not yet available)',
            transform=ax.transAxes, ha='right', va='top',
            fontsize=10, style='italic', color='red',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(plots_dir / 'chromatic_signature_expected.png', dpi=150)
    print(f"Saved: {plots_dir / 'chromatic_signature_expected.png'}")
    plt.close()


if __name__ == '__main__':
    main()
