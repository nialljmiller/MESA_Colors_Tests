#!/usr/bin/env python3
"""
Batch runner for starspots CMD grid.
Recreates MESA VI Figure 15 parameter space as a color-magnitude diagram.

Grid: 5 masses x 4 fspot values = 20 models
"""

import os
import re
import shutil
import subprocess

# Grid parameters matching MESA VI Figure 15
MASSES = [0.3, 0.6, 0.7, 0.9, 1.1]  # Solar masses
FSPOTS = [0.2, 0.4, 0.6, 0.8]       # Spot filling factors
XSPOT = 0.85                         # Fixed temperature contrast

# Paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INLIST_TEMPLATE = os.path.join(SCRIPT_DIR, "inlist_starspots_template")
INLIST_ACTIVE = os.path.join(SCRIPT_DIR, "inlist_starspots")


def create_inlist(mass, fspot):
    """Generate inlist for specific mass and fspot combination."""
    logdir = f"LOGS_M{mass}_f{fspot}"
    
    with open(INLIST_TEMPLATE, 'r') as f:
        content = f.read()
    
    # Use regex to replace values (handles varying whitespace and formats)
    content = re.sub(
        r'(initial_mass\s*=\s*)[\d.]+',
        f'\\g<1>{mass}',
        content
    )
    content = re.sub(
        r'(fspot\s*=\s*)[\d.d+-]+',
        f'\\g<1>{fspot}d0',
        content
    )
    content = re.sub(
        r"(log_directory\s*=\s*')[^']*(')",
        f"\\g<1>{logdir}\\g<2>",
        content
    )
    
    with open(INLIST_ACTIVE, 'w') as f:
        f.write(content)
    
    return logdir


def run_mesa(recompile=False):
    """Execute MESA in the work directory."""
    os.chdir(SCRIPT_DIR)
    
    if recompile:
        result = subprocess.run(["make","clean"], capture_output=True, text=True)
        result = subprocess.run(["make"], capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  ERROR: mk failed")
            print(result.stderr)
            return False
    
    # Run
    result = subprocess.run(["./rn"], capture_output=True, text=True)
    
    # Check for success
    if "termination code:" in result.stdout:
        return True
    else:
        print(f"  ERROR: MESA run failed")
        print(result.stdout[-500:] if len(result.stdout) > 500 else result.stdout)
        return False


def main():
    """Run the full grid."""
    total = len(MASSES) * len(FSPOTS)
    completed = 0
    failed = []
    
    print(f"Starting starspots CMD grid: {total} models")
    print(f"Masses: {MASSES}")
    print(f"fspot values: {FSPOTS}")
    print(f"xspot (fixed): {XSPOT}")
    print("=" * 50)
    
    first_run = True
    for mass in MASSES:
        for fspot in FSPOTS:
            run_id = f"M={mass}, fspot={fspot}"
            print(f"\n[{completed+1}/{total}] Running {run_id}...")
            
            logdir = create_inlist(mass, fspot)
            
            # Create output directory
            os.makedirs(os.path.join(SCRIPT_DIR, logdir), exist_ok=True)
            
            # Only recompile on first run
            success = run_mesa(recompile=first_run)
            first_run = False
            
            if success:
                print(f"  ✓ Completed -> {logdir}/")
                completed += 1
            else:
                print(f"  ✗ Failed")
                failed.append(run_id)
    
    # Summary
    print("\n" + "=" * 50)
    print(f"Grid complete: {completed}/{total} successful")
    if failed:
        print(f"Failed runs: {failed}")
    print("\nNext: cd python_analysis && python plot_cmd.py")


if __name__ == "__main__":
    main()