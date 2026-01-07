#!/bin/bash

# run_all.sh - Automated workflow for starspots demonstration
#
# Runs all three stages:
#   1. Evolve to main sequence
#   2. Compute spotted photometry
#   3. Compute unspotted photometry
#   4. Generate comparison plots

set -e  # Exit on error

echo "=========================================="
echo "MESA Colors: Starspots Demonstration"
echo "=========================================="

# Helper function to update inlist pointer
update_inlist() {
    local target=$1
    sed -i "s/extra_star_job_inlist_name(1) = '.*'/extra_star_job_inlist_name(1) = '${target}'/" inlist
    sed -i "s/extra_eos_inlist_name(1) = '.*'/extra_eos_inlist_name(1) = '${target}'/" inlist
    sed -i "s/extra_kap_inlist_name(1) = '.*'/extra_kap_inlist_name(1) = '${target}'/" inlist
    sed -i "s/extra_controls_inlist_name(1) = '.*'/extra_controls_inlist_name(1) = '${target}'/" inlist
    sed -i "s/extra_pgstar_inlist_name(1) = '.*'/extra_pgstar_inlist_name(1) = '${target}'/" inlist
    sed -i "s/extra_colors_inlist_name(1) = '.*'/extra_colors_inlist_name(1) = '${target}'/" inlist
}

# Stage 1: Evolve to main sequence
echo ""
echo "[Stage 1/3] Evolving 0.7 Msun star to main sequence..."
echo ""
update_inlist "inlist_evolve"
#./rn

if [ ! -f "ms_model.mod" ]; then
    echo "ERROR: ms_model.mod not created. Stage 1 failed."
    exit 1
fi
echo "✓ Main sequence model saved: ms_model.mod"

# Stage 2: Spotted photometry
echo ""
echo "[Stage 2/3] Computing photometry WITH starspots..."
echo ""
update_inlist "inlist_spotted"
./rn
echo "✓ Spotted photometry complete"

# Stage 3: Unspotted photometry
echo ""
echo "[Stage 3/3] Computing photometry WITHOUT starspots..."
echo ""
update_inlist "inlist_unspotted"
./rn
echo "✓ Unspotted photometry complete"

# Stage 4: Analysis
echo ""
echo "[Stage 4/4] Generating comparison plots..."
echo ""
mkdir -p plots
python3 python_analysis/compare_photometry.py

echo ""
echo "=========================================="
echo "COMPLETE"
echo "=========================================="
echo ""
echo "Output directories:"
echo "  LOGS_evolve/     - Evolution history"
echo "  LOGS_spotted/    - Spotted run history"
echo "  LOGS_unspotted/  - Unspotted run history"
echo "  COLORS_spotted/  - Spotted photometry"
echo "  COLORS_unspotted/- Unspotted photometry"
echo "  plots/           - Analysis figures"
echo ""
