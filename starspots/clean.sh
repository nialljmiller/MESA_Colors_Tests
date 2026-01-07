#!/bin/bash
# Clean all output files

rm -rf LOGS_* COLORS_* pgstar_out photos plots
rm -f ms_model.mod *.mod
rm -f *.png
echo "Cleaned output directories and model files"
