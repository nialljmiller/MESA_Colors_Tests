#!/bin/bash
# Clean all output from starspots CMD grid

rm -rf LOGS_M*_f*/
rm -f *.mod
rm -f inlist_template
rm -rf photos*/
rm -f python_analysis/*.pdf python_analysis/*.png

echo "Cleaned all grid outputs"
