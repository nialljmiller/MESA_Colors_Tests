Starspots CMD Test Suite
=========================

Recreates MESA VI Figure 15 as a color-magnitude diagram using MESA Colors.

Grid Parameters (matching Figure 15)
------------------------------------

==================  ============================
Parameter           Values
==================  ============================
Masses              0.3, 0.6, 0.7, 0.9, 1.1 M☉
fspot               0.2, 0.4, 0.6, 0.8
xspot               0.85 (fixed)
Metallicity         Z = 0.014
Total models        20
==================  ============================

Photometric output: LSST ugrizy → CMD axes: g-r vs M_g

File Structure
--------------

::

    starspots_cmd/
    ├── inlist                  # Main inlist (template with <<PLACEHOLDERS>>)
    ├── history_columns.list    # LSST magnitude columns
    ├── src/
    │   └── run_star_extras.f90 # Starspot gradr_factor implementation
    ├── run_grid.py             # Batch runner (modifies inlist, runs 20 models)
    ├── clean.sh                # Remove all outputs
    └── python_analysis/
        └── plot_cmd.py         # Generate CMD figure

Running
-------

.. code-block:: bash

    # 1. Compile MESA
    ./mk
    
    # 2. Run grid (20 models)
    python run_grid.py
    
    # 3. Plot CMD
    cd python_analysis
    python plot_cmd.py

Output: ``LOGS_M{mass}_f{fspot}/history.data`` for each model

Physics
-------

Starspot implementation via ``other_gradr_factor``:

- ``x_ctrl(1)`` = fspot (filling factor)
- ``x_ctrl(2)`` = xspot (temperature contrast)
- gradr_factor = 1 / [fspot × xspot⁴ + (1 - fspot)]

References: Jermyn et al. 2023 (MESA VI §7.1), Somers et al. 2020

