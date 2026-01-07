.. _starspots:

******************
starspots
******************

This test case demonstrates how MESA's starspot implementation affects 
synthetic photometry, highlighting the chromatic signature of surface 
magnetic activity.

Scientific Context
==================

Starspots on magnetically active cool stars modulate stellar brightness 
in a wavelength-dependent manner. Because spots are cooler than the 
surrounding photosphere, they emit less flux at all wavelengths, but 
the suppression is stronger in blue bands than in red. This chromatic 
signature is critical for:

- Characterizing stellar activity in exoplanet host stars
- Understanding photometric variability in M dwarfs  
- Interpreting LSST and PLATO time-series photometry

MESA implements the YREC SPOTS formalism 
(`Somers et al. 2015 <https://ui.adsabs.harvard.edu/abs/2015ApJ...807..174S>`__)
as an atmospheric boundary condition modification. Two parameters control 
the spot properties:

- **fspot**: Coverage fraction (spot filling factor)
- **xspot**: Temperature contrast, T_spot/T_photosphere

Setup
=====

This demonstration evolves a 0.7 |Msun| solar-metallicity star to main 
sequence equilibrium (5 Gyr), then compares synthetic photometry with 
and without spot coverage.

Spot parameters follow observational constraints from 
`Cao et al. 2022 <https://ui.adsabs.harvard.edu/abs/2022ApJ...924...84C>`__:

- ``fspot = 0.34`` (34% spot coverage)
- ``xspot = 0.85`` (spots 15% cooler than photosphere)

Filters: LSST ugrizy (AB magnitude system)

Running the Test
================

The workflow consists of three stages:

1. **Evolve to main sequence**::

      ./rn

   This runs ``inlist_evolve``, creating ``ms_model.mod``.

2. **Run spotted case**::

   Edit ``inlist`` to point to ``inlist_spotted``, then::

      ./rn

3. **Run unspotted case**::

   Edit ``inlist`` to point to ``inlist_unspotted``, then::

      ./rn

4. **Analyze results**::

      python python_analysis/compare_photometry.py

Alternatively, use the automated script::

   ./run_all.sh

Expected Results
================

The spotted model shows:

- **Fainter overall**: Reduced flux due to cool spot coverage
- **Chromatic signature**: Larger magnitude difference in blue (u, g) 
  than red (i, z) bands
- **Color shift**: Spotted star appears redder (higher g-r, r-i values)

The magnitude difference Δm = m_spotted - m_unspotted should increase 
toward bluer wavelengths, demonstrating how spots affect blue bands 
more strongly than red.

Output Files
============

- ``LOGS_evolve/``: History during MS evolution
- ``LOGS_spotted/``: History with spots enabled
- ``LOGS_unspotted/``: History without spots
- ``COLORS_spotted/``: Synthetic photometry with spots
- ``COLORS_unspotted/``: Synthetic photometry without spots
- ``plots/``: Analysis figures

References
==========

- Somers, G., et al. 2015, ApJ, 807, 174 (SPOTS formalism)
- Cao, L., et al. 2022, ApJ, 924, 84 (Spot parameters from observations)
- MESA VI, Jermyn et al. 2023, ApJS, 265, 15 (MESA implementation)

Last-Updated: Jan 2025
