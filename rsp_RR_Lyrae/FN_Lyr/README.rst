.. _rsp_colors_demo:

*****************************
RSP + Colors Demonstration
*****************************

This test case demonstrates MESA Colors computing phase-resolved synthetic 
photometry during Radial Stellar Pulsation (RSP) evolution of an RR Lyrae star.

Physical Setup
==============

The model is a fundamental-mode RR Lyrae pulsator with:

- Mass: 0.65 |Msun|
- Luminosity: 60 |Lsun|
- Effective Temperature: 6500 K
- Metallicity: Z = 0.0004
- Expected Period: ~0.71 days

Key Configuration
=================

This demonstration showcases the advantage of runtime synthetic photometry
by capturing colors at every RSP timestep:

.. code-block:: fortran

   ! ~1000 timesteps per pulsation cycle
   RSP_target_steps_per_cycle = 1000
   
   ! Record colors at EVERY timestep
   history_interval = 1
   
   ! Johnson UBVRI filter system
   instrument = '/colors/data/filters/Generic/Johnson'
   mag_system = 'Vega'

With these settings, the simulation produces approximately 1000 synthetic 
magnitude measurements per 0.71-day pulsation cycle—far exceeding the 
temporal resolution achievable through post-processing approaches.

Output
======

The ``history.data`` file contains:

- Standard RSP outputs: ``rsp_phase``, ``rsp_period_in_days``, ``rsp_DeltaMag``
- Synthetic magnitudes: ``U``, ``B``, ``V``, ``R``, ``I`` (absolute magnitudes)
- Stellar parameters: ``log_Teff``, ``log_L``, ``photosphere_r``

Generating Figures
==================

After running the model, generate the paper figures with:

.. code-block:: bash

   python plot_rsp_figures.py LOGS/history.data

This produces:

1. ``rsp_lightcurve.pdf`` - Multi-band light curves showing:
   
   - Wavelength-dependent amplitude (B > V > R > I)
   - Phase offset between bands (B maximum precedes I maximum)
   - Smooth sampling of rapid rise to maximum

2. ``rsp_colorloop.pdf`` - Color-magnitude diagram showing:
   
   - Counterclockwise hysteresis loop in (B-V, V) space
   - Different colors at same magnitude on rising vs falling branches
   - Direct observable signature of pulsation dynamics

Scientific Context
==================

This demonstration is relevant for:

- **LSST**: Multi-band time-domain photometry of millions of RR Lyrae stars
- **Roman GBTDS**: Infrared light curves of bulge pulsating variables
- **PLATO**: Asteroseismic characterization of classical pulsators

The ability to generate synthetic photometry at native RSP cadence enables
direct forward-modeling of observed light curves without temporal interpolation,
facilitating period-luminosity-color calibrations and variable star population
synthesis.

Reference
=========

See Section X.X of Miller, Joyce, Mocz et al. (2025), "MESA Colors: Synthetic 
Photometry During Stellar Evolution with MESA"

Last Updated: 2025
