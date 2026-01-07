.. _2M_TP_AGB:

*************************************
2M_TP_AGB (Thermal Pulse AGB Colors)
*************************************

This test case demonstrates the MESA Colors module's ability to capture rapid 
photometric evolution during helium shell flashes (thermal pulses) on the 
Asymptotic Giant Branch.

Scientific Motivation
=====================

Thermal pulses represent one of the most dramatic periodic phenomena in 
stellar evolution. During each pulse cycle:

1. **Quiescent phase** (~50,000-100,000 years): H-shell burning dominates, 
   slowly adding mass to the He layer.

2. **He-shell flash** (~100-1000 years): When the He layer reaches critical 
   mass, a thermonuclear runaway produces luminosity increases of factors 
   ~2-10 at the surface.

3. **Dredge-up** (if it occurs): Convection penetrates into processed material, 
   altering surface abundances and spectral properties.

4. **Recovery**: The star returns to quiescent H-shell burning.

These phases produce measurable photometric variations even at the surface, 
despite occurring deep in the stellar interior. A 2 M☉ AGB star experiences:

- Luminosity variations: Δlog L ~ 0.3-0.5 during pulse
- Temperature variations: ΔTeff ~ 100-300 K
- Color variations: Δ(B-V) ~ 0.1-0.3 mag

High-cadence surveys (LSST, PLATO) will detect these variations in thousands 
of AGB stars, requiring theoretical predictions at matching temporal resolution.

Why This Case Works Well for Colors
===================================

1. **Temperature range**: AGB photospheres (Teff ~ 2500-4000 K) are well 
   within the Kurucz atmosphere grid coverage.

2. **Observable timescales**: Thermal pulses produce surface variations over 
   ~100 years—long enough for MESA to resolve with reasonable timesteps, 
   short enough to be scientifically interesting.

3. **Clean physics**: Unlike novae or supernovae, no shock-dominated or 
   super-Eddington phases that would invalidate static atmosphere assumptions.

4. **TDC relevance**: The rapid He-shell flash develops faster than the 
   convective turnover time, demonstrating Time-Dependent Convection effects.

Test Case Structure
===================

**Stage 1** (``inlist_to_TP_AGB``): Evolves a 2 M☉ star from pre-main-sequence 
through core hydrogen burning, the red giant branch, core helium burning, 
and the early AGB until reaching the thermally-pulsing AGB phase. Colors is 
disabled during this long evolution (~1 Gyr) to save computation time.

**Stage 2** (``inlist_TP_AGB``): Loads the model at the onset of thermal pulses 
and continues with Colors enabled at every timestep. The simulation captures 
5-10 complete thermal pulse cycles over ~500,000 years of evolution.

Colors Configuration
====================

.. code-block:: fortran

   &colors
      use_colors = .true.
      instrument = '/colors/data/filters/Generic/Johnson'
      stellar_atm = '/colors/data/stellar_models/Kurucz2003all/'
      distance = 3.0857d19  ! 10 pc (absolute magnitudes)
      mag_system = 'Vega'
   /

The Johnson BVRI system is ideal for cool AGB stars where emission peaks 
in the red/near-IR. The Kurucz grid provides coverage down to ~3500 K; 
the coolest AGB phases may involve slight extrapolation.

Key Physics Settings
====================

.. code-block:: fortran

   ! Time-Dependent Convection for He-shell flash
   MLT_option = 'TDC'
   
   ! AGB mass loss
   AGB_wind_scheme = 'Blocker'
   Blocker_scaling_factor = 0.05d0
   
   ! Tight timestep control during rapid phases
   delta_lgL_limit = 0.05d0
   delta_lgTeff_limit = 0.01d0

Expected Output
===============

The ``LOGS/history.data`` file will contain:

**Standard columns**:
   - ``star_age``: Time since start of TP-AGB phase
   - ``log_L``, ``log_Teff``, ``log_R``: Surface properties
   - ``log_LHe``: He-burning luminosity (spikes during flash)
   - ``he_core_mass``: Growing core mass

**Colors columns** (added automatically):
   - ``Mag_bol``: Bolometric magnitude
   - ``B``, ``V``, ``R``, ``I``: Johnson magnitudes
   - ``Interp_rad``: Atmosphere interpolation diagnostic

During thermal pulses, you should see:

- ``log_LHe`` spike by 2-4 orders of magnitude
- ``log_L`` increase by 0.3-0.5 dex at surface
- ``V`` brighten by 0.5-1.5 mag
- ``B-V`` redden then recover

Running the Test
================

.. code-block:: bash

   ./mk        # Compile (if needed)
   ./rn        # Run both stages

Stage 1 takes ~30-60 minutes (evolving 1 Gyr of stellar evolution).
Stage 2 takes ~10-30 minutes (500,000 years with fine timesteps).

Analysis
========

The included ``plot.py`` generates:

1. **Multi-band light curves** showing magnitude variation through pulses
2. **Color-magnitude diagram** showing the color loop during each pulse
3. **HR diagram track** colored by time
4. **Pulse identification** marking individual thermal pulse events

.. code-block:: bash

   python plot.py

Key Demonstration Points
========================

1. **Thermal pulse photometry**: Colors captures the full photometric 
   evolution during He-shell flashes that would be undersampled by 
   post-processing at coarser intervals.

2. **Color evolution**: The (B-V) color tracks temperature variations, 
   showing the characteristic reddening during luminosity maximum.

3. **TDC effects**: Time-Dependent Convection produces different surface 
   response than equilibrium MLT during rapid flash phases.

4. **Survey-ready output**: Johnson magnitudes enable direct comparison 
   with ground-based monitoring of AGB variables.

Atmosphere Grid Considerations
==============================

The Kurucz2003all grid covers Teff ≥ 3500 K. Cool AGB stars can reach 
Teff ~ 2500-3000 K during thermal pulse maxima (when radius is largest). 
For these phases, the Colors module extrapolates the atmosphere grid.

For more accurate cool-AGB photometry, consider using:

- MARCS model atmospheres (better coverage at low Teff)
- Phoenix model atmospheres (extends to ~2000 K)

These would require adding the appropriate atmosphere grids to MESA's 
colors data directory.

References
==========

- MESA V: Thermal pulses and AGB evolution (Jermyn et al. 2023)
- AGB evolution review: Herwig (2005), ARA&A
- AGB mass loss: Blocker (1995), A&A

Last-Updated: 07Jan2025 (MESA Colors demonstration)
