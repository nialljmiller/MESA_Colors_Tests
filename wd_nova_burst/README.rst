.. _wd_nova_burst_colors:

*****************************
wd_nova_burst (Colors Module)
*****************************

This test case demonstrates the MESA Colors module's ability to capture rapid 
photometric evolution during a classical nova outburst—a thermonuclear runaway 
on the surface of an accreting white dwarf.

Scientific Motivation
=====================

Classical novae represent one of the most dramatic examples of rapid stellar 
evolution accessible to time-domain photometry. The thermonuclear runaway 
produces luminosity changes of several orders of magnitude over timescales 
of hours to days, accompanied by significant color evolution as the expanding 
photosphere cools from ~10^5 K to ~10^4 K.

High-cadence transient surveys (LSST, Roman) will detect thousands of novae 
annually, requiring theoretical infrastructure capable of predicting multi-band 
light curves at observational cadences. Traditional post-processing approaches 
risk under-sampling the rapid rise and peak phases where timesteps may be 
sub-second during the thermonuclear runaway.

The MESA Colors module computes synthetic photometry at every MESA timestep, 
ensuring that no rapid photometric evolution is lost to temporal interpolation.

Physics Ingredients
===================

This test case employs several physics components relevant to nova evolution:

* **Time-Dependent Convection (TDC)**: The ``MLT_option = 'TDC'`` setting 
  captures the non-equilibrium convective response during rapid phases when 
  the convective turnover time exceeds the evolutionary timescale.

* **High-cadence output**: ``history_interval = 1`` records photometry at 
  every MESA timestep, capturing sub-day (often sub-hour) evolution during 
  the thermonuclear runaway.

* **Super-Eddington winds**: Mass loss driven by luminosities exceeding the 
  Eddington limit affects the photospheric evolution.

Test Case Structure
===================

This test case has 2 parts:

**Part 1** (``inlist_setup``) loads ``1.1M_lgTc_7.7.mod``, a prebuilt 1.1 M☉ 
carbon-oxygen white dwarf, and evolves for a few steps to prepare for accretion. 
Colors output is disabled during this phase as the hot WD surface lies outside 
typical atmosphere grid coverage.

**Part 2** (``inlist_wd_nova_burst``) continues the evolution with accretion 
of a hydrogen-helium mixture at 10^-9 M☉/yr. The Colors module is enabled 
with LSST filters (ugrizy) and AB magnitudes. When hydrogen burning in the 
accreted envelope triggers a thermonuclear runaway, the luminosity exceeds 
10^4 L☉ and the model captures the entire outburst photometrically.

Colors Configuration
====================

The nova burst phase uses the following Colors settings:

.. code-block:: fortran

   &colors
      use_colors = .true.
      instrument = '/colors/data/filters/LSST/LSST'
      stellar_atm = '/colors/data/stellar_models/Kurucz2003all/'
      distance = 3.0857d19  ! 10 pc (absolute magnitudes)
      mag_system = 'AB'
   /

The LSST filter set was chosen for direct comparison with the Vera C. Rubin 
Observatory's transient detection pipeline. Alternative filter configurations 
(e.g., Johnson UBVRI, Roman WFI) can be substituted by changing the 
``instrument`` path.

Expected Output
===============

The history file will contain the following Colors-generated columns:

* ``Mag_bol``: Bolometric magnitude
* ``Flux_bol``: Bolometric flux (erg/s/cm²)
* ``u``, ``g``, ``r``, ``i``, ``z``, ``y``: LSST band magnitudes (AB system)

During the nova outburst:

* The bolometric magnitude will brighten by ~8-10 magnitudes
* Color evolution (e.g., g-r) will track the photospheric temperature change
* The rise to maximum will be sampled at sub-day cadence
* The decline phase will show the characteristic color reddening

Running the Test Case
=====================

.. code-block:: bash

   ./mk        # Compile
   ./rn        # Run both parts sequentially

The run produces:

* ``LOGS/history.data``: Contains all photometric columns
* ``final.mod``: Post-outburst model
* ``pgstar_out/``: Visualization frames (if pgstar enabled)

Analysis
========

The included ``plot.py`` script generates:

1. Multi-band light curves (magnitude vs. time)
2. Color-magnitude diagram (g-r vs. r)
3. Comparison of MESA timesteps with photometric evolution rate

.. code-block:: bash

   python plot.py

Key Demonstration Points
========================

1. **Rapid evolution capture**: The Colors module samples photometry at MESA's 
   native timestep resolution, which during the thermonuclear runaway can be 
   seconds to minutes—far finer than typical survey cadences.

2. **Color evolution**: The dramatic temperature changes during nova outburst 
   produce measurable color changes that would be smoothed over by coarse 
   post-processing sampling.

3. **Survey compatibility**: Output in LSST filters with AB magnitudes enables 
   direct forward-modeling comparison with transient survey data.

4. **Sub-timestep physics**: The TDC treatment ensures that convective energy 
   transport is handled appropriately even when evolution is faster than the 
   convective turnover time.

References
==========

* MESA V: Time-Dependent Convection (Jermyn et al. 2023)
* Classical nova reviews: Starrfield et al. (2016)
* LSST transient science: LSST Science Collaboration (2009)

Last-Updated: 07Jan2025 (MESA Colors demonstration)
