.. _5M_cepheid_blue_loop:

***********************************
MESA Colors Demo: Blue Loop
***********************************

Overview
========

This test case demonstrates MESA Colors' ability to track large-amplitude chromatic changes
during a blue loop evolutionary phase. A 5 |Msun| intermediate-mass star undergoes a ~2000 K
effective temperature excursion at nearly constant luminosity during core helium burning,
producing dramatic color evolution across optical and near-infrared passbands.

Science Case
============

Blue loops in intermediate-mass stars (3–9 |Msun|) represent one of the most dramatic
chromatic evolutionary phases outside of explosive events. During core helium burning, stars
with appropriate envelope structure execute excursions toward higher effective temperature
before returning to the red giant branch. This ~2000 K temperature swing occurs at nearly
constant luminosity, making blue loops an ideal testbed for atmosphere interpolation across
a wide Teff range.

The blue loop crosses the classical Cepheid instability strip multiple times, producing:

* Large-amplitude color changes (Δ(g−r) ~ 1 mag)
* CMD tracks that traverse the instability strip
* Multi-band magnitude evolution on ~10^5 yr timescales

Mission Relevance
-----------------

* **Rubin/LSST**: Stars crossing instability strips appear as rare transients
* **Roman**: Stellar population synthesis in nearby galaxies
* **Gaia**: CMD calibration for intermediate-age populations

Colors Module Demonstration
===========================

This test exercises:

* Atmosphere interpolation across wide Teff range (4000–8000 K)
* Multi-band synthetic photometry in LSST ugrizy
* Time-resolved color evolution during evolutionary transitions
* Dense temporal sampling to capture rapid CMD tracks

Test Structure
==============

This test has two parts:

**Part 1** (``inlist_start``): Creates a pre-main-sequence model and evolves through
hydrogen exhaustion and red giant branch ascent to the onset of core helium burning.

**Part 2** (``inlist_cepheid_blue_loop``): Continues evolution through the core helium
burning phase. The model executes a blue loop, crossing the instability strip boundaries
shown in the HR diagram.

Configuration
=============

Stellar Parameters
------------------

* Mass: 5 |Msun|
* Metallicity: Z = 0.008 (sub-solar, LMC-like)
* Helium: Y = 0.256
* Overshoot: Exponential, f = 0.012 (core), f = 0.022 (envelope)

Colors Configuration
--------------------

* Atmosphere grid: Kurucz2003all
* Filter set: LSST ugrizy
* Magnitude system: AB
* Distance: 10 pc (absolute magnitudes)

Output Products
===============

The ``history.data`` file contains:

* Standard stellar evolution quantities (L, Teff, log g, composition)
* MESA Colors photometry: u, g, r, i, z, y magnitudes
* Bolometric magnitude and correction
* Interpolation diagnostic (Interp_rad)

Analysis Scripts
================

``python_helpers/plot_blue_loop.py``
    Generates publication figures:
    
    * CMD tracks (g−r vs g) showing blue loop trajectory
    * Multi-band light curves (magnitude vs time)
    * Color evolution (g−r, r−i vs time)
    * Teff evolution overlaid with instability strip boundaries

Running the Test
================

.. code-block:: bash

   # Build and run
   ./mk
   ./rn

   # Or run with restart capability
   ./rn1

   # Generate analysis figures
   cd python_helpers
   python plot_blue_loop.py

Expected Results
================

The model should:

1. Complete blue loop evolution in ~6000 models
2. Cross Teff ~ 6000 K (instability strip) at least twice
3. Show Δ(g−r) > 0.8 mag color change
4. Produce continuous CMD tracks in LSST filters

Figures for Paper
=================

This demonstration produces figures for Section 4.X of the MESA Colors paper:

* **Figure Xa**: CMD track showing blue loop trajectory in g vs (g−r)
* **Figure Xb**: Multi-band magnitude evolution vs stellar age
* **Figure Xc**: Color evolution (g−r) vs age with instability strip crossing marked

Last Updated
============

January 2025 (MESA Colors paper development)

.. |Msun| replace:: M\ :sub:`☉`
