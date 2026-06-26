.. _wd_cooling_colors:

*****************************
WD Cooling Sequence with Colors
*****************************

This test suite demonstrates the MESA ``colors`` module applied to white dwarf cooling sequences, showcasing synthetic photometry in the high-gravity, low-temperature regime relevant to modern surveys.

Scientific Motivation
=====================

White dwarfs represent the final evolutionary state of low- and intermediate-mass stars. Their cooling sequences provide powerful chronometers for stellar populations and galactic archaeology. However, connecting MESA's internal structure models to observables requires:

1. **Appropriate atmosphere models**: Standard Kurucz grids do not extend to WD surface gravities (log g ~ 7–9). This test uses Koester DA atmosphere models specifically computed for hydrogen-atmosphere white dwarfs.

2. **Exotic equation of state**: WD interiors are dominated by electron degeneracy. This test employs the Skye EOS, which properly handles the degenerate regime.

3. **Survey-relevant filters**: We compute synthetic photometry in Rubin/LSST filters (ugrizy) for direct comparison with upcoming wide-field surveys.

Physical Setup
==============

**Progenitor**: 1.0 M☉ star evolved through the AGB and planetary nebula phases

**White Dwarf Properties**:
   - Final mass: ~0.54 M☉ (typical for 1 M☉ progenitor)
   - Composition: C/O core with thin He/H envelope
   - log g: ~8.0

**Cooling Range**: From hot post-AGB (~100,000 K) down to ~5,000 K

**Atmosphere Models**: Koester DA (pure hydrogen) grid
   - Teff: 5,000 – 80,000 K
   - log g: 6.5 – 9.5
   - Reference: Koester 2010, Mem.S.A.It. 81, 921

**Filters**: Rubin/LSST ugrizy system

Two-Stage Execution
===================

This test suite uses a two-stage approach:

**Stage 1** (``inlist_to_wd``): Evolves a 1 M☉ star from pre-main sequence through the thermal pulse AGB, planetary nebula ejection, and onto the white dwarf cooling track. Saves a model when the star reaches the WD cooling sequence.

**Stage 2** (``inlist_cool_wd``): Loads the saved WD model and continues the cooling sequence with the ``colors`` module enabled, computing synthetic photometry at every timestep.

Running the Test Suite
======================

To run both stages sequentially::

    ./rn

To run stages individually::

    ./rn1   # Stage 1: Create WD
    ./rn2   # Stage 2: Cool WD with colors

Stage 1 is computationally intensive (full stellar lifetime). If you want to skip directly to the colors demonstration, a pre-computed WD model (``wd_start.mod``) can be downloaded from [location TBD].

Output Products
===============

**History columns** (automatically added by colors module):
   - ``Mag_bol``: Bolometric magnitude
   - ``Flux_bol``: Bolometric flux at 10 pc
   - ``u``, ``g``, ``r``, ``i``, ``z``, ``y``: LSST filter magnitudes

**Key diagnostics**:
   - ``log_Teff``: Effective temperature
   - ``log_L``: Luminosity  
   - ``log_g``: Surface gravity
   - ``star_age``: Total age
   - ``cooling_time``: Time since WD formation

Analysis Scripts
================

``plot_wd_cooling.py``
----------------------

Generates publication-quality figures:

1. **WD Cooling Track CMD**: g vs (g-r) color-magnitude diagram showing the cooling sequence
2. **Multi-band Light Curves**: Magnitude evolution in all LSST bands vs cooling time
3. **Color Evolution**: (u-g), (g-r), (r-i) color indices vs Teff
4. **SED Evolution**: Spectral energy distributions at selected cooling ages

Usage::

    python plot_wd_cooling.py

Mission Relevance
=================

**Rubin/LSST**: Will detect billions of WDs across the Milky Way and nearby galaxies. Synthetic photometry enables:
   - Population synthesis predictions
   - Age dating via cooling sequences
   - Identification of unusual WD types

**Roman Space Telescope**: Deep infrared observations will probe cool WDs invisible to optical surveys.

**Gaia**: Parallaxes convert apparent to absolute magnitudes, enabling direct CMD comparisons.

Inlist Configuration
====================

Key ``&colors`` settings for WD applications::

    &colors
        use_colors = .true.
        instrument = '/colors/data/filters/LSST/LSST'
        stellar_atm = '/colors/data/stellar_models/Koester_DA/'
        distance = 3.0857d19  ! 10 pc for absolute magnitudes
        mag_system = 'AB'     ! AB system for LSST compatibility
    /

Key ``&eos`` settings::

    &eos
        use_skye = .true.  ! Skye EOS for degenerate matter
    /

References
==========

- Koester, D. 2010, Mem.S.A.It. 81, 921 (WD atmosphere models)
- Jermyn et al. 2021, ApJ 913, 72 (Skye EOS)
- LSST Science Collaboration 2009, arXiv:0912.0201

Data Requirements
=================

Before running, ensure you have:

1. **Koester DA atmosphere grid** processed via SED_Tools::

       stellar_atm = '/colors/data/stellar_models/Koester_DA/'

   Download from SVO: http://svo2.cab.inta-csic.es/theory/newov2/index.php?models=koester2

2. **LSST filter curves**::

       instrument = '/colors/data/filters/LSST/LSST'

   Available via SED_Tools from SVO Filter Profile Service.
