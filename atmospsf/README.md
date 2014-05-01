# Atmospheric PSF anisotropy power spectra
==========================================

This directory contains scripts and data files related to the atmospheric PSF
anisotropy power spectra.  The key issue here is that the LSST imSim simulations
(and others) suggest a particular functional form for the atmospheric PSF
anisotropy correlation function, but GalSim can only draw shears according to a
power spectrum.  So here we explore the imSim-based PSF anisotropy correlation
functions, tabulate power spectra corresponding to that correlation function for
various sets of parameters, and test the use of GalSim to generate the PSF
ellipticities to make sure that the output correlation functions are as
expected.

More specifically, the imSim-based xi_+ (defined for ellipticities corresponding
to distortions) tends to have the form
xi_+ = A / (1 + theta/theta_0)
for some amplitude A and scale length theta_0.  xi_- tends to be zero,
suggesting comparable E and B-mode power in the PSF anisotropies.

Typical values for A and theta_0 are as follows:

- The amplitude A tends to be in the range 1e-4 to 8e-4 for a 20s exposure and a
  telescope diameter of 6.7m, scaling as 1/(t_exp * diam).

- The scale length ranges from 0.1-2 degrees, with roughly uniform probability.

The files that explore this model are as follows:

1. `pk.math` is a Mathematica script that numerically integrates the correlation
functions (for various theta_0) with the appropriate Bessel function to get a
shear power spectrum, and tabulates the results in a series of files.

2. `pk_math/` is a directory containing the tabulated PSF anisotropy power spectra
from `pk.math`.

3. `plot_math_pk.py` is a script that was used to check the outputs of the
Mathematica script.  More specifically, it takes the output power spectra, tells
GalSim to generate a random shear field according to the power spectra, computes
the shear correlation function xi_+, and compares it to the original one given
above to ensure that they are consistent.

4. `great3_atmos_psf.py` is a script that randomly draws an atmospheric PSF
anisotropy field using the outputs in `pk_math/`, and makes some whisker plots and
other sanity checks.  Some of this code was later used in the GREAT3 simulation
scripts.

5. `corr2.params` is a file that was used to set parameters for corr2, the 2-point
correlation function code used by `plot_math_pk.py`.
