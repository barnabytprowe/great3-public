# Cosmic shear power spectra
============================

This directory contains scripts and data files related to the cosmic shear power
spectra for GREAT3.

Our approach was to get shear power spectra from iCosmo.  We used the following
cosmological parameters: w_0=-1, w_a=0, sigma_8=0.8, n_s=0.96, H_0=70h km/s/Mpc,
Omega_m=0.25.  We used the default iCosmo N(z) to get three values of
z_med=0.75, 1.00, 1.25.  Then, we added nuisance functions on top of those to
get something roughly but not exactly cosmological.

The files that are included are as follows:

- `getpk.pro`: This is the IDL script that was used to get the iCosmo outputs.

- tables/ includes the outputs of `getpk.pro`, which are used as inputs to the
  GREAT3 simulation scripts.

- `nuisance.py`: This script illustrates the nuisance functions we use to modify
  the shear power spectra.
