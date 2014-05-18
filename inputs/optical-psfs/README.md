# GREAT3 optical PSF models
===========================

This directory contains files used to generate the optical PSFs for GREAT3.

The contents are as follows:

1. Data files are `afta_wfirst_example_psf_exaggerated.fields_and_coefs.fits`
(for the space-based model) and are stored in the
`ground_optical_psf_zernike_coefficients_41x41` directory for the ground-based
model.

2. `GREAT3_COEFFICIENTS.ZPL` is an example ZEMAX script that was used to get the
space-based optical PSF model in the format we need for use by our scripts.

3. Ground-based PSF model: the classes are defined in `ground_optical_psf.py`,
and simple usage is demonstrated in `ground_optical_psf_example.py`.

4. Space-based PSF model: the classes are defined in `space_optical_psf.py`, and
simple usage is demonstrated in `space_optical_psf_example.py`.

The docstrings in the scripts includes usage information.  Further examples of
usage are in the great3sims package, `psf.py`.
