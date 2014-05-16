# Copyright (c) 2014, the GREAT3 executive committee (http://www.great3challenge.info/?q=contacts)
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted
# provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions
# and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of
# conditions and the following disclaimer in the documentation and/or other materials provided with
# the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to
# endorse or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
This is a simple driver script for `training_galaxy_props.py` that constructs a PSF object and pixel
scale based on command-line inputs, and then calls `training_galaxy_props.py` with that PSF.  For
GREAT3, we ran this script 6 times as follows:

    python run_props.py Euclid foo 0.05 real_galaxy_catalog_23.5_real_props_Euclid_0.05.fits
    python run_props.py Euclid foo 0.10 real_galaxy_catalog_23.5_real_props_Euclid_0.10.fits
    python run_props.py Kolmogorov 0.5 0.2 real_galaxy_catalog_23.5_real_props_Kolm_0.5_0.2.fits
    python run_props.py Kolmogorov 0.65 0.2 real_galaxy_catalog_23.5_real_props_Kolm_0.65_0.2.fits
    python run_props.py Kolmogorov 0.8 0.2 real_galaxy_catalog_23.5_real_props_Kolm_0.8_0.2.fits
    python run_props.py Kolmogorov 0.95 0.2 real_galaxy_catalog_23.5_real_props_Kolm_0.95_0.2.fits

For the two space-based cases (the first two runs), we considered two pixel scales, one of which is
relevant for single-epoch sims (0.05") and the other for multi-epoch sims (0.10").  While the PSF
was constructed to be a simplified Euclid-esque PSF, it does not matter too much whether we use
Euclid or WFIRST-AFTA parameters, since most galaxies pass all the cuts in the space-based sims
anyway.

For the ground-based sims, we consider four values of PSF FWHM (0.5", 0.65", 0.8", and 0.95"), and
when we place cuts on galaxies in the sims, we actually interpolate between those four outputs.

Note that these PSFs are simplified and, in particular, circular, so our use of them to make cuts
means that the selection of galaxies in the simulations does not depend on the direction of the PSF
ellipticity in the sims themselves.
"""
import sys
import galsim
import math
from training_galaxy_props import *

# We need to read the command-line arguments, which are:
# - PSF type: either Kolmogorov or Euclid
# - FWHM in arcsec, only used for Kolmogorov
# - pixel scale
# - output filename

if len(sys.argv) != 5:
    raise ValueError("Wrong number of command-line arguments!")

psf_type = sys.argv[1]
allowed_types = ['Kolmogorov', 'Euclid']
if psf_type not in allowed_types:
    raise ValueError("PSF type is not allowed!")
# Make the PSF object.
if psf_type == 'Kolmogorov':
    psf = galsim.Kolmogorov(fwhm = float(sys.argv[2]))
else:
    lambda_m = 800.e-9
    primary_diam_m = 1.2
    secondary_diam_m = 0.4
    jitter_rms = 0.02                  # arcsec
    lam_over_diam = (lambda_m / primary_diam_m) * 3600. * 180. / math.pi # arcsec
    obscuration = secondary_diam_m / primary_diam_m # linear obscuration ratio
    euclid_psf = galsim.OpticalPSF(
        lam_over_diam=lam_over_diam, obscuration=obscuration, coma1=-0.04, defocus=0.09,
        astig2=-0.03, astig1=0.01, oversampling=2.5)
    jitter_psf = galsim.Gaussian(sigma=jitter_rms)
    psf = galsim.Convolve(euclid_psf, jitter_psf)

# Now call the training galaxy properties script.
training_galaxy_props(psf,
                      out_filename = sys.argv[4],
                      pix_scale = float(sys.argv[3]),
                      size_factor = 0.6,
                      ps_size = 48)
