import sys
import galsim
import math
from training_galaxy_props import *

# read command-line arguments, which are:
#    PSF type: either Kolmogorov or Euclid
#    FWHM in arcsec, only used for Kolmogorov
#    pixel scale
#    output filename

if len(sys.argv) != 5:
    raise ValueError("Wrong number of command-line arguments!")

psf_type = sys.argv[1]
allowed_types = ['Kolmogorov', 'Euclid']
if psf_type not in allowed_types:
    raise ValueError("PSF type is not allowed!")
if psf_type is 'Kolmogorov':
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

training_galaxy_props(psf,
                      out_filename = sys.argv[4],
                      pix_scale = float(sys.argv[3]),
                      size_factor = 0.6,
                      ps_size = 48)
