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
This script determines the effective shapes of the 2-component galaxy models.  For galaxies for
which we use a single Sersic profile, the isophotes are elliptical and we can just use the fit
parameters for that Sersic profile to determine the shape.  For 2-component galaxies we have to make
the images and measure their shapes directly using adaptive moments.  The purpose is to get shapes
we can use when defining the B-mode shape noise field.
"""
import galsim
import pyfits
import os
import numpy as np

# Define file and directory names, etc.
in_dir = '.'
in_filename = 'real_galaxy_catalog_23.5_fits.fits'
out_dir = '.'
out_filename = 'real_galaxy_23.5_shapes.fits'
pix_scale = 0.03

# Define PSF.  We use something small, just to ensure Nyquist sampling.  What is below corresponds
# to a diffraction-limited PSF for a wavelength of 800 nm and a 4m telescope.
airy = galsim.Airy(lam_over_diam=206265.*800.e-9/4.)
pix = galsim.Pixel(pix_scale)
psf = galsim.Convolve(airy, pix)
im_psf = psf.draw(dx=pix_scale)

# Read in catalog of COSMOS galaxies with fit parameters.
infile = os.path.join(in_dir, in_filename)
dat = pyfits.getdata(infile)
n = len(dat)
print "Read in ",n," from ",infile

# Create output arrays for the following quantities: shape (e1 and e2), bulge-to-total flux ratio
# B/T (bt), and flags do_meas (did we measure the shape for a 2-component fit?) and use_bulgefit
# (should we use 2-component fit, or 1-component Sersic fit?).
e1 = np.zeros(n)
e2 = np.zeros(n)
bt = np.zeros(n)
do_meas = np.zeros(n)
use_bulgefit = np.ones(n)
# Loop over objects.
for i in range(n):
    if i % 1000 == 0:
        print "...",i

    # Get fit parameters.
    params = dat[i].field('bulgefit')
    sparams = dat[i].field('sersicfit')
    # Get the status flag for the fits.
    bstat = dat[i].field('fit_status')[0]
    sstat = dat[i].field('fit_status')[4]
    # Get the precomputed bulge-to-total ratio for the 2-component fits.
    dvc_btt = dat[i].field('fit_dvc_btt')
    # Get the precomputed median absolute deviation for the 1- and 2-component fits.
    bmad = dat[i].field('fit_mad_b')
    smad = dat[i].field('fit_mad_s')

    # First decide if we can / should use bulgefit, otherwise sersicfit.
    # This decision process depends on: the status flags for the fits, the bulge-to-total ratios (if
    # near 0 or 1, just use single component fits), the sizes for the bulge and disk (if <=0 then
    # use single component fits), the axis ratios for the bulge and disk (if <0.051 then use single
    # component fits), and a comparison of the median absolute deviations to see which is better.
    if bstat<1 or bstat>4 or dvc_btt<0.1 or dvc_btt>0.9 or np.isnan(dvc_btt) or params[9]<=0 or params[1]<=0 or params[11]<0.051 or params[3]<0.051 or smad<bmad:
        use_bulgefit[i] = 0
        # Then check if sersicfit is viable; if not, this object is a total failure:
        if sstat<1 or sstat>4 or sparams[1]<=0 or sparams[0]<=0:
            use_bulgefit[i] = -1
            do_meas[i] = -1
            e1[i] = -10.
            e2[i] = -10.
            bt[i] = -10.
            continue

    # Now, if we're supposed to use the 2-component fits, get all the parameters.
    if use_bulgefit[i]:
        bulge_q = params[11]
        bulge_beta = params[15]*galsim.radians
        bulge_hlr = 0.03*np.sqrt(bulge_q)*params[9]
        bulge_flux = 2.0*np.pi*3.607*(bulge_hlr**2)*params[8]

        disk_q = params[3]
        disk_beta = params[7]*galsim.radians
        disk_hlr = 0.03*np.sqrt(disk_q)*params[1]
        disk_flux = 2.0*np.pi*1.901*(disk_hlr**2)*params[0]
        
        bfrac = bulge_flux/(bulge_flux+disk_flux)

        # Make sure the bulge flux fraction is not nonsense.  If it is, give the outputs to indicate
        # that this is the case.
        if bfrac < 0 or bfrac > 1 or np.isnan(bfrac):
            e1[i] = -10.
            e2[i] = -10.
            do_meas[i] = -1
            use_bulgefit[i] = -1
            continue

        # Make the GalSim objects.
        bt[i] = bfrac
        bulge = galsim.Sersic(4.0, half_light_radius = bulge_hlr, flux = bulge_flux)
        disk = galsim.Exponential(half_light_radius = disk_hlr, flux = disk_flux)
        if bulge_q > 0.:
            bulge.applyShear(q = bulge_q, beta = bulge_beta)
        if disk_q > 0.:
            disk.applyShear(q = disk_q, beta = disk_beta)
        gal = bulge+disk
        obj = galsim.Convolve(gal, psf, gsparams = galsim.GSParams(maximum_fft_size=15000))
        # Try to estimate the shear using the galsim.hsm routines.  Set failure flags if it doesn't
        # work.
        try:
            im = obj.draw(dx = pix_scale)
            res = galsim.hsm.EstimateShear(im, im_psf, guess_sig_gal=20)
            e1[i] = res.corrected_e1
            e2[i] = res.corrected_e2
            do_meas[i] = 1.
        except RuntimeError:
            e1[i] = -10.0
            e2[i] = -10.0
            do_meas[i] = -1.
    # For single-component fits, life is much easier.  Just take the shape from the quoted axis
    # ratio and position angle.
    else:
        e = galsim.Shear(q = sparams[3], beta = sparams[7]*galsim.radians)
        e1[i] = e.e1
        e2[i] = e.e2

# Save results to file:
# - Bulge-to-total flux ratio
# - e1
# - e2
# - Flag indicating whether shape measurement was possible for 2-component fits
# - Flag indicating whether 2-component fits should be used
tbhdu = pyfits.new_table(pyfits.ColDefs([pyfits.Column(name='IDENT',
                                                       format='J',
                                                       array=dat.field('IDENT')),
                                         pyfits.Column(name='bulge_tot',
                                                       format='D',
                                                       array=bt),
                                         pyfits.Column(name='e1',
                                                       format='D',
                                                       array=e1),
                                         pyfits.Column(name='e2',
                                                       format='D',
                                                       array=e2),
                                         pyfits.Column(name='do_meas',
                                                       format='J',
                                                       array=do_meas),
                                         pyfits.Column(name='use_bulgefit',
                                                       format='J',
                                                       array=use_bulgefit)]
                                        ))

# Check for number of failures.
fail_ind = np.where(do_meas < -0.5)[0]
print len(fail_ind),' failures'

# Write outputs.
outfile = os.path.join(out_dir, out_filename)
print "Writing to file ",outfile
tbhdu.writeto(outfile, clobber=True)

