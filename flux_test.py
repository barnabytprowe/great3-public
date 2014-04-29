#!/sw64/bin/python2.7
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
"""This script provided an important test of the galaxy sample used for GREAT3: we compare the
apparent magnitudes in F814W from the COSMOS catalog, with the magnitudes inferred from the fluxes
in the parametric fits.  Significant discrepancies could be due to blending or masking failures,
so we use the outputs of this script as a way to flag problematic objects so they won't be used.

Note that we do not expect users of that dataset to have to run this script, because its outputs are
distributed with the sets of catalogs used for GREAT3.  We provide it merely for informational
purposes."""

import pyfits
import numpy as np
import math
import galsim
import os
import matplotlib.pyplot as plt
import scipy.stats as stats

# Define directory where the GREAT3 catalogs are to be found, two pre-existing catalogs, and the
# file name into which the outputs of this script should be placed (`out_catalog`).
dir = './great3-data'
fit_catalog = 'real_galaxy_catalog_23.5_fits.fits'
sel_catalog = 'real_galaxy_selection_info.fits'
out_catalog = 'real_galaxy_deltamag_info.fits'

# Read in the fit and selection catalogs.
fit_data = pyfits.getdata(os.path.join(dir, fit_catalog))
sel_data = pyfits.getdata(os.path.join(dir, sel_catalog))
use_bulgefit = sel_data['use_bulgefit'][:,0]
n = len(fit_data)
if len(sel_data) != n:
    raise RuntimeError("Catalog lengths do not match!")
print "Read in ",n," from fit and selection catalogs"

# For each object, calculate the flux (according to whether we use two- or single-component fits).
fit_flux = np.zeros(n)
print "Getting fluxes..."
for i in range(n):
    if i % 1000 == 0:
        print "    Galaxy ",i
    if use_bulgefit[i] == 1.:
        params = fit_data[i].field('bulgefit')

        # For documentation on the flux calculation here, see the great3sims/galaxies.py comments.
        # The short summary is that the files for the parametric fits list the surface brightness at
        # the effective radius measured along the major axis, so the conversion to flux is
        # non-trivial.
        bulge_q = params[11]
        bulge_hlr = 0.03*np.sqrt(bulge_q)*params[9] # arcsec 
        bulge_flux = 2.0*np.pi*3.607*(bulge_hlr**2)*params[8]/(0.03**2)

        disk_q = params[3]
        disk_hlr = 0.03*np.sqrt(disk_q)*params[1] # arcsec
        disk_flux = 2.0*np.pi*1.901*(disk_hlr**2)*params[0]/(0.03**2)

        fit_flux[i] = bulge_flux + disk_flux
    else:
        # Make a single Sersic model instead
        params = fit_data[i].field('sersicfit')
        gal_n = params[2]
        # Fudge the Sersic n if it is at the edge of the allowed n values in GalSim.
        gal_q = params[3]
        gal_hlr = 0.03*np.sqrt(gal_q)*params[1]
        if gal_n < 0.3: gal_n = 0.3
        tmp_ser = galsim.Sersic(gal_n, half_light_radius=1.)
        gal_bn = (1./tmp_ser.getScaleRadius())**(1./gal_n)
        prefactor = gal_n * math.gamma(2.*gal_n) * math.exp(gal_bn) / (gal_bn**(2.*gal_n))
        fit_flux[i] = 2.*np.pi*prefactor*(gal_hlr**2)*params[0]/0.03**2
print len(fit_flux[fit_flux<=0])," have zero or negative fluxes; fudging!"
# Fudge anything negative (which will not be selected due to various flags, anyway).
fit_flux[fit_flux<=0] = 1e-5
fit_mag = -2.5*np.log10(fit_flux)
catalog_mag = fit_data[0:n].field('mag_auto')

# Do a linear regression between the magnitudes from the parametric fits and from the caftalog.
# This should basically work if there aren't many outliers, which is indeed the case.
(a, b, r, tt, stderr) = stats.linregress(catalog_mag, fit_mag)
print "Linear regression params: ",a,b
best_fit = a*catalog_mag + b
# Now do an outlier rejection with delta mag ~ 0.5, just for some additional stability.
fit_diff = fit_mag - best_fit
to_use = np.abs(fit_diff)<0.5
print "After outlier rejection, using ",np.sum(to_use.astype(int))," of ",len(to_use)
(a, b, r, tt, stderr) = stats.linregress(catalog_mag[to_use], fit_mag[to_use])
print "Slope and intercept are ",a,b
best_fit_line = a*catalog_mag + b

# Scatter plot: Compare -2.5 log10(flux) from the parametric fits (`fit_mag`) with quoted magnitude
# from the COSMOS catalog (`catalog_mag`).
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(catalog_mag, fit_mag, label='Data')
ax.plot(catalog_mag, best_fit_line, label='Linear fit, a=%f, b=%f'%(a,b))
ax.set_xlabel('Catalog mag_auto')
ax.set_ylabel('Fit magnitude')
plt.legend(loc='upper left')
plt.savefig('fit_mags_scatter.png')
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(catalog_mag, fit_mag, label='Data')
ax.plot(catalog_mag, best_fit_line, label='Linear fit, a=%f, b=%f'%(a,b))
ax.set_xlabel('Catalog mag_auto')
ax.set_ylabel('Fit magnitude')
ax.set_xlim(20,24)
ax.set_ylim(-8,0)
plt.legend(loc='upper left')
plt.savefig('fit_mags_scatter_zoom.png')

# Histogram of differences between the magnitude from parametric fits vs. the expected value given
# the catalog magnitude.
fit_diff = fit_mag - best_fit_line
fig = plt.figure()
ax = fig.add_subplot(111)
n, bins, patches = ax.hist(fit_diff, 100, range=(-2,2))
ax.set_xlabel('Magnitude difference, Fit - best_fit_line')
plt.savefig('fit_mag_histdiff.png')
print "Number with difference >1 mag (either direction):",len(fit_diff[np.abs(fit_diff)>=1])
print "Number with difference >0.5 mag (either direction):",len(fit_diff[np.abs(fit_diff)>=0.5])
print "Number that are >1 mag fainter:",len(fit_diff[fit_diff>=1])

tbhdu = pyfits.new_table(pyfits.ColDefs([pyfits.Column(name='IDENT',
                                                       format='J',
                                                       array=fit_data.field('ident')),
                                         pyfits.Column(name='delta_mag',
                                                       format='D',
                                                       array=fit_diff)]
                                        ))
outfile = os.path.join(dir,out_catalog)
print "Writing to file ",outfile
tbhdu.writeto(outfile, clobber=True)
