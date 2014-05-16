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
This is a utility to check the catalogs for a training dataset, including the following questions:

* check the contents (what fields are included)
* make sure we can use it with GalSim
* pick a random entry and check all info to make sure that it makes sense.  This should include an
  image of the galaxy and PSF in addition to entries in catalog.
* histograms of various quantities
* 2d scatter plots of various quantities
* failure patterns for fits

Users will have to modify it to point to the location of the COSMOS catalogs on their system.
"""
import pyfits
import numpy as np
import galsim
import os
import matplotlib.pyplot as plt

# File name definitions, etc.  Users must change the first one to point to the proper location!
in_dir = '/Users/rmandelb/great3/data-23.5'
cat_file = 'real_galaxy_catalog_23.5.fits'
fit_file = 'real_galaxy_catalog_23.5_fits.fits'
out_dir = 'diagnostics'
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

# Check contents of catalog file.
filename = os.path.join(in_dir, cat_file)
print "Information about catalog file :"
print pyfits.info(filename)
filename = os.path.join(in_dir, fit_file)
print "Information about file with galaxy fits:"
print pyfits.info(filename)

# Make sure we can read in the catalogs / images into GalSim.
real_galaxy_catalog = galsim.RealGalaxyCatalog(cat_file, dir=in_dir)
gal = galsim.RealGalaxy(real_galaxy_catalog, random = True)
print "Read in RealGalaxyCatalog with ",real_galaxy_catalog.nobjects
print "index used: ",gal.index
print "corresponding IDENT: ",real_galaxy_catalog.ident[gal.index]

# Make sure we can read in the catalog with info about Claire's fits.
# Check correspondence in number of objects.
fit_catalog = pyfits.getdata(filename)
print "Read in catalog of ",len(fit_catalog)," fits from file ",filename
print "Field names: ",fit_catalog.names
if len(fit_catalog) != real_galaxy_catalog.nobjects:
    raise RuntimeError("Number of objects does not correspond!")

# Spew information about some random entry, to test correspondence
# with original files.  Use the random entry that we already chose.
print "Fit info for this randomly selected object:"
print fit_catalog[gal.index]
print "Real galaxy catalog info for this randomly selected object:"
print real_galaxy_catalog.ident[gal.index], \
    real_galaxy_catalog.gal_file_name[gal.index], \
    real_galaxy_catalog.PSF_file_name[gal.index], \
    real_galaxy_catalog.gal_hdu[gal.index], \
    real_galaxy_catalog.PSF_hdu[gal.index], \
    real_galaxy_catalog.pixel_scale[gal.index], \
    real_galaxy_catalog.variance[gal.index], \
    real_galaxy_catalog.mag[gal.index], \
    real_galaxy_catalog.weight[gal.index]
test_im = gal.original_image.draw(dx=real_galaxy_catalog.pixel_scale[gal.index])
outfile = os.path.join(out_dir,out_pref+'test_im.fits')
print "Writing original image for that object to ",outfile
test_im.write(outfile)

# 1d histograms of various things: mag_auto, flux_radius, zphot,
# weight, sersic n
fig = plt.figure()
ax1 = fig.add_subplot(321)
ax1.hist(fit_catalog.field('mag_auto'), bins=0.25*np.arange(26)+18)
ax1.set_xlabel('F814W magnitude')
ax2 = fig.add_subplot(322)
ax2.hist(0.03*fit_catalog.field('flux_radius'), bins=0.05*np.arange(40))
ax2.set_xlabel('Flux radius')
ax3 = fig.add_subplot(323)
ax3.hist(fit_catalog.field('zphot'), bins=0.1*np.arange(41))
ax3.set_xlabel('Photo-z')
ax4 = fig.add_subplot(324)
ax4.hist(real_galaxy_catalog.weight, bins=0.2*np.arange(20))
ax4.set_xlabel('weight')
sersic_fit = fit_catalog.field('sersicfit')
sersic_n = sersic_fit[:,2]
print sersic_n.shape
ax5 = fig.add_subplot(325)
ax5.hist(sersic_n, bins=0.25*np.arange(40))
ax5.set_xlabel('Sersic n')
outfile = os.path.join(out_dir,out_pref+'1d_hist.png')
print "Writing 1d histograms to file ",outfile
fig.tight_layout()
plt.savefig(outfile)

# 2d scatter plots: mag_auto from fits vs. real_galaxy_catalog; size
# vs. mag, mag_auto vs. zphot, size vs. weight, sersic size vs. flux_rad
fig = plt.figure()
ax1 = fig.add_subplot(321)
ax1.scatter(fit_catalog.field('mag_auto'), real_galaxy_catalog.mag)
ax1.set_xlabel('Mag in fit catalog')
ax1.set_ylabel('Mag in RGC')
ax2 = fig.add_subplot(322)
ax2.scatter(fit_catalog.field('mag_auto'), 0.03*fit_catalog.field('flux_radius'))
ax2.set_xlabel('Mag')
ax2.set_ylabel('Flux radius')
ax3 = fig.add_subplot(323)
ax3.scatter(fit_catalog.field('mag_auto'), fit_catalog.field('zphot'))
ax3.set_xlabel('Mag')
ax3.set_ylabel('Photo-z')
ax4 = fig.add_subplot(324)
ax4.scatter(0.03*fit_catalog.field('flux_radius'),real_galaxy_catalog.weight)
ax4.set_xlabel('Flux radius')
ax4.set_ylabel('weight')
sersic_q = sersic_fit[:,3]
sersic_hlr = 0.03*np.sqrt(sersic_q)*sersic_fit[:,1]
ax5 = fig.add_subplot(325)
ax5.scatter(0.03*fit_catalog.field('flux_radius'), sersic_hlr)
ax5.set_ylim((0.,3.))
ax5.set_xlim((0.,3.))
ax5.set_xlabel('Flux radius')
ax5.set_ylabel('Sersic HLR')
outfile = os.path.join(out_dir,out_pref+'2d_scatter.png')
print "Writing 2d scatter plots to file ",outfile
fig.tight_layout()
plt.savefig(outfile)

# Failure rates:
sersicfit_status = fit_catalog.field('fit_status')[:,4]
bulgefit_status = fit_catalog.field('fit_status')[:,0]
s_fail_ind = np.where((sersicfit_status == 0) | (sersicfit_status == 5))[0]
print "Number of Sersic fit failures: ",len(s_fail_ind)
b_fail_ind = np.where((bulgefit_status == 0) | (bulgefit_status == 5))[0]
print "Number of 2 component fit failures: ",len(b_fail_ind)
other_fail_ind = np.where((sersic_n < 0) & (sersicfit_status != 0) & (sersicfit_status != 5))[0]
print "Number of Sersic fit failures that were not properly flagged: ",len(other_fail_ind)
# Let's check distribution of size, mag for failures vs. for everything.
fig = plt.figure()
ax1 = fig.add_subplot(211)
ax1.hist(fit_catalog.field('mag_auto'), normed=True, facecolor='blue', alpha=0.5, label='All')
ax1.hist(fit_catalog.field('mag_auto')[s_fail_ind], normed=True, facecolor='red', alpha=0.5, label='Sersic failures')
ax1.hist(fit_catalog.field('mag_auto')[b_fail_ind], normed=True, facecolor='green', alpha=0.5, label='2 component failures')
ax1.set_xlabel('F814W magnitude')
ax2 = fig.add_subplot(212)
ax2.hist(fit_catalog.field('flux_radius'), normed=True, facecolor='blue', alpha=0.5, label='All')
ax2.hist(fit_catalog.field('flux_radius')[s_fail_ind], normed=True, facecolor='red', alpha=0.5, label='All')
ax2.hist(fit_catalog.field('flux_radius')[b_fail_ind], normed=True, facecolor='green', alpha=0.5, label='All')
ax2.set_xlabel('Flux radius')
plt.legend()
outfile = os.path.join(out_dir,out_pref+'1d_hist_failures.png')
print "Writing 1d histograms for failure cases to file ",outfile
fig.tight_layout()
plt.savefig(outfile)
