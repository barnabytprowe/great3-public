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
This script investigates the 3% of objects for which measurement in `shape_2comp.py` failed.
"""
import pyfits
import numpy as np
import matplotlib.pyplot as plt

# Define filenames etc.
fit_cat = '/Users/rmandelb/great3/data-23.5/real_galaxy_catalog_23.5_fits.fits'
shape_cat = '/Users/rmandelb/great3/data-23.5/real_galaxy_23.5_shapes.fits'
outfile = 'diagnostics/shape_failures.png'

# Read in fit information, etc.
fit_dat = pyfits.getdata(fit_cat)
n_fit = len(fit_dat)
print 'Read in ',n_fit,' from ',fit_cat

# Read in measurement information
shape_dat = pyfits.getdata(shape_cat)
n_shape = len(shape_dat)
print 'Read in ',n_shape,' from ',shape_cat
if n_fit != n_shape:
    raise RuntimeError("Inconsistent catalog sizes")
do_meas = shape_dat.field('do_meas')
do_meas_fail = do_meas[do_meas<-0.5]
print len(do_meas_fail),' failures'

# Make histograms of magnitude for the entire sample, and for the failing subset.
mag = fit_dat.field('mag_auto')
fig = plt.figure()
ax1 = fig.add_subplot(221)
ax1.hist(mag, normed=True, facecolor='blue', alpha=0.5, bins=0.5*np.arange(11)+18, label='All')
ax1.hist(mag[do_meas<-0.5], normed=True, facecolor='red', alpha=0.5, bins=0.25*np.arange(26)+18, label='Failures')
ax1.set_xlabel('F814W magnitude')

# Make histograms of size for the entire sample, and for the failing subset.
size = 0.03*fit_dat.field('flux_radius')
ax2 = fig.add_subplot(222)
ax2.hist(size, normed=True, facecolor='blue', alpha=0.5, bins=0.05*np.arange(40), label='All')
ax2.hist(size[do_meas<-0.5], normed=True, facecolor='red', alpha=0.5, bins=0.05*np.arange(40), label='Failures')
ax2.set_xlabel('flux_radius')

# Make histograms of bulge-to-total flux ratio for the entire sample, and for the failing subset.
bt = shape_dat.field('bulge_tot')
ax3 = fig.add_subplot(223)
ax3.hist(bt, normed=True, facecolor='blue', alpha=0.5, bins=0.05*np.arange(21), label='All')
ax3.hist(bt[do_meas<-0.5], normed=True, facecolor='red', alpha=0.5, bins=0.05*np.arange(21), label='Failures')
ax3.set_xlabel('B/T')

# Make histograms of Sersic n for the entire sample, and for the failing subset.
n = fit_dat.field('sersicfit')[:,2]
ax3 = fig.add_subplot(224)
ax3.hist(n, normed=True, facecolor='blue', alpha=0.5, bins=0.5*np.arange(17), label='All')
ax3.hist(n[do_meas<-0.5], normed=True, facecolor='red', alpha=0.5, bins=0.5*np.arange(17), label='Failures')
ax3.set_xlabel('n')

# Some overall formatting.
plt.legend(loc='upper right')
fig.tight_layout()

# Save to file.
plt.savefig(outfile)
