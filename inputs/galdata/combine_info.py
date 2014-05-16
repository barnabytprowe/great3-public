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
This script takes the outputs of `run_props.py` and combines them into a single file that is used by
the GREAT3 simulation scripts to select galaxies.
"""
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import os

# Define filenames, etc. (including a system-dependent directory containing the data files and
# another for the outputs).
in_dir = "/Users/rmandelb/great3/data-23.5"
out_dir = './'
cat_file_name = os.path.join(in_dir,'real_galaxy_catalog_23.5.fits')
fit_file_name = os.path.join(in_dir,'real_galaxy_catalog_23.5_fits.fits')
shape_file_name = os.path.join(in_dir,'real_galaxy_23.5_shapes.fits')
property_files = ['real_galaxy_catalog_23.5_props_Euclid_0.05.fits',
                  'real_galaxy_catalog_23.5_props_Euclid_0.10.fits',
                  'real_galaxy_catalog_23.5_props_Kolm_0.5_0.2.fits',
                  'real_galaxy_catalog_23.5_props_Kolm_0.65_0.2.fits',
                  'real_galaxy_catalog_23.5_props_Kolm_0.8_0.2.fits',
                  'real_galaxy_catalog_23.5_props_Kolm_0.95_0.2.fits']
exclude_file = []
exclude_file.append(1) # exclude file with this index, since answers are nearly the same as for the
                       # previous one (except the noise variance which simply differs by a factor of 4)
n_files = len(property_files)
out_filename = 'real_galaxy_selection_info.fits'

# Read in the basic catalog.
catalog = pyfits.getdata(cat_file_name)
n_gal = len(catalog)
print 'Read in ',n_gal,' from catalog: ',cat_file_name

# Read in the catalog with parametric fit information.
fit_catalog = pyfits.getdata(fit_file_name)
print 'Read in fits from ',fit_file_name

# Read in the catalog that has the shape information already precomputed, from `shape_2comp.py`.
shape_catalog = pyfits.getdata(shape_file_name)
print 'Read in shapes from ',shape_file_name

# Read in the list of properties from `training_galaxy_props.py` for all the different options for
# PSF and pixel scale.
props_list = []
for filename in property_files:
    props_cat = pyfits.getdata(filename)
    print 'Read in properties from ',filename
    props_list.append(props_cat)

# Eliminate failures by looping over property files.
# We do this by initializing a mask and then sequentially eliminating those that should be used from
# each of the property files.
to_use = np.ones(n_gal, dtype=np.bool)
# First, we get rid of those for which the shear couldn't be estimated, since we need that to get a
# resolution factor to use for the size cut.
for i_file in range(n_files):
    to_use = to_use & (props_list[i_file]['do_meas'] == 1)
# Also use flag for fit failures: should make sure use_bulgefit==0 or 1 for all values of seeing,
# otherwise it had some difficulty constructing the object and measuring basic properties.
for i_file in range(n_files):
    use_bulgefit = props_list[i_file]['use_bulgefit']
    to_use = to_use & ((use_bulgefit == 0) | (use_bulgefit == 1))
# Use a flag for failures when getting shapes needed for shape noise cancellation.
do_meas = shape_catalog.field('do_meas')
to_use = to_use & ((do_meas == 0) | (do_meas == 1))
# Also check for silly values of flux_frac, the fraction of flux included in the postage stamp.
# Obviously this should lie in the range 0-1 unless something extremely weird has happened.
for i_file in range(n_files):
    ff = props_list[i_file]['flux_frac']
    to_use = to_use & ((ff >= 0.) & (ff <= 1.))

# Check number of galaxies remaining.
print 'After excluding failures, using ',sum(to_use)

# Collect the info to save for all galaxies:
#    ident, to_use, arrays with resolution, max_var corrected values, flux_frac
# Also at this stage we include the correction for the incorrect flux normalization in
# `training_galaxy_props.py`, which gives the factor of 1.05e-7 as described in that script.
ident = catalog.field('ident')
resolution = np.zeros((n_gal,n_files-len(exclude_file)))
flux_frac = np.zeros((n_gal,n_files-len(exclude_file)))
max_var = np.zeros((n_gal,n_files-len(exclude_file)))
use_bulgefit = np.zeros((n_gal,n_files-len(exclude_file)))
i_file_use = 0
for i_file in range(n_files):
    if i_file not in exclude_file:
        resolution[:,i_file_use] = props_list[i_file]['resolution']
        flux_frac[:,i_file_use] = props_list[i_file]['flux_frac']
        max_var[:,i_file_use] = props_list[i_file]['noise_var_snr_20']/1.05e-7
        use_bulgefit[:,i_file_use] = props_list[i_file]['use_bulgefit']
        i_file_use += 1

# Save results to file.
tbhdu = pyfits.new_table(pyfits.ColDefs([pyfits.Column(name='IDENT',
                                                       format='J',
                                                       array=ident),
                                         pyfits.Column(name='to_use',
                                                       format='J',
                                                       array=to_use),
                                         pyfits.Column(name='resolution',
                                                       format='5D',
                                                       array=resolution),
                                         pyfits.Column(name='flux_frac',
                                                       format='5D',
                                                       array=flux_frac),
                                         pyfits.Column(name='max_var',
                                                       format='5D',
                                                       array=max_var),
                                         pyfits.Column(name='use_bulgefit',
                                                       format='5D',
                                                       array=use_bulgefit)]
                                        ))

# Write outputs.
outfile = os.path.join(out_dir, out_filename)
print "Writing to file ",outfile
tbhdu.writeto(outfile, clobber=True)
