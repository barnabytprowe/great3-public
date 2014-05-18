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
This script takes the outputs of `training_galaxy_props_real.py` for 6 different PSFs defined in
`run_props.py`, and combines them into a single file that is used by the GREAT3 simulation scripts
to select galaxies for which the use of real galaxy images is not problematic in some way.
"""
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import os

# Define filenames, etc.
property_files = ['real_galaxy_catalog_23.5_real_props_Euclid_0.05.fits',
                  'real_galaxy_catalog_23.5_real_props_Euclid_0.10.fits',
                  'real_galaxy_catalog_23.5_real_props_Kolm_0.5_0.2.fits',
                  'real_galaxy_catalog_23.5_real_props_Kolm_0.65_0.2.fits',
                  'real_galaxy_catalog_23.5_real_props_Kolm_0.8_0.2.fits',
                  'real_galaxy_catalog_23.5_real_props_Kolm_0.95_0.2.fits']
sn_index = 3 # index of file on list to use to get S/N in original image: should be 2-5, since did
             # not calculate for 0-1.  (It's the original image, so we only really needed to do it
             # once.)
exclude_file = []
n_files = len(property_files)
out_dir = './'
out_filename = 'real_galaxy_image_selection_info.fits'

# Read in the list of properties from `training_galaxy_props_real.py` for all the different options
# for PSF and pixel scale.
props_list = []
for filename in property_files:
    props_cat = pyfits.getdata(filename)
    n_gal = len(props_cat)
    print 'Read in properties from ',filename,' for ',n_gal
    props_list.append(props_cat)

# Info to save for all galaxies:
sn_ellip_gauss = props_list[sn_index]['sn_ellip_gauss']
min_var_white = np.zeros((n_gal,n_files-len(exclude_file)))
i_file_use = 0
for i_file in range(n_files):
    if i_file not in exclude_file:
        min_var_white[:,i_file_use] = props_list[i_file]['min_var_white']
        i_file_use += 1

# Save results to file.
tbhdu = pyfits.new_table(pyfits.ColDefs([pyfits.Column(name='sn_ellip_gauss',
                                                       format='D',
                                                       array=sn_ellip_gauss),
                                         pyfits.Column(name='min_var_white',
                                                       format='6D',
                                                       array=min_var_white)]
                                        ))
# Write outputs.
outfile = os.path.join(out_dir, out_filename)
print "Writing to file ",outfile
tbhdu.writeto(outfile, clobber=True)
