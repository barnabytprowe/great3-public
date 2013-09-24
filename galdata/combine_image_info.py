import pyfits
import numpy as np
import matplotlib.pyplot as plt
import os

# define filenames, etc.
property_files = ['real_galaxy_catalog_23.5_real_props_Euclid_0.05.fits',
                  'real_galaxy_catalog_23.5_real_props_Euclid_0.10.fits',
                  'real_galaxy_catalog_23.5_real_props_Kolm_0.5_0.2.fits',
                  'real_galaxy_catalog_23.5_real_props_Kolm_0.65_0.2.fits',
                  'real_galaxy_catalog_23.5_real_props_Kolm_0.8_0.2.fits',
                  'real_galaxy_catalog_23.5_real_props_Kolm_0.95_0.2.fits']
sn_index = 3 # index of file on list to use to get S/N in original image: should be 2-5, since did not calculate for 0-1
exclude_file = []
n_files = len(property_files)
out_dir = './'
out_filename = 'real_galaxy_image_selection_info.fits'

# read in data
props_list = []
for filename in property_files:
    props_cat = pyfits.getdata(filename)
    n_gal = len(props_cat)
    print 'Read in properties from ',filename,' for ',n_gal
    props_list.append(props_cat)

# info to save for all:
sn_ellip_gauss = props_list[sn_index]['sn_ellip_gauss']
min_var_white = np.zeros((n_gal,n_files-len(exclude_file)))
i_file_use = 0
for i_file in range(n_files):
    if i_file not in exclude_file:
        min_var_white[:,i_file_use] = props_list[i_file]['min_var_white']
        i_file_use += 1

# save results to file
tbhdu = pyfits.new_table(pyfits.ColDefs([pyfits.Column(name='sn_ellip_gauss',
                                                       format='D',
                                                       array=sn_ellip_gauss),
                                         pyfits.Column(name='min_var_white',
                                                       format='6D',
                                                       array=min_var_white)]
                                        ))
# write outputs
outfile = os.path.join(out_dir, out_filename)
print "Writing to file ",outfile
tbhdu.writeto(outfile, clobber=True)
