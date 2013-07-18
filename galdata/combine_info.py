import pyfits
import numpy as np
import matplotlib.pyplot as plt
import os

# define filenames, etc.
cat_file_name = '/Users/rmandelb/great3/data-23.5/real_galaxy_catalog_23.5.fits'
fit_file_name = '/Users/rmandelb/great3/data-23.5/real_galaxy_catalog_23.5_fits.fits'
shape_file_name = 'real_galaxy_23.5_shapes.fits'
property_files = ['real_galaxy_catalog_23.5_props_Euclid_0.05.fits',
                  'real_galaxy_catalog_23.5_props_Euclid_0.10.fits',
                  'real_galaxy_catalog_23.5_props_Kolm_0.5_0.2.fits',
                  'real_galaxy_catalog_23.5_props_Kolm_0.65_0.2.fits',
                  'real_galaxy_catalog_23.5_props_Kolm_0.8_0.2.fits',
                  'real_galaxy_catalog_23.5_props_Kolm_0.95_0.2.fits']
exclude_file = []
exclude_file.append(1) # exclude file with this limit
n_files = len(property_files)
out_dir = './'
out_filename = 'real_galaxy_selection_info.fits'

# read in data
catalog = pyfits.getdata(cat_file_name)
n_gal = len(catalog)
print 'Read in ',n_gal,' from catalog: ',cat_file_name

fit_catalog = pyfits.getdata(fit_file_name)
print 'Read in fits from ',fit_file_name

shape_catalog = pyfits.getdata(shape_file_name)
print 'Read in shapes from ',shape_file_name

props_list = []
for filename in property_files:
    props_cat = pyfits.getdata(filename)
    print 'Read in properties from ',filename
    props_list.append(props_cat)

# eliminate failures by looping over property files
# initialize a sort of mask and then sequentially decide which to use
to_use = np.ones(n_gal, dtype=np.bool)
for i_file in range(n_files):
    to_use = to_use & (props_list[i_file]['do_meas'] == 1)
# also use flag for fit failures
bulgefit_status = fit_catalog.field('fit_status')[:,0]
to_use = to_use & ((bulgefit_status != 0) & (bulgefit_status != 5))
# use flag for failures when getting shapes needed for shape noise cancellation
do_meas = shape_catalog.field('do_meas')
to_use = to_use & ((do_meas == 0) | (do_meas == 1))
# also check for silly values of flux_frac
for i_file in range(n_files):
    ff = props_list[i_file]['flux_frac']
    to_use = to_use & ((ff >= 0.) & (ff <= 1.))

# check number left
print 'After excluding failures, using ',sum(to_use)

# info to save for all:
# ident, to_use, arrays with resolution, max_var corrected values, flux_frac
ident = catalog.field('ident')
resolution = np.zeros((n_gal,n_files-len(exclude_file)))
flux_frac = np.zeros((n_gal,n_files-len(exclude_file)))
max_var = np.zeros((n_gal,n_files-len(exclude_file)))
i_file_use = 0
for i_file in range(n_files):
    if i_file not in exclude_file:
        resolution[:,i_file_use] = props_list[i_file]['resolution']
        flux_frac[:,i_file_use] = props_list[i_file]['flux_frac']
        max_var[:,i_file_use] = props_list[i_file]['noise_var_snr_20']/1.05e-7
        i_file_use += 1

# save results to file
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
                                                       array=max_var)]
                                        ))

# write outputs
outfile = os.path.join(out_dir, out_filename)
print "Writing to file ",outfile
tbhdu.writeto(outfile, clobber=True)
