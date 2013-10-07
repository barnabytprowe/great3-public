# A quick utility to successively impose some cuts based on what's in great3sims/galaxies.py, and
# see what happens to various properties of the galaxies.
import numpy as np
import matplotlib.pyplot as plt
import galsim
import pyfits
import os

# Define the set of quantities to plot, directories, etc..
plot_vals = ["e", "sn", "hlr", "mag", "bt", "zphot", "flux_radius"]
plot_min_val = [0, 0, 0, 18, 0, 0, 0]
plot_max_val = [1, 6, 2, 23.5, 1, 2, 2]
obs_type = "space"
n_vals = len(plot_vals)
gal_dir = "./great3_data"
rgc_file = 'real_galaxy_catalog_23.5.fits'
rgc_fits_file = 'real_galaxy_catalog_23.5_fits.fits'
rgc_im_sel_file = 'real_galaxy_image_selection_info.fits'
rgc_sel_file = 'real_galaxy_selection_info.fits'
rgc_shapes_file = 'real_galaxy_23.5_shapes.fits'
rgc_dmag_file = 'real_galaxy_deltamag_info.fits'
rgc_mask_file = 'real_galaxy_mask_info.fits'

# Read in all the innumerable catalogs needed to make cuts.
rgc = galsim.RealGalaxyCatalog(rgc_file, dir=gal_dir, preload=False)
fit_catalog = pyfits.getdata(os.path.join(gal_dir, rgc_fits_file))
shapes_catalog = pyfits.getdata(os.path.join(gal_dir, rgc_shapes_file))
selection_catalog = pyfits.getdata(os.path.join(gal_dir, rgc_sel_file))
use_bulgefit = selection_catalog.field('use_bulgefit')[:,0]
im_selection_catalog = pyfits.getdata(os.path.join(gal_dir, rgc_im_sel_file))
original_sn = im_selection_catalog.field('sn_ellip_gauss')
noise_fail_val = 1.e-10 # number to assign for negative noise variance values, then discard
if obs_type == "ground":
    ind = -2
else:
    ind = 0
noise_min_var = im_selection_catalog.field('min_var_white')[:,ind]
noise_max_var = selection_catalog.field('max_var')[:,ind]
noise_max_var[noise_max_var < noise_fail_val] = noise_fail_val
flux_frac = selection_catalog.field('flux_frac')[:,ind]
resolution = selection_catalog.field('resolution')[:,ind]
dmag_catalog = pyfits.getdata(os.path.join(gal_dir, rgc_dmag_file))
dmag = dmag_catalog.field('delta_mag')
mask_catalog = pyfits.getdata(os.path.join(gal_dir, rgc_mask_file))
average_mask_adjacent_pixel_count = mask_catalog['average_mask_adjacent_pixel_count']
peak_image_pixel_count = mask_catalog['peak_image_pixel_count']
peak_image_pixel_count[peak_image_pixel_count == 0.] = 1.e-4
min_mask_dist_pixels = mask_catalog['min_mask_dist_pixels']

# Define the plot quantities.
e = np.sqrt(shapes_catalog["e1"]**2 + shapes_catalog["e2"]**2)
params = fit_catalog['sersicfit']
sn = params[:,2]
gal_q = params[:,3]
hlr = 0.03*np.sqrt(gal_q)*params[:,1]
mag = fit_catalog["mag_auto"]
bt = shapes_catalog["bulge_tot"]
zphot = fit_catalog["zphot"]
flux_radius = 0.03*fit_catalog["flux_radius"]

# Define the numbers that will be used for cuts.
min_flux_frac = 0.99
min_resolution = 1./3
approx_sn_gal = np.zeros(rgc.nobjects)
if obs_type == "ground":
    variance = 0.0046
else:
    variance = 0.00138
approx_sn_gal[noise_max_var > noise_fail_val] = \
    20.0*np.sqrt(noise_max_var[noise_max_var > noise_fail_val] / variance)
sn_min = 17.0
sn_max = 100.0
mask_cond = np.logical_or.reduce(
    [min_mask_dist_pixels > 11,
     average_mask_adjacent_pixel_count/peak_image_pixel_count < 0.2
     ])
cuts = [selection_catalog.field('to_use') == 1,
        np.abs(dmag) < 0.8,
        flux_frac >= min_flux_frac,
        resolution >= min_resolution,
        approx_sn_gal >= sn_min,
        approx_sn_gal <= sn_max,
        noise_max_var > noise_fail_val,
        shapes_catalog.field('do_meas') > -0.5,
        e < 1.,
        original_sn >= 20.,
        noise_min_var <= 0.96*variance,
        mask_cond
        ]
cut_names = ["to_use",
             "delta_mag_fit",
             "flux fraction",
             "resolution",
             "sn_lower_lim",
             "sn_upper_lim",
             "max_noise_var",
             "shape meas",
             "e < 1.",
             "original S/N >= 20",
             "noise_min_var",
             "masking cuts"
        ]
n_cuts = len(cuts)
print "Testing ",n_cuts," cuts applied in sequence"

# Set up plots for plot_vals, and loop over cut imposition to show what happens to distributions.
for plot_val_ind in range(len(plot_vals)):
    fig = plt.figure()
    ax = plt.subplot(111)
    n, bins, patches = ax.hist(eval(plot_vals[plot_val_ind]), bins=25, label='Full sample: %d'%rgc.nobjects,
                               normed=True, range=(plot_min_val[plot_val_ind],
                                                   plot_max_val[plot_val_ind]),
                               alpha = 0.5
                               )
    ax.set_xlabel(plot_vals[plot_val_ind])
    to_use = np.ones(rgc.nobjects).astype(bool)
    for cut_ind in range(n_cuts):
        to_use = np.logical_and.reduce([to_use, cuts[cut_ind]])
        if cut_ind == n_cuts - 1:
            ax.hist(
                eval(plot_vals[plot_val_ind])[to_use], bins=bins,
                label='%s: %d'%(cut_names[cut_ind], len(eval(plot_vals[plot_val_ind])[to_use])),
                normed=True,
                alpha=0.5
                )
    plt.legend()
    outfile = "galaxy_props_"+obs_type+"_"+plot_vals[plot_val_ind]+".png"
    print "Saving figure to file ",outfile
    plt.savefig(outfile)
