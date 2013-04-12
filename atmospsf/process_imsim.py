# This script can be used to process imSim simulations that were produced with atmospheric (but not
# optics) PSF, to study the statistical properties of the PSFs.

import numpy as np
import matplotlib.pyplot as plt
import os, sys
import scipy.optimize as opti
import galsim
import pylab

# First, get access to the power spectrum estimator.  It's not yet in GalSim or installed on my
# system, so have to add its dir to path.  This should work on other people's machines if they are
# sitting in the directory containing this file (process_imsim.py) in the great3-private repo.
pse_path = os.path.abspath('../../metrics')
sys.path.append(pse_path)
import pse

# Helper function: read in a single file (`filename') and the corresponding FITS cube ('cubename').
# `verbose' determines how much information gets outputted to stdout, 'do_plot' determines whether
# or not to display plots, 'n_ell' determines how many logarithmic ell bins to use for the PS
# estimation, and 'use_hsm' determines whether we use HSM moments or sextractor ones.
def process_file(filename, cubename, n_ell, verbose=True, do_plot=False, use_hsm=False):
    # Read in the data
    data = np.loadtxt(filename).transpose()

    # Separate into the fields we want
    x = 0.2*data[1] # arcsec (convert from 0.2" LSST pixels)
    y = 0.2*data[2] # arcsec (convert from 0.2" LSST pixels)
    e1 = data[3]
    e2 = data[4]
    fwhm = 3600.*data[6] # to get arcsec

    # Also, remeasure shapes from adaptive moments using the FITS cube
    imlist = galsim.fits.readCube(cubename)
    if len(imlist) != len(x):
        raise RuntimeError("Image list is not the same size as the catalog!")
    e1_hsm = np.zeros_like(e1)
    e2_hsm = np.zeros_like(e2)
    fwhm_hsm = np.zeros_like(fwhm)
    for imindex in range(len(imlist)):
        res = imlist[imindex].FindAdaptiveMom()
        e1_hsm[imindex] = res.observed_shape.e1
        e2_hsm[imindex] = res.observed_shape.e2
        fwhm_hsm[imindex] = 2.35482*res.moments_sigma*0.2 # arcsec
    if verbose:
        print "Comparing my shapes with those from files from Chihway"
    if do_plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(e1, e1_hsm, marker='+', linestyle='None')
        ax.plot(e2, e2_hsm, marker='o', linestyle='None')
        ax.set_xlabel('Sextractor e')
        ax.set_ylabel('HSM e')
        plt.savefig('diagnostics/compare_e.jpg')
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(fwhm, fwhm_hsm, marker='+', linestyle='None')
        ax.set_xlabel('Sextractor FWHM')
        ax.set_ylabel('HSM FWHM')
        plt.savefig('diagnostics/compare_fwhm.jpg')

    if use_hsm:
        # Process FWHM to get fluctuations around mean.
        dfwhm = fwhm_hsm/np.mean(fwhm_hsm)-1.0
        # Remove the large-scale correlations due to atmospheric dispersion, which are roughly like
        # a constant non-zero mean ellipticity.  For the PS, that will look like something peaked at
        # k=0, which we do not represent with our lensing engine.  Hence we remove it here, and when
        # simulating atmospheric PSFs, we can simply draw a random value of e (from within a range
        # given by the sims) to add after the fact.
        de1 = e1_hsm-np.mean(e1_hsm)
        de2 = e2_hsm-np.mean(e2_hsm)
    else:
        dfwhm = fwhm/np.mean(fwhm)-1.0
        de1 = e1-np.mean(e1)
        de2 = e2-np.mean(e2)

    # Get some basic statistics:
    if verbose:
        print "Some statistics of file ",filename
        print "Number of stars in input grid (expected 8820):",data.shape[1]
        print "Min, max, mean x (arcsec):",np.min(x), np.max(x), np.mean(x)
        print "Min, max, mean y (arcsec):",np.min(y), np.max(y), np.mean(y)
        print "Mean FWHM and std dev of fluctuations around mean:",np.mean(fwhm), np.std(dfwhm)
        print "Mean e1, std dev of e1:",np.mean(e1), np.std(e1)
        print "Mean e2, std dev of e2:",np.mean(e2), np.std(e2) 
        print "Mean HSM sigma and std dev of fluctations around mean:",np.mean(fwhm_hsm), np.std(fwhm_hsm)
        print "Corr coeff between FWHM from catalog and HSM sigma:",np.sum(fwhm_hsm*fwhm)/np.sqrt(np.sum(fwhm_hsm**2)*np.sum(fwhm**2))
        print "Mean e1, std dev of e1 from HSM:",np.mean(e1_hsm), np.std(e2_hsm)
        print "Mean e2, std dev of e2 from HSM:",np.mean(e2_hsm), np.std(e2_hsm)
        print "Corr coeff between e1 from catalog and HSM:",np.sum(e1*e1_hsm)/np.sqrt(np.sum(e1**2)*np.sum(e1_hsm**2))
        print "Showing grid:"
    if do_plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(x, y, marker='+', linestyle='None')
        ax.set_xlabel('x (arcsec)')
        ax.set_ylabel('y (arcsec)')
        plt.xlim((0,1000))
        plt.ylim((0,1000))
        ax.set_title('Locations of stars')
        plt.savefig('diagnostics/all_stars.jpg')

    # Now define expected grid locations
    x_grid_min = 200.
    y_grid_min = 48.
    dgrid = 360. # in arcsec: 10 degrees / 100
    ngrid = 21 # because this grid is only covering 2 degrees
    x_grid_index = (x-x_grid_min)/dgrid
    round_x_grid_index = np.round(x_grid_index)
    y_grid_index = (y-y_grid_min)/dgrid
    round_y_grid_index = np.round(y_grid_index)
    vecuse = (np.abs(round_x_grid_index-x_grid_index)<0.05) & (np.abs(round_y_grid_index-y_grid_index)<0.05)
    x_use = round_x_grid_index[vecuse].astype(np.int)
    y_use = round_y_grid_index[vecuse].astype(np.int)
    de1_use = de1[vecuse].astype(np.float32)
    de2_use = de2[vecuse].astype(np.float32)
    dfwhm_use = dfwhm[vecuse].astype(np.float32)
    if verbose:
        print "After grid selection:"
        print "Number of grid points:",len(x_use)
        print "Min,max x,y:",np.min(x_use), np.max(x_use), np.min(y_use), np.max(y_use)
    if do_plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(x, y, marker='o', linestyle='None')
        ax.plot(x_use, y_use, marker='+', linestyle='None')
        ax.set_xlabel('x index')
        ax.set_ylabel('y index')
        plt.savefig('diagnostics/gridstars.jpg')
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(x[vecuse], x_use, marker='o', linestyle='None')
        ax.plot(y[vecuse], y_use, marker='+', linestyle='None')
        ax.set_xlabel('x or y value')
        ax.set_ylabel('x or y index')
        plt.savefig('diagnostics/gridstars_ind.jpg')

    # Make grids of the interesting parameters.
    de1_grid = np.zeros((ngrid, ngrid))
    de2_grid = np.zeros((ngrid, ngrid))
    dfwhm_grid = np.zeros((ngrid, ngrid))
    de1_grid[y_use, x_use] = de1_use
    de2_grid[y_use, x_use] = de2_use
    dfwhm_grid[y_use, x_use] = dfwhm_use
    testvec = de1_grid.flatten()
    if do_plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.hist(dfwhm_grid,bins=20,histtype='step',normed=True)
        ax.set_xlabel('dfwhm')
        ax.set_title('histogram of dfwhm values from grid')
        plt.savefig('diagnostics/histdfwhm.jpg')
        pylab.figure()
        pylab.pcolor(dfwhm_grid)
        pylab.colorbar()
        pylab.savefig('diagnostics/dfwhm_grid_color.jpg')
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if use_hsm:
            ax.scatter(x, y, marker='o', c=fwhm_hsm, s=50*(30*(fwhm_hsm/np.mean(fwhm_hsm)-1)+1))
        else:
            ax.scatter(x, y, marker='o', c=fwhm, s=50*(30*(fwhm/np.mean(fwhm)-1)+1))
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        plt.xlim((0,1000))
        plt.ylim((0,1000))
        ax.set_title('marker color and size reflects FWHM')
        plt.savefig('diagnostics/fwhm_all.jpg')

    # PS estimation based on ellipticities
    theta = ngrid*dgrid/3600. # total size of grid in degrees, required for PSE
    atmos_e_pse = pse.PowerSpectrumEstimator(ngrid, theta, n_ell)
    ell, cee, cbb, ceb = atmos_e_pse.estimate(de1_grid, de2_grid)

    # Check best-fit power-law on large k (small scales).
    subell = ell[(ell>=500) & (ell<= 2500)]
    sublgell = np.log10(ell[(ell>=500) & (ell<= 2500)])-np.log10(1000.)
    sublgcee = np.log10(subell**2 * cee[(ell>=500) & (ell<= 2500)]/(2.*np.pi))
    sublgcbb = np.log10(subell**2 * cbb[(ell>=500) & (ell<= 2500)]/(2.*np.pi))
    A = np.vander(sublgell, 2)
    (coeffs_ee, residuals, rank, sing_vals) = np.linalg.lstsq(A, sublgcee)
    (coeffs_bb, residuals, rank, sing_vals) = np.linalg.lstsq(A, sublgcbb)
    if verbose:
        print "ell^2 * P_EE fit: amplitude at ell=1000 and power-law index are",10**coeffs_ee[1], coeffs_ee[0]
        print "ell^2 * P_EE fit: amplitude at ell=1000 and power-law index are",10**coeffs_bb[1], coeffs_bb[0]
    f_ee = 10**(coeffs_ee[0]*sublgell + coeffs_ee[1])
    f_bb = 10**(coeffs_bb[0]*sublgell + coeffs_bb[1])
    f_kol = 10**((1./3.)*sublgell + 0.5*(coeffs_ee[1]+coeffs_bb[1]))
    if do_plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(ell, ell*(ell+1)*cee/(2.*np.pi), marker='+', linestyle='None', 
                label='Observed EE power')
        ax.plot(ell, ell*(ell+1)*cbb/(2.*np.pi), marker='+', linestyle='None', 
                label='Observed BB power')
        ax.plot(ell, ell*(ell+1)*ceb/(2.*np.pi), marker='+', linestyle='None',
                label='Observed EB power')
        ax.plot(subell, f_ee, label='Best-fit power law (EE)')
        ax.plot(subell, f_bb, label='Best-fit power law (BB)')
        ax.plot(subell, f_kol, label='Kolmogorov expectation')
        lim = (1.e-6, 1.e-3)
        plt.ylim(lim)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('ell')
        ax.set_ylabel('ell*(ell+1)*C_ell/2pi')
        ax.set_title('Power spectrum from atmospheric PSF shapes')
        plt.legend(loc=4)
        min_ell = 2.*np.pi/(theta*np.pi/180.)
        max_ell = np.sqrt(2.)*ngrid*min_ell/2.
        ax.plot(np.array((min_ell,min_ell)), np.array(lim), color='black')
        ax.plot(np.array((max_ell,max_ell)), np.array(lim), color='black')
        plt.savefig('diagnostics/ps.jpg')

    # Compare with atmospheric power spectra from papers, etc.

    # Reconstruct "kappa" field based on e1, e2, and compare with dfwhm.  Check statistical
    # properties, etc.

    # Save / output / return information, incl. mean ellipticities
    # Average over all the many sub-grids?

# Main function: for now, just specify a single file directly.
expstring = '3min'
imnumstring = '92'
infile = '/Users/rmandelb/great3/chihway-sims/catalogs/'+expstring+'/focalplane_'+imnumstring+'.txt'
incube = '/Users/rmandelb/great3/chihway-sims/images/'+expstring+'/cube_'+imnumstring+'.fits.gz'
n_ell = 6
do_plot = True # make plots
verbose = True # Spew loads of diagnostics
use_hsm = True # use HSM moments instead of sextractor ones
process_file(infile, incube, n_ell, verbose=verbose, do_plot=do_plot, use_hsm=use_hsm)

# Later, will need to write code that crawls through all directories and processes everything,
# stores and postprocesses results.

