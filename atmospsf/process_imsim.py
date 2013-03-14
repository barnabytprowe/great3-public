# This script can be used to process imSim simulations that were produced with atmospheric (but not
# optics) PSF, to study the statistical properties of the PSFs.

import numpy as np
import matplotlib.pyplot as plt

# Helper function: read in a single file
def process_file(filename, verbose):
    # Read in the data
    data = np.loadtxt(filename).transpose()

    # Separate into the fields we want
    x = 0.2*data[1] # arcsec (convert from 0.2" LSST pixels)
    y = 0.2*data[2] # arcsec (convert from 0.2" LSST pixels)
    e1 = data[3]
    e2 = data[4]
    fwhm = 3600.*data[6] # to get arcsec

    # Process FWHM to get fluctuations around mean
    dfwhm = fwhm/np.mean(fwhm)-1.0

    # Get some basic statistics:
    if verbose:
        print "Some statistics of file ",filename
        print "Number of stars in input grid (expected 8820):",data.shape[1]
        print "Min, max, mean x (arcsec):",np.min(x), np.max(x), np.mean(x)
        print "Min, max, mean y (arcsec):",np.min(y), np.max(y), np.mean(y)
        print "Mean FWHM and std dev of fluctuations around mean:",np.mean(fwhm), np.std(dfwhm)
        print "Mean e1, std dev of e1:",np.mean(e1), np.std(e1)
        print "Mean e2, std dev of e2:",np.mean(e2), np.std(e2) 
        print "Showing grid:"
        #fig = plt.figure()
        #ax = fig.add_subplot(111)
        #ax.plot(x, y, marker='+')
        #ax.set_xlabel('x (arcsec)')
        #ax.set_ylabel('y (arcsec)')
        #plt.xlim((0,1000))
        #plt.ylim((0,1000))
        #ax.set_title('Locations of stars')
        #plt.show()

    # Now define expected grid locations
    x_grid_min = 200.
    y_grid_min = 48.
    dgrid = 360. # 10 degrees / 100
    ngrid = 21 # because this grid is only covering 2 degrees
    x_grid_index = (x-x_grid_min)/dgrid
    round_x_grid_index = np.round(x_grid_index)
    y_grid_index = (y-y_grid_min)/dgrid
    round_y_grid_index = np.round(y_grid_index)
    x_use = round_x_grid_index[(np.abs(round_x_grid_index-x_grid_index)<0.05) & (np.abs(round_y_grid_index-y_grid_index)<0.05)].astype(np.int)
    y_use = round_y_grid_index[(np.abs(round_x_grid_index-x_grid_index)<0.05) & (np.abs(round_y_grid_index-y_grid_index)<0.05)].astype(np.int)
    e1_use = e1[(np.abs(round_x_grid_index-x_grid_index)<0.05) & (np.abs(round_y_grid_index-y_grid_index)<0.05)].astype(np.int)
    e2_use = e2[(np.abs(round_x_grid_index-x_grid_index)<0.05) & (np.abs(round_y_grid_index-y_grid_index)<0.05)].astype(np.int)
    dfwhm_use = dfwhm[(np.abs(round_x_grid_index-x_grid_index)<0.05) & (np.abs(round_y_grid_index-y_grid_index)<0.05)].astype(np.int)
    if verbose:
        print "After grid selection:"
        print "Number of grid points:",len(x_use)
        print "Min,max x,y:",np.min(x_use), np.max(x_use), np.min(y_use), np.max(y_use)
        #fig = plt.figure()
        #ax = fig.add_subplot(111)
        #ax.plot(x_use, y_use, marker='+')
        #ax.set_xlabel('x index')
        #ax.set_ylabel('y index')
        #plt.show()

    # Make grids of the interesting parameters.
    e1_grid = np.zeros((ngrid, ngrid))
    e2_grid = np.zeros((ngrid, ngrid))
    dfwhm_grid = np.zeros((ngrid, ngrid))
    e1_grid[x_use, y_use] = e1_use
    e2_grid[x_use, y_use] = e2_use
    dfwhm_grid[x_use, y_use] = dfwhm_use

# Main function: for now, just specify a single file directly.
infile = '/Users/rmandelb/great3/chihway-sims/catalogs/1min/focalplane_47.txt'
verbose = 1
process_file(infile, verbose)

# Later, will need to write code that crawls through all directories and processes everything.
