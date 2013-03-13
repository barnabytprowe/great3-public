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
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(x, y, marker='+')
        ax.set_xlabel('x (arcsec)')
        ax.set_ylabel('y (arcsec)')
        ax.set_title('Locations of stars')
        plt.show()

# Main function: for now, just specify a single file directly.
infile = '/Users/rmandelb/great3/chihway-sims/catalogs/1min/focalplane_47.txt'
verbose = 1
process_file(infile, verbose)

# Later, will need to write code that crawls through all directories and processes everything.
