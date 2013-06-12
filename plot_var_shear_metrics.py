import os
import numpy as np
import matplotlib.pyplot as plt
import galsim
import g3metrics

NIMS = 200
NGRID = 100        # Each image contains a grid of NGRID x NGRID galaxies
DX_GRID = 0.1      # Grid spacing (must be in degrees)
NOISE_SIGMA = 0.05 # Expected noise on each shear after shape noise pushed largely into B-mode
NBINS_ANGULAR = 3  # Number of angular bins for correlation function metric
MIN_SEP = 0.5
MAX_SEP = 10.

NTRUESETS = 20      # Don't necessarily need to have NIMS input shears. But easiest if
                    # NTRUESETS is an integral fraction of NIMS..

CFID = 1.e-4 # Fiducial, "target" m and c values
MFID = 1.e-3 #

PLOT = False # Plot?
# Plotting ranges of interest
CMIN = CFID
CMAX = 1.e-2
MMIN = MFID
MMAX = 1.e-1
NBINS_TEST = 5      # Number of bins to plot in the ranges above
NMONTE = 30         # Number of montecarlo samples
NOISE_SIGMA = 0.05  # Noise due to pixel shot noist on a shear estimate, per galaxy

#GALSIM_DIR=os.path.join("/Path", "To", "Your", "Repo")
GALSIM_DIR=os.path.join("/Users", "browe", "great3", "galsim")

if __name__ == "__main__":

    reference_ps = g3metrics.read_ps(galsim_dir=GALSIM_DIR)

    # Make the truth catalogues (a list of 2D, NGRIDxNGRID numpy arrays), reusing the reference_ps
    # each time for simplicity
    g1true_list, g2true_list = g3metrics.make_var_truth_catalogs(
        NTRUESETS, NIMS, [reference_ps] * NTRUESETS, ngrid=NGRID, dx_grid=DX_GRID,
        grid_units=galsim.degrees)

    # Generate arrays of values for test values of c and m
    cvals = CMIN * (CMAX / CMIN)**(np.arange(NBINS_TEST) / float(NBINS_TEST - 1.)) # geom. series
    mvals = MMIN * (MMAX / MMIN)**(np.arange(NBINS_TEST) / float(NBINS_TEST - 1.))
    cgrid, mgrid = np.meshgrid(cvals, mvals) # 2D arrays covering full space

    # Create empty storage arrays
    qCF1 = np.empty((NBINS_TEST, NBINS_TEST, NMONTE))
    mCF1 = np.empty((NBINS_TEST, NBINS_TEST, NMONTE))
    cCF1 = np.empty((NBINS_TEST, NBINS_TEST, NMONTE))

    # Then generate submissions, and truth submissions
    for c, i in zip(cvals, range(NBINS_TEST)):

        print "Calculating CF metrics with c_i = "+str(c)
        for m, j in zip(mvals, range(NBINS_TEST)):

            print "Calculating CF metrics with m_i = "+str(m)
            for krepeat in range(NMONTE):

                # Make a fake submission
                theta, mapEsubs, mapBsubs, maperrsubs, mapEtrues, mapBtrues = \
                    g3metrics.make_submission_var_shear_CF(
                        c1=c, c2=c, m1=m, m2=m, g1true_list=g1true_list, g2true_list=g2true_list,
                        noise_sigma=NOISE_SIGMA, dx_grid=DX_GRID, nbins=NBINS_ANGULAR,
                        min_sep=MIN_SEP, max_sep=MAX_SEP)
                qCF1_tmp, c_tmp, m_tmp = g3metrics.metricMapCF_var_shear(
                    mapEsubs, maperrsubs, mapEtrues, NTRUESETS, nbins=NBINS_ANGULAR,
                    min_sep=MIN_SEP, max_sep=MAX_SEP, plot=PLOT)
                qCF1[i, j, krepeat] = qCF1_tmp
                mCF1[i, j, krepeat] = m_tmp
                cCF1[i, j, krepeat] = c_tmp
                print "Completed "+str(krepeat + 1)+"/"+str(NMONTE)+" Monte Carlo realizations"

    print "Saving output"
    np.save('qCF1.npy', qCF1)
    np.save('mCF1.npy', mCF1)
    np.save('cCF1.npy', cCF1)
    

