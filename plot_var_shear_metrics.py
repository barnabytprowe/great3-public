import os
import numpy as np
import matplotlib.pyplot as plt
import galsim
import g3metrics

NIMS = 8
NGRID = 500         # Each image contains a grid of NGRID x NGRID galaxies
DX_GRID = 0.02      # Grid spacing (must be in degrees)
NOISE_SIGMA = 0.05  # Expected noise on each shear after shape noise pushed largely into B-mode
NBINS_ANGULAR = 15  # Number of angular bins for correlation function metric
MIN_SEP = DX_GRID
MAX_SEP = 10.

USE_ERRORS = True
CORRECT_B = False
OUTFILE_SUFFIX = "fine_useerrors"

NTRUESETS = 8      # Don't necessarily need to have NIMS input shears. But easiest if
                   # NTRUESETS is an integral fraction of NIMS..

CFID = 1.e-4 # Fiducial, "target" m and c values
MFID = 1.e-3 #

PLOT = False # Plot while calculating?
# Plotting ranges of interest
CMIN = CFID
CMAX = 1.e-2
MMIN = MFID
MMAX = 1.e-1
NBINS_TEST = 3      # Number of bins to plot in the ranges above
NMONTE = 50         # Number of montecarlo samples

#GALSIM_DIR=os.path.join("/Path", "To", "Your", "Repo")
GALSIM_DIR=os.path.join("/Users", "browe", "great3", "galsim")

QFILE = os.path.join('results', 'qCF1'+OUTFILE_SUFFIX+'.npy')
MFILE = os.path.join('results', 'mCF1'+OUTFILE_SUFFIX+'.npy')
CFILE = os.path.join('results', 'cCF1'+OUTFILE_SUFFIX+'.npy')

if __name__ == "__main__":

    if os.path.isfile(QFILE) and os.path.isfile(MFILE) and os.path.isfile(CFILE):
        qCF1 = np.load(QFILE)
        mCF1 = np.load(MFILE)
        cCF1 = np.load(CFILE)
        # Generate arrays of values for test values of c and m
        cvals = CMIN * (CMAX / CMIN)**(np.arange(NBINS_TEST) / float(NBINS_TEST - 1.)) # geo series
        mvals = MMIN * (MMAX / MMIN)**(np.arange(NBINS_TEST) / float(NBINS_TEST - 1.))
        cgrid, mgrid = np.meshgrid(cvals, mvals) # 2D arrays covering full space
 
    else:
        reference_ps = g3metrics.read_ps(galsim_dir=GALSIM_DIR)
        # Make the truth catalogues (a list of 2D, NGRIDxNGRID numpy arrays), reusing the
        # reference_ps each time for simplicity
        g1true_list, g2true_list = g3metrics.make_var_truth_catalogs(
            NTRUESETS, NIMS, [reference_ps] * NTRUESETS, ngrid=NGRID, dx_grid=DX_GRID,
            grid_units=galsim.degrees)
        # Generate arrays of values for test values of c and m
        cvals = CMIN * (CMAX / CMIN)**(np.arange(NBINS_TEST) / float(NBINS_TEST - 1.)) # geo series
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
                            c1=c, c2=c, m1=m, m2=m, g1true_list=g1true_list,
                            g2true_list=g2true_list, noise_sigma=NOISE_SIGMA, dx_grid=DX_GRID,
                            nbins=NBINS_ANGULAR, min_sep=MIN_SEP, max_sep=MAX_SEP)
                    #import pdb; pdb.set_trace()
                    qCF1_tmp, c_tmp, m_tmp = g3metrics.metricMapCF_var_shear_mc(
                        mapEsubs, maperrsubs, mapEtrues, mapBtrues, NTRUESETS, nbins=NBINS_ANGULAR,
                        min_sep=MIN_SEP, max_sep=MAX_SEP, plot=PLOT, correct_B_theory=CORRECT_B,
                        use_errors=USE_ERRORS)
                    qCF1[i, j, krepeat] = qCF1_tmp
                    mCF1[i, j, krepeat] = m_tmp
                    cCF1[i, j, krepeat] = c_tmp
                    print "Completed "+str(krepeat + 1)+"/"+str(NMONTE)+" Monte Carlo realizations"
        print "Saving output"
        np.save(QFILE, qCF1)
        np.save(MFILE, mCF1)
        np.save(CFILE, cCF1)

    # Get basic statistics
    qmean = np.mean(qCF1, axis=2)
    mmean = np.mean(mCF1, axis=2)
    cmean = np.mean(cCF1, axis=2)
    qstds = np.std(qCF1, axis=2)
    mstds = np.std(mCF1, axis=2)
    cstds = np.std(cCF1, axis=2)
    qerrs = qstds / np.sqrt(NMONTE)
    merrs = mstds / np.sqrt(NMONTE)
    cerrs = cstds / np.sqrt(NMONTE)

    # Plot the bestfitting m values and the standard deviation of the NMONTE results
    plt.clf()
    plt.axhline(ls='--', color='k')
    for i in range(len(cvals)):
        plt.errorbar(
            mvals * (1. + 0.03 * i), mmean[i, :], yerr=mstds[i, :], fmt='+',
            label='c = %.2e'%cvals[i])
    plt.plot(mvals, mvals, 'k')
    plt.xscale('log')
    plt.xlim(3.e-4, 3.e0)
    plt.legend()
    plt.title('Best fitting m values and standard deviation of test popn.')
    plt.ylabel('Best fitting m')
    plt.xlabel('Input m')
    plt.savefig(os.path.join('plots', 'mvals_stds_CF1_'+OUTFILE_SUFFIX+'.png'))
    # Plot the bestfitting m values and the standard errors of the NMONTE results
    plt.clf()
    plt.axhline(ls='--', color='k')
    for i in range(len(cvals)):
        plt.errorbar(
            mvals * (1. + 0.03 * i), mmean[i, :], yerr=merrs[i, :], fmt='+',
            label='c = %.2e'%cvals[i])
    plt.plot(mvals, mvals, 'k')
    plt.xscale('log')
    plt.xlim(3.e-4, 3.e0)
    plt.legend()
    plt.title('Best fitting m values and standard error on mean')
    plt.ylabel('Best fitting m')
    plt.xlabel('Input m')
    plt.savefig(os.path.join('plots', 'mvals_errs_CF1_'+OUTFILE_SUFFIX+'.png'))

    # Plot the bestfitting c^2 values and the standard deviation of the NMONTE results
    plt.clf()
    plt.axes([0.175, 0.125, 0.775, 0.775])
    plt.axhline(ls='--', color='k')
    for i in range(len(mvals)):
        plt.errorbar(
            cvals**2 * (1. + 0.03 * i), cmean[:, i], yerr=cstds[:, i], fmt='+',
            label='m = %.2e'%mvals[i])
    plt.plot(cvals**2, cvals**2, 'k')
    plt.xscale('log')
    plt.xlim(3.e-9, 1.e-1)
    plt.legend()
    plt.title(r'Best fitting c$^2$ values and standard deviation of test popn.')
    plt.ylabel(r'Best fitting c$^2$')
    plt.xlabel(r'Input c$^2$')
    plt.savefig(os.path.join('plots', 'cvals_stds_CF1_'+OUTFILE_SUFFIX+'.png'))
    # Plot the bestfitting c^2 values and the standard errors of the NMONTE results
    plt.clf()
    plt.axes([0.175, 0.125, 0.775, 0.775])
    plt.axhline(ls='--', color='k')
    for i in range(len(mvals)):
        plt.errorbar(
            cvals**2 * (1. + 0.03 * i), cmean[:, i], yerr=cerrs[:, i], fmt='+',
            label='m = %.2e'%mvals[i])
    plt.plot(cvals**2, cvals**2, 'k')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(3.e-9, 1.e-1)
    plt.legend()
    plt.title(r'Best fitting c$^2$ values and standard error on mean')
    plt.ylabel(r'Best fitting c$^2$')
    plt.xlabel(r'Input c$^2$')
    plt.savefig(os.path.join('plots', 'cvals_errs_CF1_'+OUTFILE_SUFFIX+'.png'))

    # Plot the Q metric values and the standard deviation of the NMONTE results
    plt.clf()
    plt.axhline(ls='--', color='k')
    for i in range(len(mvals)):
        plt.errorbar(
            cvals**2 * (1. + 0.03 * i), qmean[:, i], yerr=qstds[:, i], fmt='+',
            label='m = %.2e'%mvals[i])
    plt.xscale('log')
    plt.xlim(1.e-3, 1.e-10)
    plt.legend()
    plt.title('Mean Q metric and standard deviation of test popn.')
    plt.ylabel('QCF1')
    plt.xlabel(r'Input c$^2$')
    plt.savefig(os.path.join('plots', 'Qs_stds_CF1_'+OUTFILE_SUFFIX+'.png'))
    # Plot Q values and the standard errors of the NMONTE results
    plt.clf()
    plt.axhline(ls='--', color='k')
    for i in range(len(mvals)):
        plt.errorbar(
            cvals**2 * (1. + 0.03 * i), qmean[:, i], yerr=qerrs[:, i], fmt='+',
            label='m = %.2e'%mvals[i])
    plt.xscale('log')
    plt.xlim(1.e-3, 1.e-10)
    plt.legend()
    plt.title('Mean Q metric and standard error on mean')
    plt.ylabel('QCF1')
    plt.xlabel(r'Input c$^2$')
    plt.savefig(os.path.join('plots', 'Qs_errs_CF1_'+OUTFILE_SUFFIX+'.png'))
