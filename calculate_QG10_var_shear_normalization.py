import os
import cPickle
import matplotlib.pyplot as plt
import numpy as np
import galsim
import g3metrics

NIMS = int(200)
NGRID = 100        # Each image contains a grid of NGRID x NGRID galaxies
DX_GRID = 0.1      # Grid spacing (must be in degrees)
NOISE_SIGMA = 0.05 # Expected noise on each shear after shape noise pushed largely into B-mode
NBINS=8            # Number of bins of PS for metric
Q10_SCALING = 1

CTEST = [1.e-4, 3.e-4, 1.e-3, 3.e-3, 1.e-2, 3.e-2]
MTEST = [1.e-3, 3.e-3, 1.e-2, 3.e-2, 1.e-1, 3.e-1]
NREPEAT = 1

#GALSIM_DIR=os.path.join("/Path", "To", "Your", "Repo")
GALSIM_DIR=os.path.join("/Users", "browe", "great3", "galsim")

reference_ps = g3metrics.read_ps(galsim_dir=GALSIM_DIR)

# Make the truth catalogues (a list of 2D, NGRIDxNGRID numpy arrays), reusing the reference_ps
# each time for simplicity
g1true_list, g2true_list = g3metrics.make_var_truth_catalogs(
    1, NIMS, [reference_ps,], ngrid=NGRID, grid_units=galsim.degrees)

# Define some empty storage arrays
qG10unnorm = np.empty((len(CTEST), len(MTEST), NREPEAT))
qG10norm = np.empty((len(CTEST), len(MTEST), NREPEAT))
qQuadPSunnorm = np.empty((len(CTEST), len(MTEST), NREPEAT))
qQuadPSnorm = np.empty((len(CTEST), len(MTEST), NREPEAT))

plt.clf()
# Then generate submissions, and truth submissions
for c, i in zip(CTEST, range(len(CTEST))):

    for m, j in zip(MTEST, range(len(MTEST))):

        for krepeat in range(NREPEAT):

            # Make a fake submission
            ksub, pEsubs, pBsubs, pEtrues, pBtrues = g3metrics.make_submission_var_shear(
                c1=c, c2=c, m1=m, m2=m, g1true_list=g1true_list, g2true_list=g2true_list,
                noise_sigma=NOISE_SIGMA, dx_grid=DX_GRID, nbins=NBINS)
            # Calculate the G10 metric for this realization
            qG10, mean_pEestG10, mean_pEtrueG10, mean_diffG10 = g3metrics.metricQuadPS_var_shear(
                ksub, pEsubs, [NOISE_SIGMA**2] * NIMS, pEtrues, scaling=Q10_SCALING,
                dx_grid=DX_GRID)
            # Calculate the QuadPS metric for this realization
            qQuadPS, mean_pEestQuadPS, mean_pEtrueQuadPS, mean_diffQuadPS = \
                g3metrics.metricQuadPS_var_shear(
                    ksub, pEsubs, [NOISE_SIGMA**2] * NIMS, pEtrues, scaling=Q10_SCALING,
                    dx_grid=DX_GRID)
            # Save the results
            qG10unnorm[i, j, krepeat] = qG10
            qQuadPSunnorm[i, j, krepeat] = qQuadPS
            # Plot the first run, just to see what this looks like...
            if i == 0 and krepeat == 0:
                plt.loglog(ksub, mean_pEtrueG10 * ksub**2)
                plt.loglog(
                    ksub, mean_pEestG10 * ksub**2,
                    label=r"m$_i$="+str(m)+", c$_i$="+str(c)+"    Q10 = "+str(qG10))
            print c, m, krepeat, qG10, qQuadPS

plt.ylabel(r'k$^2$P(k)')
plt.xlabel(r'k')
plt.ylim(1.e-5, 1.e-1)
plt.legend()
plt.savefig('plots/example_k2PkQG10_'+str(NREPEAT)+'.png')
plt.show()
