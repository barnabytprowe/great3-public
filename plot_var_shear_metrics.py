import os
import numpy as np
import galsim
import g3metrics

NIMS = 200
NGRID = 100   # Each image contains a grid of NGRID x NGRID galaxies
DX_GRID = 0.1 # Grid spacing (must be in degrees)
NOISE_SIGMA = 0.05 # Expected noise on each shear after shape noise pushed largely into B-mode
NBINS=8       # Number of bins of PS for metric

CTEST=1.e-4
MTEST=[1.e-3, 3.e-3, 1.e-2, 3.e-2]

GALSIM_DIR=os.path.join("Users", "browe", "great3", "galsim")

referense_ps = g3metrics.read_ps(galsim_dir=GALSIM_DIR)

# Make the truth catalogues (a list of 2D, NGRIDxNGRID numpy arrays), reusing the reference_ps
# each time for simplicity
g1true_list, g2true_list = g3metrics.make_var_truth_catalogs(
    1, NIMS, [reference_ps,], ngrid=NGRID, grid_units=galsim.degrees)

# Then generate submissions, and truth submissions
for m in MTEST:
    ksub, pEsubs, pBsubs, pEtrues, pBtrues = g3metrics.make_submission_var_shear(
        c1=CTEST, c2=CTEST, m1=m, m2=m, g1true_list=g1true_list, g2true_list=g2true_list,
        noise_sigma=NOISE_SIGMA, dx_grid=DX_GRID, nbins=NBINS)
    
