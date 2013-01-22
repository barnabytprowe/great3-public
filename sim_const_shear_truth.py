#!/usr/bin/env python

import os
import numpy as np

NIMS = 200               # Number of images per set, always 200 for G10
SIGMA_TRUE = 0.04        # Standard deviation of true input shears for normal distribution
RANGE_TRUE = 0.08        # Range of true input shears for a uniform distribution
NTRUE = 50               # Don't necessarily need to have NIMS input shears. But easiest if NTRUE an
                         # integral fraction of NIMS...

def make_truth_normal_dist():
    """Generate truth catalogues with a N(0, SIGMA_TRUE) distribution of input truth values in
    g1 and g2.
    """
    if not os.path.isdir('g3truth'):
        os.mkdir('g3truth')

    imagen = np.arange(NIMS, dtype=int) + 1
    g1true = np.repeat(np.random.randn(NTRUE) * SIGMA_TRUE, NIMS/NTRUE)  # these should have NIMS
    g2true = np.repeat(np.random.randn(NTRUE) * SIGMA_TRUE, NIMS/NTRUE)  # elements, with repetition
    np.savetxt(
        './g3truth/g3_const_shear_truth.dat', np.array((imagen, g1true, g2true)).T, 
        fmt=('%d', '%14.7f', '%14.7f'))
    return

def make_truth_uniform_dist():
    """Generate truth catalogues with a U(-RANGE_TRUE, RANGE_TRUE) distribution of input truth
    values in g1 and g2.
    """
    if not os.path.isdir('g3truth'):
        os.mkdir('g3truth')

    imagen = np.arange(NIMS, dtype=int) + 1
    g1true = np.repeat((2. * np.random.rand(NTRUE) - 1.) * SIGMA_TRUE, NIMS/NTRUE)
    g2true = np.repeat((2. * np.random.rand(NTRUE) - 1.) * SIGMA_TRUE, NIMS/NTRUE)
    np.savetxt(
        './g3truth/g3_const_shear_truth.dat', np.array((imagen, g1true, g2true)).T, 
        fmt=('%d', '%14.7f', '%14.7f'))
    return

if __name__ is "__main__":
    make_truth_normal_dist()
#    make_truth_uniform_dist()
