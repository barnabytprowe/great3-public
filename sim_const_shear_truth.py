#!/usr/bin/env python

from sys import argv
import os
import numpy as np

NIMS = 200                # Number of images per set, always 200 for G10
SIGMA_TRUE = 0.03         # Standard deviation of true input shears
NTRUE = 20                # Don't necessarily need to have NIMS input shears, but easiest if it's an
                          # integral fraction...

if not os.path.isdir('g3truth'):
    os.mkdir('g3truth')

imagen = np.arange(NIMS, dtype=int) + 1
g1true = np.repeat(np.random.randn(NTRUE) * SIGMA_TRUE, NIMS/NTRUE)  # these should have NIMS
g2true = np.repeat(np.random.randn(NTRUE) * SIGMA_TRUE, NIMS/NTRUE)  # elements, with repetition...

np.savetxt(
    './g3truth/g3_const_shear_truth.dat', np.array((imagen, g1true, g2true)).T, 
    fmt=('%d', '%14.7f', '%14.7f'))
