#!/usr/bin/env python

from sys import argv
import os
import numpy as np

NIMS = 200                # Number of images per set, was 200 for G10
NGALS_PER_IM = 1e4        # Number of galaxies per image, was 100x100 for G10
NOISE_SIGMA = 0.05        # The effective standard deviation on a single shear due to pixel noise

if len(argv) != 6:
    print ""
    print "usage: ./sim_const_shear_submission.py label m1 m2 c1 c2"
    print ""
    print "Generate fake GREAT3 submissions and save them to "
    print " ./g3subs/g3_const_shear_sub.<label>.dat"
    print ""
    exit(1)
else:
    label = argv[1]
    m1 = float(argv[2])
    m2 = float(argv[3])
    c1 = float(argv[4])
    c2 = float(argv[5])

# Read in the truth catalog (assuming that sim_const_shear_truth.py has been run)
truth = np.loadtxt('./g3truth/g3_const_shear_truth.dat')
if truth.shape[0] != NIMS:
    raise ValueError("Truth table in ./g3truth/g3_const_shear_truth.dat does not have NIMS="+
                     str(NIMS)+" entries as expected.")
g1true = truth[:, 1]
g2true = truth[:, 2]

# Then ready some arrays for the output submission
g1sub = np.empty(NIMS)
g2sub = np.empty(NIMS)

# Loop generating subs (doing this the long way - could use central limit theorem but this is super safe!)
for i in range(NIMS):
    g1gals = (1. + m1) * g1true[i] + c1 + np.random.randn(NGALS_PER_IM) * NOISE_SIGMA
    g2gals = (1. + m2) * g2true[i] + c2 + np.random.randn(NGALS_PER_IM) * NOISE_SIGMA
    g1sub[i] = np.mean(g1gals)
    g2sub[i] = np.mean(g2gals)

# Save output
if not os.path.isdir('g3subs'):
    os.mkdir('g3subs')
np.savetxt(
    './g3subs/g3_const_shear_sub.'+label+'.dat', np.array((truth[:, 0], g1sub, g2sub)).T,
    fmt=('%d', '%14.7f', '%14.7f'))

