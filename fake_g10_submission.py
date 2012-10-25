#!/usr/bin/env python

import sys
import os
import numpy as np

NIMS = 200                # Number of images per set, always 200 for G10
ISETS = range(27)[1:]     # List of sets, e.g. [1, 2, 3, ..., 26]
    
true_dir = 'g10truth'    # Place where true shear cats (one per set) are stored
output_dir = 'g10fakes'  # Directory in which to output fake submissions

if len(sys.argv) != 5:
    print "usage: ./fake_g10_submission.py <label> <m [dbl]> <c [dbl]> <sigma_g [dbl]>"
    sys.exit(1)
else:
    label, m, c, sigma_g = (sys.argv[1], np.float(sys.argv[2]), np.float(sys.argv[3]),
                            np.float(sys.argv[4]))

# Loop through sets and make NIMS fake submissions for each set
for iset in ISETS:
    setfile = os.path.join(true_dir, 'g10_truth_xyg1g2_set'+str(iset)+'.dat')  # shears the same
                                                                               # within sets
    if not os.path.isfile(setfile) and os.path.isfile(setfile+'.gz'):
        from subprocess import check_call
        check_call(['gunzip', setfile+'.gz'])
    data = np.loadtxt(setfile)
    g1 = data[:, 2] # get shears from truth cat
    g2 = data[:, 3]
    for jim in xrange(NIMS):
        # Here we add shear bias and noise, g_obs = (1 + m) * g_true + c + noise
        g1_noisy = c + (1. + m) * g1 + sigma_g * np.random.normal(loc=0.0, scale=1.0, size=len(g1))
        g2_noisy = c + (1. + m) * g2 + sigma_g * np.random.normal(loc=0.0, scale=1.0, size=len(g2))
        # Write the output
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        outfile = os.path.join(
            output_dir, 'g10_v1_objs_'+str(iset)+'_'+str(jim + 1)+'.'+label+'.dat')
        output = open(outfile, 'w')
        print 'Writing '+outfile
        for k in xrange(len(g1)):
            output.write('%e  %e  %e  %e \n' % (data[k, 0], data[k, 1], g1_noisy[k], g2_noisy[k]))
        output.close()
