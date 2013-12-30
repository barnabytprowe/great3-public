"""@file tabulate_const_shear_metric_rev1.py

Tabulate values of the new constant shear metrics (revision 1: Dec 13/Jan 14) for given biases, for
inclusion in the GREAT3 Handbook.
"""

import sys
import os
import numpy as np
path, module = os.path.split(__file__)
sys.path.append(os.path.join(path, "..", "server", "great3")) # Appends the folder
                                                              # great3-private/server/great3 to
                                                              # sys.path
import evaluate

NTEST = 1000
NOISE_SIGMA = {"ground": 0.15, "space": 0.10}
CVALS = evaluate.CFID * 10.**(.5 * np.arange(5))
MVALS = evaluate.MFID * 10.**(.5 * np.arange(5))

truth_dir = "/Users/browe/great3/beta/truth" # Modify to wherever truth is unpacked
# Temporarily write this into evaluate.TRUTH_DIR - HACKY!
evaluate.TRUTH_DIR = truth_dir

EXPERIMENT = "control" # Use the control truth values


if __name__ == "__main__":

    # Dicts containing arrays for storing Q_c values versus m and c, for ground and space
    qc = {"ground": np.empty((NTEST, len(CVALS))), "space": np.empty((NTEST, len(CVALS)))}
    qm = {"ground": np.empty((NTEST, len(MVALS))), "space": np.empty((NTEST, len(MVALS)))}
    for obs_type in ("ground", "space"):

        # First we build the truth table
        print "Building truth tables for control-"+obs_type+"-constant"
        _, g1true, g2true = evaluate.get_generate_const_truth(EXPERIMENT, obs_type)
    	print "Calculating Q_c values versus c for control-"+obs_type+"-constant data in GREAT3"
    	for jc, cval in enumerate(CVALS):

    		for itest in xrange(NTEST):

    			qc[obs_type][itest, jc] = 


