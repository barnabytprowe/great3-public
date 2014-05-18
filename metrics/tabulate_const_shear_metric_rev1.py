# Copyright (c) 2014, the GREAT3 executive committee (http://www.great3challenge.info/?q=contacts)
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted
# provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions
# and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of
# conditions and the following disclaimer in the documentation and/or other materials provided with
# the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to
# endorse or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""@file tabulate_const_shear_metric_rev1.py

Tabulate values of the new constant shear metrics (revision 1: Feb 14) for given biases.

See https://github.com/barnabytprowe/great3-public/wiki/The-revised-GREAT3-metrics,-February-2014
"""

import sys
import os
import tempfile
import numpy as np
import g3metrics
path, module = os.path.split(__file__)
sys.path.append(os.path.join(path, "..", "server", "great3")) # Appends the folder
                                                              # great3-private/server/great3 to
                                                              # sys.path
import evaluate

NTEST = 1000
NGALS_PER_IMAGE = 10000
NOISE_SIGMA = {"ground": 0.15, "space": 0.10}
CVALS = evaluate.CFID * 10.**(.5 * np.arange(5))
MVALS = evaluate.MFID * 10.**(.5 * np.arange(5))

TRUTH_DIR = "/Users/browe/great3/beta/truth" # Modify to wherever truth is unpacked

EXPERIMENT = "control" # Use the control truth values


if __name__ == "__main__":

    import cPickle
    # Dicts containing arrays for storing Q_c values versus m and c, for ground and space
    qc = {"ground": np.empty((NTEST, len(CVALS))), "space": np.empty((NTEST, len(CVALS)))}
    qm = {"ground": np.empty((NTEST, len(MVALS))), "space": np.empty((NTEST, len(MVALS)))}
    coutfile = os.path.join("results", "tabulated_const_Q_c_versus_c_norm2.pkl")
    moutfile = os.path.join("results", "tabulated_const_Q_c_versus_m_norm2.pkl")

    if not os.path.isfile(coutfile):
        for obs_type in ("ground", "space",):

            print "Calculating Q_c values versus c for control-"+obs_type+"-constant data in GREAT3"
            print "NOISE_SIGMA = "+str(NOISE_SIGMA[obs_type])
            # First we build the truth table
            print "Building truth tables for control-"+obs_type+"-constant"
            subfield_index, g1true, g2true = evaluate.get_generate_const_truth(
                EXPERIMENT, obs_type, truth_dir=TRUTH_DIR)
            rotations = evaluate.get_generate_const_rotations(
                EXPERIMENT, obs_type, truth_dir=TRUTH_DIR)
            for jc, cval in enumerate(CVALS):

                for itest in xrange(NTEST):

                    # Build the submission
                    g1sub, g2sub = g3metrics.make_submission_const_shear(
                        cval, 0., evaluate.MFID, evaluate.MFID, g1true, g2true, NGALS_PER_IMAGE,
                        NOISE_SIGMA[obs_type], rotate_cs=rotations)
                    fdsub, subfile = tempfile.mkstemp(suffix=".dat")
                    with os.fdopen(fdsub, "wb") as fsub:
                        np.savetxt(
                            fsub, np.array((subfield_index, g1sub, g2sub)).T, fmt="%d %e %e")
                    # Then evaluate Q_c
                    qc[obs_type][itest, jc] = evaluate.q_constant(
                        subfile, EXPERIMENT, obs_type, just_q=True, truth_dir=TRUTH_DIR)
                    os.remove(subfile)
                    print "Test %4d / %4d (c+ = %.3e) Q_c = %.4f" % (
                        itest+1, NTEST, cval, qc[obs_type][itest, jc])

                print "Mean Q_c = "+str(qc[obs_type][:, jc].mean())+" for "+str(NTEST)+\
                        " sims (with c+ = "+str(cval)+", obs_type = "+str(obs_type)+")"
                print

        print "Saving pickled Q_c versus c dict to "+coutfile
        with open(coutfile, "wb") as fout:
            cPickle.dump(qc, fout)
        print
    else:
        with open(coutfile, "rb") as fin:
            qc = cPickle.load(fin)

    if not os.path.isfile(moutfile):
        for obs_type in ("ground", "space",):

            print "Calculating Q_c values versus m for control-"+obs_type+"-constant data in GREAT3"
            print "NOISE_SIGMA = "+str(NOISE_SIGMA[obs_type])
            # First we build the truth table
            print "Building truth tables for control-"+obs_type+"-constant"
            subfield_index, g1true, g2true = evaluate.get_generate_const_truth(
                EXPERIMENT, obs_type, truth_dir=TRUTH_DIR)
            rotations = evaluate.get_generate_const_rotations(
                EXPERIMENT, obs_type, truth_dir=TRUTH_DIR)
            for jm, mval in enumerate(MVALS):

                for itest in xrange(NTEST):

                    # Build the submission
                    g1sub, g2sub = g3metrics.make_submission_const_shear(
                        evaluate.CFID, 0., mval, mval, g1true, g2true, NGALS_PER_IMAGE,
                        NOISE_SIGMA[obs_type], rotate_cs=rotations)
                    fdsub, subfile = tempfile.mkstemp(suffix=".dat")
                    with os.fdopen(fdsub, "wb") as fsub:
                        np.savetxt(fsub, np.array((subfield_index, g1sub, g2sub)).T, fmt="%d %e %e")
                    # Then evaluate Q_c
                    qm[obs_type][itest, jm] = evaluate.q_constant(
                        subfile, EXPERIMENT, obs_type, just_q=True, truth_dir=TRUTH_DIR)
                    os.remove(subfile)
                    print "Test %4d / %4d (m = %.3e) Q_c = %.4f" % (
                        itest+1, NTEST, mval, qm[obs_type][itest, jm])

                print "Mean Q_c = "+str(qm[obs_type][:, jm].mean())+" for "+str(NTEST)+\
                    " sims (with m = "+str(mval)+", obs_type = "+str(obs_type)+")"
                print

        print "Saving pickled Q_c versus m dict to "+moutfile
        with open(moutfile, "wb") as fout:
            cPickle.dump(qm, fout)
        print
    else:
        with open(moutfile, "rb") as fin:
            qm = cPickle.load(fin)

    print ""
    print "Table of Q_c (ground sims) at constant m = mfid = "+str(evaluate.MFID)
    print "    c+       Q   "
    for c, Q in zip(CVALS, np.mean(qc["ground"], axis=0)):

        print "{:8f} {:8.3f}".format(c, Q)

    print ""
    print "Table of Q_c (space sims) at constant m = mfid = "+str(evaluate.MFID)
    print "    c+       Q   "
    for c, Q in zip(CVALS, np.mean(qc["space"], axis=0)):

        print "{:8f} {:8.3f}".format(c, Q)

    print ""
    print "Table of Q_c (ground sims) at constant c = cfid = "+str(evaluate.CFID)
    print "    m        Q   "
    for m, Q in zip(MVALS, np.mean(qm["ground"], axis=0)):

        print "{:8f} {:8.3f}".format(m, Q)

    print ""
    print "Table of Q_c (space sims) at constant c = cfid = "+str(evaluate.CFID)
    print "    m        Q   "
    for m, Q in zip(MVALS, np.mean(qm["space"], axis=0)):

        print "{:8f} {:8.3f}".format(m, Q)

    



