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
"""@file tabulate_correlated_const_shear_metric_rev1.py

Tabulate values of the new constant shear metrics (revision 1: Feb 14) for given biases, including
intra branch correlations.

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
import tabulate_const_shear_metric_rev1


# Set the following module-scope constants equal to the uncorrelated tabulator code values, they
# are unlikely to change
NTEST = tabulate_const_shear_metric_rev1.NTEST
NGALS_PER_IMAGE = tabulate_const_shear_metric_rev1.NGALS_PER_IMAGE
NOISE_SIGMA = tabulate_const_shear_metric_rev1.NOISE_SIGMA
TRUTH_DIR = tabulate_const_shear_metric_rev1.TRUTH_DIR
EXPERIMENT = tabulate_const_shear_metric_rev1.EXPERIMENT

# These values might change in tests, set them here (will be obvious if they are not kept in sync!)
CVALS = evaluate.CFID * 10.**(.5 * np.arange(5))
MVALS = evaluate.MFID * 10.**(.5 * np.arange(5))

# Inter-submission subfield correlation coefficient (estimated from tests using HSM REGAUSS and
# im3shape, see calculate_im3shape_regauss_correlation.py)
RHO = 0.6

NRUNS = 100 # Number of runs to do from which to estimate the average std(Q_c), see comment below.

# Filenames for final mean standard deviation results averaged over NRUNS runs:
CSTDOUTFILE = os.path.join(
    "results","tabulated_correlated_std_const_Q_c_versus_c_NRUNS"+str(NRUNS)+".pkl")
MSTDOUTFILE = os.path.join(
    "results","tabulated_correlated_std_const_Q_c_versus_m_NRUNS"+str(NRUNS)+".pkl")

if __name__ == "__main__":

    import cPickle
    # Note I've updated this code because I noticed that the standard deviation of correlated Q_c
    # values, which is the main thing of interest that gets calculated by this script, is a strong
    # function of what the *mean* Q_c is per bin.  As these are all correlated, this changes quite
    # a bit from run to run!  Therefore I am going to repeat the calculations below NRUN times to
    # sample from a range of these results and take the average std(Q_c) to get something more
    # representative
    qclist = []
    qmlist = []
    for krun in range(NRUNS):

        # Dicts containing arrays for storing Q_c values versus m and c, for ground and space
        qc = {"ground": np.empty((NTEST, len(CVALS))), "space": np.empty((NTEST, len(CVALS)))}
        qm = {"ground": np.empty((NTEST, len(MVALS))), "space": np.empty((NTEST, len(MVALS)))}
        coutfile = os.path.join(
            "results", "tabulated_correlated_const_Q_c_versus_c_norm"+str(krun)+".pkl")
        moutfile = os.path.join(
            "results", "tabulated_correlated_const_Q_c_versus_m_norm"+str(krun)+".pkl")
        if not os.path.isfile(coutfile):
            for obs_type in ("ground", "space",):

                print "Calculating Q_c values versus c for control-"+obs_type+\
                    "-constant data in GREAT3"
                print "NOISE_SIGMA = "+str(NOISE_SIGMA[obs_type])
                print "RHO = "+str(RHO)
                # First we build the truth table
                print "Building truth tables for control-"+obs_type+"-constant"
                subfield_index, g1true, g2true = evaluate.get_generate_const_truth(
                    EXPERIMENT, obs_type, truth_dir=TRUTH_DIR)
                rotations = evaluate.get_generate_const_rotations(
                    EXPERIMENT, obs_type, truth_dir=TRUTH_DIR)
                for jc, cval in enumerate(CVALS):

                    # Build the submissions
                    g1sub, g2sub = g3metrics.make_multiple_submissions_const_shear(
                        cval, 0., evaluate.MFID, evaluate.MFID, g1true, g2true, NGALS_PER_IMAGE,
                        NOISE_SIGMA[obs_type], rotate_cs=rotations, nsubmissions=NTEST, rho=RHO)
                    # Loop over submissions evaluating metric
                    for itest in xrange(NTEST):

                        fdsub, subfile = tempfile.mkstemp(suffix=".dat")
                        with os.fdopen(fdsub, "wb") as fsub:
                            np.savetxt(
                                fsub, np.array((subfield_index, g1sub[:, itest],
                                g2sub[:, itest])).T, fmt="%d %e %e")
                        qc[obs_type][itest, jc] = evaluate.q_constant(
                            subfile, EXPERIMENT, obs_type, just_q=True, truth_dir=TRUTH_DIR)
                        os.remove(subfile)

                    print "mean(Q_c), std(Q_c) = "+str(qc[obs_type][:, jc].mean())+", "+\
                        str(qc[obs_type][:, jc].std())+" for "+str(NTEST)+" sims (with c+ = "+\
                        str(cval)+", obs_type = "+str(obs_type)+")"
                    print

            print "Saving pickled Q_c versus c dict to "+coutfile
            with open(coutfile, "wb") as fout: cPickle.dump(qc, fout)
            print
        else:
            with open(coutfile, "rb") as fin: qc = cPickle.load(fin)
        if not os.path.isfile(moutfile):
            for obs_type in ("ground", "space",):

                print "Calculating Q_c values versus m for control-"+obs_type+\
                    "-constant data in GREAT3"
                print "NOISE_SIGMA = "+str(NOISE_SIGMA[obs_type])
                print "RHO = "+str(RHO)
                # First we build the truth table
                print "Building truth tables for control-"+obs_type+"-constant"
                subfield_index, g1true, g2true = evaluate.get_generate_const_truth(
                    EXPERIMENT, obs_type, truth_dir=TRUTH_DIR)
                rotations = evaluate.get_generate_const_rotations(
                    EXPERIMENT, obs_type, truth_dir=TRUTH_DIR)
                for jm, mval in enumerate(MVALS):

                    # Build the submissions
                    g1sub, g2sub = g3metrics.make_multiple_submissions_const_shear(
                        evaluate.CFID, 0., mval, mval, g1true, g2true, NGALS_PER_IMAGE,
                        NOISE_SIGMA[obs_type], rotate_cs=rotations, nsubmissions=NTEST, rho=RHO)
                    # Loop over submissions evaluating the metric
                    for itest in xrange(NTEST):

                        fdsub, subfile = tempfile.mkstemp(suffix=".dat")
                        with os.fdopen(fdsub, "wb") as fsub:
                            np.savetxt(
                                fsub, np.array((subfield_index, g1sub[:, itest],
                                g2sub[:, itest],)).T, fmt="%d %e %e")
                        qm[obs_type][itest, jm] = evaluate.q_constant(
                            subfile, EXPERIMENT, obs_type, just_q=True, truth_dir=TRUTH_DIR)
                        os.remove(subfile)

                    print "mean(Q_c), std(Q_c) = "+str(qm[obs_type][:, jm].mean())+", "+\
                        str(qm[obs_type][:, jm].std())+" for "+str(NTEST)+" sims (with m = "+\
                        str(mval)+", obs_type = "+str(obs_type)+")"
                    print

            print "Saving pickled Q_c versus m dict to "+moutfile
            with open(moutfile, "wb") as fout: cPickle.dump(qm, fout)
            print
        else:
            with open(moutfile, "rb") as fin: qm = cPickle.load(fin)
        print
        print "Table of Q_c (ground sims) at constant m = mfid = "+str(evaluate.MFID)
        print "    c+       Q_c    std(Q_c)"
        for c, Q, dQ in zip(CVALS, np.mean(qc["ground"], axis=0), np.std(qc["ground"], axis=0)):

            print "{:8f} {:8.3f} {:8.3f}".format(c, Q, dQ)

        print
        print "Table of Q_c (space sims) at constant m = mfid = "+str(evaluate.MFID)
        print "    c+       Q_c    std(Q_c)"
        for c, Q, dQ in zip(CVALS, np.mean(qc["space"], axis=0), np.std(qc["space"], axis=0)):

            print "{:8f} {:8.3f} {:8.3f}".format(c, Q, dQ)

        print
        print "Table of Q_c (ground sims) at constant c = cfid = "+str(evaluate.CFID)
        print "    m        Q_c    std(Q_c)"
        for m, Q, dQ in zip(MVALS, np.mean(qm["ground"], axis=0), np.std(qm["ground"], axis=0)):

            print "{:8f} {:8.3f} {:8.3f}".format(m, Q, dQ)

        print
        print "Table of Q_c (space sims) at constant c = cfid = "+str(evaluate.CFID)
        print "    m        Q_c    std(Q_c)"
        for m, Q, dQ in zip(MVALS, np.mean(qm["space"], axis=0), np.std(qm["space"], axis=0)):

            print "{:8f} {:8.3f} {:8.3f}".format(m, Q, dQ)

        qclist.append(qc)
        qmlist.append(qm)

    # Then tabulate the average statistics across these NRUNS runs
    qcmean = {"ground": np.zeros(len(CVALS)), "space": np.zeros(len(CVALS))}
    sqcmean = {"ground": np.zeros(len(CVALS)), "space": np.zeros(len(CVALS))}
    qmmean = {"ground": np.zeros(len(MVALS)), "space": np.zeros(len(MVALS))}
    sqmmean = {"ground": np.zeros(len(MVALS)), "space": np.zeros(len(MVALS))}
    for obs_type in ("ground", "space"):

        for qc, qm in zip(qclist, qmlist):

            qcmean[obs_type] += np.mean(qc[obs_type], axis=0) / float(NRUNS)
            qmmean[obs_type] += np.mean(qm[obs_type], axis=0) / float(NRUNS)
            sqcmean[obs_type] += np.std(qc[obs_type], axis=0) / float(NRUNS)
            sqmmean[obs_type] += np.std(qm[obs_type], axis=0) / float(NRUNS)

        print
        print
        print "OVERALL table of average Q_c and std(Q_c) from NRUNS = "+str(NRUNS)+" tests"
        print "At constant m1 = m2 = mfid = "+str(evaluate.MFID)+", "+obs_type+" sims"
        print "    c        Q_c    std(Q_c)"
        for c, Q, dQ in zip(CVALS, qcmean[obs_type], sqcmean[obs_type]):

            print "{:8f} {:8.3f} {:8.3f}".format(c, Q, dQ)

        print "At constant c1 = cfid = "+str(evaluate.CFID)+", "+obs_type+" sims"
        print "    m        Q_c    std(Q_c)"
        for m, Q, dQ in zip(MVALS, qmmean[obs_type], sqmmean[obs_type]):

            print "{:8f} {:8.3f} {:8.3f}".format(m, Q, dQ)

    with open(CSTDOUTFILE, "wb") as fout: cPickle.dump(sqcmean, fout)
    with open(MSTDOUTFILE, "wb") as fout: cPickle.dump(sqmmean, fout)
