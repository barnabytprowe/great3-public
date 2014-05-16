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
"""@file tabulate_correlated_variable_shear_metric_rev1.py

Tabulate values of the new constant shear metrics (revision 1: Feb 14) for given biases, including
intra branch correlations.

See https://github.com/barnabytprowe/great3-public/wiki/The-revised-GREAT3-metrics,-February-2014
"""

import sys
import os
import numpy as np
path, module = os.path.split(__file__)
sys.path.append(os.path.join(path, "..", "server", "great3")) # Appends the folder
                                                              # great3-private/server/great3 to
                                                              # sys.path
import evaluate
import tabulate_variable_shear_metric_rev1

# Set the following module-scope constants equal to the uncorrelated tabulator code values, they
# are unlikely to change
TRUTH_DIR = tabulate_variable_shear_metric_rev1.TRUTH_DIR
EXPERIMENT = tabulate_variable_shear_metric_rev1.EXPERIMENT

# These values might change in tests, set them here (will be obvious if they are not kept in sync!)
CVALS = evaluate.CFID * 10.**(.5 * np.arange(5))
MVALS = evaluate.MFID * 10.**(.5 * np.arange(5))

# Inter-submission subfield correlation coefficient (estimated from tests using HSM REGAUSS and
# im3shape, see calculate_variable_i3_regauss_correlation.py)
RHO = 0.45

NTEST = tabulate_variable_shear_metric_rev1.NTEST
NRUNS = 100 # Number of runs to do from which to estimate the average std(Q_v), see comment below.

# Filenames for final mean standard deviation results averaged over NRUNS runs:
CSTDOUTFILE = os.path.join(
    "results","tabulated_correlated_std_variable_Q_v_versus_c_NRUNS"+str(NRUNS)+".pkl")
MSTDOUTFILE = os.path.join(
    "results","tabulated_correlated_std_variable_Q_v_versus_m_NRUNS"+str(NRUNS)+".pkl")


def make_multiple_variable_submissions(nsubmissions, map_E_ref, map_E_unitc, c, m, cholesky):
    """Make a fake submission based on map_E_true, map_E_unitc.
    """
    # Basic checks first
    if cholesky.shape[0] != cholesky.shape[1]:
        raise ValueError("Input Cholesky matrix must be square")
    if cholesky.shape[0] != nsubmissions * evaluate.NBINS_THETA:
        raise ValueError(
            "Input Cholesky matrix must be square with dimensions "+
            "(nsubmissions * evaluate.NBINS_THETA)")
    # Define the array for storing all the submissions (note that this will be flattened as required
    # later)
    map_E_submissions = np.empty((evaluate.NBINS_THETA, nsubmissions, evaluate.NFIELDS))
    for ifield in range(evaluate.NFIELDS):

        map_E_ref_field = map_E_ref[
            ifield * evaluate.NBINS_THETA: (ifield + 1) * evaluate.NBINS_THETA]
        map_E_unitc_field = map_E_unitc[
            ifield * evaluate.NBINS_THETA: (ifield + 1) * evaluate.NBINS_THETA]
        map_E_subs_field = (
            (1. + 2. * m + m**2) * np.tile(map_E_ref_field, nsubmissions) +
            c**2 * np.tile(map_E_unitc_field, nsubmissions) + 
            np.dot(cholesky, np.random.randn(evaluate.NBINS_THETA * nsubmissions)))
        map_E_submissions[:, :, ifield] = np.reshape(
            map_E_subs_field, (evaluate.NBINS_THETA, nsubmissions), order='F')
        #import pdb; pdb.set_trace()

    # Reshape and return
    # DEBUGGING
    if False:
        ret = np.reshape(
            map_E_submissions, (evaluate.NBINS_THETA * evaluate.NFIELDS, nsubmissions), order='F')
    else:
        ret = np.reshape(
            np.transpose(map_E_submissions, (0, 2, 1)),
            (evaluate.NBINS_THETA * evaluate.NFIELDS, nsubmissions), order='F')
    #import pdb; pdb.set_trace()
    return ret


def write_submission(map_E_submission, outfile="test_submission.asc"):
    """Make a fake ASCII submission given an input vector of map_E values, and save it to outfile
    """
    field = np.arange(evaluate.NBINS_THETA * evaluate.NFIELDS) / evaluate.NBINS_THETA
    theta = evaluate.EXPECTED_THETA
    map_B = np.zeros(evaluate.NBINS_THETA * evaluate.NFIELDS)
    maperr = np.zeros(evaluate.NBINS_THETA * evaluate.NFIELDS)
    # Save in ASCII format
    with open(outfile, "wb") as fout:
        fout.write("# field_index  theta [deg]  map_E  map_B  maperr\n")
        np.savetxt(
            fout, np.array((field, theta, map_E_submission, map_B, maperr)).T,
            fmt=" %2d %.18e %.18e %.18e %.18e")
    # Job done, return
    return


if __name__ == "__main__":

    import cPickle
    import tempfile
    import pyfits
    import g3metrics
    import calculate_variable_cholesky
    import test_evaluate
    # Note I noticed that the standard deviation of correlated Q_v values, which is the main thing
    # of interest that gets calculated by this script, is a strong function of what the *mean* Q_v
    # is per bin.  As these are all correlated, this changes quite a bit from run to run!
    # Therefore I am going to repeat the calculations below NRUN times to sample from a range of
    # these results and take the average std(Q_v) to get something more representative
    qclist = []
    qmlist = []
    for krun in range(NRUNS):

        # Dicts containing arrays for storing Q_v values versus m and c, for ground and space
        qc = {"ground": np.empty((NTEST, len(CVALS))), "space": np.empty((NTEST, len(CVALS)))}
        qm = {"ground": np.empty((NTEST, len(MVALS))), "space": np.empty((NTEST, len(MVALS)))}
        coutfile = os.path.join(
            "results", "tabulated_correlated_variable_Q_v_versus_c_norm"+str(krun)+".pkl")
        moutfile = os.path.join(
            "results", "tabulated_correlated_variable_Q_v_versus_m_norm"+str(krun)+".pkl")
        # Test to see if coutfile is already made, if not make it, if so load it
        if not os.path.isfile(coutfile):
            for obs_type in ("ground", "space"):

                print "Calculating Q_v values versus c for control-"+obs_type+\
                    "-constant data in GREAT3"
                print "Loading Cholesky decomposition matrix from "+\
                    calculate_variable_cholesky.CHOLESKY_OUTFILE[obs_type]
                cholesky = pyfits.getdata(calculate_variable_cholesky.CHOLESKY_OUTFILE[obs_type])
                print "RHO = "+str(RHO)
                # First we build the truth table
                print "Getting/generating truth tables for control-"+obs_type+"-constant"
                field_ref, theta_ref, map_E_ref, _, maperr_ref = \
                    evaluate.get_generate_variable_truth(
                        EXPERIMENT, obs_type, truth_dir=TRUTH_DIR,
                        mape_file_prefix=evaluate.MAPEOBS_FILE_PREFIX,
                        file_prefixes=("galaxy_catalog", "galaxy_catalog"),
                        suffixes=("_intrinsic", ""), make_plots=False)
                # Get the unitc term
                map_E_unitc = 2. * test_evaluate.make_unitc(
                    EXPERIMENT, obs_type, truth_dir=TRUTH_DIR)
                # Loop over c values
                for jc, cval in enumerate(CVALS):

                    # Build the submissions (includes inter-method and inter-bin correlations)
                    map_E_field_subs = make_multiple_variable_submissions(
                        NTEST, map_E_ref, map_E_unitc, cval, evaluate.MFID, cholesky)
                    # Loop over submissions evaluating metric
                    for itest in xrange(NTEST):

                        fdsub, subfile = tempfile.mkstemp(suffix=".dat")
                        os.close(fdsub)
                        write_submission(map_E_field_subs[:, itest], outfile=subfile)
                        qc[obs_type][itest, jc] = evaluate.q_variable(
                            subfile, EXPERIMENT, obs_type, truth_dir=TRUTH_DIR)
                        os.remove(subfile)

                    print "mean(Q_v), std(Q_v) = "+str(qc[obs_type][:, jc].mean())+", "+\
                        str(qc[obs_type][:, jc].std())+" for "+str(NTEST)+" sims (with c = "+\
                        str(cval)+", obs_type = "+str(obs_type)+")"
                    print

            print "Saving pickled Q_v versus c dict to "+coutfile
            with open(coutfile, "wb") as fout: cPickle.dump(qc, fout)
            print
        else:
            with open(coutfile, "rb") as fin: qc = cPickle.load(fin)
        # Test to see if moutfile is already made, if not make it, if so load it
        if not os.path.isfile(moutfile):
            for obs_type in ("ground", "space"):

                print "Calculating Q_v values versus m for control-"+obs_type+\
                    "-constant data in GREAT3"
                print "Loading Cholesky decomposition matrix from "+\
                    calculate_variable_cholesky.CHOLESKY_OUTFILE[obs_type]
                cholesky = pyfits.getdata(calculate_variable_cholesky.CHOLESKY_OUTFILE[obs_type])
                print "RHO = "+str(RHO)
                # First we build the truth table
                print "Getting/generating truth tables for control-"+obs_type+"-constant"
                field_ref, theta_ref, map_E_ref, _, maperr_ref = \
                    evaluate.get_generate_variable_truth(
                        EXPERIMENT, obs_type, truth_dir=TRUTH_DIR,
                        mape_file_prefix=evaluate.MAPEOBS_FILE_PREFIX,
                        file_prefixes=("galaxy_catalog", "galaxy_catalog"),
                        suffixes=("_intrinsic", ""), make_plots=False)
                # Get the unitc term
                map_E_unitc = 2. * test_evaluate.make_unitc(
                    EXPERIMENT, obs_type, truth_dir=TRUTH_DIR)
                # Loop over m values
                for jm, mval in enumerate(MVALS):

                    # Build the submissions (includes inter-method and inter-bin correlations)
                    map_E_field_subs = make_multiple_variable_submissions(
                        NTEST, map_E_ref, map_E_unitc, evaluate.CFID, mval, cholesky)
                    # Loop over submissions evaluating metric
                    for itest in xrange(NTEST):

                        fdsub, subfile = tempfile.mkstemp(suffix=".dat")
                        os.close(fdsub)
                        write_submission(map_E_field_subs[:, itest], outfile=subfile)
                        qm[obs_type][itest, jm] = evaluate.q_variable(
                            subfile, EXPERIMENT, obs_type, truth_dir=TRUTH_DIR)
                        os.remove(subfile)

                    print "mean(Q_v), std(Q_v) = "+str(qm[obs_type][:, jm].mean())+", "+\
                        str(qm[obs_type][:, jm].std())+" for "+str(NTEST)+" sims (with m = "+\
                        str(mval)+", obs_type = "+str(obs_type)+")"
                    print

            print "Saving pickled Q_v versus m dict to "+moutfile
            with open(moutfile, "wb") as fout: cPickle.dump(qm, fout)
            print
        else:
            with open(moutfile, "rb") as fin: qm = cPickle.load(fin)

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
