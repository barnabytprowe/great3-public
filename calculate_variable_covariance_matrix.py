"""@file calculate_variable_submission_covariance_matrix.py

For both ground and space, calculate the variance per bin of the aperture mass dispersion statistic
around the truth values, and the inter-bin covariances.
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
import test_evaluate

NTEST = 1000 # Should do AT LEAST 1000
NGALS_PER_IMAGE = 10000
NOISE_SIGMA = {"ground": 0.15, "space": 0.10}
TRUTH_DIR = "/Users/browe/great3/beta/truth" # Modify to wherever truth is unpacked
EXPERIMENT = "control" # Use the control truth values

def decompose_cov(cov, dumb=False):
    """Take in a covariance matrix, and return the tuple `(diag, correlation_matrix)`.

    I think this is a much better way to visualize covariance matrices, by splitting them into a
    dimensionless matrix specifying correlations and a dimensional vector of diagonal variances.

    What I call the `correlation_matrix` above is simply the matrix of product moment correlation
    coefficients (PMCCs) relating to an input covariance matrix `cov`, defined as:

        correlation_matrix[i, j] = cov[i, j] / math.sqrt(cov[i, i] * cov[j, j])

    i.e. just a square matrix with unity on the diagonals and the relevant PMCC on off-diagonals.

    To preserve the full information content of `cov` this code also returns

        diag = numpy.diag(cov)

    These are diagonal variance values of the input covariance matrix.  `diag` together with
    `correlation_matrix` fully specifies `cov`.

    If the kwarg `dumb` is set then the code does this calculation the dumbest and least efficient
    way possible in Python, an optional behaviour used for checking results.
    """
    if cov.shape[0] != cov.shape[1]:
        raise ValueError("Input covariance matrix must be square!")
    if not dumb:
        raise NotImplementedError("Sorry decompose_cov only does dumb at the moment!")
    else:
        dim = cov.shape[0]
        diag = np.diag(cov)
        sqrt_diag = np.sqrt(diag)
        correlation_matrix = np.identity(dim)
        for i in xrange(cov.shape[0]):

            for j in xrange(cov.shape[1]):

                if i == j:
                    pass
                else:
                    correlation_matrix[i, j] = cov[i, j] / (sqrt_diag[i] * sqrt_diag[j])

    # Return diag values and correlation matrix
    return diag, correlation_matrix

if __name__ == "__main__":

    import cPickle
    # Dicts containing arrays for mapE values for ground and space
    mapE = {
        "ground": np.empty((evaluate.NBINS_THETA * evaluate.NFIELDS, NTEST)),
        "space":  np.empty((evaluate.NBINS_THETA * evaluate.NFIELDS, NTEST))}
    mapEoutfile = os.path.join("results", "tabulated_map_E_N"+str(NTEST)+".pkl")
    # Start the loop over ground/space
    if not os.path.isfile(mapEoutfile):
        for obs_type in ("ground", "space",):

            print "Calculating map_E values for control-"+obs_type+"-variable data in GREAT3"
            print "NOISE_SIGMA = "+str(NOISE_SIGMA[obs_type])
            # First we build the truth table
            print "Building truth tables for control-"+obs_type+"-variable"
            # Get the x,y, true intrinsic ellips and shears for making fake submissions
            _, x, y, g1true, g2true = test_evaluate.get_variable_gtrue(
                EXPERIMENT, obs_type, truth_dir=TRUTH_DIR)
            _, _, _, g1int, g2int = test_evaluate.get_variable_gsuffix(
                EXPERIMENT, obs_type, truth_dir=TRUTH_DIR)
            # Then loop over specified number of tests         
            for itest in xrange(NTEST):

                # Build the submission
                fdsub, subfile = tempfile.mkstemp(suffix=".dat")
                result = test_evaluate.make_variable_submission(
                    x, y, g1true, g2true, g1int, g2int, 0., 0., 0., 0., outfile=subfile,
                    noise_sigma=NOISE_SIGMA[obs_type])
                os.close(fdsub)
                data = np.loadtxt(subfile)
                mapE[obs_type][:, itest] = data[:, 2]
                print "Completed test %4d / %4d" % (itest+1, NTEST)
                os.remove(subfile)
        
        print "Saving pickled map_E tables to "+mapEoutfile
        with open(mapEoutfile, "wb") as fout:
            cPickle.dump(mapE, fout)
        print
    else:
        with open(mapEoutfile, "rb") as fin:
            mapE = cPickle.load(fin)

    # Then calculate the covariance matrix, per field
    cov = {
        "ground": np.empty((evaluate.NBINS_THETA, evaluate.NBINS_THETA, evaluate.NFIELDS)),
        "space":  np.empty((evaluate.NBINS_THETA, evaluate.NBINS_THETA, evaluate.NFIELDS))}
    # Dicts for storing the mean (across fields) covariance, diagonal variances and correlation
    # matrix
    cov_mean = {}
    diag_mean = {}
    corr_mean = {}
    for obs_type in ("ground", "space"):

        # Calculate each field's covariance matrix
        for ifield in range(evaluate.NFIELDS):

            mapE_field = mapE[obs_type][
                ifield * evaluate.NBINS_THETA: (ifield + 1) * evaluate.NBINS_THETA, :]
            cov[obs_type][:, :, ifield] = np.cov(mapE_field)
            # Inspect field to field variation in covariance matrix
            #import matplotlib.pyplot as plt
            #plt.clf()
            #plt.imshow(np.log10(np.abs(cov[obs_type][:, :, ifield])), interpolation="nearest")
            #plt.colorbar()
            #plt.show()
            # Results: no appearance of strong systematic differences between fields, using mean
            # for sims is fine I think...


        # Plot the mean covariance matrix across all fields, log10(abs(C)), the correlations and the
        # diagonal errors
        cov_mean[obs_type] = np.mean(cov[obs_type], axis=2)
        diag_mean[obs_type], corr_mean[obs_type] = decompose_cov(cov_mean[obs_type], dumb=True)

        print "Plotting results for "+obs_type
        import matplotlib.pyplot as plt
        plt.clf()
        plt.imshow(np.log10(np.abs(np.mean(cov[obs_type], axis=2))), interpolation="nearest")
        plt.colorbar()
        plt.title(
            r"log$_{10}$|<C$_{ij}>$| for "+obs_type+" sims (NTEST="+str(NTEST)+", all fields)")
        plt.savefig(os.path.join("plots", "cov_mean_"+obs_type+"_N"+str(NTEST)+".png"))

        plt.clf()
        plt.imshow(corr_mean[obs_type], interpolation="nearest")
        plt.colorbar()
        plt.title(
            r"$\rho_{ij}$ for "+obs_type+" sims (NTEST="+str(NTEST)+", all fields)")
        plt.savefig(os.path.join("plots", "corr_mean_"+obs_type+"_N"+str(NTEST)+".png"))

        plt.clf()
        plt.loglog(evaluate.EXPECTED_THETA[:evaluate.NBINS_THETA], np.sqrt(diag_mean[obs_type]),
            label="NTEST = "+str(NTEST))
        plt.xlabel(r"$\theta$ [degrees]", fontsize="large")
        plt.ylabel(r"$\sigma_{M(\theta_i)}$", fontsize="x-large")
        plt.title(r"Uncertainties (diagonal only) on $M(\theta_i)$ for "+obs_type+" data")
        plt.legend()
        plt.savefig(os.path.join("plots", "diag_cov_mean_"+obs_type+"_N"+str(NTEST)+".png"))







