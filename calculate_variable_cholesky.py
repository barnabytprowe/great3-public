"""@file calculate_variable_cholesky.py

For both ground and space, construct a large covariance matrix for NTEST submissions and then
calculate and save (for future use) its Cholesky decomposition.
"""

import sys
import os
import numpy as np
path, module = os.path.split(__file__)
sys.path.append(os.path.join(path, "..", "server", "great3")) # Appends the folder
                                                              # great3-private/server/great3 to
                                                              # sys.path
import calculate_variable_covariance_matrix
import tabulate_variable_shear_metric_rev1

# The number of submissions to use: should match the tabulated (non-correlated results) for
# consistency
NTEST = tabulate_variable_shear_metric_rev1.NTEST
RHO = 0.45

def construct_full_covariance(cov, ntest, rho=0.):
    """Take an input covariance, and make an ntest * cov.shape dimensioned-covariance matrix.
    """
    if cov.shape[0] != cov.shape[1]:
        raise ValueError("Input covariance matrix must be square")
    print "Building full "+str(cov.shape[0] * ntest)+" x "+str(cov.shape[1] * ntest)+\
        " covariance matrix"
    full_cov = np.empty((cov.shape[0] * ntest, cov.shape[1] * ntest))
    for i in xrange(ntest):

        for j in xrange(ntest):

            if i == j:
                full_cov[
                    i * cov.shape[0]: (i + 1) * cov.shape[0],
                    j * cov.shape[1]: (j + 1) * cov.shape[1]] = cov
            else:
                full_cov[
                    i * cov.shape[0]: (i + 1) * cov.shape[0],
                    j * cov.shape[1]: (j + 1) * cov.shape[1]] = rho * cov
    # Return
    return full_cov

def get_full_covs(ntest, rho,
                  inground=calculate_variable_covariance_matrix.COV_MEAN_OUTFILE["ground"],
                  inspace=calculate_variable_covariance_matrix.COV_MEAN_OUTFILE["space"],
                  outground=os.path.join("results", "cholesky_ground.npy"),
                  outspace=os.path.join("results", "cholesky_space.npy")):
    """Get and save the full size Cholesky matrices for both ground and space.
    """
    # First do ground
    print "Loading covariance matrix from "+inground+"..."
    cov_ground = np.load(inground)
    print "Constructing full covariance matrix..."
    full_cov_ground = construct_full_covariance(cov_ground, ntest, rho=rho)
    # DEBUG: Write out for checking
    #import pyfits
    #pyfits.writeto("full_cov_ground.fits", full_cov_ground, clobber=True)
    print "Performing Cholesky decomposition..."
    cholesky_ground = np.linalg.cholesky(full_cov_ground)
    np.save(outground, cholesky_ground)
    print "Saved Cholesky decomposition to "+outground
    # Then do space
    print "Loading covariance matrix from "+inspace+"..."
    cov_space = np.load(inspace)
    print "Constructing full covariance matrix..."
    full_cov_space = construct_full_covariance(cov_space, ntest, rho=rho)
    print "Performing Cholesky decomposition..."
    cholesky_space = np.linalg.cholesky(full_cov_space)
    np.save(outspace, cholesky_space)
    print "Saved Cholesky decomposition to "+outspace
    return


if __name__ == "__main__":

    get_full_covs(NTEST, RHO)