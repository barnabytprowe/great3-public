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
"""@file calculate_variable_cholesky.py

For both ground and space, construct a large covariance matrix for NTEST submissions and then
calculate and save (for future use) its Cholesky decomposition.  Used in
https://github.com/barnabytprowe/great3-public/wiki/The-revised-GREAT3-metrics,-February-2014
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

CHOLESKY_OUTFILE = {
    "ground": os.path.join("results", "cholesky_ground.fits"),
    "space": os.path.join("results", "cholesky_space.fits")}


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

def get_full_covs(ntest=NTEST, rho=RHO,
                  inground=calculate_variable_covariance_matrix.COV_MEAN_OUTFILE["ground"],
                  inspace=calculate_variable_covariance_matrix.COV_MEAN_OUTFILE["space"],
                  outground=CHOLESKY_OUTFILE["ground"], outspace=CHOLESKY_OUTFILE["space"]):
    """Get and save the full size Cholesky matrices for both ground and space.
    """
    import pyfits
    # First do ground
    print "Loading covariance matrix from "+inground+"..."
    cov_ground = np.load(inground)
    print "Constructing full covariance matrix..."
    full_cov_ground = construct_full_covariance(cov_ground, ntest, rho=rho)
    # DEBUG: Write out for checking
    #pyfits.writeto(inground+".full.fits", full_cov_ground, clobber=True)
    print "Performing Cholesky decomposition..."
    cholesky_ground = np.linalg.cholesky(full_cov_ground)
    pyfits.writeto(outground, cholesky_ground)
    print "Saved Cholesky decomposition to "+outground
    # Then do space
    print "Loading covariance matrix from "+inspace+"..."
    cov_space = np.load(inspace)
    print "Constructing full covariance matrix..."
    full_cov_space = construct_full_covariance(cov_space, ntest, rho=rho)
    # DEBUG: Write out for checking
    #pyfits.writeto(inspace+".full.fits", full_cov_space, clobber=True)
    print "Performing Cholesky decomposition..."
    cholesky_space = np.linalg.cholesky(full_cov_space)
    pyfits.writeto(outspace, cholesky_space)
    print "Saved Cholesky decomposition to "+outspace
    return


if __name__ == "__main__":

    get_full_covs(ntest=NTEST, rho=RHO)
