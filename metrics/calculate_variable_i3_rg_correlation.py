#!/usr/bin/env python

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
"""@file calculate_variable_i3_rg_correlation.py
Script used to calculate the product moment correlation coeffient for variable shear submissions
from im3shape and the great3_ec_regauss method.

See https://github.com/barnabytprowe/great3-public/wiki/The-revised-GREAT3-metrics,-February-2014
"""
import os
import sys
import numpy as np
import numpy.linalg as linalg
sys.path.append(os.path.join("..", "server", "great3"))
import evaluate
import test_evaluate

EXPERIMENT = "control"
OBS_TYPE = "ground"
SHEAR_TYPE = "variable"

# Paths to control-ground-variable filenames (im3shape score ~XXX, regauss score ~YYY; TODO: update
# these!), make sure these are present!
I3FILE = os.path.join("results", "im3shape-great3-beta_control-ground-variable.asc")
RGFILE = os.path.join(
    "results", "great3_ec-regauss-example-example_1_cgv-2013-12-27T15:32:20.450921+00:00.g3")

TRUTH_DIR = "/Users/browe/great3/beta/truth" # Modify to wherever truth is unpacked
CORR2_EXEC = "corr2" # Mike's corr2 executable
MAPESHEAR_FILE_PREFIX = "mapEshear_" # \
MAPEINT_FILE_PREFIX = "mapEint_"     #  } Prefixes for the different types of variable truth cats
MAPEOBS_FILE_PREFIX = "mapEobs_"     # /

USE_ALL_BINS = True  # If set, does not simply stick to the evaluate.USEBINS bin selection for
                     # estimating the product moment correlation coefficient rho

if __name__ == "__main__":

    # Load up relevant data
    i3data = np.loadtxt(I3FILE)
    rgdata = np.loadtxt(RGFILE)
    # Extract the salient parts of the submissions from data
    field_i3 = i3data[:, 0].astype(int)
    theta_i3 = i3data[:, 1]
    map_E_i3 = i3data[:, 2]
    field_rg = rgdata[:, 0].astype(int)
    theta_rg = rgdata[:, 1]
    map_E_rg = rgdata[:, 2]
    # Load/generate the truth shear signal
    field_shear, theta_shear, map_E_shear, _, maperr_shear = evaluate.get_generate_variable_truth(
        EXPERIMENT, OBS_TYPE, truth_dir=TRUTH_DIR, logger=None, corr2_exec=CORR2_EXEC,
        mape_file_prefix=MAPESHEAR_FILE_PREFIX, suffixes=("",), make_plots=False)
    # Then generate the intrinsic only map_E, useful for examinging plots, including the maperr
    # (a good estimate of the relative Poisson errors per bin) which we will use to provide a weight
    field_int, theta_int, map_E_int, _, maperr_int = evaluate.get_generate_variable_truth(
        EXPERIMENT, OBS_TYPE, truth_dir=TRUTH_DIR, logger=None, corr2_exec=CORR2_EXEC,
        mape_file_prefix=MAPEINT_FILE_PREFIX, suffixes=("_intrinsic",), make_plots=False)
    # Then generate the theory observed = int + shear combined map signals - these are our reference
    # Note this uses the new functionality of get_generate_variable_truth for adding shears
    field_ref, theta_ref, map_E_ref, _, maperr_ref = evaluate.get_generate_variable_truth(
        EXPERIMENT, OBS_TYPE, truth_dir=TRUTH_DIR, logger=None, corr2_exec=CORR2_EXEC,
        mape_file_prefix=MAPEOBS_FILE_PREFIX, file_prefixes=("galaxy_catalog", "galaxy_catalog"),
        suffixes=("_intrinsic", ""), make_plots=False)

    # We can use this information to construct the linear model
    # c**2 + (1 + 2 * m + m**2) * <map_E_ref> = a + b * <map_E_ref>
    # ...and fit it to the submissions

    # Hmmm actually a linear model doesn't work because we do want to constrain c**2 to be
    # positive... Try using optimize instead

    if USE_ALL_BINS:
        # HACKY: Temporarily set evaluate.USEBINS to be all true to use *all* angular bins for
        # calculating the product moment correlation coeff
        evaluate.USEBINS = np.ones(evaluate.USEBINS.shape, dtype=bool)

    # Get truth as a full catalog
    _, x, y, g1true, g2true = test_evaluate.get_variable_gtrue(
        EXPERIMENT, OBS_TYPE, truth_dir=TRUTH_DIR)
    # Then make a variable submission to learn what pure c1 and c2 look like in map^2
    import tempfile
    fdtmp, tmpfile = tempfile.mkstemp(suffix=".dat")
    result = test_evaluate.make_variable_submission(
        x, y, np.zeros_like(g1true), np.zeros_like(g2true), np.zeros_like(g1true),
        np.zeros_like(g2true), 1., 0., 0., 0., outfile=tmpfile, noise_sigma=0.)
    os.close(fdtmp)
    map_E_c1 = np.loadtxt(tmpfile)[:, 2]
    os.remove(tmpfile)
    fdtmp, tmpfile = tempfile.mkstemp(suffix=".dat")
    result = test_evaluate.make_variable_submission(
        x, y, np.zeros_like(g1true), np.zeros_like(g2true), np.zeros_like(g1true),
        np.zeros_like(g2true), 0., 1., 0., 0., outfile=tmpfile, noise_sigma=0.)
    os.close(fdtmp)
    map_E_c2 = np.loadtxt(tmpfile)[:, 2]
    os.remove(tmpfile)
    # Then average these (I checked they are very similar as you would expect) to get our "unit c"
    # term for the modelling.  Note that previously I summed these, but I think this was an error as
    # that actually corresponds approximately to a shear field with g1 = g2 = 1, i.e. |g| = sqrt(2)
    map_E_unitc = .5 * (map_E_c1 + map_E_c2)

    # Then we quickly define the error normalized difference vector (essentially gets used as chisq
    # later by optimize).  Note this defines the model we are making of the submissions:
    #
    #  <map_E_model> = <map_E_unitc> * c^2 + <map_E_ref> * (1 + 2 * m + m^2)
    def map_diff_func(cm_array, mapEsub, maperrsub, mapEref, mapEunitc):
        """Difference of an m-c model of the aperture mass statistic and the submission.
        """
        retval = (
                mapEunitc * cm_array[0]**2 + mapEref * (1. + 2. * cm_array[1] + cm_array[1]**2)
                - mapEsub
            ) / maperrsub
        return retval

    # Use optimize to determmine the (non-linear) model parameters and covariance matrix using the
    # variance of the model residuals
    import scipy.optimize
    # First im3shape
    results_i3 = scipy.optimize.leastsq(
        map_diff_func, np.array([0., 0.]),
        args=(
            map_E_i3[evaluate.USEBINS],
            maperr_ref[evaluate.USEBINS],
            map_E_ref[evaluate.USEBINS], # Note we use the reference errors, this will appropriately
                                         # weight different bins and is not noisy
            map_E_unitc[evaluate.USEBINS]), full_output=True)
    ci3 = results_i3[0][0]
    mi3 = results_i3[0][1]
    map_E_i3_model = map_E_unitc * ci3**2 + map_E_ref * (1. + 2. * mi3 + mi3**2)
    residual_variance_i3 = np.var(((map_E_i3 - map_E_i3_model) / maperr_ref)[evaluate.USEBINS])
    covcm_i3 = results_i3[1] * residual_variance_i3
    sigc_i3 = np.sqrt(covcm_i3[0, 0])
    sigm_i3 = np.sqrt(covcm_i3[1, 1])
    # Then REGAUSS
    results_rg = scipy.optimize.leastsq(
        map_diff_func, np.array([0., 0.]),
        args=(
            map_E_rg[evaluate.USEBINS],
            maperr_ref[evaluate.USEBINS],
            map_E_ref[evaluate.USEBINS], # Note we use the reference errors, this will appropriately
                                         # weight different bins and is not noisy
            map_E_unitc[evaluate.USEBINS]), full_output=True)
    crg = results_rg[0][0]
    mrg = results_rg[0][1]
    map_E_rg_model = map_E_unitc * crg**2 + map_E_ref * (1. + 2. * mrg + mrg**2)
    residual_variance_rg = np.var(((map_E_rg - map_E_rg_model) / maperr_ref)[evaluate.USEBINS])
    covcm_rg = results_rg[1] * residual_variance_rg
    sigc_rg = np.sqrt(covcm_rg[0, 0])
    sigm_rg = np.sqrt(covcm_rg[1, 1])

    # Print a summary of the results
    print "c_i3 = %+.5f +/- %.5f" % (ci3, sigc_i3)
    print "c_rg = %+.5f +/- %.5f" % (crg, sigc_rg)
    print "m_i3 = %+.5f +/- %.5f" % (mi3, sigm_i3)
    print "m_rg = %+.5f +/- %.5f" % (mrg, sigm_rg)
    
    # Get the model corrected differences ()
    dEi3 = (map_E_i3 - map_E_i3_model)[evaluate.USEBINS]
    dErg = (map_E_rg - map_E_rg_model)[evaluate.USEBINS]
    
    # Number of bins we're going to use (those for which evaluate.USEBINS == True)
    ntouse = np.sum(evaluate.USEBINS) / evaluate.NFIELDS

    # Define cov_theta, covariance as a function of angular bin
    cov_theta = np.empty((2, 2, ntouse))
    rho_theta = np.empty(ntouse)
    for i in range(ntouse):

        cov_theta[:, :, i] = np.cov(dEi3[i::ntouse], dErg[i::ntouse])
        rho_theta[i] = cov_theta[1, 0, i] / np.sqrt(cov_theta[0, 0, i] * cov_theta[1, 1, i])
    
    # Make a plot of these results: use the overall standard deviation as the error estimate
    import matplotlib.pyplot as plt
    plt.clf()
    plt.errorbar(
        evaluate.EXPECTED_THETA[evaluate.USEBINS][:ntouse], rho_theta,
        yerr=np.ones(ntouse) * np.sqrt(1. / (evaluate.NFIELDS - 3.)), # Fisher approx derived from
                                                                      # Wikipedia page on PMCC
        fmt='+')
    plt.xscale("log")
    plt.axhline(ls='--', color='k')
    plt.xlabel(r"$\theta$ [deg]", size="large")
    plt.ylabel(r"Product moment correlation coefficient $\rho$", size="large")
    plotfile = "./plots/rho_vs_theta_variable"
    plt.ylim(-0.5, 1.2)
    if USE_ALL_BINS:
        plt.title("Per bin PMCC for all angular bins")
        plotfile += "_allbins"
    else:
        plt.title("Per bin PMCC for only angular bins used by metric")
        plotfile += "_usebins"
    plt.savefig(plotfile+".png")
    plt.show()

    # Calculate covariance matrix for all data
    cov = np.cov(dEi3, dErg)
    # Calculate correlation coefficient
    rho = cov[1, 0] / np.sqrt(cov[0, 0] * cov[1, 1])
    print "PMCC rho = %.3f +/- %.3f" % (rho, np.sqrt(1. / (evaluate.NFIELDS * ntouse - 3.)))
