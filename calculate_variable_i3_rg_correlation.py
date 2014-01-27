#!/usr/bin/env python
import os
import sys
import numpy as np
sys.path.append(os.path.join("..", "server", "great3"))
import evaluate

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

    # Now, previous experience tells me that, without taking inter-bin correlations into account, it
    # is difficult to get an unbiased simultaneous estimate of m and c from correlation function
    # results.  Given the current setup it *might* be worth retrying, as there have been many
    # changes since the last time.  But for now I am going to not correct, and see how we do...
    
    # Get the differences
    dEi3 = (map_E_i3 - map_E_ref)
    dErg = (map_E_rg - map_E_ref)
    
    # Define cov_theta
    cov_theta = np.empty(evaluate.NBINS_THETA)
    for i in range(evaluate.NBINS_THETA):
        
    # Calculate covariance matrices
    cov = np.cov(dEi3, dErg)
    #cov2 = np.cov(d2i3, d2rg)
    # Calculate correlation coefficient
    rho = cov[1, 0] / np.sqrt(cov[0, 0] * cov[1, 1])
    #rho2 = cov2[1, 0] / np.sqrt(cov2[0, 0] * cov2[1, 1])
    #print "Correlation coefficient (g1) rho = "+str(rho1) 
    #print "Correlation coefficient (g2) rho = "+str(rho2) 
