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
"""@file calculate_QG3S_var_shear_normalization.py

Script for calculating G3S (both quadratic and absolute normed) metric results for a range of biases
(normalization comes from the ensemble mean value at fiducial bias).
"""

import sys
import os
import cPickle
import matplotlib.pyplot as plt
import numpy as np
import galsim
import g3metrics


NIMS = 8
NGRID = 500         # Each image contains a grid of NGRID x NGRID galaxies
DX_GRID = 0.02      # Grid spacing (must be in degrees)
NOISE_SIGMA = 0.05  # Expected noise on each shear after shape noise pushed largely into B-mode
NBINS_ANGULAR = 15  # Number of angular bins for correlation function metric
MIN_SEP = DX_GRID
MAX_SEP = 10.

NFIELDS = 8      # Don't necessarily need to have NIMS input shear fields. Easiest if
                 # NFIELDS is an integral fraction of NIMS..

CFID = 2.e-4 # Fiducial, "target" m and c values
MFID = 2.e-3 #

PLOT = False # Plot while calculating?
# Plotting ranges of interest
CMIN = CFID
CMAX = 2.e-2
MMIN = MFID
MMAX = 2.e-1
NBINS_TEST = 5      # Number of bins to plot in the ranges above
NMONTE = 100        # Number of montecarlo samples

SCALEPS = 2.        # Factor by which to multiply the PS in order to increase m sensitivity

# Generate arrays of values for test values of c and m
CVALS = CMIN * (CMAX / CMIN)**(np.arange(NBINS_TEST) / float(NBINS_TEST - 1)) # geo series
MVALS = MMIN * (MMAX / MMIN)**(np.arange(NBINS_TEST) / float(NBINS_TEST - 1))
CGRID, MGRID = np.meshgrid(CVALS, MVALS) # 2D arrays covering full space

#GALSIM_DIR=os.path.join("/Path", "To", "Your", "Repo")
GALSIM_DIR=os.path.join("/Users", "browe", "great3", "GalSim")
#GALSIM_DIR=os.path.join("/home", "browe", "great3", "64", "GalSim")

OUTFILE = os.path.join(
    'results', 'normalization_G3S_N'+str(NMONTE)+'_noise_sigma'+str(NOISE_SIGMA)+'.pkl')

if __name__ == "__main__":

    # Load up the reference PS
    reference_ps = g3metrics.read_ps(galsim_dir=GALSIM_DIR, scale=SCALEPS)
    # Define some empty storage arrays
    qG3S_AMD_unnorm = np.empty((len(CVALS), len(MVALS), NMONTE))
    qG3S_QMD_unnorm = np.empty((len(CVALS), len(MVALS), NMONTE))

    # Then loop
    for krepeat in range(NMONTE):

        # Make the truth catalogues (a list of 2D, NGRIDxNGRID numpy arrays), reusing the reference
        # ps each time for simplicity
        print "Generating truth tables for Monte Carlo realization = "+str(krepeat + 1)+"/"+\
            str(NMONTE)
        g1true_list, g2true_list = g3metrics.make_var_truth_catalogs(
            NFIELDS, NIMS, [reference_ps,] * NFIELDS, ngrid=NGRID, dx_grid=DX_GRID,
            grid_units=galsim.degrees)

        # Then generate submissions, and truth submissions
        for c, i in zip(CVALS, range(NBINS_TEST)):
            print "Calculating CF metrics with c_i = "+str(c)
            for m, j in zip(MVALS, range(NBINS_TEST)):
                print "Calculating CF metrics with m_i = "+str(m)
                # Make a fake submission
                theta, mapEsubs, mapBsubs, maperrsubs, mapEtrues, mapBtrues = \
                    g3metrics.make_submission_var_shear_CF(
                        c1=c, c2=c, m1=m, m2=m, g1true_list=g1true_list,
                        g2true_list=g2true_list, noise_sigma=NOISE_SIGMA, dx_grid=DX_GRID,
                        nbins=NBINS_ANGULAR, min_sep=MIN_SEP, max_sep=MAX_SEP)
                # Calculate the metrics
                qG3S_AMD_unnorm[i, j, krepeat] = g3metrics.metricG3S_AMD(
                    mapEsubs, maperrsubs, mapEtrues, mapBtrues, NFIELDS, normalization=1.)
                qG3S_QMD_unnorm[i, j, krepeat] = g3metrics.metricG3S_QMD(
                    mapEsubs, maperrsubs, mapEtrues, mapBtrues, NFIELDS, normalization=1.)

        print str(__file__)+": Completed "+str(krepeat + 1)+"/"+str(NMONTE)+\
                    " Monte Carlo realizations"
        sys.stdout.write('\n')

    print ""
    print "Writing results to "+OUTFILE
    cPickle.dump((qG3S_AMD_unnorm, qG3S_QMD_unnorm), open(OUTFILE, 'wb'))
    print ""
