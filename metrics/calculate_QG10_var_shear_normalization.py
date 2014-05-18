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
"""@file calculate_QG10_var_shear_normalization.py

Script for calculating G10-style metric results for a range of biases (normalization comes from the
ensemble mean value at fiducial bias).
"""

import sys
import os
import cPickle
import matplotlib.pyplot as plt
import numpy as np
import galsim
import g3metrics

NIMS = int(200)
NGRID = 100        # Each image contains a grid of NGRID x NGRID galaxies
DX_GRID = 0.1      # Grid spacing (must be in degrees)
NOISE_SIGMA = 0.05 # Expected noise on each shear after shape noise pushed largely into B-mode
NBINS=8            # Number of bins of PS for metric
Q10_SCALING = 1

CTEST = [1.e-4, 3.e-4, 1.e-3, 3.e-3, 1.e-2]
MTEST = [1.e-3, 3.e-3, 1.e-2, 3.e-2, 1.e-1]
NREPEAT = 1000

#GALSIM_DIR=os.path.join("/Path", "To", "Your", "Repo")
GALSIM_DIR=os.path.join("/Users", "browe", "great3", "galsim")

OUTFILE = os.path.join(
    'results', 'normalizationv3_G10_QuadPS_N'+str(NREPEAT)+'_noise_sigma'+str(NOISE_SIGMA)+'.pkl')

if __name__ == "__main__":

    reference_ps = g3metrics.read_ps(galsim_dir=GALSIM_DIR)

    # Make the truth catalogues (a list of 2D, NGRIDxNGRID numpy arrays), reusing the reference_ps
    # each time for simplicity
    g1true_list, g2true_list = g3metrics.make_var_truth_catalogs(
        1, NIMS, [reference_ps,], ngrid=NGRID, grid_units=galsim.degrees)
    # Define some empty storage arrays
    qG10unnorm = np.empty((len(CTEST), len(MTEST), NREPEAT))
    qQuadPSunnorm = np.empty((len(CTEST), len(MTEST), NREPEAT))

    # TEMP: Make all the PS realizations (of the truth) the same to see if this alters noise props
    g1true_list = [g1true_list[0],] * len(g1true_list)
    g2true_list = [g2true_list[0],] * len(g2true_list)
    
    # Then generate submissions, and truth submissions
    for c, i in zip(CTEST, range(len(CTEST))):

        print "Calculating PS metrics with c_i = "+str(c)
        for m, j in zip(MTEST, range(len(MTEST))):

            print "Calculating PS metrics with m_i = "+str(m)
            for krepeat in range(NREPEAT):

                # Make a fake submission
                ksub, pEsubs, pBsubs, pEtrues, pBtrues = g3metrics.make_submission_var_shear(
                    c1=c, c2=c, m1=m, m2=m, g1true_list=g1true_list, g2true_list=g2true_list,
                    noise_sigma=NOISE_SIGMA, dx_grid=DX_GRID, nbins=NBINS)
                # Calculate the G10 metric for this realization
                qG10, mean_pEestG10, mean_pEtrueG10, mean_diffG10 = g3metrics.metricG10_var_shear(
                    ksub, pEsubs, [NOISE_SIGMA**2] * NIMS, pEtrues, scaling=Q10_SCALING,
                    dx_grid=DX_GRID)
                # Calculate the QuadPS metric for this realization
                qQuadPS, mean_pEestQuadPS, mean_pEtrueQuadPS, mean_diffQuadPS = \
                    g3metrics.metricQuadPS_var_shear(
                        ksub, pEsubs, [NOISE_SIGMA**2] * NIMS, pEtrues, scaling=Q10_SCALING,
                        dx_grid=DX_GRID)
                # Save the results
                qG10unnorm[i, j, krepeat] = qG10
                qQuadPSunnorm[i, j, krepeat] = qQuadPS
                # Plot something to stdout so the user knows progress is being made
                sys.stdout.write('.')
                sys.stdout.flush()

            sys.stdout.write('\n')

    print ""
    print "Writing results to "+OUTFILE
    cPickle.dump((qG10unnorm, qQuadPSunnorm), open(OUTFILE, 'wb'))
    print ""

