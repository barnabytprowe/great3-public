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
"""@file plot_variable_submission.py

Handy command line executable script for plotting up a GREAT3 variable shear submission.
"""

# Constants
NFIELDS = 10
NBINS_THETA = 15
YLIM_EMODE = 2.e-5
YLIM_BMODE = 2.e-5


def plot(submission_filename, output_filename, nfields=NFIELDS, nbins_theta=NBINS_THETA,
         ylim_emode=YLIM_EMODE, ylim_bmode=YLIM_BMODE):
    """Plot a submission.
    """ 
    import numpy as np
    import matplotlib.pyplot as plt

    # Load the data from the input submission
    data = np.loadtxt(submission_filename)
    field, theta, map_E, map_B, maperr = (
        data[:, 0].astype(int), data[:, 1], data[:, 2], data[:, 3], data[:, 4])

    # Then plot (largely borrowed from the code in server/great3/evaluate.py)
    plt.figure(figsize=(10, 8))
    plt.subplot(211)
    for ifield in range(nfields):
        plt.semilogx(
            theta[ifield * nbins_theta: (ifield + 1) * nbins_theta],
            map_E[ifield * nbins_theta: (ifield + 1) * nbins_theta], label="Field "+str(ifield))

    plt.ylim(-ylim_emode, ylim_emode)
    plt.title(submission_filename+" E-mode")
    plt.ylabel("Ap. Mass Dispersion")
    plt.axhline(ls="--", color="k")
    plt.legend()
    plt.subplot(212)
    for ifield in range(nfields):
        plt.semilogx(
            theta[ifield * nbins_theta: (ifield + 1) * nbins_theta],
            map_B[ifield * nbins_theta: (ifield + 1) * nbins_theta], label="Field "+str(ifield))

    plt.ylim(-ylim_bmode, ylim_bmode)
    plt.title(submission_filename+" B-mode")
    plt.xlabel("Theta [degrees]")
    plt.ylabel("Ap. Mass Dispersion")
    plt.axhline(ls="--", color="k")
    plt.legend()
    plt.savefig(output_filename)
    return


if __name__ == "__main__":

    import sys
    # Get the input and output filenames from the command line
    if len(sys.argv) != 3:
        print "plot_variable_submission.py"
        print "usage: ./plot_variable_submission.py input_submission output_filename"
        sys.exit(1)
    submission_filename = sys.argv[1]
    output_filename = sys.argv[2]
    plot(submission_filename, output_filename)
 
