#!/usr/bin/env python
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

NFIELDS = 10
NBINS_THETA = 15

if __name__ == "__main__":

    # Get the input and output filenames from the command line
    if len(sys.argv) != 3:
        print "plot_variable_submission.py"
        print "usage: ./plot_variable_submission.py input_submission output_filename"
        sys.exit(1)
    submission_filename = sys.argv[1]
    output_filename = sys.argv[2]

    # Load the data from the input submission
    data = np.loadtxt(submission_filename)
    field, theta, map_E, map_B, maperr = (
        data[:, 0].astype(int), data[:, 1], data[:, 2], data[:, 3], data[:, 4])

    # Then plot (largely borrowed from the code in server/great3/evaluate.py)
    plt.figure(figsize=(10, 8))
    plt.subplot(211)
    for ifield in range(NFIELDS):
        plt.semilogx(
            theta[ifield * NBINS_THETA: (ifield + 1) * NBINS_THETA],
            map_E[ifield * NBINS_THETA: (ifield + 1) * NBINS_THETA], label="Field "+str(ifield))

    plt.ylim(-2.e-5, 2.e-5)
    plt.title(submission_filename+" E-mode")
    plt.ylabel("Ap. Mass Dispersion")
    plt.axhline(ls="--", color="k")
    plt.legend()
    plt.subplot(212)
    for ifield in range(NFIELDS):
        plt.semilogx(
            theta[ifield * NBINS_THETA: (ifield + 1) * NBINS_THETA],
            map_B[ifield * NBINS_THETA: (ifield + 1) * NBINS_THETA], label="Field "+str(ifield))

    plt.ylim(-2.e-5, 2.e-5)
    plt.title(submission_filename+" B-mode")
    plt.xlabel("Theta [degrees]")
    plt.ylabel("Ap. Mass Dispersion")
    plt.axhline(ls="--", color="k")
    plt.legend()
    plt.savefig(output_filename) 
