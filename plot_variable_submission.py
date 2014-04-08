#!/usr/bin/env python

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
 
