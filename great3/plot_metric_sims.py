#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

CVALS = (2.e-4, 2.e-3, 2.e-2)
MVALS = (2.e-3, 2.e-2, 2.e-1)

def plot_vs_m(data, cvals=CVALS, mvals=MVALS):
    """Return a plot of mean and std Q values versus m along x
    """
    fig = plt.figure()
    for ic in range(len(cvals)):
        mean_data = np.mean(data[:, ic, :], axis=0)
        std_data = np.std(data[:, ic, :], axis=0)
        plt.errorbar(
            np.asarray(mvals) * (1.05**ic), mean_data, yerr=std_data,
            label='Input c = %.2e'%cvals[ic])
    plt.xscale('log')
    plt.xlim(4.e-1, 1.e-3)
    plt.legend(loc=2)
    plt.ylabel(r'Q$_V$')
    plt.xlabel('Input m')
    return fig

def plot_vs_c(data, cvals=CVALS, mvals=MVALS):
    """Return a plot of mean and std Q values versus c along x
    """
    fig = plt.figure()
    for im in range(len(mvals)):
        mean_data = np.mean(data[:, :, im], axis=0)
        std_data = np.std(data[:, :, im], axis=0)
        plt.errorbar(
            np.asarray(cvals) * (1.05**im), mean_data, yerr=std_data,
            label='Input m = %.2e'%mvals[im])
    plt.xscale('log')
    plt.xlim(4.e-2, 1.e-4)
    plt.legend(loc=2)
    plt.ylabel(r'Q$_V$')
    plt.xlabel('Input c')
    return fig


if __name__ == "__main__":

    import sys
    if len(sys.argv) != 4:
        print "usage: ./plot_metric_sims.py INFILE CPLOTFILE MPLOTFILE"
    data = np.load(sys.argv[1])
    figc = plot_vs_c(data)
    figm = plot_vs_m(data) 
    figc.savefig(sys.argv[2])
    figm.savefig(sys.argv[3])
