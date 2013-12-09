#!/usr/bin/env python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import count_objects
sys.path.append("..")
import great3sims.mapper

if __name__ == "__main__":

    if len(sys.argv) == 1 or len(sys.argv) > 4:
        print "usage: ./plot_detect_galaxy_snr.py FITSFILE1 [... FITSFILE2] OUTFILE"
        sys.exit(1)
    infile1 = sys.argv[1]
    data1 = count_objects.file_count(infile1, sex_exec="sex")
    snr1 = data1["FLUX_AUTO"] / data1["FLUXERR_AUTO"]
    infilesplit1 = infile1.split("/")
    label1 = infilesplit1[-4]+"/"+infilesplit1[-3]+"/"+infilesplit1[-2]+"/"+infilesplit1[-1]
    plt.hist(snr1, bins=100, range=(0, 100), label=label1, histtype="step")
    if len(sys.argv) == 4:
        infile2 = sys.argv[2]
        outfile = sys.argv[3]
        data2 = count_objects.file_count(infile2, sex_exec="sex")
        snr2 = data2["FLUX_AUTO"] / data2["FLUXERR_AUTO"]
        infilesplit2 = infile2.split("/")
        label2 = infilesplit2[-4]+"/"+infilesplit2[-3]+"/"+infilesplit2[-2]+"/"+infilesplit2[-1]
        plt.hist(snr2, bins=100, range=(0, 100), label=label2, histtype="step")
    else:
        outfile = sys.argv[2]
    plt.legend()
    plt.xlabel("FLUX_AUTO / FLUXERR_AUTO")
    plt.ylabel("Counts")
    print "Saving plot to "+outfile
    plt.savefig(outfile)
    plt.show()
    
