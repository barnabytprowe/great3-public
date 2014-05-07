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
"""@file plot_detect_galaxy_snr.py

Quick and dirty command line tool for making plots of either a single FITS image file's SExtractor
SNR distribution, or two distributions for two separate FITS files, overlaid.  This latter
functionality was useful for validating that the multiepoch branch image SNRs were smaller than
the control etc. SNRs by an appropriate factor.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import count_objects
sys.path.append("..")
import great3sims.mapper

if __name__ == "__main__":

    if len(sys.argv) <= 2 or len(sys.argv) > 4:
        print "usage: ./plot_detect_galaxy_snr.py FITSFILE1 [... FITSFILE2] OUTFILE"
        sys.exit(1)
    infile1 = sys.argv[1]
    data1 = count_objects.file_count(infile1, sex_exec="sex")
    snr1 = data1["FLUX_AUTO"] / data1["FLUXERR_AUTO"]
    infilesplit1 = infile1.split("/")
    label1 = infilesplit1[-4]+"/"+infilesplit1[-3]+"/"+infilesplit1[-2]+"/"+infilesplit1[-1]
    plt.hist(snr1, bins=100, range=(0, 100), label=label1, histtype="step")
    #import itertools
    #for x, y, snr in itertools.izip(data1['X_IMAGE'], data1['Y_IMAGE'], snr1):
    #    if snr > 100.: print "%.4f %.4f %.4f" % (x, y, snr)
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
    
