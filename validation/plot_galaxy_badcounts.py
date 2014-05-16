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
"""@file plot_galaxy_badcounts.py

Convenience script for finding and plotting the histogram of object SNRs for the GREAT3 galaxy image
with an object number count furthest from the 10000 required.  Uses a "badcount" dict as output by 
`count_objects.py`
"""

# usage: ./plot_galaxy_badcounts.py DICTFILE
import os
from sys import argv
import cPickle
import pyfits
import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":

    if len(argv) != 2:
        print "usage: ./plot_galaxy_badcounts.py DICTFILE"
        import sys
        sys.exit(1)
    # List for storing object counts 
    nobs = []
    # Setup a couple of variables to get the largest nobs
    biggest_nobs = 10000
    biggest_nobs_field = ""
    baddict = cPickle.load(open(argv[1]))
    for key, subdict in baddict.iteritems():

        if "starfield_image" not in os.path.split(key)[-1]:
            nobs.append(subdict["nobs"])
            if np.abs(nobs[-1] - 10000) > np.abs(biggest_nobs - 10000):
                biggest_nobs = nobs[-1]
                biggest_nobs_field = key

    # Then plot the hist
    plt.hist(
        nobs, bins=50,
        label="Worst offender: "+os.path.split(os.path.split(biggest_nobs_field)[0])[-1]+"/"+
        os.path.split(biggest_nobs_field)[-1]+" with "+str(biggest_nobs))
    plt.xlabel("N objects")
    plt.legend()
    splitname = ((argv[1]).rsplit("_"))
    plt.title("Fields in "+str(splitname[1])+" with N objects != 10000")
    outfile = "histogram_"+((argv[1]).rstrip(".p"))+".png"
    print "Saving plot to "+outfile
    plt.savefig(outfile)
    plt.show()
     

