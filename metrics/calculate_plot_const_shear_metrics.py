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
"""@file calculate_plot_const_shear_metrics.py

Simulate constant shear metrics for a range of biases, and plot.
"""
import os
import cPickle
import numpy as np
import g3metrics

NIMS = 200               # Number of images per set, always 200 for G10
NGALS_PER_IM = 10000     # In GREAT08/GREAT10 there were 10000 galaxies per image
TRUE_MIN = 0.01          # Range of true input |shear| for random selection from an annulus around
TRUE_MAX = 0.05          # the origin
NFIELDS = 10             # Don't necessarily need to have NIMS input shears, but rather NFIELDS
                         # where NFIELDS * NSUBFIELDS = NIMS

CFID = 2.e-4 # Fiducial, "target" m and c values
MFID = 2.e-3 #

# Plotting ranges of interest
CMIN = CFID
CMAX = 2.e-2
MMIN = MFID
MMAX = 2.e-1

NBINS = 5     # Number of bins to plot in the ranges above
NMONTE = 100  # Number of montecarlo realizations
NOISE_SIGMA = 0.05  # Noise due to pixel shot noist on a shear estimate, per galaxy

# Generate arrays of values
cvals = CMIN * (CMAX / CMIN)**(np.arange(NBINS) / float(NBINS - 1.)) # geometric series
mvals = MMIN * (MMAX / MMIN)**(np.arange(NBINS) / float(NBINS - 1.))
cgrid, mgrid = np.meshgrid(cvals, mvals) # 2D arrays covering full space if needed

# Create empty storage arrays
Q08 = np.empty((NBINS, NBINS, NMONTE))
QZ1 = np.empty((NBINS, NBINS, NMONTE))
QZ2 = np.empty((NBINS, NBINS, NMONTE))

# File for storing pickled output
OUTFILE = os.path.join('results', 'const_shear_metrics.pkl')

if __name__ == "__main__":

    # Load up the outfile for plotting, or calculate and save before plotting
    if os.path.isfile(OUTFILE):
        Q08, QZ1, QZ2 = cPickle.load(open(OUTFILE, 'rb'))

    else:
        for krepeat in range(NMONTE):

            # Generate the truth tables for this realization
            g1true, g2true = g3metrics.make_const_truth_uniform_annular(
                NFIELDS, NIMS, range_min=TRUE_MIN, range_max=TRUE_MAX)    
            # Loop over each c, m combination
            for i in range(NBINS):
                for j in range(NBINS):
                    # Make the submissions
                    g1sub, g2sub = g3metrics.make_submission_const_shear(
                        cvals[i], cvals[i], mvals[j], mvals[j], g1true, g2true,
                        ngals_per_im=NGALS_PER_IM, noise_sigma=NOISE_SIGMA)
                    # Calculate the metrics and store
                    Q08[i, j, krepeat] = g3metrics.metricQ08_const_shear(
                        g1sub, g2sub, g1true, g2true, nfields=NFIELDS)
                    QZ1[i, j, krepeat] = g3metrics.metricQZ1_const_shear(
                        g1sub, g2sub, g1true, g2true, cfid=CFID, mfid=MFID)[0]
                    QZ2[i, j, krepeat] = g3metrics.metricQZ2_const_shear(
                        g1sub, g2sub, g1true, g2true, cfid=CFID, mfid=MFID)[0]

            print "Calculated const shear metrics for "+str(krepeat + 1)+"/"+\
                str(NMONTE)+" realizations"

        # Save the results as a tuple of NumPy arrays
        cPickle.dump((Q08, QZ1, QZ2), open(OUTFILE, 'wb'))

    # Get basic statistics
    q08mean = np.mean(Q08, axis=2)
    qZ1mean = np.mean(QZ1, axis=2)
    qZ2mean = np.mean(QZ2, axis=2)
    q08stds = np.std(Q08, axis=2)
    qZ1stds = np.std(QZ1, axis=2)
    qZ2stds = np.std(QZ2, axis=2)
    q08errs = q08stds / np.sqrt(NMONTE)
    qZ1errs = qZ1stds / np.sqrt(NMONTE)
    qZ2errs = qZ2stds / np.sqrt(NMONTE)

    # Normalise QZs to to Q=1000 for fiducial case (should be approx OK, but might be biased)
    qZ1norm = 1000. / qZ1mean[0, 0]
    qZ2norm = 1000. / qZ2mean[0, 0]
    qZ1norm_err = qZ1norm * (qZ1errs[0, 0] / qZ1mean[0, 0])
    qZ2norm_err = qZ2norm * (qZ2errs[0, 0] / qZ2mean[0, 0])
    print ("Normalization factor for QZ1 = %.3f +/- %.3f" % (qZ1norm, qZ1norm_err))
    print ("Normalization factor for QZ2 = %.3f +/- %.3f" % (qZ2norm, qZ2norm_err))
    qZ1mean *= qZ1norm
    qZ2mean *= qZ2norm
    
    # Plot the Q metric and the standard deviation of the NMONTE results, versus c first
    import matplotlib.pyplot as plt
    # QZ1
    plt.clf()
    plt.axhline(ls='--', color='k')
    for i in range(len(mvals)):
        plt.errorbar(
            cvals * (1.00**i), qZ1mean[:, i], yerr=qZ1stds[:, i], label='m = %.0e' % mvals[i])
    plt.xscale('log')
    plt.xlim(3.e-2, 2.e-5)
    plt.legend()
    plt.title('Mean QZ1 metric and standard deviation of '+str(NMONTE)+' realizations')
    plt.ylabel('Q')
    plt.xlabel('c')
    plt.savefig(os.path.join('plots', 'const_QZ1_vs_c_N'+str(NMONTE)+'.png'))

    # QZ2
    plt.clf()
    plt.axhline(ls='--', color='k')
    for i in range(len(mvals)):
        plt.errorbar(
            cvals * (1.00**i), qZ2mean[:, i], yerr=qZ2stds[:, i], label='m = %.0e' % mvals[i])
    plt.xscale('log')
    plt.xlim(3.e-2, 2.e-5)
    plt.legend()
    plt.title('Mean QZ2 metric and standard deviation of '+str(NMONTE)+' realizations')
    plt.ylabel('Q')
    plt.xlabel('c')
    plt.savefig(os.path.join('plots', 'const_QZ2_vs_c_N'+str(NMONTE)+'.png'))

    # Plot the Q metric and the standard deviation of the NMONTE results, versus m
    # QZ1
    plt.clf()
    plt.axhline(ls='--', color='k')
    for i in range(len(cvals)):
        plt.errorbar(
            mvals * (1.00**i), qZ1mean[i, :], yerr=qZ1stds[i, :], label='c = %.0e' % cvals[i])
    plt.xscale('log')
    plt.xlim(3.e-1, 2.e-4)
    plt.legend()
    plt.title('Mean QZ1 metric and standard deviation of '+str(NMONTE)+' realizations')
    plt.ylabel('Q')
    plt.xlabel('m')
    plt.savefig(os.path.join('plots', 'const_QZ1_vs_m_N'+str(NMONTE)+'.png'))

    # QZ2
    plt.clf()
    plt.axhline(ls='--', color='k')
    for i in range(len(cvals)):
        plt.errorbar(
            mvals * (1.00**i), qZ2mean[i, :], yerr=qZ2stds[i, :], label='c = %.0e' % cvals[i])
    plt.xscale('log')
    plt.xlim(3.e-1, 2.e-4)
    plt.legend()
    plt.title('Mean QZ2 metric and standard deviation of '+str(NMONTE)+' realizations')
    plt.ylabel('Q')
    plt.xlabel('m')
    plt.savefig(os.path.join('plots', 'const_QZ2_vs_m_N'+str(NMONTE)+'.png'))

    print ""
    print "Table of constant shear Q at constant c = cfid = "+str(CFID)
    print "    m        Q   "
    for m, Q in zip(mvals, qZ1mean[0, :]):
        print "{:8f} {:8.3f}".format(m, Q)

    print ""
    print "Table of constant shear Q at constant m = mfid = "+str(MFID)
    print "    c        Q   "
    for c, Q in zip(cvals, qZ1mean[:, 0]):
        print "{:8f} {:8.3f}".format(c, Q)

