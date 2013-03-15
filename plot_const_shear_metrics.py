import os
import numpy as np
import g3metrics

NIMS = 300               # Number of images per set, always 200 for G10
NGALS_PER_IM = 10000     # In GREAT08/GREAT10 there were 10000 galaxies per image
TRUE_SIGMA = 0.04        # Standard deviation of true input shears for normal distribution
TRUE_RANGE = 0.08        # Range of true input shears for a uniform distribution
NTRUESETS = 50           # Don't necessarily need to have NIMS input shears. But easiest if
                         # NTRUESETS is an integral fraction of NIMS..

CFID = 2.e-4 # Fiducial, "target" m and c values
MFID = 2.e-3 #

# Plotting ranges of interest
CMIN = CFID
CMAX = 1.e-1
MMIN = MFID
MMAX = 1.e0

NBINS = 7 # Number of bins to plot in the ranges above
NMONTE = 30  # Number of montecarlo samples
NOISE_SIGMA = 0.10  # Noise due to pixel shot noist on a shear estimate, per galaxy

# Generate arrays of values
cvals = CMIN * (CMAX / CMIN)**(np.arange(NBINS) / float(NBINS - 1.)) # geometric series
mvals = MMIN * (MMAX / MMIN)**(np.arange(NBINS) / float(NBINS - 1.))

cgrid, mgrid = np.meshgrid(cvals, mvals) # 2D arrays covering full space

# Generate the truth tables
g1true, g2true = g3metrics.make_truth_normal_dist(NTRUESETS, NIMS, true_sigma=TRUE_SIGMA)

# Create empty storage arrays
Q08_vs_m = np.empty(NBINS)
QZ1_vs_m = np.empty(NBINS)
QZ2_vs_m = np.empty(NBINS)
Q08_vs_c = np.empty(NBINS)
QZ1_vs_c = np.empty(NBINS)
QZ2_vs_c = np.empty(NBINS)

import matplotlib.pyplot as plt

for j in range(NMONTE):

    # Loop over mvalues making indepented submissions at each c, m combination
    for i in range(NBINS):
        g1sub, g2sub = g3metrics.make_submission_const_shear(CFID, CFID, mvals[i], mvals[i],
                                                           g1true, g2true,
                                                           ngals_per_im=NGALS_PER_IM,
                                                           noise_sigma=NOISE_SIGMA)
        Q08_vs_m[i] = g3metrics.metricQ08_const_shear(g1sub, g2sub, g1true, g2true,
                                                    ntruesets=NTRUESETS)
        QZ1_vs_m[i] = g3metrics.metricQZ1_const_shear(g1sub, g2sub, g1true, g2true,
                                                    cfid=CFID, mfid=MFID)[0]
        QZ2_vs_m[i] = g3metrics.metricQZ2_const_shear(g1sub, g2sub, g1true, g2true,
                                                    cfid=CFID, mfid=MFID)[0]

    if j == 0:
        plt.loglog(mvals, QZ1_vs_m, 'r', label='QZ1')
        plt.loglog(mvals, QZ2_vs_m, 'g', label='QZ2')
        plt.loglog(mvals, Q08_vs_m, 'b', label='Q08')
    else:
        plt.loglog(mvals, QZ1_vs_m, 'r')
        plt.loglog(mvals, QZ2_vs_m, 'g')
        plt.loglog(mvals, Q08_vs_m, 'b')

plt.title('Q vs m: all c = '+str(CFID)+'\n '+
          'NIMS = '+str(NIMS)+', NOISE_SIGMA = '+str(NOISE_SIGMA)+', TRUE_SIGMA = '+str(TRUE_SIGMA))
plt.xlabel('m = m1 = m2')
plt.ylabel('Q')
plt.legend()
if not os.path.isdir('./plots'):
    os.mkdir('./plots')
outfile = './plots/Q_vs_m_const_shear_nims'+str(NIMS)+'_nsig'+str(NOISE_SIGMA)+\
    '_truesig'+str(TRUE_SIGMA)+'.png'
plt.savefig(outfile)

plt.figure()
for j in range(NMONTE):

    # Loop over mvalues making indepented submissions at each c, m combination
    for i in range(NBINS):
        g1sub, g2sub = g3metrics.make_submission_const_shear(cvals[i], cvals[i], MFID, MFID,
                                                           g1true, g2true,
                                                           ngals_per_im=NGALS_PER_IM,
                                                           noise_sigma=NOISE_SIGMA)
        Q08_vs_c[i] = g3metrics.metricQ08_const_shear(g1sub, g2sub, g1true, g2true,
                                                    ntruesets=NTRUESETS)
        QZ1_vs_c[i] = g3metrics.metricQZ1_const_shear(g1sub, g2sub, g1true, g2true,
                                                    cfid=CFID, mfid=MFID)[0]
        QZ2_vs_c[i] = g3metrics.metricQZ2_const_shear(g1sub, g2sub, g1true, g2true,
                                                    cfid=CFID, mfid=MFID)[0]

    if j == 0:
        plt.loglog(cvals, QZ1_vs_c, 'r', label='QZ1')
        plt.loglog(cvals, QZ2_vs_c, 'g', label='QZ2')
        plt.loglog(cvals, Q08_vs_c, 'b', label='Q08')
    else:
        plt.loglog(cvals, QZ1_vs_c, 'r')
        plt.loglog(cvals, QZ2_vs_c, 'g')
        plt.loglog(cvals, Q08_vs_c, 'b')

plt.title('Q vs c: all m = '+str(MFID)+'\n '+
          'NIMS = '+str(NIMS)+', NOISE_SIGMA = '+str(NOISE_SIGMA)+', TRUE_SIGMA = '+str(TRUE_SIGMA))
plt.xlabel('c = c1 = c2')
plt.ylabel('Q')
plt.legend()
outfile = './plots/Q_vs_c_const_shear_nims'+str(NIMS)+'_nsig'+str(NOISE_SIGMA)+'_truesig'+\
    str(TRUE_SIGMA)+'.png'
plt.savefig(outfile)
