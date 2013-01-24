
import numpy as np
import g3utils
import sim_const_shear_submission as sim_sub
import sim_const_shear_truth as sim_truth

CMIN = 1.e-4
CMAX = 1.e-2

MMIN = 1.e-3
MMAX = 1.e-1

NBINS = 5

NOISE_SIGMA = 0.05

# Generate arrays of values
cvals = CMIN * (CMAX / CMIN)**(np.arange(NBINS) / float(NBINS - 1.)) # geometric series
mvals = MMIN * (MMAX / MMIN)**(np.arange(NBINS) / float(NBINS - 1.))

cgrid, mgrid = np.meshgrid(cvals, mvals) # 2D arrays covering full space

# Generate the truth tables
g1true, g2true = sim_truth.make_truth_normal_dist()

Q08_vs_m = np.empty(NBINS)
QZ1_vs_m = np.empty(NBINS)
QZ2_vs_m = np.empty(NBINS)
Q08_vs_c = np.empty(NBINS)
QZ1_vs_c = np.empty(NBINS)
QZ2_vs_c = np.empty(NBINS)
for i in range(NBINS):
    g1sub, g2sub = sim_sub.make_submission(cvals[0], cvals[0], mvals[i], mvals[i],
                                           noise_sigma=NOISE_SIGMA)
    Q08_vs_m[i] = g3utils.metricQ08_const_shear(g1sub, g2sub, g1true, g2true)
    QZ1_vs_m[i] = g3utils.metricQZ1_const_shear(g1sub, g2sub, g1true, g2true)[0]
    QZ2_vs_m[i] = g3utils.metricQZ2_const_shear(g1sub, g2sub, g1true, g2true)[0]
    g1sub, g2sub = sim_sub.make_submission(cvals[i], cvals[i], mvals[0], mvals[0],
                                           noise_sigma=NOISE_SIGMA)
    Q08_vs_c[i] = g3utils.metricQ08_const_shear(g1sub, g2sub, g1true, g2true)
    QZ1_vs_c[i] = g3utils.metricQZ1_const_shear(g1sub, g2sub, g1true, g2true)[0]
    QZ2_vs_c[i] = g3utils.metricQZ2_const_shear(g1sub, g2sub, g1true, g2true)[0]


import matplotlib.pyplot as plt

plt.loglog(mvals, Q08_vs_m, label='Q08')
plt.loglog(mvals, QZ1_vs_m, label='QZ1')
plt.loglog(mvals, QZ2_vs_m, label='QZ2')
plt.title('Metrics versus m [c = 1.e-4]')
plt.xlabel('m = m1 = m2')
plt.ylabel('Q')
plt.legend()

plt.figure()
plt.loglog(cvals, Q08_vs_c, label='Q08')
plt.loglog(cvals, QZ1_vs_c, label='QZ1')
plt.loglog(cvals, QZ2_vs_c, label='QZ2')
plt.title('Metrics versus c [m = 1.e-3]')
plt.xlabel('c = c1 = c2')
plt.ylabel('Q')
plt.legend()

plt.show()






