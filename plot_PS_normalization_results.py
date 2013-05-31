import os
import cPickle
import numpy as np
import matplotlib.pyplot as plt
import calculate_QG10_var_shear_normalization as calcnorm

NORMFILE = os.path.join('results', 'normalizationv2_G10_QuadPS_N1000_noise_sigma0.05.pkl')
# Load up the tuple of results
norm_data = cPickle.load(open(NORMFILE, 'rb'))

qG10unnorm = norm_data[0]
qQuadPSunnorm = norm_data[1]

ctest = calcnorm.CTEST
mtest = calcnorm.MTEST

plt.clf()
if not os.path.isdir('plots'):
    os.mkdir('plots')

# First plot the histograms of QG10 for each of the ctest values for m = mtest[0] (fiducial)
plotfile1 = os.path.join('plots', 'normalization_QG10_vs_c.png')
for c, i in zip(ctest, range(len(ctest))):
    plt.hist(qG10unnorm[i, 0, :].flatten(), bins=30, range=(0., 6.e7), label=r'c$_i$='+str(c))
plt.xlabel(r'QG10 (unnormalized fiducial m$_i$='+str(mtest[0])+')')
plt.legend()
plt.savefig(plotfile1)

# Then plot the histograms of QQuadPS for each of the ctest values for m = mtest[0] (fiducial)
plotfile2 = os.path.join('plots', 'normalization_QQuadPS_vs_c.png')
plt.clf()
for c, i in zip(ctest, range(len(ctest))):
    plt.hist(qQuadPSunnorm[i, 0, :].flatten(), bins=30, range=(0., 6.e7), label=r'c$_i$='+str(c))
plt.xlabel(r'QQuadPS (unnormalized fiducial m$_i$='+str(mtest[0])+')')
plt.legend()
plt.savefig(plotfile2)

# Then calculate the normalization factors for the QG10 and QQuadPS metrics
normQG10 = np.mean(qG10unnorm[:, 0, :].flatten())
normQQuadPS = np.mean(qQuadPSunnorm[:, 0, :].flatten())
# Use this to make normalized arrays of test case metric values
qG10 = 1000. * qG10unnorm / normQG10
qQuadPS = 1000. * qQuadPSunnorm / normQQuadPS

# Then plot histograms of normalized values as a function of m
# QG10
plt.clf()
plotfile3 = os.path.join('plots', 'hists_QG10_vs_m.png')
for m, j in zip(mtest, range(len(mtest))):
    plt.hist(qG10[:, j, :].flatten(), bins=75, range=(0., 2500.), label=r'm$_i$='+str(m))
plt.legend()
plt.xlabel('QG10')
plt.savefig(plotfile3)

# QQuadPS
plt.clf()
plotfile4 = os.path.join('plots', 'hists_QQuadPS_vs_m.png')
for m, j in zip(mtest, range(len(mtest))):
    plt.hist(qQuadPS[:, j, :].flatten(), bins=75, range=(0., 2500.), label=r'm$_i$='+str(m))
plt.legend()
plt.xlabel('QQuadPS')
plt.savefig(plotfile4)

# Then plot means and standard deviations of these histograms
# QG10
plt.clf()
plotfile5 = os.path.join('plots', 'QG10_vs_m.png')
meanQG10 = []
stdQG10 = []
for m, j in zip(mtest, range(len(mtest))):
    meanQG10.append(np.mean(qG10[:, j, :].flatten()))
    stdQG10.append(np.std(qG10[:, j, :].flatten()))
plt.errorbar(np.log10(mtest), meanQG10, yerr=stdQG10)
plt.xlim(-.5, -3.5)
plt.xlabel(r'log$_{10}$(m$_i$)')
plt.ylabel('QG10')
plt.savefig(plotfile5)

# QQuadPS
plt.clf()
plotfile6 = os.path.join('plots', 'QQuadPS_vs_m.png')
meanQQuadPS = []
stdQQuadPS = []
for m, j in zip(mtest, range(len(mtest))):
    meanQQuadPS.append(np.mean(qQuadPS[:, j, :].flatten()))
    stdQQuadPS.append(np.std(qQuadPS[:, j, :].flatten()))
plt.errorbar(np.log10(mtest), meanQQuadPS, yerr=stdQQuadPS)
plt.xlim(-.5, -3.5)
plt.xlabel(r'log$_{10}$(m$_i$)')
plt.ylabel('QQuadPS')
plt.savefig(plotfile6)



