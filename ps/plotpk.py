import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.optimize as opti

print "Getting Cell data from iCosmo outputs..."
file1 = 'tables/cosmo-fid.zmed0.75.out'
file2 = 'tables/cosmo-fid.zmed1.00.out'
file3 = 'tables/cosmo-fid.zmed1.25.out'

data1 = np.loadtxt(file1)
data2 = np.loadtxt(file2)
data3 = np.loadtxt(file3)

ell = data1[:,0]
cell1 = data1[:,1]
cell2 = data2[:,1]
cell3 = data3[:,1]

# I think that we have to do P = ell^2 C_ell/(2 pi)
p1 = cell1*(ell**2)/(2.*np.pi)
p2 = cell2*(ell**2)/(2.*np.pi)
p3 = cell3*(ell**2)/(2.*np.pi)

print "Making figure of C_ell vs. ell"
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(ell,cell1,color='blue',label='zmed=0.75')
ax.plot(ell,cell2,color='red',label='zmed=1.00')
ax.plot(ell,cell3,color='green',label='zmed=1.25')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('ell')
ax.set_ylabel('C_ell')
ax.set_title('Fiducial cosmology, 3 median redshifts')
plt.legend()
ax.plot(np.array((36.,36.)),np.array((1.e-13,1.e-7)),color='black')
ax.plot(np.array((3600.,3600.)),np.array((1.e-13,1.e-7)),color='black')
plt.savefig('plots/fiducial_cell.eps')

print "Fitting Power to power-law functions..."
sublgell = np.log10(ell[(ell >= 36) & (ell <= 3600)])
sublgp1 = np.log10(p1[(ell >= 36) & (ell <= 3600)])
sublgp2 = np.log10(p2[(ell >= 36) & (ell <= 3600)])
sublgp3 = np.log10(p3[(ell >= 36) & (ell <= 3600)])
powerlawfit = lambda p, x, y: (y - (p[0] + p[1]*x))
pinit = [-7.0, 1.0]
out1, success = opti.leastsq(powerlawfit, pinit, args=(sublgell, sublgp1))
print out1[0], out1[1]
out2, success = opti.leastsq(powerlawfit, pinit, args=(sublgell, sublgp2))
print out2[0], out2[1]
out3, success = opti.leastsq(powerlawfit, pinit, args=(sublgell, sublgp3))
print out3[0], out3[1]
subell = 10**sublgell
fitp1 = 10**(out1[0] + out1[1]*sublgell)
fitp2 = 10**(out2[0] + out2[1]*sublgell)
fitp3 = 10**(out3[0] + out3[1]*sublgell)
realp1 = 10**sublgp1
realp2 = 10**sublgp2
realp3 = 10**sublgp3
ratio1 = realp1/fitp1
ratio2 = realp2/fitp2
ratio3 = realp3/fitp3

print "Plotting power..."
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(ell,p1,color='blue',label='zmed=0.75')
ax.plot(subell,fitp1,color='blue',linestyle='dotted',label='(best power law)')
ax.plot(ell,p2,color='red',label='zmed=1.00')
ax.plot(subell,fitp2,color='red',linestyle='dotted',label='(best power law)')
ax.plot(ell,p3,color='green',label='zmed=1.25')
ax.plot(subell,fitp3,color='green',linestyle='dotted',label='(best power law)')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('ell')
ax.set_ylabel('P')
ax.set_title('Fiducial cosmology, 3 median redshifts')
plt.legend(loc=4)
ax.plot(np.array((36.,36.)),np.array((1.e-7,1.e-3)),linestyle='dashed',color='black')
ax.plot(np.array((3600.,3600.)),np.array((1.e-7,1.e-3)),linestyle='dashed',color='black')
plt.savefig('plots/fiducial_p.eps')

print "Plotting ratio of fit to real power..."
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(subell,ratio1,color='blue')
ax.plot(subell,ratio2,color='red')
ax.plot(subell,ratio3,color='green')
ax.set_xscale('log')
ax.set_xlabel('ell')
ax.set_ylabel('ratio of fitted to real power')
plt.savefig('plots/p_fit_ratio.eps')
