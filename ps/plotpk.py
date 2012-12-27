import numpy as np
import matplotlib.pyplot as plt
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
A = np.vander(sublgell, 5)
(coeffs1, residuals, rank, sing_vals) = np.linalg.lstsq(A, sublgp1)
(coeffs2, residuals, rank, sing_vals) = np.linalg.lstsq(A, sublgp2)
(coeffs3, residuals, rank, sing_vals) = np.linalg.lstsq(A, sublgp3)
print coeffs1
print coeffs2
print coeffs3
f1 = coeffs1[0]*(sublgell**4) + coeffs1[1]*(sublgell**3) + coeffs1[2]*(sublgell**2) + coeffs1[3]*sublgell + coeffs1[4]
f2 = coeffs2[0]*(sublgell**4) + coeffs2[1]*(sublgell**3) + coeffs2[2]*(sublgell**2) + coeffs2[3]*sublgell + coeffs2[4]
f3 = coeffs3[0]*(sublgell**4) + coeffs3[1]*(sublgell**3) + coeffs3[2]*(sublgell**2) + coeffs3[3]*sublgell + coeffs3[4]
subell = 10**sublgell
fitp1 = 10**f1
fitp2 = 10**f2
fitp3 = 10**f3
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
ax.plot(subell,fitp1,color='blue',linestyle='dotted',label='(4th order polynomial)')
ax.plot(ell,p2,color='red',label='zmed=1.00')
ax.plot(subell,fitp2,color='red',linestyle='dotted',label='(4th order polynomial)')
ax.plot(ell,p3,color='green',label='zmed=1.25')
ax.plot(subell,fitp3,color='green',linestyle='dotted',label='(4th order polynomial)')
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
