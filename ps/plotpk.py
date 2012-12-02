import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

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

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(ell,p1,color='blue',label='zmed=0.75')
ax.plot(ell,p2,color='red',label='zmed=1.00')
ax.plot(ell,p3,color='green',label='zmed=1.25')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('ell')
ax.set_ylabel('P')
ax.set_title('Fiducial cosmology, 3 median redshifts')
plt.legend()
ax.plot(np.array((36.,36.)),np.array((1.e-7,1.e-3)),linestyle='dashed',color='black')
ax.plot(np.array((3600.,3600.)),np.array((1.e-7,1.e-3)),linestyle='dashed',color='black')
plt.savefig('plots/fiducial_p.eps')
