import numpy as np
import scipy.special as spec
import math
import matplotlib.pyplot as plt

# define constants up here
ell_min = 36.
ell_max = 3600.
ell_piv = 360.
ellfile = 'tables/cosmo-fid.zmed1.00.out'

# Get dimensionless shapelets functions at some numpy array of position x.
# These are defined as B_n(x, beta) = (1/sqrt(beta)) phi_n(x/beta).
def bn(n, x, beta):
    phi_n_x = phin(n, x/beta)
    return phi_n_x/np.sqrt(beta)

# Get shapelets basis functions at some numpy array of position x.
def phin(n, x):
    # get the H_n function from scipy
    hn = spec.hermite(n)
    # evaluate it at our values of x
    hn_x = hn(x)
    # now put in the exp factor
    phi_n_x = hn_x * np.exp(-(x**2)/2.)
    # finally, normalize it
    phi_n_x /= np.sqrt((2.**n) * np.sqrt(np.pi) * math.factorial(n))
    return phi_n_x

# make an array of ell values that we might want to use-- take from PS files
print 'Making ell array'
elldata = np.loadtxt(ellfile)
ell = elldata[1:,0]
ell_use = ell[(ell >= ell_min) & (ell <= ell_max)]
lgell_use = np.log10(ell_use)

# define our x vector = log10(ell/ell_piv)
x = lgell_use - np.log10(ell_piv)

# make the bn for a few values of n, beta
print 'Generating shapelets as a function of ell'
b0_beta1 = bn(0.0, x, 0.3)
b0_beta2 = bn(0.0, x, 0.6)
b0_beta3 = bn(0.0, x, 0.9)
b1_beta1 = bn(1.0, x, 0.3)
b1_beta2 = bn(1.0, x, 0.6)
b1_beta3 = bn(1.0, x, 0.9)
b2_beta1 = bn(2.0, x, 0.3)
b2_beta2 = bn(2.0, x, 0.6)
b2_beta3 = bn(2.0, x, 0.9)
b3_beta1 = bn(3.0, x, 0.3)
b3_beta2 = bn(3.0, x, 0.6)
b3_beta3 = bn(3.0, x, 0.9)

# plot a bunch of bn for some value of beta, a few n
print 'Plotting...'
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(ell_use, b0_beta1, color='blue', label='n=0, beta=0.3')
ax.plot(ell_use, b1_beta1, color='red', label='n=1, beta=0.3')
ax.plot(ell_use, b1_beta2, color='red', linestyle='dashed', label='n=1, beta=0.6')
ax.plot(ell_use, b2_beta1, color='green', label='n=2, beta=0.3')
ax.plot(ell_use, b3_beta1, color='magenta', label='n=3, beta=0.3')
ax.set_xscale('log')
ax.set_xlabel('ell')
ax.set_ylabel('shapelets basis function values')
ax.set_title('Shapelets functions to modify power spectra')
plt.legend(loc=4)
plt.savefig('plots/shapelets_vary_n.eps')

# plot a bunch of bn for 2 values of n, multiple different beta's
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(ell_use, b0_beta1, color='blue', label='n=0, beta=0.3')
ax.plot(ell_use, b0_beta2, color='red', label='n=0, beta=0.6')
ax.plot(ell_use, b0_beta3, color='green', label='n=0, beta=0.9')
ax.plot(ell_use, b3_beta1, color='magenta', label='n=3, beta=0.3')
ax.plot(ell_use, b3_beta2, color='cyan', label='n=3, beta=0.6')
ax.plot(ell_use, b3_beta3, color='black', label='n=3, beta=0.9')
ax.set_xscale('log')
ax.set_xlabel('ell')
ax.set_ylabel('shapelets basis function values')
ax.set_title('Shapelets functions to modify power spectra')
plt.legend(loc=4)
plt.savefig('plots/shapelets_vary_beta.eps')
