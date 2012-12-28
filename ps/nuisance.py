import numpy as np
import scipy.special as spec
import math
import matplotlib.pyplot as plt
import galsim

# define constants up here
ell_min = 36.
ell_max = 3600.
ell_piv = 360.
ellfile = 'tables/cosmo-fid.zmed1.00.out'
n_trials = 6
max_order = 5
beta = 0.3
max_amp = 0.1

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
n_ell = len(lgell_use)
print 'Using ',n_ell,' ell values'

# define our x vector = log10(ell/ell_piv)
x = lgell_use - np.log10(ell_piv)

# initialize the RNG to generate the random numbers to use as coefficients
ud = galsim.UniformDeviate()

# get the bn for our chosen value of beta
print 'Getting the basis functions at our ell values'
bvals = np.zeros((max_order+1, n_ell))
for order in range(0, max_order):
    bvals[order, :] = bn(order, x, beta)

# make the nuisance function for however many trials we want
print 'Generating ',n_trials,' nuisance functions'
nuisance = np.zeros((n_trials, n_ell))
for i in range(n_trials):
    for j in range(0, max_order):
        this_rand = -max_amp + 2*max_amp*ud()
        nuisance[i,:] += this_rand*bvals[j,:]

# plot a bunch of bn for some value of beta, a few n
print 'Plotting...'
fig = plt.figure()
ax = fig.add_subplot(111)
for i in range(n_trials):
    ax.plot(ell_use, nuisance[i,:])
ax.set_xscale('log')
ax.set_xlabel('ell')
ax.set_ylabel('nuisance function')
ax.set_title('randomly generated shapelets-based nuisance functions')
plt.savefig('plots/shapelets_nuisance.eps')
