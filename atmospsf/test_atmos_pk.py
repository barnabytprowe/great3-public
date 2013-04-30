import numpy as np
import matplotlib.pyplot as plt
import galsim
from scipy.special import jv

# This class defines the integrand for when we are trying to get
#   P(k) = (1/2pi) \int theta dtheta J_0(k theta) xi_+(theta).
# We have a special functional form for xi(theta), which is
#   xi_+(theta) = amp / (1 + theta/theta_0).
# This means we can store amp and theta, and evaluate the integrand for arbitrary theta, without
#   having to use a lookup table.
class p_integrand:
    def __init__(self, theta_0, k):
        self.theta_0 = theta_0
        self.k = k
    def __call__(self, theta):
        return theta * jv(0, self.k*theta)/(1+theta/self.theta_0)

# This function takes a vector of k values for which the user wants a power spectrum P(k), and
# integrates our special form for xi_+(theta) in order to get P(k).
def corrfunc_to_ps(amp, theta_0, k, max_r):
    # For now, assume all units for angles are arcsec, all units for k are arcsec^-1.
    # We are trying to get P(k) for our atmospheric PSF, which has the special property that xi_-=0,
    # implying P_E(k) = P_B(k).  Thus, there is some P(k) = P_E(k) + P_B(k) which relates to xi_+
    # via
    # P(k) = (1/2pi) \int theta dtheta J_0(k theta) xi_+(theta).
    # We have a special functional form for xi(theta), which is
    #   xi_+(theta) = amp / (1 + theta/theta_0).
    # Later, we will have to take this P(k), and when using the galsim lensing engine, assign P_E =
    # P_B = 0.5 P(k) that comes from this function.
    # We cannot really integrate to infinity, so the user can input some large value of max_r to use
    # as the upper limit in the integral.
    p = np.zeros_like(k)
    for k_ind in range(len(k)):
        print k_ind, k[k_ind]
        integrand = p_integrand(theta_0, k[k_ind])
        p[k_ind] = galsim.integ.int1d(integrand, 0., max_r,
                                      rel_err=1.e-4, abs_err=1.e-8)
    p *= amp/(2.*np.pi)
    return p

# Set up our grid (use GREAT3 standard values), since that determines the values of k that we
# actually care about.  Except we kind of want to oversample in order to not screw up the shear
# power on small scales, so we will oversample the grid.  And, we want to subsample so as to be able
# to offset our grids and include small-scale power (probably by 10, but for now just do 5 for
# computational efficiency).
oversample = 5
subsample = 5
ngrid = 100*oversample*subsample
theta_arcsec = 10.*3600*oversample
dtheta_arcsec = theta_arcsec / ngrid
# Set up the k values corresponding to this grid.
k_min = 2.*np.pi / theta_arcsec
k_max = np.sqrt(2.)*np.pi/dtheta_arcsec
# Now make a vector of k values to use to tabulate P(k)
test_k = np.logspace(np.log10(k_min), np.log10(k_max), num=20)
amp = 1.e-4
print "To PS, 1"
p1 = corrfunc_to_ps(amp, 100., test_k, 2*theta_arcsec)
print "To PS, 2"
p2 = corrfunc_to_ps(amp, 1000., test_k, 6*theta_arcsec)
print "To PS, 3"
p3 = corrfunc_to_ps(amp, 3000., test_k, 9.*theta_arcsec)
print "Plotting"
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(test_k, p1, label='theta0=100 arcsec')
ax.plot(test_k, p2, label='theta0=1000 arcsec')
ax.plot(test_k, p3, label='theta0=3000 arcsec')
# and just for fun, let's plot a cosmological PS.  For that, we have to convert the iCosmo one to
# arcsec instead of radians.
pk_file = '/Users/rmandelb/great3/GalSim/examples/data/cosmo-fid.zmed1.00.out'
pk_data = np.loadtxt(pk_file).transpose()
rad_to_arcsec = 3600.*180./np.pi
k_arcsec = pk_data[0]/rad_to_arcsec
pk_arcsec = pk_data[1]*(rad_to_arcsec**2)
ax.plot(k_arcsec, pk_arcsec, label='cosmic shear, zmed=1')
ell_20 = 20./rad_to_arcsec
ell_2000 = 2000./rad_to_arcsec
ax.plot((ell_20, ell_20), (1.e-5,1.e5), label='ell=20')
ax.plot((ell_2000, ell_2000), (1.e-5,1.e5), label='ell=2000')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('k [1/arcsec]')
ax.set_ylabel('P(k) [arcsec^2]')
plt.legend()
plt.savefig('diagnostics/test_atmos_pk.jpg')

# Take the outputs and run Mike's corr2 program to check that the outputs correspond to
# inputs.
