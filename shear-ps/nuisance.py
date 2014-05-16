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
"""This file contains definitions for some nuisance functions (based on shapelets) that were used to
modify the cosmological shear power spectra.  All it does is make some plots of modified shear power
spectra to show that they are still cosmological-ish and basically not crazy."""
import numpy as np
import scipy.special as spec
import math
import matplotlib.pyplot as plt
import galsim

# Define constants up here: ell ranges, and shapelets order to use for modifying the power
# spectrum.
ell_min = 36.
ell_max = 3600.
ell_piv = 360.
ellfile = 'tables/cosmo-fid.zmed1.00.out'
n_trials = 6
max_order = 5
beta = 0.3
max_amp = 0.1

def bn(n, x, beta):
    """Routine to get dimensionless shapelets functions at some NumPy array of position `x`.

    These are defined as
        B_n(x, beta) = (1/sqrt(beta)) phi_n(x/beta).
    """
    phi_n_x = phin(n, x/beta)
    return phi_n_x/np.sqrt(beta)

def phin(n, x):
    """Routine to get shapelets basis functions at some NumPy array of position `x`."""
    # Get the H_n function from scipy.
    hn = spec.hermite(n)
    # Evaluate it at our values of x.
    hn_x = hn(x)
    # Now put in the exp factor.
    phi_n_x = hn_x * np.exp(-(x**2)/2.)
    # Finally, normalize it.
    phi_n_x /= np.sqrt((2.**n) * np.sqrt(np.pi) * math.factorial(n))
    return phi_n_x

# Make an array of ell values that we might want to use-- take from PS files.
print 'Making ell array'
elldata = np.loadtxt(ellfile)
ell = elldata[1:,0]
ell_use = ell[(ell >= ell_min) & (ell <= ell_max)]
lgell_use = np.log10(ell_use)
n_ell = len(lgell_use)
print 'Using ',n_ell,' ell values'

# Define our x vector = log10(ell/ell_piv).
x = lgell_use - np.log10(ell_piv)

# Initialize the RNG to generate the random numbers to use as coefficients.
ud = galsim.UniformDeviate()

# Get the bn for our chosen value of beta.
print 'Getting the basis functions at our ell values'
bvals = np.zeros((max_order+1, n_ell))
for order in range(0, max_order):
    bvals[order, :] = bn(order, x, beta)

# Make the nuisance function for however many trials we want.
print 'Generating ',n_trials,' nuisance functions'
nuisance = np.zeros((n_trials, n_ell))
for i in range(n_trials):
    for j in range(0, max_order):
        this_rand = -max_amp + 2*max_amp*ud()
        nuisance[i,:] += this_rand*bvals[j,:]

# Plot a bunch of nuisance functions for some value of beta, a few n.
# Note: we basically are going to multiply the power spectra by 1+these functions.
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
