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
"""This script was used to check the outputs of the Mathematica script pk.math.  More specifically,
it takes the output power spectra, tells GalSim to generate a random shear field according to the
power spectra, computes the shear correlation function xi_+, and compares it to the original one
given above to ensure that they are consistent.
"""
import numpy as np
import matplotlib.pyplot as plt
import galsim
import subprocess
import os
import sys

def run_corr2_ascii(x,y,e1,e2,min_sep,max_sep,nbins):
    """This is a simple utility to run corr2 (the two-point correlation function code) on an ascii
    file."""
    f = open('temp.cat','w')
    for (i,j), value in np.ndenumerate(x):
        xij = value
        yij = y[i,j]
        g1ij = e1[i,j]
        g2ij = e2[i,j]
        f.write('%e  %e  %e  %e\n'%(xij,yij,g1ij,g2ij))
    f.close()
    subprocess.Popen(['corr2','corr2.params',
                      'min_sep=%f'%min_sep,'max_sep=%f'%max_sep,'nbins=%f'%nbins]).wait()
    results = np.loadtxt('temp.e2')
    os.remove('temp.cat')
    os.remove('temp.e2')
    return results

fig = plt.figure()
ax = fig.add_subplot(111)
# Define values for theta_0, the scale length in the atmospheric PSF correlation function (which had
# been pre-determined for mathematica calculations and used to make filenames, meaning this part of
# the code can't be changed without changing those other pieces of code and rerunning Mathematica).
min_theta0 = 360.  # 0.1 degree
max_theta0 = 7200. # 2 degrees
ntheta0 = 11
theta0 = np.logspace(np.log10(min_theta0),np.log10(max_theta0),ntheta0)

# Loop over the values of theta_0, read in the P(k) from the precomputed files, and plot.
for t0 in theta0:
    strval = str(int(round(t0)))
    infile = 'pk_math/Pk'+strval+'.dat'
    pk_dat = np.loadtxt(infile).transpose()
    k = pk_dat[0]
    pk = 1.e-4*2.*np.pi*pk_dat[1] # Put in the proper normalization here.
    ax.plot(k, pk, label='theta0='+strval)
# And just for fun, let's plot a cosmological shear power spectrum as well.  For that, we have to
# convert the iCosmo outputs to arcsec instead of radians.
pk_file = '../ps/tables/cosmo-fid.zmed1.00.out'
pk_data = np.loadtxt(pk_file).transpose()
rad_to_arcsec = 3600.*180./np.pi
k_arcsec = pk_data[0]/rad_to_arcsec
pk_arcsec = pk_data[1]*(rad_to_arcsec**2)
ax.plot(k_arcsec, pk_arcsec, label='<gg>, z=1')
# Note some of the interesting scales.
ell_20 = 20./rad_to_arcsec
ell_2000 = 2000./rad_to_arcsec
ax.plot((ell_20, ell_20), (1.e-6,1.e5))
ax.plot((ell_2000, ell_2000), (1.e-6,1.e5))
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('k [1/arcsec]')
ax.set_ylabel('P(k) [arcsec^2]')
plt.legend()
print "Saving plot to file"
plt.savefig('diagnostics/test_pk_math.jpg')
# The outputs are quite interesting: for the lowest k (or ell) in the range probed by our grid, the
# atmospheric PSF anisotropy power spectrum is above the cosmological one.  However, for the highest
# k the situation is reversed, and the cosmological shear power spectrum is higher.  This is roughly
# what we expect, but it's good to confirm. 

# Next we want to explicitly check what happens when we use GalSim to draw shears according to those
# power spectra.  So we start by setting up our grid (using GREAT3 standard values for the sub-grids
# used for PSF determination for ground-based simulations, i.e., 2x2 degrees), since that determines
# the values of k that we actually care about.  We want to oversample the grid in order to not screw
# up the shear power on small scales.  And, we want to subsample so as to be able to offset our
# grids and include small-scale power (probably by 10, but for now, 5).
oversample = 20
subsample = 5
ngrid = 20*oversample*subsample
theta_deg = 2.*oversample
theta_arcsec = theta_deg*3600
dtheta_arcsec = theta_arcsec / ngrid
# Set up the k values corresponding to this grid.
k_min = 2.*np.pi / theta_arcsec
k_max = np.sqrt(2.)*np.pi/dtheta_arcsec
# The shear generation is expensive, so let's not do it every time, particularly if I just want to
# mess with the P(k) calculation.
do_shear = True
# Also we'll make an option to just use some of the grid, to save time.
use_subgrid = True
if do_shear:
    for t0 in theta0:
        print "Getting shears using GalSim for theta0=",t0
        strval = str(int(round(t0)))
        infile = 'pk_math/Pk'+strval+'.dat'
        pk_dat = np.loadtxt(infile).transpose()
        k = pk_dat[0]
        pk = 1.e-4*2.*np.pi*pk_dat[1] # put in the normalization here

        tab_pk = galsim.LookupTable(k, 0.5*pk, x_log=True, f_log=True)
        # Take the outputs and run corr2 to check that the outputs correspond to inputs.  Use 0.5*P
        # but give that to both P_E and P_B.
        ps = galsim.PowerSpectrum(tab_pk, tab_pk, units=galsim.arcsec)

        # Note, typically the buildGrid() method returns shear, not distortion e.  For nearly round
        # things, e~2*shear.  But for the atmospheric PSF stuff, I gave the code amplitudes
        # corresponding to some ellipticity variance, so the results should correspond to e1/e2 in
        # amplitude as well.  The alternative approach would have been to use shear variances
        # (factor of 4 lower) and then to take the e that come out and multiply by 2.  This seemed
        # silly so I didn't do it.  However, it does mean that the kappa fluctuations are too high
        # by a factor of 4, so I must reduce them.
        e1, e2, kappa = ps.buildGrid(grid_spacing=dtheta_arcsec,
                                     ngrid=ngrid,
                                     get_convergence=True)
        kappa /= 4.
        grid_range = (dtheta_arcsec) * np.arange(ngrid)
        x, y = np.meshgrid(grid_range, grid_range)
        if not use_subgrid:
            results = run_corr2_ascii(x,y,e1,e2,dtheta_arcsec,np.sqrt(2)*dtheta_arcsec*ngrid,ngrid)
            theta_corr2 = results[:,1]
            xip_corr2 = results[:,2]
            xim_corr2 = results[:,3]    
        else:
            n_iter = 0
            for x_ind in range(oversample):
                x_start = x_ind*ngrid/(oversample)
                x_end = (x_ind+1)*ngrid/(oversample)
                for y_ind in range(oversample):
                    print "  Subsample ",n_iter
                    y_start = y_ind*ngrid/(oversample)
                    y_end = (y_ind+1)*ngrid/(oversample)
                    e1_sub = e1[x_start:x_end, y_start:y_end]
                    e2_sub = e2[x_start:x_end, y_start:y_end]
                    x_sub = x[x_start:x_end, y_start:y_end]
                    y_sub = y[x_start:x_end, y_start:y_end]
                    results = run_corr2_ascii(x_sub,y_sub,e1_sub,e2_sub,dtheta_arcsec,np.sqrt(2)*dtheta_arcsec*ngrid/oversample,ngrid/(oversample))
                    if n_iter == 0:
                        theta_corr2 = results[:,1]
                        xip_corr2 = results[:,2]
                        xim_corr2 = results[:,3]
                    else:
                        xip_corr2 += results[:,2]
                        xim_corr2 += results[:,3]
                    n_iter += 1
            xip_corr2 /= oversample**2
            xim_corr2 /= oversample**2

        print "Plotting ellipticity correlation function."
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(theta_corr2, xip_corr2, label='observed xi+')
        ax.plot(theta_corr2, xim_corr2, label='observed xi-')
        theory = 1.e-4/(1+theta_corr2/t0)
        ax.plot(theta_corr2, 1.e-4/(1+theta_corr2/t0), label='input xi+')
        ratio = xip_corr2/theory
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('theta [arcsec]')
        ax.set_ylabel('xi')
        ax.set_title('theta0='+strval)
        plt.legend(loc=3)
        print "Saving plot to file"
        plt.savefig('diagnostics/test_pk_math_xi_'+strval+'.jpg')

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(theta_corr2, ratio)
        plt.ylim((0.5,1.5))
        ax.set_xscale('log')
        ax.set_xlabel('theta [arcsec]')
        ax.set_ylabel('observed/theory xi_+')
        plt.savefig('diagnostics/test_pk_math_xi_'+strval+'_ratio.jpg')
