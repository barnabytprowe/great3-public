import galsim
import numpy as np
import matplotlib.pyplot as plt
import subprocess
from scipy.special import jv

# This uses Mike Jarvis's corr2 program for calculating the correlation function.
# It is available at https://code.google.com/p/mjarvis/ and needs to be installed separately.

### set up basic parameters ###
# how many times to do this (to beat down noise)
n_realization = 100
# file containing theoretical P(k), with fake values added above ell=2000
pkfile = 'test_pse/ps.wmap7lcdm.2000.dat'
theory_tab = galsim.LookupTable(file=pkfile, interpolant='linear')
# prefix for figures
figpref = 'test_corr/galsim_test_corr2_'
# number of ell bins used for the estimation
n_ell = 15
# N for our grid used for estimating shears
grid_nx = 100
# length of grid in one dimension (degrees)
theta = 10. # degrees
dtheta = theta/grid_nx

full_theory_integral = True  # Set this to False to try to account for the fact that the FFT
                             # generation of the shear field ignores large scale power.
                             # (i.e. power in k modes < 2pi/theta)

extra_res = 1       # Extra resolution factor for g1,g2 grid.
                    # This means we use more of the large scale (small k) power in building
                    # the shear field.  Setting extra_res = 10 works reasonably well.

fair_grid = False   # Set to True if you don't want to keep the shear grid the same size when
                    # extra_res > 1.  We only need the extra_res for building the shapes to avoid 
                    # biases.  We don't need all the extra values for the measurements.
                    # However, they don't hurt -- they help beat down the noise, so if this 
                    # is False, you can get away with far fewer realizations.

# parameters for corr2:
min_sep = dtheta
max_sep = grid_nx * np.sqrt(2) * dtheta
nbins = 100

# parameters about how to extrapolate to small r
fit_n = grid_nx/4    # How many points to use in the fit
fit_quad = False     # Use quadratic term in logr?

min_ell = 2.*np.pi/(theta * np.pi/180)
max_ell = min_ell * grid_nx/np.sqrt(2.)

class xi_integrand:
    def __init__(self, pk, r, n):
        self.pk = pk
        self.r = r
        self.n = n
    def __call__(self, k):
        return k * self.pk(k) * jv(self.n, self.r*k)
        
def calculate_xi(r, pk):
    """Calculate xi+(r) or xi-(r) from a power spectrum.
    """
    #print 'Start calculate_xi'
    # xi+/-(r) = 1/2pi int(dk k P(k) J0/4(kr), k=0..inf)

    if full_theory_integral:
        # These are closer to the true correlation function.  But they don't match the observations.
        int_min = pk.x_min
        int_max = pk.x_max
    else:
        # These manage to get reasonably close agreement with the theory
        int_min = 2.*np.pi/(theta * extra_res * np.pi/180. * np.sqrt(2.))
        int_max = 2.*np.pi/(dtheta * np.pi/180.)
        #print 'int range = ',int_min,int_max

    rrad = r * np.pi/180.  # Convert to radians

    xip = np.zeros_like(r)
    xim = np.zeros_like(r)
    for (xi, n) in [ (xip, 0), (xim, 4) ]:
        for i in range(len(r)):
            integrand = xi_integrand(pk, rrad[i], n)
            xi[i] = galsim.integ.int1d(integrand, int_min, int_max,
                                       rel_err=1.e-6, abs_err=1.e-12)
        xi /= 2. * np.pi
    return xip, xim

class extended_xi:
    def __init__(self, r, xi, w):
        #print 'start extended_xi'
        #print 'r = ',r
        #print 'xi = ',xi
        #print 'w = ',w
        nonzero = (xi != 0.)
        #print 'r[nonzero] = ',r[nonzero]
        #print 'xi[nonzero] = ',xi[nonzero]
        #print 'w[nonzero] = ',w[nonzero]
        self.tab = galsim.LookupTable(r[nonzero],xi[nonzero])
        if np.all(xi[nonzero][0:fit_n] > 0.):
            self.power = True
            # Fit a power law to the small-r points to use as extrapolation to smaller r values.
            logr = np.log(r[nonzero][0:fit_n])
            logxi = np.log(xi[nonzero][0:fit_n])
            ww = w[nonzero][0:fit_n]
            #print 'logr = ',logr
            #print 'logxi = ',logxi
            if fit_quad:
                A = np.vstack([ww*logr**2, ww*logr, ww*np.ones(fit_n)]).T
                self.a, self.b, self.c = np.linalg.lstsq(A,ww*logxi)[0]
            else:
                A = np.vstack([ww*logr, ww*np.ones(fit_n)]).T
                self.b, self.c = np.linalg.lstsq(A,ww*logxi)[0]
                self.a = 0
            #print 'a, b, c = ',self.a, self.b, self.c
        else:
            self.power = False

    def __call__(self, r):
        if (r < self.tab.x_min): 
            logr = np.log(r)
            if self.power: return np.exp((self.a*logr+self.b)*logr+self.c)
            else: return 0.
        elif (r > self.tab.x_max): return 0.
        else: return self.tab(r)

def calculate_ps(ell, r, xip, xim, xix, w):
    """Calculate EE, BB, EB power spectra from xi+, xi-, xix.
    """
    #print 'Start calculate_ps'
    # Same integrals as above, so can still use the xi_integrand class
    # P_EE = 1/2 2pi int(dr r (xi+(r) J0(kr) + xi-(r) J4(kr)), r=0..inf)
    # P_BB = 1/2 2pi int(dr r (xi+(r) J0(kr) - xi-(r) J4(kr)), r=0..inf)
    # P_EB = 2pi int(dr r (xix(r) J0(kr)), r=0..inf) 
    #                              ^ Not sure if this is right.  Might be 2.  Or 4. ??

    rrad = r * np.pi/180.  # Convert to radians

    int_min = 0.
    int_max = rrad[xip!=0.][-1]
    #print 'int range = ',int_min,int_max

    xipint = np.zeros_like(ell)
    ximint = np.zeros_like(ell)
    xixint = np.zeros_like(ell)
    for (integ, xi,n) in [ (xipint, xip,0), (ximint, xim, 4), (xixint, xix, 0) ]:
        xifunc = extended_xi(rrad,xi,w)
        for i in range(len(ell)):
            integrand = xi_integrand(xifunc, ell[i], n)
            integ[i] = galsim.integ.int1d(integrand, int_min, int_max,
                                          rel_err=1.e-6, abs_err=1.e-12)
        integ *= 2. * np.pi
    # Convert to Delta^2
    xipint *= ell * (ell+1)/(2.*np.pi)
    ximint *= ell * (ell+1)/(2.*np.pi)
    xixint *= ell * (ell+1)/(2.*np.pi)

    return (xipint + ximint)*0.5, (xipint - ximint)*0.5, xixint

# Utility for plotting results
def doplot(ell, t, e, b, eb, pref, string, title, rat=None, lim=(1e-7,1e-4), bin_theory=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(ell, t, label='theory')
    if bin_theory is not None:
        ax.plot(ell, bin_theory, label='binned theory')
    ax.plot(ell, e, label='Observed EE power')
    ax.plot(ell, b, label='Observed BB power')
    ax.plot(ell, eb, label='Observed EB power')
    plt.ylim(lim)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('ell')
    ax.set_ylabel('ell*(ell+1)*C_ell/2pi')
    ax.set_title(title)
    plt.legend(loc=4)
    min_ell = 2.*np.pi/(theta*np.pi/180.)
    max_ell = np.sqrt(2.)*grid_nx*min_ell/2.
    ax.plot(np.array((min_ell,min_ell)), np.array(lim), color='black')
    ax.plot(np.array((max_ell,max_ell)), np.array(lim), color='black')
    figfile = pref + string + '.jpg'
    plt.savefig(figfile)
    print 'Wrote to file ',figfile

    #if rat is not None:
    if False:
        ratio = rat/t
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(ell, ratio, label='vs. theory')
        if bin_theory is not None:
            bin_ratio = rat/bin_theory
            ax.plot(ell, bin_ratio, label='vs. binned theory')
        plt.ylim(0.5,1.5)
        ax.set_xscale('log')
        ax.set_xlabel('ell')
        ax.set_ylabel('Observed / theory')
        ax.set_title(title)
        plt.legend(loc=4)
        ax.plot(np.array((min_ell,min_ell)), np.array((0.5,1.5)), color='black')
        ax.plot(np.array((max_ell,max_ell)), np.array((0.5,1.5)), color='black')
        ax.plot(np.array((0.5*min_ell,1.5*max_ell)),np.array((1.,1.)),color='black')
        figfile = pref + string + '.ratio.jpg'
        plt.savefig(figfile)
        print 'Wrote to file ',figfile

def doplot_corr(r, t_xip, t_xim, xip, xim, xix, pref, string, title, lim=(1e-8,2e-5)):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    nonzero = (xip != 0.)
    ax.plot(r, t_xip, 'black', label='theory xi+')
    ax.plot(r, t_xim, 'grey', label='theory xi-')
    ax.plot(r[nonzero], xip[nonzero], 'blue', label='Observed xi+')
    ax.plot(r[nonzero], xim[nonzero], 'green', label='Observed xi-')
    ax.plot(r[nonzero], xix[nonzero], 'red', label='Observed xix')
    ax.plot(r, -t_xip, 'black', ls='dashed')
    ax.plot(r, -t_xim, 'grey', ls='dashed')
    ax.plot(r[nonzero], -xip[nonzero], 'blue', ls='dashed')
    ax.plot(r[nonzero], -xim[nonzero], 'green', ls='dashed')
    ax.plot(r[nonzero], -xix[nonzero], 'red', ls='dashed')
    plt.ylim(lim)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('r')
    ax.set_ylabel('xi')
    ax.set_title(title)
    plt.legend(loc='upper right')
    ax.plot(np.array((min_sep,min_sep)), np.array(lim), color='black')
    ax.plot(np.array((max_sep,max_sep)), np.array(lim), color='black')
    figfile = pref + string + '.jpg'
    plt.savefig(figfile)
    print 'Wrote to file ',figfile

    if False:
        ratio = xip[nonzero]/t_xip[nonzero]
        fig = plt.figure()
        if (t_xim != 0.).any():
            ax = fig.add_subplot(211)
        else:
            ax = fig.add_subplot(111)
        ax.plot(r[nonzero], ratio, label='xi+ vs. theory')
        plt.ylim(0.5,1.5)
        ax.set_xscale('log')
        if (t_xim != 0.).any():
            ax.set_xlabel('r')
        ax.set_ylabel('Observed xi+ / theory')
        ax.set_title(title)
        plt.legend(loc='lower right')
        ax.plot(np.array((min_sep,min_sep)), np.array((0.5,1.5)), color='black')
        ax.plot(np.array((max_sep,max_sep)), np.array((0.5,1.5)), color='black')
        ax.plot(np.array((0.5*min_sep,1.5*max_sep)),np.array((1.,1.)),color='black')

        if (t_xim != 0.).any():
            ratio = xim[nonzero]/t_xim[nonzero]
            ax = fig.add_subplot(212)
            ax.plot(r[nonzero], ratio, label='xi- vs. theory')
            plt.ylim(0.5,1.5)
            ax.set_xscale('log')
            ax.set_xlabel('r')
            ax.set_ylabel('Observed xi- / theory')
            plt.legend(loc='lower right')
            ax.plot(np.array((min_sep,min_sep)), np.array((0.5,1.5)), color='black')
            ax.plot(np.array((max_sep,max_sep)), np.array((0.5,1.5)), color='black')
            ax.plot(np.array((0.5*min_sep,1.5*max_sep)),np.array((1.,1.)),color='black')

        figfile = pref + string + '.ratio.jpg'
        plt.savefig(figfile)
        print 'Wrote to file ',figfile

def run_corr2_fits(x,y,g1,g2):
    import pyfits
    import os
    # Use fits binary table for faster I/O. (Converting to/from strings is slow.)
    assert x.shape == y.shape
    assert x.shape == g1.shape
    assert x.shape == g2.shape
    x_col = pyfits.Column(name='x', format='1D', array=x.flatten() )
    y_col = pyfits.Column(name='y', format='1D', array=y.flatten() )
    g1_col = pyfits.Column(name='g1', format='1D', array=g1.flatten() )
    g2_col = pyfits.Column(name='g2', format='1D', array=g2.flatten() )
    cols = pyfits.ColDefs([x_col, y_col, g1_col, g2_col])
    table = pyfits.new_table(cols)
    phdu = pyfits.PrimaryHDU()
    hdus = pyfits.HDUList([phdu,table])
    hdus.writeto('temp.fits',clobber=True)
    subprocess.Popen(['corr2','corr2.params','file_name=temp.fits',
                      'x_col=x', 'y_col=y', 'g1_col=g1', 'g2_col=g2',
                      'min_sep=%f'%min_sep,'max_sep=%f'%max_sep,'nbins=%f'%nbins]).wait()
    results = np.loadtxt('temp.e2')
    os.remove('temp.fits')
    os.remove('temp.e2')
    return results


def run_corr2_ascii(x,y,g1,g2):
    f = open('temp.cat','w')
    for (i,j), value in np.ndenumerate(x):
        xij = value
        yij = y[i,j]
        g1ij = g1[i,j]
        g2ij = g2[i,j]
        f.write('%e  %e  %e  %e\n'%(xij,yij,g1ij,g2ij))
    f.close()
    subprocess.Popen(['corr2','corr2.params','file_name=temp.cat',
                      'x_col=1', 'y_col=2', 'g1_col=3', 'g2_col=4',
                      'min_sep=%f'%min_sep,'max_sep=%f'%max_sep,'nbins=%f'%nbins]).wait()
    results = np.loadtxt('temp.e2')
    os.remove('temp.cat')
    os.remove('temp.e2')
    return results
 

# define the galsim PowerSpectrum objects for the case of only E power, only B, and E+B
test_ps_e=galsim.PowerSpectrum(e_power_function = theory_tab, units='radians')
test_ps_b=galsim.PowerSpectrum(b_power_function = theory_tab, units='radians')
test_ps_eb=galsim.PowerSpectrum(e_power_function = theory_tab,
                                b_power_function = theory_tab, units='radians')

# Set up arrays to store results.
xip_e = np.zeros(nbins)
xim_e = np.zeros(nbins)
xix_e = np.zeros(nbins)
xip_b = np.zeros(nbins)
xim_b = np.zeros(nbins)
xix_b = np.zeros(nbins)
xip_eb = np.zeros(nbins)
xim_eb = np.zeros(nbins)
xix_eb = np.zeros(nbins)
#r = np.logspace(np.log10(min_sep),np.log10(max_sep),nbins)
r = np.zeros(nbins)
w = np.zeros(nbins)
ell = np.logspace(np.log10(min_ell),np.log10(max_ell),n_ell)

print "Averaging measured power spectra over realizations: ",n_realization
for ireal in range(n_realization):
    print 'ireal = ',ireal
    g1, g2 = test_ps_e.buildGrid(grid_spacing=dtheta, ngrid=grid_nx*extra_res,
                                 units=galsim.degrees)
    grid_range = dtheta * np.arange(grid_nx*extra_res)
    x, y = np.meshgrid(grid_range, grid_range)
    if fair_grid:
        g1 = g1[0:grid_nx,0:grid_nx]
        g2 = g2[0:grid_nx,0:grid_nx]
        x = x[0:grid_nx,0:grid_nx]
        y = y[0:grid_nx,0:grid_nx]
    results = run_corr2_fits(x,y,g1,g2)
    r += results[:,1]
    xip_e += results[:,2]
    xim_e += results[:,3]
    xix_e += results[:,5]
    w += results[:,7]

    g1, g2 = test_ps_b.buildGrid(grid_spacing=dtheta, ngrid=grid_nx*extra_res,
                                 units=galsim.degrees)
    if fair_grid:
        g1 = g1[0:grid_nx,0:grid_nx]
        g2 = g2[0:grid_nx,0:grid_nx]
    results = run_corr2_fits(x,y,g1,g2)
    r += results[:,1]
    xip_b += results[:,2]
    xim_b += results[:,3]
    xix_b += results[:,5]
    w += results[:,7]

    g1, g2 = test_ps_eb.buildGrid(grid_spacing=dtheta, ngrid=grid_nx*extra_res, 
                                  units=galsim.degrees)
    if fair_grid:
        g1 = g1[0:grid_nx,0:grid_nx]
        g2 = g2[0:grid_nx,0:grid_nx]
    results = run_corr2_fits(x,y,g1,g2)
    r += results[:,1]
    xip_eb += results[:,2]
    xim_eb += results[:,3]
    xix_eb += results[:,5]
    w += results[:,7]

# get averages
r /= 3*n_realization
xip_e /= n_realization
xim_e /= n_realization
xix_e /= n_realization
xip_b /= n_realization
xim_b /= n_realization
xix_b /= n_realization
xip_eb /= n_realization
xim_eb /= n_realization
xix_eb /= n_realization

#print 'average r = ',r
#print 'average xip_e = ',xip_e
#print 'average xip_b = ',xip_b
#print 'average xip_eb = ',xip_eb

print "Convert between corr and ps"
theory_p = np.zeros_like(ell)
for i in range(n_ell):
    theory_p[i] = theory_tab(ell[i])*ell[i]*(ell[i]+1)/(2.*np.pi)

theory_xip, theory_xim = calculate_xi(r,theory_tab)

ee_e, bb_e, eb_e = calculate_ps(ell, r, xip_e, xim_e, xix_e, w)
ee_b, bb_b, eb_b = calculate_ps(ell, r, xip_b, xim_b, xix_b, w)
ee_eb, bb_eb, eb_eb = calculate_ps(ell, r, xip_eb, xim_eb, xix_eb, w)

print "Making figures of dimensionless power, and writing to files"
doplot(ell, theory_p, ee_e, bb_e, eb_e, 
       figpref, 'input_pe_ps', 'Input P_E only', rat=ee_e)
doplot(ell, theory_p, ee_b, bb_b, eb_b,
       figpref, 'input_pb_ps', 'Input P_B only', rat=bb_b)
doplot(ell, theory_p, ee_eb, bb_eb, eb_eb, 
       figpref, 'input_peb_ps', 'Input P_EB', rat=(ee_eb+bb_eb)/2)

doplot_corr(r, theory_xip, theory_xim, xip_e, xim_e, xix_e, 
            figpref, 'input_pe_corr', 'Input P_E only')
doplot_corr(r, theory_xip, -theory_xim, xip_b, xim_b, xix_b, 
            figpref, 'input_pb_corr', 'Input P_B only')
doplot_corr(r, 2*theory_xip, np.zeros_like(r), xip_eb, xim_eb, xix_eb, 
            figpref, 'input_peb_corr', 'Input P_EB')

