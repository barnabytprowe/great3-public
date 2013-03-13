import galsim
import numpy as np
import pse
import matplotlib.pyplot as plt

### set up basic parameters ###
# how many times to do this (to beat down noise)
n_realization = 100
# file containing theoretical P(k), with fake values added above ell=2000
pkfile = 'test_pse/ps.wmap7lcdm.2000.dat'
theory_tab = galsim.LookupTable(file=pkfile, interpolant='linear')
# prefix for figures
figpref = 'test_pse/galsim_test_pse_'
# number of ell bins used for the estimation
n_ell = 15
# N for our grid used for estimating shears
grid_nx = 100
# length of grid in one dimension (degrees)
theta = 10. # degrees
dtheta = theta/grid_nx
# set verbosity
verbose = False


def doplot(ell, t, e, b, eb, pref, string, title, rat=None, lim=(1e-7,1e-4),
           bin_theory=None):
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

    if rat is not None:
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

test_ps_e=galsim.PowerSpectrum(e_power_function = theory_tab, units='radians')
test_ps_b=galsim.PowerSpectrum(b_power_function = theory_tab, units='radians')
test_ps_eb=galsim.PowerSpectrum(e_power_function = theory_tab,
                                b_power_function = theory_tab, units='radians')

e_p_e = np.zeros((n_ell, n_realization))
e_p_b = np.zeros((n_ell, n_realization))
e_p_eb = np.zeros((n_ell, n_realization))
b_p_e = np.zeros((n_ell, n_realization))
b_p_b = np.zeros((n_ell, n_realization))
b_p_eb = np.zeros((n_ell, n_realization))
eb_p_e = np.zeros((n_ell, n_realization))
eb_p_b = np.zeros((n_ell, n_realization))
eb_p_eb = np.zeros((n_ell, n_realization))

print "Averaging measured power spectra over realizations: ",n_realization
for ireal in range(n_realization):
    if verbose:
        print "Getting shears on a grid with E power only"
    g1, g2 = test_ps_e.buildGriddedShears(grid_spacing=dtheta, ngrid=grid_nx, units=galsim.degrees)
    g2 = -1.*g2 # this might not be necessary, though it was for the GREAT10 PS estimation code
    pse_e = pse.PowerSpectrumEstimator(grid_nx, theta, n_ell)
    if ireal == 0:
        ell, cee_e, cbb_e, ceb_e, c_binned_theory_nowt = pse_e.estimate(g1, g2, theory_func = theory_tab)
        ell, cee_e, cbb_e, ceb_e, c_binned_theory = pse_e.estimate(g1, g2, weight_EE=True,
                                                                   theory_func = theory_tab)
    else:
        ell, cee_e, cbb_e, ceb_e = pse_e.estimate(g1, g2, weight_EE=True)

    if verbose:
        print "Getting shears on a grid with B power only"
    g1, g2 = test_ps_b.buildGriddedShears(grid_spacing=dtheta, ngrid=grid_nx, units=galsim.degrees)
    g2 = -1.*g2 # this might not be necessary, though it was for the GREAT10 PS estimation code
    pse_b = pse.PowerSpectrumEstimator(grid_nx, theta, n_ell)
    ell, cee_b, cbb_b, ceb_b = pse_b.estimate(g1, g2)

    if verbose:
        print "Getting shears on a grid with E and B power"
    g1, g2 = test_ps_eb.buildGriddedShears(grid_spacing=dtheta, ngrid=grid_nx, units=galsim.degrees)
    g2 = -1.*g2 # this might not be necessary, though it was for the GREAT10 PS estimation code
    pse_eb = pse.PowerSpectrumEstimator(grid_nx, theta, n_ell)
    ell, cee_eb, cbb_eb, ceb_eb = pse_eb.estimate(g1, g2, weight_EE=True)

    e_p_e[:,ireal] = cee_e
    e_p_b[:,ireal] = cbb_e
    e_p_eb[:,ireal] = ceb_e
    b_p_e[:,ireal] = cee_b
    b_p_b[:,ireal] = cbb_b
    b_p_eb[:,ireal] = ceb_b
    eb_p_e[:,ireal] = cee_eb
    eb_p_b[:,ireal] = cbb_eb
    eb_p_eb[:,ireal] = ceb_eb

# get average and theoretical powers in a dimensionless way, ell(ell+1) C_ell/(2pi)
e_avgp_e = np.mean(e_p_e,1)*ell*(ell+1)/(2.*np.pi)
e_avgp_b = np.mean(e_p_b,1)*ell*(ell+1)/(2.*np.pi)
e_avgp_eb = np.mean(e_p_eb,1)*ell*(ell+1)/(2.*np.pi)
b_avgp_e = np.mean(b_p_e,1)*ell*(ell+1)/(2.*np.pi)
b_avgp_b = np.mean(b_p_b,1)*ell*(ell+1)/(2.*np.pi)
b_avgp_eb = np.mean(b_p_eb,1)*ell*(ell+1)/(2.*np.pi)
eb_avgp_e = np.mean(eb_p_e,1)*ell*(ell+1)/(2.*np.pi)
eb_avgp_b = np.mean(eb_p_b,1)*ell*(ell+1)/(2.*np.pi)
eb_avgp_eb = np.mean(eb_p_eb,1)*ell*(ell+1)/(2.*np.pi)

theory_p = np.zeros_like(ell)
theory_p_binned = c_binned_theory*ell*(ell+1.)/(2.*np.pi)
theory_p_binned_nowt = c_binned_theory_nowt*ell*(ell+1.)/(2.*np.pi)
for i in range(n_ell):
    theory_p[i] = theory_tab(ell[i])*ell[i]*(ell[i]+1)/(2.*np.pi)

print "Making figures of dimensionless power, and writing to files"
doplot(ell, theory_p, e_avgp_e, e_avgp_b, e_avgp_eb, figpref, 'input_pe',
       'Input P_E only', rat=e_avgp_e, bin_theory=theory_p_binned)
doplot(ell, theory_p, b_avgp_e, b_avgp_b, b_avgp_eb, figpref, 'input_pb',
       'Input P_B only', rat=b_avgp_b, bin_theory=theory_p_binned_nowt)
doplot(ell, theory_p, eb_avgp_e, eb_avgp_b, eb_avgp_eb, figpref, 'input_peb',
       'Input P_EB')
