import galsim
import numpy as np
import pse

### set up basic parameters ###
# how many times to do this (to beat down noise)
n_realization = 50
# file containing theoretical P(k), with fake values added above ell=2000
pkfile = 'test_pse/ps.wmap7lcdm.2000.dat'
theory_tab = galsim.LookupTable(file=pkfile)
# name of text file to store results
outfile = 'test_pse/galsim_test_pse.out'
# number of ell bins used for the estimation
n_ell = 20
# N for our grid used for estimating shears
grid_nx = 100
# length of grid in one dimension (degrees)
theta = 10. # degrees
dtheta = theta/grid_nx

test_ps_e=galsim.PowerSpectrum(e_power_function = pkfile, units='radians')
test_ps_b=galsim.PowerSpectrum(b_power_function = pkfile, units='radians')
test_ps_eb=galsim.PowerSpectrum(e_power_function = pkfile, b_power_function = pkfile, units='radians')

e_p_e = np.zeros((n_ell, n_realization))
e_p_b = np.zeros((n_ell, n_realization))
e_p_eb = np.zeros((n_ell, n_realization))
b_p_e = np.zeros((n_ell, n_realization))
b_p_b = np.zeros((n_ell, n_realization))
b_p_eb = np.zeros((n_ell, n_realization))
eb_p_e = np.zeros((n_ell, n_realization))
eb_p_b = np.zeros((n_ell, n_realization))
eb_p_eb = np.zeros((n_ell, n_realization))

for ireal in range(n_realization):
    print "Iteration ",ireal

    print "Getting shears on a grid with E power only"
    g1, g2 = test_ps_e.buildGriddedShears(grid_spacing=dtheta, ngrid=grid_nx, units=galsim.degrees)
    g2 = -1.*g2 # this might not be necessary, though it was for the GREAT10 PS estimation code
    pse_e = pse.PowerSpectrumEstimator(grid_nx, theta, n_ell)
    ell, cee_e, cbb_e, ceb_e = pse_e.estimate(g1, g2)

    print "Getting shears on a grid with B power only"
    g1, g2 = test_ps_b.buildGriddedShears(grid_spacing=dtheta, ngrid=grid_nx, units=galsim.degrees)
    g2 = -1.*g2 # this might not be necessary, though it was for the GREAT10 PS estimation code
    pse_b = pse.PowerSpectrumEstimator(grid_nx, theta, n_ell)
    ell, cee_b, cbb_b, ceb_b = pse_b.estimate(g1, g2)

    print "Getting shears on a grid with E and B power"
    g1, g2 = test_ps_eb.buildGriddedShears(grid_spacing=dtheta, ngrid=grid_nx, units=galsim.degrees)
    g2 = -1.*g2 # this might not be necessary, though it was for the GREAT10 PS estimation code
    pse_eb = pse.PowerSpectrumEstimator(grid_nx, theta, n_ell)
    ell, cee_eb, cbb_eb, ceb_eb = pse_eb.estimate(g1, g2)

    e_p_e[:,ireal] = cee_e
    e_p_b[:,ireal] = cbb_e
    e_p_eb[:,ireal] = ceb_e
    b_p_e[:,ireal] = cee_b
    b_p_b[:,ireal] = cbb_b
    b_p_eb[:,ireal] = ceb_b
    eb_p_e[:,ireal] = cee_eb
    eb_p_b[:,ireal] = cbb_eb
    eb_p_eb[:,ireal] = ceb_eb

e_avgp_e = np.mean(e_p_e,1)
e_avgp_b = np.mean(e_p_b,1)
e_avgp_eb = np.mean(e_p_eb,1)
b_avgp_e = np.mean(b_p_e,1)
b_avgp_b = np.mean(b_p_b,1)
b_avgp_eb = np.mean(b_p_eb,1)
eb_avgp_e = np.mean(eb_p_e,1)
eb_avgp_b = np.mean(eb_p_b,1)
eb_avgp_eb = np.mean(eb_p_eb,1)

theory_p = np.zeros_like(ell)
for i in range(n_ell):
    theory_p[i] = theory_tab(ell[i])

data_all = (ell, theory_p, e_avgp_e, e_avgp_b, e_avgp_eb, b_avgp_e, b_avgp_b, b_avgp_eb, eb_avgp_e, eb_avgp_b, eb_avgp_eb)
data = np.column_stack(data_all)
np.savetxt(outfile, data)

