"""@file test_bmode_noise.py

This module started out as the script used to demonstrate that B-mode shape cancellation would work
even for the COSMOS galaxies which don't have the p(|e|) you would get from a purely Gaussian random
B-mode field.  Here we tested that the ranked-ordering prescription described in the handbook would
actually work.

In later tests for the validation of the GREAT3 catalogues, there were some additions to this script
to create output for the aperture mass dispersion split into E/B modes as a function of the
kmin_factor and kmax_factor input to the PowerSpectrum instance .buildGrid() method.  This turned
out to be important as it leads to power leaking from the pure B-mode noise into the E-mode! 

As a result, the metric calculation for Q_v was required to take into account the small quantity of
leaked B-mode noise from the galaxy intrinsic (pre-lensing) ellipticities when comparing submissions
to the reference "truth" E-mode aperture mass dispersion.
"""

import os
import pyfits
import numpy as np
import galsim

# Please modify this filepath to point to the correct location for real_galaxy_23.5_shapes.fits if
# necessary
SHAPECAT = os.path.join("..", "inputs", "galdata", "real_galaxy_23.5_shapes.fits")
RANDOM_SEED = 1335133
NGRID = 500
SIGMA_NOISE = 0.05

def generate_bmode_shears(var, ngrid, rng=None, kmax_factor=16, kmin_factor=1):
    """Generate b-mode shears, returns: g1, g2, gmag, gphi
    """
    ps = galsim.PowerSpectrum(
        b_power_function=lambda karr : var * np.ones_like(karr) / float(kmax_factor)**2)
    g1, g2 = ps.buildGrid(
        grid_spacing=1., ngrid=ngrid, rng=rng, kmax_factor=kmax_factor, kmin_factor=kmin_factor)
    # These g1, g2 *shears* do not have to be on the unit disc, so we have to convert them to a |g|
    # like ellipticity via the result for a circle following such a shear, using Schneider (2006)
    # eq. 12.  Note, this itself will mean that the ellips are not pure B-mode!!
    gmag = np.sqrt(g1 * g1 + g2 * g2)
    for i in range(g1.shape[0]): # ugly but my arm is broken
        for j in range(g2.shape[1]):
            if gmag[i, j] > 1.:
                g1[i, j] = g1[i, j] / gmag[i, j]**2
                g2[i, j] = g2[i, j] / gmag[i, j]**2
                gmag[i, j] = np.sqrt(g1[i, j]**2 + g2[i, j]**2) 
            else:
                pass
    gphi = .5 * np.arctan2(g2, g1)
    return g1, g2, gmag, gphi

def select_from_catalog(gmag, ngrid):
    """Returns an ngrid x ngrid array with values selected at random from the input array gmag.
    """
    try:
        selection = np.random.choice(gmag, (ngrid, ngrid))
    except AttributeError:
        import random
        gmag_indices = range(len(gmag))
        selection = [
            gmag[index] for index in [random.choice(gmag_indices) for i in range(ngrid**2)]]
        selection = np.reshape(np.asarray(selection), (ngrid, ngrid))
    return selection

def match_gmags_by_rank(gmagt, gmag, ngrid):
    """Takes an array of target gmags (gmagt) and identifies the corresponding rank-sorted gmag
    for each element, then returns these as an nside x nside array.
    """
    isorted_gmagt = np.argsort(gmagt.flatten())
    gmag_sample_sorted = np.sort(gmag.flatten())
    sorted_gmag_dict = {}
    for i in range(ngrid * ngrid):
        sorted_gmag_dict[isorted_gmagt[i]] = gmag_sample_sorted[i]
    gmags = np.reshape([sorted_gmag_dict[i] for i in range(ngrid * ngrid)], (ngrid, ngrid))
    return gmags


if __name__ == "__main__":

    print "Reading shape data from "+SHAPECAT
    data = pyfits.getdata(SHAPECAT)
    # Select well-measured e1 and e2
    do_meas = data.field('do_meas')
    e1 = data.field('e1')[do_meas > -0.5]
    e2 = data.field('e2')[do_meas > -0.5]
    emag = np.sqrt(e1 * e1 + e2 * e2)
    ephi = .5 * np.arctan2(e2, e1)
    # Convert these to g1, g2
    gmag = emag / (1. + np.sqrt(1. - e1 * e1 - e2 * e2))
    gphi = ephi
    # Catch the 16 nans (|e| > 1) and discard
    gmag = gmag[~np.isnan(gmag)]
    gphi = gphi[~np.isnan(gmag)]
    g1 = gmag * np.cos(2. * gphi)
    g2 = gmag * np.sin(2. * gphi)
    # Estimate total var from each component (assume isotropic) to use when generating B-mode
    # noise as a Gaussian field
    gvar = g1.var() + g2.var()
    
    # Get some target 'pure' B-mode shears
    g1t, g2t, gmagt, gphit = generate_bmode_shears(gvar, NGRID)

    # Set up a power spectrum estimator
    my_pse = galsim.pse.PowerSpectrumEstimator(NGRID, NGRID * 180. / np.pi)
    # Get the e-mode leakage of the 'target' intrinsic shape noise, present due to the unit disc
    # upper limit for |g| <= 1
    ell, eet, bbt, ebt = my_pse.estimate(g1t, g2t)

    if not os.path.isdir('plots'): os.mkdir('plots')

    import matplotlib.pyplot as plt

    plt.clf()
    plt.loglog(ell, eet, color='k', label='Best case intrinsic E')
    plt.loglog(ell, bbt, color='r', label='Best case intrinsic B')
    plt.ylabel('Power')
    plt.ylim(1.e-6, 1.e1)
    plt.legend()
    plt.axhline(SIGMA_NOISE**2, ls='--', color='k')
    plt.savefig('./plots/bestcase_bmode_shapenoise.png')

    plt.clf()
    gmag_selected = select_from_catalog(gmag, NGRID)
    g1crude = gmag_selected * np.cos(2. * gphit)
    g2crude = gmag_selected * np.sin(2. * gphit)
    ell, eec, bbc, ebc = my_pse.estimate(g1crude, g2crude)
    plt.loglog(ell, eec, color='k', label='Crude (unmatched |g|) intrinsic E')
    plt.loglog(ell, bbc, color='r', label='Crude (unmatched |g|) intrinsic B')
    plt.ylabel('Power')
    plt.ylim(1.e-6, 1.e1)
    plt.legend()
    plt.axhline(SIGMA_NOISE**2, ls='--', color='k')
    plt.savefig('./plots/crude_bmode_shapenoise.png')

    plt.clf()
    gmags = match_gmags_by_rank(gmagt, gmag_selected, NGRID)
    g1s = gmags * np.cos(2. * gphit)
    g2s = gmags * np.sin(2. * gphit)
    ell, ees, bbs, ebs = my_pse.estimate(g1s, g2s)
    plt.loglog(ell, ees, color='k', label='Rank matched |g| intrinsic E')
    plt.loglog(ell, bbs, color='r', label='Rank matched |g| intrinsic B')
    plt.ylabel('Power')
    plt.ylim(1.e-6, 1.e1)
    plt.legend()
    plt.axhline(SIGMA_NOISE**2, ls='--', color='k')
    plt.savefig('./plots/ranked_bmode_shapenoise.png')

    print "Mean leaked E = "+str(ees.mean())
    print "SIGMA_NOISE**2 = "+str(SIGMA_NOISE**2)
    print "Ratio = "+str(SIGMA_NOISE**2 / ees.mean())

    plt.clf()
    plt.subplot(211)
    plt.title('COSMOS fits P(|g|)')
    plt.hist(gmag_selected.flatten(), range=(0, 1), bins=30)
    plt.ylim(0, 22500)
    plt.subplot(212)
    plt.title('Gaussian field B-mode P(|g|)')
    plt.hist(gmagt.flatten(), range=(0, 1), bins=30, color='r')
    plt.ylim(0, 22500)
    plt.xlabel('P(|g|)')
    plt.savefig('./plots/pe_hist.png')

    # Trying sending this to corr2
    x, y = np.meshgrid(np.arange(NGRID) * 10. / float(NGRID), np.arange(NGRID) * 10. / float(NGRID))
    import sys
    sys.path.append(os.path.join("..", "presubmission_script"))
    import presubmission
    for g1, g2, typestring in zip((g1t, g1s), (g2t, g2s), ("Pure", "Ranked")):

        results = presubmission.run_corr2(
            x.flatten(), y.flatten(), g1, g2, np.ones(NGRID * NGRID))
        ylim = 2.e-5
        theta = []
        mapE = []
        mapB = []
        mmxa = []
        mmxb = []
        maperr = []
        for line in results:

            theta.append(float(line[0]))
            mapE.append(float(line[1]))
            mapB.append(float(line[2]))
            mmxa.append(float(line[3]))
            mmxb.append(float(line[4]))
            maperr.append(float(line[5]))

        theta = np.asarray(theta)
        mapE = np.asarray(mapE)
        mapB = np.asarray(mapB)
        mmxa = np.asarray(mmxa)
        mmxb = np.asarray(mmxb)
        maperr = np.asarray(maperr)
        plt.clf()
        plt.errorbar(theta, mapE, fmt="k-", yerr=maperr, lw=2., label=r"Map$^2$")
        plt.xscale("log")
        plt.errorbar(theta, mapB, fmt="r-", yerr=maperr, label=r"Mx$^2$", lw=2.)
        plt.errorbar(theta, mmxa, fmt="g-", yerr=maperr, label="MMx(a)", lw=2.)
        plt.errorbar(theta, mmxb, fmt="b-", yerr=maperr, label="MMx(b)", lw=2.)
        plt.axhline(color="k", ls="--")
        plt.legend()
        #print mapE
        #print mapB
        plt.ylim(-ylim, ylim)
        plt.xlabel("R [degrees]")
        plt.title(typestring+" B-mode shapes")
        plt.savefig("./plots/map_bmode_only_"+typestring+".png")
        print mapE.mean(), np.abs(mapE).max()
        print mapB.mean(), np.abs(mapB).max()
        plt.show()


