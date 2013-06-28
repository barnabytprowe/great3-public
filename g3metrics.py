"""Module containing GREAT3 metric calculation utilities.
"""
import numpy as np
import galsim

def make_const_truth_normal_dist(ntrue, nims, true_sigma=0.03,
                                 saveto='./g3truth/g3_const_shear_truth.dat'):
    """Generate truth catalogues with a N(0, TRUE_SIGMA) distribution of input truth values in
    g1 and g2.

    Set saveto=None to disable hardcopy output.  nims should be an integer multiple of ntrue.
    """
    if (nims % ntrue) != 0:
        raise ValueError("nims must be divisible by ntrue.")
    imagen = np.arange(nims, dtype=int) + 1
    g1true = np.repeat(np.random.randn(ntrue) * true_sigma, nims/ntrue)  # these should have NIMS
    g2true = np.repeat(np.random.randn(ntrue) * true_sigma, nims/ntrue)  # elements, with repetition
    if saveto is not None:
        import os
        if not os.path.isdir(os.path.split(saveto)[0]):
            os.mkdir(os.path.split(saveto)[0])
        np.savetxt(saveto, np.array((imagen, g1true, g2true)).T, fmt=('%d', '%14.7f', '%14.7f'))
    return g1true, g2true

def make_const_truth_uniform_dist(ntrue, nims, true_range=0.03,
                                  saveto='./g3truth/g3_const_shear_truth.dat'):
    """Generate truth catalogues with a U(-TRUE_RANGE, TRUE_RANGE) distribution of input truth
    values in g1 and g2.

    Set saveto=None to disable hardcopy output. nims should be an integer multiple of ntrue.
    """
    if (nims % ntrue) != 0:
        raise ValueError("nims must be divisible by ntrue.")
    if not os.path.isdir('g3truth'):
        os.mkdir('g3truth')
    imagen = np.arange(nims, dtype=int) + 1
    g1true = np.repeat((2. * np.random.rand(ntrue) - 1.) * true_range, nims/ntrue)
    g2true = np.repeat((2. * np.random.rand(ntrue) - 1.) * true_range, nims/ntrue)
    if saveto is not None:
        import os
        if not os.path.isdir(os.path.split(saveto)[0]):
            os.mkdir(os.path.split(saveto)[0])
        np.savetxt(saveto, np.array((imagen, g1true, g2true)).T, fmt=('%d', '%14.7f', '%14.7f'))
    return g1true, g2true

def read_ps(galsim_dir=None):
    """Read in a Power Spectrum stored in the GalSim repository.

    Returns a galsim.PowerSpectrum object.
    """
    import os
    if galsim_dir is None:
        raise ValueError(
            "You must supply a directory for your GalSim install via the `galsim_dir` kwarg.")
    tab_ps = galsim.LookupTable(
        file=os.path.join(galsim_dir, 'examples', 'data', 'cosmo-fid.zmed1.00_smoothed.out'),
        interpolant='linear')
    # Put this table into an E-mode power spectrum and return
    ret = galsim.PowerSpectrum(tab_ps, None, units=galsim.radians)
    return ret

def make_var_truth_catalogs(ntrue, nims, ps_list, ngrid=100, dx_grid=0.1, grid_units=galsim.degrees,
                            rng=None):
    """Generate truth catalogues of g1, g2 in a grid of galaxy locations for each of nims images.

    Inputs
    ------
    ntrue      Number of underlying true power spectra to realize, must be the same as
               `len(ps_list)`.  Note this parameter is mostly retained for syntactical similarity
               with const truth table generation scripts, aiding checking.
    nims       The number of images across which to make realizations of the power spectra in
               `ps_list`. nims should be an integer multiple of ntrue.
    ps_list    List of galsim.PowerSpectrum objects of length `ntrue`, used to select the
               power spectra for use in the realizations of g1, g2.
    ngrid      Number of galaxies along a side of the (assumed square) image grid
    dx_grid    Image grid spacing in units of `grid_units`.
    grid_units Units of grid spacing, must be a galsim.Angle instance.
    rng        A galsim.BaseDeviate instance used to generate random numbers, if None one will be
               initialized internally.

    Returns the tuple `(g1true_list, g2true_list)`, each element of which is a list of `nims` 2D
    NumPy arrays of dimensions ngrid x ngrid, containg the realizations of the shear field g1 and g2
    components at the grid points specified.
    """
    if (nims % ntrue) != 0:
        raise ValueError("nims must be divisible by ntrue.")
    if len(ps_list) != ntrue:
        raise ValueError("Input ps_list does not have ntrue elements.")
    if rng is None:
        rng = galsim.BaseDeviate()
    if not isinstance(rng, galsim.BaseDeviate):
        raise TypeError("Input rng not a galsimBaseDeviate or derived class instance.") 
   # Number of sets of ntrue images
    nsets = nims / ntrue
    # Setup empty lists for storing the g1, g2 NumPy arrays
    g1true_list = []
    g2true_list = []
    # Loop through the power spectra and build the gridded shear realizations for each image
    for ps in ps_list:
        g1, g2 = ps.buildGrid(grid_spacing=dx_grid, ngrid=ngrid, units=grid_units, rng=rng)
        for i in range(nsets):
            g1true_list.append(g1)
            g2true_list.append(g2)
    return g1true_list, g2true_list

def make_submission_var_shear_PS(c1, c2, m1, m2, g1true_list, g2true_list, noise_sigma, dx_grid=0.1,
                                 nbins=8, label=None, calculate_truth=True):
    """Make a fake var shear submission in the form of a Power Spectrum.

    BARNEY NOTE: In the real data we should maybe do this in the (x, y) coordinate frame determined
    by the primary direction of the PSF ellipticity.

    Arguments
    ---------
    * Provided c1, c2, m1, m2 shear estimation bias values.
    * Two lists of truth tables g1true_grid_list, g2true_grid_list, list of 2D NumPy arrays
      containing the variable shears at each grid point.
    * Image grid dx_grid spacing in units of degrees.

    Outputs a set of tables of k, P_E(k),... etc.
    """
    # Get the number of images from the truth table
    nims = len(g1true_list)
    if len(g2true_list) != nims:
        raise ValueError("Supplied g1true, g2true not matching length.")
    # Then ready an empty list (will store arrays) for the output submission
    if calculate_truth:
        pEtrue_list = []
        pBtrue_list = []
    pEsub_list = []
    pBsub_list = []
    for i in range(nims):
        g1gals = (1. + m1) * g1true_list[i] + c1 + noise_sigma * np.random.randn(
            *g1true_list[i].shape) # pleasantly magic asterisk *args functionality
        g2gals = (1. + m2) * g2true_list[i] + c2 + noise_sigma * np.random.randn(
            *g2true_list[i].shape)
        # Setup the power spectrum estimator
        pse = galsim.pse.PowerSpectrumEstimator(
            g1gals.shape[0], dx_grid * g1gals.shape[0], nbins)
        # Then estimate the true signals if asked, and the noisy submission ones
        if calculate_truth:
            k, pEtrue_tmp, pBtrue_tmp, pEBtrue_tmp = pse.estimate(
                g1true_list[i], g2true_list[i])
            pEtrue_list.append(pEtrue_tmp)
            pBtrue_list.append(pBtrue_tmp)
        ell, pEsub_tmp, pBsub_tmp, pEBsub_tmp = pse.estimate(g1gals, g2gals)
        pEsub_list.append(pEsub_tmp)
        pBsub_list.append(pBsub_tmp)
    if calculate_truth:
        ret = k, pEsub_list, pBsub_list, pEtrue_list, pBtrue_list
    else:
        ret = k, pEsub_list, pBsub_list
    return ret

def make_submission_var_shear_CF(c1, c2, m1, m2, g1true_list, g2true_list, noise_sigma,
                                 dx_grid=0.1, nbins=8, label=None, calculate_truth=True,
                                 min_sep=0.1, max_sep=10., verbose=False):
    """Make a fake var shear submission in a Correlation Function.

    BARNEY NOTE: In the real data we should maybe do this in the (x, y) coordinate frame determined
    by the primary direction of the PSF ellipticity.

    Arguments
    ---------
    * Provided c1, c2, m1, m2 shear estimation bias values.
    * Two lists of truth tables g1true_grid_list, g2true_grid_list, list of 2D NumPy arrays
      containing the variable shears at each grid point.
    * Image grid dx_grid spacing in units of degrees.

    Outputs a set of tables of
        theta, ApMass_E(theta), Ap_Mass_B(theta), Ap_Mass_err(theta), ...
    """
    # Get the number of images from the truth table
    nims = len(g1true_list)
    if len(g2true_list) != nims:
        raise ValueError("Supplied g1true, g2true not matching length.")
    # Test for the size of a (square assumed) grid
    ngrid = g1true_list[0].shape[0]
    if (g1true_list[0].shape[1] != ngrid or g2true_list[0].shape[0] != ngrid or
        g2true_list[1].shape[1] != ngrid):
        raise ValueError("Input g1 and g2 true should be square and the same shape!")
    # Then ready an empty list (will store arrays) for the output submission
    if calculate_truth:
        mEtrue_list = []
        mBtrue_list = []
    mEsub_list = []
    mBsub_list = []
    merrsub_list = []
    for i in range(nims):
        if verbose: print "Generating aperture mass submission for image "+str(i + 1)+"/"+str(nims)
        # Build the x, y grid
        x, y = np.meshgrid(np.arange(ngrid) * dx_grid, np.arange(ngrid) * dx_grid)
        g1gals = (1. + m1) * g1true_list[i] + c1 + noise_sigma * np.random.randn(
            *g1true_list[i].shape) # pleasantly magic asterisk *args functionality
        g2gals = (1. + m2) * g2true_list[i] + c2 + noise_sigma * np.random.randn(
            *g2true_list[i].shape)
        # Then estimate the true signals if asked, and the noisy submission ones
        if calculate_truth:
            results_truth = run_corr2_ascii(
                x, y, g1true_list[i], g2true_list[i], min_sep=min_sep, max_sep=max_sep, nbins=nbins,
                temp_cat='temp.cat', params_file='corr2.params', m2_file_name='temp.m2',
                xy_units='degrees', sep_units='degrees')
            mEtrue_list.append(results_truth[:, 1])
            mBtrue_list.append(results_truth[:, 2])
        results = run_corr2_ascii(
                x, y, g1gals, g2gals, min_sep=min_sep, max_sep=max_sep, nbins=nbins,
                temp_cat='temp.cat', params_file='corr2.params', m2_file_name='temp.m2',
                xy_units='degrees', sep_units='degrees')
        theta = results[:, 0]
        mEsub_list.append(results[:, 1])
        mBsub_list.append(results[:, 2])
        merrsub_list.append(results[:, 4])
    if calculate_truth:
        ret = theta, mEsub_list, mBsub_list, merrsub_list, mEtrue_list, mBtrue_list
    else:
        ret = theta, mEsub_list, mBsub_list, merrsub_list
    return ret

def make_submission_const_shear(c1, c2, m1, m2, g1true, g2true, ngals_per_im, noise_sigma,
                                label=None):
    """Make a fake const shear submission.

    BARNEY NOTE: In the real data we should do this in the (x, y) coordinate frame determined
    by the primary direction of the PSF ellipticity.

    Arguments
    ---------
    * Provided c1, c2, m1, m2 shear estimation bias values
    * Truth tables g1true, g2true
    * The number of galaxies per image, ngals_per_im, represented by each truth table element
    * The noise_sigma standard deviation on the shear estimate for each galaxy due to noise

    Saves output to ./g3subs/g3_const_shear_sub.<label>.dat if label is not `None`
    """
    # Get the number of images from the truth table
    nims = len(g1true)
    if len(g2true) != nims:
        raise ValueError("Supplied g1true, g2true not matching length.")

    # Then ready some arrays for the output submission
    g1sub = np.empty(nims)
    g2sub = np.empty(nims)

    # Loop generating subs (doing this the long way - could use central limit theorem but this is
    # super safe!)
    for i in range(nims):
        g1gals = (1. + m1) * g1true[i] + c1 + np.random.randn(ngals_per_im) * noise_sigma
        g2gals = (1. + m2) * g2true[i] + c2 + np.random.randn(ngals_per_im) * noise_sigma
        g1sub[i] = np.mean(g1gals)
        g2sub[i] = np.mean(g2gals)
    # Save output if requested
    if label is not None:
        import os
        if not os.path.isdir('g3subs'):
            os.mkdir('g3subs')
        np.savetxt(
            './g3subs/g3_const_shear_sub.'+label+'.dat', np.array((truth[:, 0], g1sub, g2sub)).T,
            fmt=('%d', '%14.7f', '%14.7f'))
    return g1sub, g2sub

def _calculateSvalues(xarr, yarr, sigma2=1.):
    """Calculates the intermediate S values required for basic linear regression.

    See, e.g., Numerical Recipes (Press et al 1992) Section 15.2.
    """
    if len(xarr) != len(yarr):
        raise ValueError("Input xarr and yarr differ in length!")
    if len(xarr) <= 1:
        raise ValueError("Input arrays must have 2 or more values elements.")

    S = len(xarr) / sigma2
    Sx = np.sum(xarr / sigma2)
    Sy = np.sum(yarr / sigma2)
    Sxx = np.sum(xarr * xarr / sigma2)
    Sxy = np.sum(xarr * yarr / sigma2)
    return (S, Sx, Sy, Sxx, Sxy)

def run_corr2_ascii(x, y, e1, e2, min_sep=0.1, max_sep=10., nbins=8, temp_cat='temp.cat',
                    params_file='corr2.params', m2_file_name='temp.m2', xy_units='degrees',
                    sep_units='degrees'):
    import os
    import subprocess
    import tempfile
    # Create temporary, unique files for I/O
    catfile = tempfile.mktemp(suffix=temp_cat)
    m2file = tempfile.mktemp(suffix=m2_file_name)
    f = open(catfile, 'wb')
    for (i, j), value in np.ndenumerate(x):
        xij = value
        yij = y[i,j]
        g1ij = e1[i,j]
        g2ij = e2[i,j]
        f.write('%e  %e  %e  %e\n'%(xij, yij, g1ij, g2ij))
    f.close()
    subprocess.Popen([
        'corr2', params_file, 'file_name='+str(catfile), 'm2_file_name='+str(m2file),
        'min_sep=%f'%min_sep, 'max_sep=%f'%max_sep, 'nbins=%f'%nbins,
        'x_units='+str(xy_units), 'y_units='+str(xy_units), 'sep_units='+str(sep_units)]).wait()
    results = np.loadtxt(m2file)
    os.remove(catfile)
    os.remove(m2file)
    return results

def fitline(xarr, yarr):
    """Fit a line y = a + b * x to input x and y arrays by least squares.

    Returns the tuple (a, b, Var(a), Cov(a, b), Var(b)), after performing an internal estimate of
    measurement errors from the best-fitting model residuals.

    See Numerical Recipes (Press et al 1992; Section 15.2) for a clear description of the details
    of this simple regression.
    """
    # Get the S values (use default sigma2, best fit a and b still valid for stationary data)
    S, Sx, Sy, Sxx, Sxy = _calculateSvalues(xarr, yarr)

    # Get the best fit a and b
    Del = S * Sxx - Sx * Sx
    a = (Sxx * Sy - Sx * Sxy) / Del
    b = (S * Sxy - Sx * Sy) / Del

    # Use these to estimate the sigma^2 by residuals from the best-fitting model
    ymodel = a + b * xarr
    sigma2 = np.mean((yarr - ymodel)**2)
    # And use this to get model parameter error estimates
    var_a  = sigma2 * Sxx / Del
    cov_ab = - sigma2 * Sx / Del
    var_b  = sigma2 * S / Del
    return (a, b, var_a, cov_ab, var_b)

def metricQ08_const_shear(g1est, g2est, g1true, g2true, ntruesets):
    """Calculate a GREAT08-style Q value for constant shear branch results.

    Each element of the input arrays g1est, g2est is assumed to be the average of the shear
    estimates in a large image, of which there will be some decent-sized number (~200).   The arrays
    g1true, g2true are the corresponding truth values.  There are ntruesets independent truth values
    in each array.
    """
    if len(g1est) % ntruesets != 0:
        raise ValueError("Number of separate truth values ntrue is not a divisor for the number "+
                         "of input array elements without remainder.")
    nset = len(g1est) / ntruesets  # Number of elements in each set of distinct truth values
    g1set = []
    g2set = []
    for i in range(ntruesets):
        g1set.append(np.mean((g1est - g1true)[i * nset: (i + 1) * nset]))
        g2set.append(np.mean((g2est - g2true)[i * nset: (i + 1) * nset]))
    return 1.e-4 / (0.5 * ((sum(g1set) / float(ntruesets))**2 + (sum(g2set) / float(ntruesets))**2))

def metricQZ1_const_shear(g1est, g2est, g1true, g2true, cfid=1.e-4, mfid=1.e-3):
    """Calculate a metric along the lines suggested by Joe Zuntz in Pittsburgh (option 1).
    """
    c1, m1, var_c1, cov_c1m1, var_m1 = fitline(g1true, g1est - g1true)
    c2, m2, var_c2, cov_c2m2, var_m2 = fitline(g2true, g2est - g2true)
    sig_c1 = np.sqrt(var_c1)
    sig_m1 = np.sqrt(var_m1)
    sig_c2 = np.sqrt(var_c2)
    sig_m2 = np.sqrt(var_m2)
    Q = 2000. / np.sqrt((c1 / cfid)**2 + (c2 / cfid)**2 + (m1 / mfid)**2 + (m2 / mfid)**2)
    return (Q, c1, m1, c2, m2, sig_c1, sig_m1, sig_c2, sig_m2)

def metricQZ2_const_shear(g1est, g2est, g1true, g2true, cfid=1.e-4, mfid=1.e-3):
    """Calculate a metric along the lines suggested by Joe Zuntz in Pittsburgh (option 2).
    """
    c1, m1, var_c1, cov_c1m1, var_m1 = fitline(g1true, g1est - g1true)
    c2, m2, var_c2, cov_c2m2, var_m2 = fitline(g2true, g2est - g2true)
    sig_c1 = np.sqrt(var_c1)
    sig_m1 = np.sqrt(var_m1)
    sig_c2 = np.sqrt(var_c2)
    sig_m2 = np.sqrt(var_m2)
    Q = 4000. / (np.abs(c1 / cfid) + np.abs(c2 / cfid) + np.abs(m1 / mfid) + np.abs(m2 / mfid))
    return (Q, c1, m1, c2, m2, sig_c1, sig_m1, sig_c2, sig_m2)

def metricQZ3_const_shear(g1est, g2est, g1true, g2true, cfid=1.e-4, mfid=1.e-3):
    """Calculate a metric along the lines suggested by Joe Zuntz in Pittsburgh (option 3).
    """
    c1, m1, var_c1, cov_c1m1, var_m1 = fitline(g1true, g1est - g1true)
    c2, m2, var_c2, cov_c2m2, var_m2 = fitline(g2true, g2est - g2true)
    sig_c1 = np.sqrt(var_c1)
    sig_m1 = np.sqrt(var_m1)
    sig_c2 = np.sqrt(var_c2)
    sig_m2 = np.sqrt(var_m2)
    Q = 500. * np.sqrt((cfid / c1)**2 + (cfid / c2)**2 + (mfid / m1)**2 + (mfid / m2)**2)
    return (Q, c1, m1, c2, m2, sig_c1, sig_m1, sig_c2, sig_m2)

def metricG10_var_shear(k, pEsub_list, varsub_list, pEtrue_list, scaling=0.001, dx_grid=0.1):
    """Crude attempt at coding up the G10 metric, with variance subtraction.
    """
    mean_pEsub = pEsub_list[0] - (varsub_list[0] * (dx_grid * np.pi / 180.)**2) # See below...
    mean_pEtrue = pEtrue_list[0]
    mean_diff = mean_pEsub - pEtrue_list[0]
    for pEsub, varsub, pEtrue in zip(pEsub_list[1:], varsub_list[1:], pEtrue_list[1:]):
        varsub *= (dx_grid * np.pi / 180.)**2 # Convert the variance per grid point into radian^2
        mean_pEsub += (pEsub - varsub)
        mean_pEtrue += pEtrue 
        mean_diff += (pEsub - varsub - pEtrue)
    mean_pEsub /= len(pEsub_list)
    mean_pEtrue /= len(pEtrue_list)
    mean_diff /= len(pEsub_list)
    I_tilde_over_k = np.abs(mean_diff) * k # comes from C11
    Q = scaling / np.sum(I_tilde_over_k)
    return Q, mean_pEsub, mean_pEtrue, mean_diff

def metricQuadPS_var_shear(k, pEsub_list, varsub_list, pEtrue_list, scaling=0.001, dx_grid=0.1):
    """Crude attempt at coding up a similar-to-G10 metric, but using a sum in quadrature, with
    variance subtraction.
    """
    mean_pEsub = pEsub_list[0] - (varsub_list[0] * (dx_grid * np.pi / 180.)**2) # See below...
    mean_pEtrue = pEtrue_list[0]
    mean_diff = mean_pEsub - pEtrue_list[0]
    for pEsub, varsub, pEtrue in zip(pEsub_list[1:], varsub_list[1:], pEtrue_list[1:]):
        varsub *= (dx_grid * np.pi / 180.)**2 # Convert the variance per grid point into radian^2
        mean_pEsub += (pEsub - varsub)
        mean_pEtrue += pEtrue 
        mean_diff += (pEsub - varsub - pEtrue)
    mean_pEsub /= len(pEsub_list)
    mean_pEtrue /= len(pEtrue_list)
    mean_diff /= len(pEsub_list)
    I_tilde_over_k = (mean_diff) * k # comes from C11
    Q = scaling / np.sqrt(np.sum(I_tilde_over_k**2))
    return Q, mean_pEsub, mean_pEtrue, mean_diff

def calculate_map_unitc(ngrid=100, dx_grid=0.1, nbins=8, min_sep=0.1, max_sep=10., plotfile=None):
    """Calculate the aperture mass statistic for this geometry due a constant input ellipticity
    c1=c2=1.
    """
    xygrid = np.arange(ngrid) * dx_grid
    x, y = np.meshgrid(xygrid, xygrid)
    e1unitc = np.ones_like(x)
    e2unitc = np.ones_like(x)
    results_unitc = run_corr2_ascii(
        x, y, e1unitc, e2unitc, min_sep=min_sep, max_sep=max_sep, nbins=nbins)
    if plotfile is not None:
        import matplotlib.pyplot as plt
        plt.plot(results_unitc[:, 0], results_unitc[:, 1])
        plt.plot(results_unitc[:, 0], results_unitc[:, 2])
        plt.xlabel('theta [degrees]')
        plt.ylabel(r'E mode <M$_{ap}$$^2$> for unit c$_1$=c$_2$=1')
        print "Saving plot of unit c1=c2=1 <Map^2> to "+plotfile
        plt.savefig(plotfile)
    return results_unitc[:, 0], results_unitc[:, 1], results_unitc[:, 2]

def map_squared_diff_func(cm_array, mapEsub, maperrsub, mapEtrue, mapEunitc):
    """Squared difference of an m-c model of the aperture mass statistic and the submission.
    """
    retval = ((
        mapEsub - (mapEunitc * cm_array[0] + mapEtrue * (1. + 2. * cm_array[1]))
        ) / maperrsub)**2
    return retval

def metricMapCF_var_shear_mc(mapEsub_list, maperrsub_list, mapEtrue_list, mapBtrue_list, ntruesets,
                             ngrid=100, dx_grid=0.1, nbins=8, cfid=1.e-4, mfid=1.e-3, min_sep=0.1,
                             max_sep=10., plot=False, select_by_B_leakage=0.,
                             correct_B_theory=True, use_errors=False):
    """The ntruesets must be an integer divisor of len(mapEsub_list)
    """
    import scipy.optimize
    # First calculate what an input unit c1=c2=1 looks like
    theta, mapE_unitc, mapB_unitc = calculate_map_unitc(
        ngrid=ngrid, dx_grid=dx_grid, nbins=nbins, min_sep=min_sep, max_sep=max_sep)
    # Calculate the number of images per set of realizations (trueset)
    nperset = len(mapEsub_list) / ntruesets
    c2s = []
    ms = []
    for iset in range(ntruesets):
        mapEsub_mean = mapEsub_list[iset * nperset]
        maperrsub_mean = maperrsub_list[iset * nperset]
        mapEtrue_mean = mapEtrue_list[iset * nperset]
        mapBtrue_mean = mapBtrue_list[iset * nperset]
        for jimage in range(nperset)[1:]:
            mapEsub_mean += mapEsub_list[iset * nperset + jimage]
            maperrsub_mean += maperrsub_list[iset * nperset + jimage]
            mapEtrue_mean += mapEtrue_list[iset * nperset + jimage]
            mapBtrue_mean += mapBtrue_list[iset * nperset + jimage]
        # Divide by nperset to get the mean mapE etc.
        mapEsub_mean /= float(nperset)
        maperrsub_mean /= (float(nperset) * np.sqrt(ntruesets))
        mapEtrue_mean /= float(nperset)
        mapBtrue_mean /= float(nperset)
        # Optionally only select the regions where |true B| < select_by_B_leakage * |true E|
        if select_by_B_leakage > 0.:
            use_for_fit = (np.abs(mapBtrue_mean) / np.abs(mapEtrue_mean) < select_by_B_leakage)
        else:
            use_for_fit = np.array([True,] * len(mapEsub_mean), dtype=bool)

        # Ready the terms to go to the fitter
        submission = mapEsub_mean[use_for_fit]
        if use_errors:
            errors = maperrsub_mean[use_for_fit]
        else:
            errors = np.ones_like(submission)
        truth = mapEtrue_mean[use_for_fit]
        unit_contamination = mapE_unitc[use_for_fit]
        if correct_B_theory:
            submission -= mapBtrue_mean[use_for_fit]
            truth -= mapBtrue_mean[use_for_fit]
            unit_contamination -= mapB_unitc[use_for_fit]
            #import pdb; pdb.set_trace()
        # Use scipy.optimize to fit m and c (note should change this to numpy.leastsq since the
        # model can be made approximately linear 
        results = scipy.optimize.leastsq(
            map_squared_diff_func, np.array([0., 0.]),
            args=(
                submission,
                errors,
                truth,
                unit_contamination))
        c2 = results[0][0]
        m = results[0][1]
        if plot:
            import os
            import matplotlib.pyplot as plt
            if not os.path.isdir('plots'): os.mkdir('plots')
            plt.errorbar(
                theta, submission, yerr=maperrsub_mean, label='Map submission', color='r')
            plt.plot(theta, truth , 'g--', label='Map truth realizations')
            plt.plot(
                theta, truth * (1. + 2. * m) + unit_contamination * c2, 'b',
                label='Best fitting linear model')
            plt.legend()
            plt.title(
                r'Best fitting m = '+str(m)+',  c$^2$ = '+str(c2)+' \n'+
                'Set '+str(iset + 1)+'/'+str(ntruesets)+' ('+str(nperset)+' images)')
            plt.xlabel('theta [degrees]')
            plt.xscale('log')
            plt.ylabel('E-mode Map')
            plt.savefig(os.path.join(
                'plots', 'aperture_mass_metric_set'+str(iset+1)+'of'+str(ntruesets)+'.png'))
            plt.show()
        c2s.append(c2)
        ms.append(m)
    c2 = np.mean(np.array(c2s))
    m = np.mean(np.array(ms))
    Q = np.sqrt(2.) * 1000. / np.sqrt(np.abs(c2 / cfid**2) + (m / mfid)**2)
    return Q, c2, m

def calc_diffs_E(mapEsub_list, maperrsub_list, mapEtrue_list, ntruesets, use_errors=False,
                 select_by_B_leakage=0.):
    """Calculate differences between all submitted and target E-mode aperture mass statistics,
    averaged over fields, optionally normalized by errors and excluding large B leakage regions.
    """
    # Calculate the number of images per set of realizations (trueset)
    nperset = len(mapEsub_list) / ntruesets
    diffs = np.zeros((len(mapEsub_list[0]), ntruesets)) # Assumes all elements of map lists same
                                                        # length, but an exception will be thrown
                                                        # later if this is not so...
    for iset in range(ntruesets):
        mapEsub_mean = mapEsub_list[iset * nperset]
        maperrsub_mean = maperrsub_list[iset * nperset]
        mapEtrue_mean = mapEtrue_list[iset * nperset]
        mapBtrue_mean = mapBtrue_list[iset * nperset]
        for jimage in range(nperset)[1:]:
            mapEsub_mean += mapEsub_list[iset * nperset + jimage]
            maperrsub_mean += maperrsub_list[iset * nperset + jimage]
            mapEtrue_mean += mapEtrue_list[iset * nperset + jimage]
            mapBtrue_mean += mapBtrue_list[iset * nperset + jimage]
        # Divide by nperset to get the mean mapE etc.
        mapEsub_mean /= float(nperset)
        maperrsub_mean /= (float(nperset) * np.sqrt(ntruesets))
        mapEtrue_mean /= float(nperset)
        mapBtrue_mean /= float(nperset)
        # Optionally only select the regions where |true B| < select_by_B_leakage * |true E|
        if select_by_B_leakage > 0.:
            use_for_fit = (np.abs(mapBtrue_mean) / np.abs(mapEtrue_mean) < select_by_B_leakage)
        else:
            use_for_fit = np.array([True,] * len(mapEsub_mean), dtype=bool)

        # Ready the terms to use for calculating the dispersion statistic
        submission = mapEsub_mean[use_for_fit]
        if use_errors:
            errors = maperrsub_mean[use_for_fit]
        else:
            errors = np.ones_like(submission)
        truth = mapEtrue_mean[use_for_fit]
        diffs[:, iset] = (submission - truth) / errors
    # And return...
    return diffs

def metricG3S_AMD(mapEsub_list, maperrsub_list, mapEtrue_list, mapBtrue_list, ntruesets,
                  plot=False, select_by_B_leakage=0., use_errors=False, normalization=1000.):
    """A metric based on the absolute differences between a submitted E-mode aperture mass (Map)
    statistic and the target.
    """
    diffs = calc_diffs_E(
        mapEsub_list, maperrsub_list, mapEtrue_list, ntruesets, use_errors=use_errors,
        select_by_B_leakage=select_by_B_leakage)
    mean_error = np.mean(np.asarray(maperrsub_list))
    return normalization * mean_error / np.mean(np.abs(diffs))

def metricG3S_QMD(mapEsub_list, maperrsub_list, mapEtrue_list, mapBtrue_list, ntruesets,
                  plot=False, select_by_B_leakage=0., use_errors=False, normalization=1000.):
    """A metric based on summed (in quadrature) differences between a submitted E-mode aperture
    mass (Map) statistic and the target.
    """
    diffs = calc_diffs_E(
        mapEsub_list, maperrsub_list, mapEtrue_list, ntruesets, use_errors=use_errors,
        select_by_B_leakage=select_by_B_leakage)
    mean_error = np.mean(np.asarray(maperrsub_list))
    return normalization * mean_error / np.sqrt(np.mean(diffs**2))
