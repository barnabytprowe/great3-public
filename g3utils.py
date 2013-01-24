"""Module containing GREAT3 metric calculation utilities.
"""

import numpy as np

def make_truth_normal_dist(ntrue, nims, true_sigma=0.03,
                           saveto='./g3truth/g3_const_shear_truth.dat'):
    """Generate truth catalogues with a N(0, TRUE_SIGMA) distribution of input truth values in
    g1 and g2.

    Set saveto=None to disable hardcopy output.
    """
    imagen = np.arange(nims, dtype=int) + 1
    g1true = np.repeat(np.random.randn(ntrue) * true_sigma, nims/ntrue)  # these should have NIMS
    g2true = np.repeat(np.random.randn(ntrue) * true_sigma, nims/ntrue)  # elements, with repetition
    if saveto is not None:
        import os
        if not os.path.isdir(os.path.split(saveto)[0]):
            os.mkdir(os.path.split(saveto)[0])
        np.savetxt(saveto, np.array((imagen, g1true, g2true)).T, fmt=('%d', '%14.7f', '%14.7f'))
    return g1true, g2true

def make_truth_uniform_dist(ntrue, nims, true_range=0.03,
                           saveto='./g3truth/g3_const_shear_truth.dat'):
    """Generate truth catalogues with a U(-TRUE_RANGE, TRUE_RANGE) distribution of input truth
    values in g1 and g2.

    Set saveto=None to disable hardcopy output.
    """
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

def make_submission_const_shear(c1, c2, m1, m2, g1true, g2true, ngals_per_im, noise_sigma,
                                label=None):
    """Make a fake const shear submission.

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
    Q = 500. * np.sqrt((cfid / c1)**2 + (cfid / c2)**2 + (mfid / m1)**2 + (mfid / m2)**2)
    return (Q, c1, m1, c2, m2, sig_c1, sig_m1, sig_c2, sig_m2)



