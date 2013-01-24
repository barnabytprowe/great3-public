"""Module containing GREAT3 metric calculation utilities.
"""

import numpy as np

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

def metricQ08_const_shear(g1est, g2est, g1true, g2true):
    """Calculate a GREAT08-style Q value for constant shear branch results.

    Each element of the input arrays g1est, g2est is assumed to be the average of the shear
    estimates in a large image, of which there will be some decent-sized number (~200).   The
    arrays g1true, g2true are the corresponding truth values.
    """
    return 1.4 / np.sqrt(np.mean((g1est - g1true)**2 + (g2est - g2true)**2))

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
