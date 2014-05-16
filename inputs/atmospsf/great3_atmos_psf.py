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
"""This script randomly draws an atmospheric PSF anisotropy field using the outputs in pk_math/, and
makes some whisker plots and other sanity checks.  Some of this code was later used in the GREAT3
simulation scripts.
"""
import galsim
import numpy as np
import pylab

def atmosPSFParams(t_exp=20, diam=6.7, rng=None, N=1):
    """A subroutine to generate random parameters for the atmospheric PSF.

    This routine generates a random value of seeing from some distribution, and random values for
    the two parameters (A, theta_0) that determine the correlation function of PSF anisotropies:

        xi_+(theta) = A / (1 + theta/theta_0)

    These values are then rescaled from the fiducial values of exposure time t_exp=20s and
    telescope diameter diam=6.7m to the input values of t_exp and diam, where xi is proportional
    to 1/t_exp and 1/diam.

    The inputs are exposure time in seconds and telescope diameter in meters.  Optionally, the user
    can input a galsim.BaseDeviate (rng) to use as the pseudo-random number generator, and can use
    the keyword N to instruct the routine to return more than one set of parameters (in which case
    the return values are NumPy arrays instead of floats).

    The numbers that are returned are A (dimensionless), theta_0 (arcsec), and seeing FWHM (arcsec),
    either as three individual floats (if N=1) or as three NumPy arrays of length N (if N>1).
    """
    # Check for sane ranges of exposure time and diameter.
    if t_exp < 0:
        raise ValueError("Invalid input for exposure time!")
    if diam < 0:
        raise ValueError("Invalid input for telescope diameter!")
    if N <= 0:
        raise ValueError("Invalid input for N!")

    # Check for input RNG or make a new one.
    if rng is None:
        uniform_deviate = galsim.UniformDeviate()
    elif isinstance(rng, galsim.BaseDeviate):
        # If it's a BaseDeviate, we can convert to UniformDeviate
        uniform_deviate = galsim.UniformDeviate(rng)
    else:
        raise TypeError("The rng provided to atmosPSFParams is not a BaseDeviate")

    # Distribution for average seeing comes from the distribution for Mauna Kea, as read off Fig 1
    # in
    # http://www2.keck.hawaii.edu/optics/ScienceCase/TechSciInstrmnts/Products_SeeingVarMaunaKea.pdf
    # (R band, at zenith).  We use the distribution for the central part of the night which is an
    # intermediate case compared to the others.  Currently there is no option for correction for
    # observations that are not at zenith or not in R band.
    fwhm_arcsec = 0.05 + 0.10*np.arange(17)
    freq = (0., 0., 0., 7.5, 19., 20., 17., 13., 9., 5., 3.5, 2., 1., 1., 0.5, 0.0, 0.0)
    dd = galsim.DistDeviate(uniform_deviate, 
                            function=galsim.LookupTable(fwhm_arcsec, freq, interpolant='linear'))
    if N > 1:
        rand_seeing = np.zeros(N)
        for i in range(N):
            rand_seeing[i] = dd()
    else:
        rand_seeing = dd()

    # Now work on the PSF anisotropies:
    # For A, take a range based on ImSim for 20s exposures: 1e-4 to 8e-4 (could debate about what
    # distribution to use, but for now, flat isn't too bad).  This has to be rescaled according to
    # the exposure time and diameter.
    min_A = 1.e-4
    max_A = 8.e-4
    if N > 1:
        rand_A = np.zeros(N)
        for i in range(N):
            rand_A[i] = min_A + (max_A-min_A)*uniform_deviate()
    else:
        rand_A = min_A + (max_A-min_A)*uniform_deviate()
    rand_A *= (20./t_exp)*(6.7/diam)

    # For theta_0, take a range based on ImSim, from 0.1-2 degree (again, flat distribution).
    min_theta_0 = 0.1
    max_theta_0 = 2.0
    if N > 1:
        rand_theta_0 = np.zeros(N)
        for i in range(N):
            rand_theta_0[i] = min_theta_0 + (max_theta_0-min_theta_0)*uniform_deviate()
    else:
        rand_theta_0 = min_theta_0 + (max_theta_0-min_theta_0)*uniform_deviate()
    rand_theta_0 *= 3600. # need arcsec

    return rand_A, rand_theta_0, rand_seeing

def makeAtmosPSFPk(A, theta0, interpolate=False):
    """A routine to reconstruct the atmospheric PSF given the input parameters.

    This routine requires an amplitude and a scale length for the atmospheric PSF ellipticity
    correlation function xi_+(theta), defined as

        xi_+(theta) = A / (1 + theta/theta_0)

    The input theta0 must be in arcsec.  The function returns two NumPy arrays containing k
    [arcsec^-1] and P(k) [arcsec^2], respectively.  The P(k) that is returned corresponds to the one
    that should be given to the GalSim lensing engine for both E and B power.

    Currently, we only have P(k) for quantized values of theta_0 from 0.1-2 degrees (11 values,
    logarithmically spaced).  For now, we simply round to the nearest one in log space and use
    that.  There is a flag to request interpolation, but that is not currently implemented.  Seems
    like overkill but could be done if people think it's important.
    """
    if interpolate:
        raise NotImplementedError("Sorry, cannot yet interpolate between gridded theta_0 values")
    if A <= 0 or theta0 <= 0:
        raise ValueError("Input A, theta0 may not be negative or zero!")

    # Set up our grid in theta0, which is defined based on how we did the mathematica calculations.
    min_theta0 = 360. # 0.1 degree
    max_theta0 = 7200. # 2 degrees
    if theta0 < min_theta0 or theta0 > max_theta0:
        raise ValueError("Input theta0 is outside of the range of values with tabulated P(k)!")
    log_min_theta0 = np.log10(min_theta0)
    log_max_theta0 = np.log10(max_theta0)
    ntheta0 = 11
    d_log_theta0 = (log_max_theta0 - log_min_theta0) / (ntheta0 - 1)
    theta0_grid = np.logspace(np.log10(min_theta0),np.log10(max_theta0),ntheta0)

    # Figure out where our theta0 value lies on the grid
    input_log_theta0 = np.log10(theta_0)
    theta0_index = int( round( (input_log_theta0 - log_min_theta0) / d_log_theta0) )
    strval = str(int(theta0_grid[theta0_index]))
    infile = 'pk_math/Pk'+strval+'.dat'
    pk_dat = np.loadtxt(infile).transpose()
    k = pk_dat[0]
    # Normalize by amplitude A, by factor of (2pi) in front of integral to get P(k), and 0.5 since
    # the power has to be split between E and B.
    pk = A*np.pi*pk_dat[1]
    return k, pk

def getAtmosPSFGrid(k, Pk, ngrid=20, dtheta_arcsec=360., kmin_factor=15,
                    subsample=10, rng=None, return_basic=True):
    """A routine to build an anisotropy and size fluctuation grid for the atmospheric PSF P(k).

    This routine takes NumPy arrays of k and P(k) values to be used for the E and B power for the
    atmospheric PSF (same P(k) for E and B), and returns grids of atmospheric e1, e2, and fractional
    fluctuations in size.

    Note that the input power spectra are for the ellipticity (distortion) fluctuations.  The
    lensing engine works in terms of shear, and for nearly round objects, shear~distortion/2.  Since
    we gave the lensing engine the PS of distortion fluctuations, the return values can be used
    directly for PSF ellipticity (distortion), but the kappa values are too high by a factor of 2.

    The routine requires k and Pk, and has a number of optional parameters that default to what will
    be used for GREAT3:

       ngrid           Number of grid points for galaxies in one dimension (default 20, for the 2x2
                       degree subfields used for PSF estimation).
       dtheta_arcsec   Spacing between grid points for the galaxies in arcsec (default 360, i.e.,
                       0.1 degree).
       kmin_factor     Factor by which to spatially extend the grids to get a smaller kmin value,
                       growing them larger by this factor in each dimension so as to properly
                       represent the large-scale shear correlations (default 20).
       subsample          Factor by which to subsample the grid, so as to get several offset
                          realizations of the same atmosphere (default 10).  [Note, this is like
                          kmax_factor for the lensing engine work on GalSim issue #377, but there,
                          the idea is to represent smaller kmax without actually getting back a more
                          densely packed grid.  Here, we actually do want a more densely packed
                          grid, so I'm calling it subsample to suggest the explicit subsampling of
                          the grid.
       rng                RNG to use for the generation of these fields.
       return_basic       Return the basic grid that is not enlarged by kmin_factor (if True), or
                          actually return the huge grid (if False).  [Default = True]

    The results are returned as three NumPy arrays for the PSF e1, PSF e2, and fractional change in
    size of the PSF.
    """
    # Set up the PowerSpectrum object.
    tab_pk = galsim.LookupTable(k, Pk, x_log=True, f_log=True)
    ps = galsim.PowerSpectrum(tab_pk, tab_pk, units=galsim.arcsec)

    # Use buildGrid() to get the grid.  Note that with the code in GalSim issue #377, this code could
    # be simplified.  It automatically takes care of the `kmin_factor` expansion of the grid.  This
    # would also lead to simplification of the code below where the return values are selected based
    # on return_basic.
    e1, e2, kappa = ps.buildGrid(grid_spacing = dtheta_arcsec/subsample,
                                 ngrid = ngrid*kmin_factor*subsample,
                                 get_convergence = True,
                                 rng = rng)
    # Redefine the kappa's
    kappa /= 2.

    # Take subset of the overly huge grid.  Since the grid is periodic, can just take one corner, don't
    # have to try to be in the middle.
    if return_basic is False:
        # make grid positions
        ntot = ngrid*kmin_factor*subsample
        grid_spacing = dtheta_arcsec/subsample/60.
        min = (-ntot/2 + 0.5) * grid_spacing
        max = (ntot/2 - 0.5) * grid_spacing
        x, y = np.meshgrid(np.arange(min,max+grid_spacing,grid_spacing),
                           np.arange(min,max+grid_spacing,grid_spacing))

        return e1, e2, kappa, x, y
    else:
        n_use = ngrid*subsample
        grid_spacing = dtheta_arcsec/subsample/60.
        min = (-n_use/2 + 0.5) * grid_spacing
        max = (n_use/2 - 0.5) * grid_spacing
        x, y = np.meshgrid(np.arange(min,max+grid_spacing,grid_spacing),
                           np.arange(min,max+grid_spacing,grid_spacing))
        return e1[0:n_use, 0:n_use], e2[0:n_use, 0:n_use], kappa[0:n_use, 0:n_use], x-min, y-min

def drawWhiskerShear(g1, g2, x, y, title=None):
    """Draw shear whisker plots from an array of g1, g2 values.
    """
    g = (g1**2 + g2**2)**0.5
    theta = 0.5*np.arctan2(g2, g1)
    gx = g*np.cos(theta)
    gy = g*np.sin(theta)
    pylab.figure()
    magnify = 50.
    res=pylab.quiver(x, y, magnify*gx, magnify*gy, scale=16, headwidth=0, pivot='middle', units='width')
    tmpstr = str(np.median(g))
    if title is not None:
        titlestr = ", median shear = {0:.4}".format(np.median(g))
        pylab.title(title+titlestr)
    pylab.xlabel('X position [arcmin]')
    pylab.ylabel('Y position [arcmin]')
    pylab.show()        

def drawSizeVariation(kappa, x, y, title=None):
    """Show variation in PSF size across the field of view, based on a grid of "kappa" values.
    """
    pylab.figure()
    pylab.pcolor(x, y, kappa)
    pylab.colorbar()
    if title is not None:
        pylab.title(title)
    pylab.xlabel('X position [arcmin]')
    pylab.ylabel('Y position [arcmin]')
    pylab.show()

def drawBoth(g1, g2, kappa, x, y, title=None, outfile=None):
    """Draw both types of plots from arrays of g1, g2, kappa values.  Write to file `outfile`.
    """
    g = (g1**2 + g2**2)**0.5
    theta = 0.5*np.arctan2(g2, g1)
    gx = g*np.cos(theta)
    gy = g*np.sin(theta)
    pylab.figure()
    magnify = 50.
    pylab.pcolor(x/60., y/60., kappa)
    pylab.colorbar()
    res=pylab.quiver(x/60., y/60., magnify*gx, magnify*gy, scale=16, headwidth=0, pivot='middle', units='width')
    tmpstr = str(np.median(g))
    if title is not None:
        titlestr = ", median shear = {0:.4}".format(np.median(g))
        pylab.title(title+titlestr)
    pylab.xlabel('X position [deg]')
    pylab.ylabel('Y position [deg]')
    if outfile is None:
        pylab.show()
    else:
        pylab.savefig(outfile)

# Main function:
#    Given an exposure time and telescope size, generate a random P(k), use it to simulate PSF
#    anisotropy and size variation across the FOV, and check the outputs in various ways.
if __name__ == "__main__":

    # Define some defaults
    t_exp = 120.0 # exposure time in seconds
    diam = 4.0    # telescope size
    handbook = False # are we doing a combined size/whisker plot for handbook?
    outfile = 'test.png'

    # Currently set up to do a random realization with seed initialized from the time.
    A, theta_0, seeing = atmosPSFParams(t_exp=t_exp, diam=diam)
    # These return values are floats.  If we had requested more than one random realization, we
    # would have gotten NumPy arrays.  Currently we require t_exp and diam to be a single number,
    # but would be trivial to generalize so the user can supply arrays of values instead.

    print "For exposure time of",t_exp,"s and a",diam,"m telescope, this random realization of"
    print "the atmosphere has A, theta_0, seeing = ",A,theta_0,seeing

    # Make the atmospheric P(k)
    print "Making the P(k)"
    k, Pk = makeAtmosPSFPk(A, theta_0)

    # Make the gridded e1, e2.
    print "Getting the gridded PSF parameters"
    # Will want in general to subsample by 10 so we can have multiple grids sampling the same PSF.
    # But for the sake of easy viewing, use subsample=3.
    e1, e2, kappa, x, y = getAtmosPSFGrid(k, Pk, subsample=3)

    if handbook:
        print "Plotting PSF shear, kappa in a single plot"
        drawBoth(e1/2., e2/2., kappa, x, y, title="PSF patterns in 2x2 degree field", outfile=outfile)
    else:
        print "Plotting the PSF shear as a whisker plot"
        drawWhiskerShear(e1/2., e2/2., x, y, title='PSF shear in 2x2 degree field')
        print "Plotting the PSF kappa as a color plot"
        drawSizeVariation(kappa, x, y, title='PSF size variation in 2x2 degree field')
