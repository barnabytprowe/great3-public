import galsim
import numpy as np
from . import constants

def makeBuilder(shear_type, obs_type, multiepoch, ps_dir):
    """Return a ShearBuilder appropriate for the given options.

    @param[in] shear_type  Shear field type: either "constant" or "variable"
    @param[in] obs_type    Observation type: either "ground" or "space"
    @param[in] multiepoch  Multiepoch?  True or False
    @param[in] ps_dir      Directory with tabulated iCosmo shear power spectra.
    """
    if shear_type == 'variable':
        return VariableShearBuilder(ps_dir=ps_dir, obs_type=obs_type, multiepoch=multiepoch)
    elif shear_type == 'constant':
        return ConstantShearBuilder(obs_type=obs_type, multiepoch=multiepoch)
    else:
        raise ValueError("Invalid shear_type: %s - must be 'constant' or 'variable'" % shear_type)

class ShearBuilder(object):

    def generateFieldParameters(self, rng, field_index):
        """Return a dict of metaparameters for the given field.
        These will be passed to generateCatalog when it is called.

        @param[in] rng          A galsim.BaseDeviate to be used for
                                any random numbers.
        @param[in] field_index  Index of the field of images being simulated
        """
        raise NotImplementedError("ShearBuilder is abstract.")

    def generateSubfieldParameters(self, rng, subfield_index, field_parameters):
        """Return a dict of metaparameters for the given subfield.
        These will be passed to generateCatalog when it is called.

        @param[in] rng                     A galsim.BaseDeviate to be used for
                                           any random numbers.
        @param[in] subfield_index          Index of patch of sky being simulated
        @param[in] field_parameters        Output from generateFieldParameters.
        """
        raise NotImplementedError("ShearBuilder is abstract.")

    def generateEpochParameters(self, rng, subfield_index, epoch_index):
        """Return a dict of metaparameters for the given epoch.
        These will be passed to generateCatalog when it is called.

        @param[in] rng          A galsim.BaseDeviate to be used for
                                any random numbers.
        @param[in] epoch_index  Index of the epoch of a given subfield being simulated.
        """
        raise NotImplementedError("ShearBuilder is abstract.")

    def generateCatalog(self, rng, catalog, parameters, offsets, subfield_index):
        """Fill the g1 and g2 columns of the given catalog with the
        lensing shear values to apply at each point.

        @param[in]     rng         A galsim.BaseDeviate to be used for any random numbers.
        @param[in,out] catalog     A structured NumPy array to fill.  The 'index', 'x',
                                   and 'y' columns will already be filled, and should be
                                   used to evaluate the shear field for spatially-variable
                                   shear.  The 'g1' and 'g2' columns should be filled with this
                                   method.  Other columns may be present as well, and should be
                                   ignored.
        @param[in]     parameters  A dict of metaparameters, as returned by the
                                   generateParameters() method.
        @param[in]     offsets     Offsets of this subfield with respect to the first in the field.
        @param[in] subfield_index  Index of this subfield within the branch, used to construct a
                                   unique object ID for each galaxy.
        """
        raise NotImplementedError("ShearBuilder is abstract.")

class ConstantShearBuilder(ShearBuilder):
    """ShearBuilder for constant shear fields.
    """

    def __init__(self, obs_type, multiepoch, min_g=1E-2, max_g=5E-2):
        self.min_g = min_g
        self.max_g = max_g
        self.obs_type = obs_type
        self.multiepoch = multiepoch

    def generateFieldParameters(self, rng, field_index):
        theta = 2.0*rng()*np.pi
        shear_rng = galsim.DistDeviate(rng, function = lambda g : g,
                                       x_min = self.min_g, x_max = self.max_g)
        g = shear_rng()
        return dict(g1=g*np.cos(2.0*theta), g2=g*np.sin(2.0*theta), mu=1.)

    def generateSubfieldParameters(self, rng, subfield_index, field_parameters):
        # For the constant shear case, the results for the subfield are the same as those for the
        # field.  So just pass the field_parameters along without doing anything else.
        return field_parameters

    def generateEpochParameters(self, rng, subfield_index, epoch_index):
        raise NotImplementedError("ConstantShearBuilder makes shear parameters at field level!")

    def generateCatalog(self, rng, catalog, parameters, offsets, subfield_index):
        # offsets is not actually relevant yet, but we put it as a placeholder for variable shear,
        # where it's necessary.
        xsize = constants.xsize[self.obs_type][self.multiepoch]
        ysize = constants.ysize[self.obs_type][self.multiepoch]
        catalog["x_field_pos"] = (catalog["x"]+1+0.5*xsize)/xsize-1
        catalog["y_field_pos"] = (catalog["y"]+1+0.5*ysize)/ysize-1
        for record in catalog:
            record["g1"] = parameters["g1"]
            record["g2"] = parameters["g2"]
            record["mu"] = parameters["mu"]
            record["ID"] = 1e6*subfield_index + 1e3*record["x_field_pos"] + record["y_field_pos"]

class VariableShearBuilder(ShearBuilder):
    """ShearBuilder for variable shear fields.
    """
    # Define some basic parameters related to our grid settings.  Note that the real ell_max is not
    # as given here, because we use multiple offset grids to access smaller scales / higher ell.
    # However, since this is only used to define nuisance functions, it's okay.  What it means is
    # that the nuisance function will not be very significant at those higher ell, i.e., the power
    # spectrum will be closer to cosmological.
    ell_min = 360./constants.image_size_deg
    ell_max = ell_min*constants.nrows/2.
    # We define a pivot/central scale for the nuisance function based on a geometric mean of the
    # min, max usable ell values.
    ell_piv = np.sqrt(ell_min*ell_max)
    # Define some things we need to construct filenames for tabulated cosmological shear power
    # spectra.
    infile_pref = 'cosmo-fid.zmed'
    zmed_str = ['0.75', '1.00', '1.25']
    nzmed = len(zmed_str)
    # Choose some values for the shapelets nuisance functions.  We go up to order 5 and allow a 
    # maximum variation in P(k) due to any given shapelets order that is ~10% (max_amp).
    max_order = 5
    beta = 0.3
    max_amp = 0.1
    # Set up empty cache for power spectrum parameters and galsim.PowerSpectrum object
    cached_ps_tab = None
    cached_ps_nuisance = None
    cached_ps = None
    # Include a multiplicative factor for the shear power spectrum to make metric more sensitive to
    # m.  This will be applied to all randomly-constructed shear power spectra.
    mult_factor = 2.

    def __init__(self, ps_dir, obs_type, multiepoch):
        self.ps_dir = ps_dir
        self.obs_type = obs_type
        self.multiepoch = multiepoch

    def generateFieldParameters(self, rng, field_index):
        # Need to decide on parameters for the cosmological part of the power spectrum first.
        # For this, we just need to choose a single random number to interpolate between our
        # cosmological P(k) that are tabulated from iCosmo.
        #
        # This will just be a random number from a uniform distribution that goes from 0 to N_P-1
        # where N_P is the number that are tabulated.
        ps_tab = (self.nzmed-1)*rng()

        # Now, need to choose parameters for the nuisance shapelets functions.  It is an array of
        # self.max_order shapelets parameters according to the specified maximum amplitude.
        ps_nuisance = np.zeros(self.max_order)
        for i_order in range(self.max_order):
            ps_nuisance[i_order] = -self.max_amp + 2*self.max_amp*rng()

        return dict(ps_tab=ps_tab, ps_nuisance=ps_nuisance)

    def generateSubfieldParameters(self, rng, subfield_index, field_parameters):
        # For the variable shear case, the power spectrum for each subfield is the same as for the
        # subfield overall.  So just pass the field_parameters along without doing anything else.
        return field_parameters

    def generateEpochParameters(self, rng, subfield_index, epoch_index):
        raise NotImplementedError("VariableShearBuilder makes shear parameters at field level!")

    def generateCatalog(self, rng, catalog, parameters, offsets, subfield_index):
        # Here's where most of the work happens.
        #
        # We need a cache for a grid of shear values covering the entire field, i.e., including all
        # possible positions in all subfields (modulo sub-pixel offsets from the subfield grid -
        # we're not trying to take those into account).  If there is nothing in the cache for this
        # field, then make a new grid and save it in the cache.
        #
        ps_tab = parameters["ps_tab"]
        ps_nuisance = parameters["ps_nuisance"]
        if ps_tab != self.cached_ps_tab:
            # If nothing is cached for this power spectrum, then first we have to define the power
            # spectrum in a way that the galsim lensing engine can use it.
            # Begin by identifying and reading in the proper files for the cosmological part of the
            # power spectrum.
            file_index = np.floor(ps_tab)
            residual = ps_tab - file_index
            import os
            infile1 = os.path.join(self.ps_dir , 
                                   self.infile_pref + self.zmed_str[int(file_index)]+'.out')
            data1 = np.loadtxt(infile1).transpose()
            ell = data1[0]
            p1 = data1[1]
            infile2 = os.path.join(self.ps_dir , 
                                   self.infile_pref + self.zmed_str[int(file_index)+1]+'.out')
            data2 = np.loadtxt(infile2).transpose()
            p2 = data2[1]
            # Now make a geometric mean to get the cosmological power spectrum.
            p_cos = (p1**(1.-residual))*(p2**residual)
            p_cos *= self.mult_factor
            # Construct the shapelets nuisance functions
            x = np.log10(ell/self.ell_piv)
            n_ell = len(ell)
            b_values = np.zeros((self.max_order, n_ell))
            for order in range(0, self.max_order):
                b_values[order,:] = self._bn(order, x, self.beta)
            nuisance_func = np.zeros(n_ell)
            for order in range(0, self.max_order):
                nuisance_func += ps_nuisance[order]*b_values[order,:]
            p_use = p_cos*(1.0+nuisance_func)
            # Note: units for ell, p_use are 1/radians and radians^2, respectively.

            # Now, we have arrays we can use to make a power spectrum object with E-mode power
            # only.  While we are at it, we cache it and its parameters.
            ps_lookup = galsim.LookupTable(ell, p_use, x_log=True, f_log=True)
            self.cached_ps = galsim.PowerSpectrum(ps_lookup, units = galsim.radians)
            self.cached_ps_tab = ps_tab
            self.cached_ps_nuisance = ps_nuisance

            # Define the grid on which we want to get shears.
            # This is a little tricky: we have a setup for subfield locations within the field that
            # is defined in builder.py function generateSubfieldOffsets.  The first subfield is
            # located at the origin, and to represent it alone, we would need a constants.nrows x
            # constants.ncols grid of shears.  But since we subsample by a parameter given as
            # constants.subfield_grid_subsampling, each grid dimension must be larger by that
            # amount.
            if constants.nrows != constants.ncols:
                raise NotImplementedError("Currently variable shear grids require nrows=ncols")
            n_grid = constants.subfield_grid_subsampling * constants.nrows
            grid_spacing = constants.image_size_deg / n_grid

            # Run buildGrid() to get the shears and convergences on this grid.  However, we also
            # want to effectively change the value of k_min that is used for the calculation, to get
            # a reasonable shear correlation function on large scales without excessive truncation.
            # TODO: check that this value of kmin_factor is adequate!  Right now it is determined by
            # the fact that our grid gives ell_min=36, and the iCosmo outputs only go to ell=10, so
            # using kmin_factor>3 would require some extrapolation.
            # We also define a grid center such that the position of the first pixel is (0,0).
            grid_center = 0.5 * (constants.image_size_deg - grid_spacing)
            self.cached_ps.buildGrid(grid_spacing = grid_spacing,
                                     ngrid = n_grid,
                                     units = galsim.degrees,
                                     rng = rng,
                                     center = (grid_center, grid_center),
                                     kmin_factor = 3)
            # Now that our cached PS has a grid of shears / convergences, we can use getLensing() to
            # get the quantities we need for a lensing measurement at any position, so this part of
            # the calculation is done.

        # Now get the shears/convergences for each galaxy position in the
        # catalog.  This is fastest if done all at once, with one call to getLensing.  And this is
        # actually slightly tricky, because we have to take into account: 
        #    (1) The position of the galaxy within the subfield.
        #    (2) The offset of the subfield with respect to the field.
        # And make sure we've gotten the units right for both of these.  We are ignoring centroid
        # shifts of order 1 pixel (max 0.2" for ground data) which can occur within an image.
        #
        # We can define object indices in x, y directions - i.e., make indices that range
        # from 0 to constants.nrows-1.
        xsize = constants.xsize[self.obs_type][self.multiepoch]
        ysize = constants.ysize[self.obs_type][self.multiepoch]
        x_ind = (catalog["x"]+1+0.5*xsize)/xsize-1
        y_ind = (catalog["y"]+1+0.5*ysize)/ysize-1
        # Turn this into (x, y) positions within the subfield, in degrees.
        x_pos = x_ind * constants.image_size_deg / constants.nrows
        y_pos = y_ind * constants.image_size_deg / constants.ncols
        # But now we have to add the subfield offset.  These are calculated as a fraction of the
        # separation between galaxies, so we have to convert to degrees.
        x_pos += offsets[0] * constants.image_size_deg / constants.nrows
        y_pos += offsets[1] * constants.image_size_deg / constants.ncols
        catalog["g1"], catalog["g2"], catalog["mu"] = \
            self.cached_ps.getLensing(pos=(x_pos, y_pos), units=galsim.degrees)
        # Previous numbers were in degrees.  But now we need to save some numbers for ID generation,
        # which have to be ints.  So we will save them in units of subfield grid spacing, i.e.,
        # within a given subfield, galaxies are spaced by constants.subfield_grid_subsampling.
        # Right now x_ind, y_ind are integers (spaced by 1) and offsets[0] and offsets[1] span the
        # range (0, 1/constants.subfield_grid_subsampling), so the line below has
        # constants.subfield_grid_subsampling multiplying both.
        catalog["x_field_pos"] = np.round(
            constants.subfield_grid_subsampling * (offsets[0] + x_ind)).astype(int)
        catalog["y_field_pos"] = np.round(
            constants.subfield_grid_subsampling * (offsets[1] + y_ind)).astype(int)
        for record in catalog:
            record["ID"] = 1e6*subfield_index + 1e3*record["x_field_pos"] + record["y_field_pos"]

    def _bn(self, n, x, beta):
        """A helper function to compute shapelets functions for a given order n, for specified x and
        width beta.
        """
        phi_n_x = self._phin(n, x/beta)
        return phi_n_x / np.sqrt(beta)

    def _hermite(self, n, x):
        try:
            import scipy.special as spec
            # get the H_n function from scipy
            hn = spec.hermite(n)
            # evaluate it at our x
            return hn(x)
        except:
            # If you don't have scipy, use a simple recursion relation:
            if n == 0:
                return 1.
            else:
                hkm1 = 1.
                hk = 2.*x
                for k in range(1,n):
                    hk, hkm1 = 2.*x*hk - 2.*k*hkm1, hk
                return hk

    def _phin(self, n, x):
        """A helper function defining shapelets basis functions at an array of positions x, for
        order n."""
        import math
        hn_x = self._hermite(n,x)
        # Put in the exponential factor, and properly normalize it.
        phi_n_x = hn_x * np.exp(-(x**2)/2.)
        phi_n_x /= np.sqrt((2.**n) * np.sqrt(np.pi) * math.factorial(n))
        return phi_n_x
