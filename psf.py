import galsim
import numpy as np
from . import constants

def makeBuilder(obs_type, variable_psf, multiepoch, shear_type, atmos_ps_dir):
    """Return a PSFBuilder appropriate for the given options.

    @param[in] obs_type     Observation type: either "ground" or "space".
    @param[in] variable_psf If True, we need a spatially-varying PSF.
    @param[in] multiepoch   If True, this is a multiepoch simulation, and we are determining the PSF
                            of multiple exposures simultaneously.  If False, the PSF should be that
                            of a single exposure.
    @param[in] shear_type   Shear type: either "constant" or "variable".  We need this for the
                            determination of how many subfields to use per field, and whether there
                            are so few PSFs that we have to force them to follow the seeing
                            distribution, not any other reason.
    @param[in] atmos_ps_dir Directory with tabulated atmospheric PSF anisotropy power spectra.
    """
    if obs_type == "space" or obs_type == "ground":
        if variable_psf:
            return VariablePSFBuilder(obs_type, multiepoch, shear_type, atmos_ps_dir)
        else:
            return ConstPSFBuilder(obs_type, multiepoch, shear_type)
    else:
        raise ValueError("Invalid obs_type: %s - must be 'ground' or 'space'" % obs_type)

class PSFBuilder(object):

    def generateFieldParameters(self, rng, field_index):
        """Return a list of dicts of metaparameters for the given field of observations.  These
        will be passed to generateCatalog when it is called.

        @param[in] rng          A galsim.UniformDeviate to be used for
                                any random numbers.
        @param[in] field_index  Index of the field of images being simulated.

        Each entry in the list corresponds to one of the epochs (i.e., a single entry for single
        epoch observations, or multiple entries for multiepoch).
        """
        raise NotImplementedError("PSFBuilder is abstract.")

    def generateEpochParameters(self, rng, subfield_index, epoch_index, field_parameters):
        """Return a dict of metaparameters for the given subfield and
        epoch.  These will be passed to generateCatalog when it
        is called.

        @param[in] rng               A galsim.UniformDeviate to be used for
                                     any random numbers.
        @param[in] subfield_index    Index of the simulated patch of sky
                                     we are "observing" with this PSF.
        @param[in] epoch_index       Index of this "epoch" among
                                     those of the same subfield.
        @param[in] field_parameters  Output from generateFieldParameters.

        A "schema" entry is required in the returned dict: a list of
        (name, type) tuples that define the catalog fields that will
        be filled by this builder.  These field names may be used
        in catalogs that also have fields from other builders, so
        care should be taken to ensure that each schema entry are unique.
        """
        raise NotImplementedError("PSFBuilder is abstract.")

    def generateCatalog(self, rng, catalog, parameters, offsets, normalized):
        """Fill columns of the given catalog with per-position PSF parameters.

        @param[in]     rng         A galsim.UniformDeviate to be used for any random numbers.
        @param[in,out] catalog     A structured NumPy array to fill.  The 'index', 'x',
                                   and 'y' columns will already be filled, and should be
                                   used as the points at which to evaluate the PSF
                                   spatial variation (if any).
                                   All columns defined by the 'schema' entry in the dict
                                   returned by generateParameters() will also be present,
                                   and should be filled.  Other columns may be present as
                                   well, and should be ignored.
        @param[in]     parameters  A dict of metaparameters, as returned by the
                                   generateParameters() method.
        @param[in]     offsets     Offsets of this subfield with respect to the first in the field.
        @param[in]     normalized  Whether star objects should be flux-normalized (as
                                   appropriate for convolving galaxies) or have a distribution
                                   of fluxes (as appropriate for making real-PSF star images).
        """
        raise NotImplementedError("PSFBuilder is abstract.")

    def makeConfigDict(self):
        """Make the 'psf' portion of the config dict.
        """
        raise NotImplementedError("PSFBuilder is abstract.")

    def makeGalSimObject(self, record, parameters):
        """Given a catalog record and a dict of metaparameters (as generated
        by generateCatalog() and generateParameters(), respectively), create
        a galsim.GSObject that represents the PSF (not including the pixel
        component).

        NOTE: if we want to use the galsim config interface to create images,
        we should replace this method with one that writes the relevant
        section of a config file.
        """
        raise NotImplementedError("PSFBuilder is abstract.")

# A quick helper function that sets any Catalog entries in a dict to use index=0
def UseZeroIndex(d):
    if 'type' in d and d['type'] == 'Catalog':
        # If this is a catlog entry, use index = 0.
        d['index'] = 0
    else:
        # Else recurse onto all sub-dicts.
        for key in d:
            if isinstance(d[key],dict): UseZeroIndex(d[key])
            elif isinstance(d[key],list): UseZeroIndex(d[key])


class ConstPSFBuilder(PSFBuilder):
    """PSFBuilder for experiments with uniform PSF across the images.

    Depending on obs_type and multiepoch, outputs are slightly different.
    """
    # Set some basic parameters here.  For optical PSFs, this should include some dicts for space
    # vs. ground using the obs_type parameter.

    # For lam/diam, Euclid has 0.1375".
    #               WFIRST has a range for the different bands, but with 2.4m mirror and 1.1 micron,
    #               it is 0.0945".
    # In terms of sampling, we use the correspondence between OpticalPSF <-> Airy 
    # and FWHM >~ 2 pix. For an Airy, FWHM is very close to lam/diam.  We want to be
    # Nyquist sampled at 0.05" pixels, so we need lam/diam > 0.1" (roughly).  For the multiepoch
    # case in space, we want to be undersampled by a factor of 2 at 0.1" pixels, so we need
    # lam/diam < 0.2".  Let's give a cushion on these numbers as well.
    #
    # And then for the ground, let's consider typical wavelengths between 500-900 nm, and diam
    # between 2-8m.  This gives a total range for lam/diam = 0.013 -> 0.083.
    min_lam_over_diam = {
        "space" : 0.12, # arcsec
        "ground" : 0.013, # arcsec
        }
    max_lam_over_diam = {
        "space" : 0.17, # arcsec
        "ground" : 0.083,
        }

    # Choose some value for obscuration by the secondary, using a uniform distribution between the
    # min/max values given here.  We use the same range for space / ground, but to enable
    # flexibility later on, this is also a dict.
    min_obscuration = {
        "space" : 0.1,
        "ground" : 0.1,
        }
    max_obscuration = {
        "space" : 0.5,
        "ground" : 0.5,
        }

    # Aberrations for space telescopes: the baseline WFIRST2.4 PSF model (without tilt,
    # misalignment, etc.) seem to have aberrations roughly in the range [-0.01, 0.01] (per
    # aberration).  This gives a RMS summed over all 8 aberrations of ~0.03 in units of wavelength.
    # I'm going to assume they are uncorrelated (not necessarily true!).  Furthermore, an additional
    # RMS of 0.07 is allowed above the design.  So for now, let's give a total RMS of 0.08.  Euclid
    # may allow larger aberrations when considered in units of wavelength.
    #
    # Ground telescopes tend to have larger aberrations.  For now I'm putting in 4x larger than
    # space, but I need to check out Hironao's PR to get a better number for this.
    rms_aberration = {
        "space" : 0.08,
        "ground" : 0.41,
        }

    # We have to define the set of aberrations that we'll use, and some names that will go into the
    # schema to specify them.  The purpose of setting it up this way is that if we decide we want
    # our constant PSF branch PSFs to be simpler, we can easily remove aberrations without
    # changing too much of the code.
    use_aber = ["defocus", "astig1", "astig2", "coma1", "coma2", "trefoil1", "trefoil2", "spher"]
    n_aber = len(use_aber)
    # Decide on the relative weight given to each aberration.  Those with higher weights contribute
    # a greater variation to the total rms_aberration.  For space, we are representing primarily the
    # aberrations on top of the design residual, for which there is no real model and so having
    # uniform weights seems like the only way to go.  For ground, based on variation of the optical
    # PSF model across the FOV and including defocus / tilt / etc., it seems that variation in
    # defocus is more important than the other aberrations, so we give defocus a substantially
    # higher weight and treat the others equally.
    aber_weights = {
        "space" : np.ones(n_aber),
        "ground" : np.array((0.36, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07))
        }
    opt_schema_pref = "opt_psf_"

    # Struts: we choose an integer number of struts from the values on the following list.
    # The angle from vertical is chosen at random.
    # For now the thickness is not a free parameter, though it could be.  The possible number of
    # struts differs for space and ground.  Because of issues with accuracy of rendering PSFs with
    # struts in the latter case (see discussion on
    # https://github.com/rmandelb/great3-private/pull/53), we do not have any struts for
    # ground-based sims.
    strut_list = {
        "space" : [0, 2, 3, 4, 6],
        "ground" : [0],
        }

    # For space-based PSFs only: include jitter and charge diffusion.  The plan as stated on #42,
    # based on discussion with Lance Miller and several WFIRST people, is to make jitter a Gaussian
    # with RMS of 0.005-0.015" per axis, e=0-0.3 (uniform distribution), and a random direction.  So
    # let's make a total sigma in the range sqrt(2)*(0.005-0.015).
    min_jitter_sigma = np.sqrt(2.)*0.005 # arcsec
    max_jitter_sigma = np.sqrt(2.)*0.005
    min_jitter_e = 0.
    max_jitter_e = 0.3
    # And for charge diffusion, the plan is a Gaussian with sigma=0.1-0.2pix (uniform
    # distribution), with e=0-0.2 (uniform distribution), always the same direction.  Let's use
    # always the same pixel scale for determining that sigma, 0.05".  For simplicity since we always
    # want the same direction of ellipicity, let's say that it's always in +e1.
    min_charge_sigma = 0.1*0.05 # arcsec
    max_charge_sigma = 0.2*0.05 # arcsec
    min_charge_e1 = 0.
    max_charge_e1 = 0.2

    # Atmospheric PSF info:
    atmos_scheme_pref = "atmos_psf_"
    # Draw seeing from a distribution:
    fwhm_arcsec = 0.05 + 0.10*np.arange(17)
    #Old version is commented out.  We use a new version that has only 0.5-0.85
    #freq = (0., 0., 0., 7.5, 19., 20., 17., 13., 9., 5., 3.5, 2., 1., 1., 0.5, 0.0, 0.0)
    freq = (0., 0., 0., 0.0, 0.0, 20., 17., 13., 9., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    # Later on we will draw from this distribution using:
    # dd = galsim.DistDeviate(uniform_deviate,
    #                        function=galsim.LookupTable(fwhm_arcsec, freq, interpolant='linear'))
    # We also need a value of atmospheric PSF anisotropy.  Based on the atmospheric PSF
    # investigations, we use the following prescription:

    # The total variance (over both components) of PSF ellipticity is in the range 1e-4 - 8e-4 for
    # Texp=20s and diam = 6.7m.  Let's say for now that typically the diam value might be 0.5-1.2
    # times the fiducial value, multiplying the variance by 0.84-1.7.  Let's also say our exposure
    # times are typically 30s-3min, multiplying the variance by 0.11-0.67.  The lowest possible
    # variance coming out of all of these factors is ~9e-6, and the highest is ~9e-4, but those
    # extreme values require some conspiracy between a few random numbers.  Also for the sake of
    # GREAT3 it is somewhat inconvenient if the PSF is *too* round, so I would like to exclude the
    # values of variance below 1e-4.  So let's choose the norm of the atmospheric PSF ellipticity
    # as sqrt(variance) for variance in the range 1e-4 -> 9e-4 (uniform).
    min_atmos_psf_e = np.sqrt(1.e-4)
    max_atmos_psf_e = np.sqrt(9.e-4)

    def __init__(self, obs_type, multiepoch, shear_type):
        self.obs_type = obs_type
        self.multiepoch = multiepoch
        self.variable_psf = False # by definition, so it doesn't get passed in as an arg
        self.shear_type = shear_type

        # If (self.variable_psf == False and self.shear_type == "variable") and self.multiepoch ==
        # False, then we only have very few PSFs in each branch.  This includes 4 branches: control
        # experiment + variable shear, real galaxy + variable shear, in both cases for space and
        # ground sims.  For the 2 of these branches that are ground-based sims, this can be
        # problematic, since we want to properly sample the distribution of atmospheric seeing, and
        # it's hard to do that with so few random numbers.  So, for this case only, we are going to
        # draw random values of seeing from the CDF.  For the sake of simplicity, for this one case
        # only, we are going to store some information we need to track these without duplication.
        if self.shear_type == "variable" and not self.multiepoch and self.obs_type == "ground":
            n_subfields_per_field = \
                constants.n_subfields_per_field[self.shear_type][self.variable_psf]
            # Only use non-deep fields for this. For deep fields, force it to be something around
            # the median seeing.  (For the case where we're not forcing anything, there will be
            # enough deep observations that there should be some near median.  But for the case
            # where we have very few PSFs, there is only 1 for the deep data, so it seems fair to
            # require it to be something decent-ish.)
            n_fields = (constants.n_subfields - constants.n_deep_subfields) / n_subfields_per_field
            self.dpercentile = (1./float(n_fields))
            self.force_distribution = True
            self.used_percentile = np.zeros(n_fields).astype(bool)
        else:
            self.force_distribution = False

    def generateFieldParameters(self, rng, field_index):
        # Define some constants that we will need throughout.
        if self.multiepoch:
            n_epochs = constants.n_epochs
        else:
            n_epochs = 1

        # Make RNG
        uniform_deviate = galsim.UniformDeviate(rng)
        gaussian_deviate = galsim.GaussianDeviate(rng)
        if self.obs_type == "ground":
            dist_deviate = galsim.DistDeviate(
                rng,
                function = galsim.LookupTable(self.fwhm_arcsec, self.freq, interpolant='linear'))

        # Get the basic parameters per epoch.  In previous versions of the code, the parameters also
        # depended on the subfield.  However, now we keep all subfields within the field the same,
        # and change the number of subfields per field if we want to change things around.  So these
        # arrays only have one dimension with length n_epochs.
        aber_dict = dict()
        for aber in self.use_aber:
            aber_dict[aber]=np.zeros(n_epochs)
        lam_over_diam = np.zeros(n_epochs)
        obscuration = np.zeros(n_epochs)
        n_struts = np.zeros(n_epochs)
        strut_angle = np.zeros(n_epochs)
        pad_factor = np.zeros(n_epochs)
        if self.obs_type == "ground":
            atmos_psf_fwhm = np.zeros(n_epochs)
            atmos_psf_e = np.zeros(n_epochs)
            atmos_psf_beta = np.zeros(n_epochs)
        else:
            jitter_sigma = np.zeros(n_epochs)
            jitter_e = np.zeros(n_epochs)
            jitter_beta = np.zeros(n_epochs)
            charge_sigma = np.zeros(n_epochs)
            charge_e1 = np.zeros(n_epochs)

        # Get numbers for various parameters that will be used to make an OpticalPSF:
        # * We choose a range of lam_over_diam motivated by upcoming telescopes and stored within
        #   the class.
        # * We draw random aberrations according to some overall RMS for the total over all
        #   aberrations, assuming they are independent.
        # * We choose a reasonable range of obscuration, using min/max defined in class.
        # * We choose a reasonable range of struts, using the list of options defined in class.
        # * We choose a random strut_angle, within the range expected given the number of struts.
        for i_epoch in range(n_epochs):
            lam_over_diam[i_epoch] = (
                uniform_deviate() * (self.max_lam_over_diam[self.obs_type] -
                                     self.min_lam_over_diam[self.obs_type]) +
                self.min_lam_over_diam[self.obs_type] )

            tmp_vec = np.zeros(self.n_aber)
            for ind_ab in range(self.n_aber):
                tmp_vec[ind_ab] = \
                    self.rms_aberration[self.obs_type] * self.aber_weights[self.obs_type][ind_ab] *\
                    gaussian_deviate() / np.sqrt(np.sum(self.aber_weights[self.obs_type]**2))
                aber_dict[self.use_aber[ind_ab]][i_epoch] = tmp_vec[ind_ab]

            obscuration[i_epoch] = (
                uniform_deviate() * (self.max_obscuration[self.obs_type] -
                                     self.min_obscuration[self.obs_type]) +
                self.min_obscuration[self.obs_type] )

            strut_rand_ind = \
                int( np.floor( uniform_deviate() * len(self.strut_list[self.obs_type]) ) )
            tmp_n_struts = self.strut_list[self.obs_type][strut_rand_ind]
            n_struts[i_epoch] = tmp_n_struts

            if tmp_n_struts > 0:
                strut_angle[i_epoch] = 360.*uniform_deviate()
            else:
                strut_angle[i_epoch] = 0.

            # This is the default.  But option here to go larger if desired.
            pad_factor[i_epoch] = 1.5

            if self.obs_type == "ground":
                # for seeing values, check whether we need to force them to follow the distribution
                # (because there are very few PSFs in the branch) or not.
                if not self.force_distribution:
                    atmos_psf_fwhm[i_epoch] = dist_deviate()
                else:
                    # First check whether we're in a deep field or not.  For deep fields, we simply
                    # force something median-ish.  Otherwise, do a more complicated calculation.
                    n_subfields_per_field = \
                        constants.n_subfields_per_field[self.shear_type][self.variable_psf]
                    n_shallow_fields = \
                        (constants.n_subfields - constants.n_deep_subfields) / n_subfields_per_field
                    if field_index < n_shallow_fields:
                        # Draw a random uniform deviate from 0-1, and choose which bin it should go
                        # into.  Check if that bin has already been taken.  If yes, then keep doing
                        # this until it's not taken.  If not, then take that bin (so no other field
                        # can), and take the value of seeing corresponding to that value of the CDF.
                        # Note, the code below will only work for single epoch sims, but we should
                        # have taken care of that logically, above, when setting the value of
                        # self.force_distribution.
                        while True:
                            test_rand = uniform_deviate()
                            test_bin = int(test_rand/self.dpercentile)
                            if not self.used_percentile[test_bin]: break
                        self.used_percentile[test_bin] = True
                        atmos_psf_fwhm[i_epoch] = dist_deviate.val(test_rand)
                    else:
                        test_rand = 0.4 + 0.2 * uniform_deviate()
                        atmos_psf_fwhm[i_epoch] = dist_deviate.val(test_rand)

                # for other atmospheric PSF parameters, just draw at random even if there are
                # only a few PSFs per branch
                atmos_psf_e[i_epoch] = \
                    uniform_deviate()*(self.max_atmos_psf_e - self.min_atmos_psf_e) \
                    + self.min_atmos_psf_e
                atmos_psf_beta[i_epoch] = uniform_deviate()*180.0
            else:
                jitter_sigma[i_epoch] = \
                    uniform_deviate() * (self.max_jitter_sigma - \
                                             self.min_jitter_sigma) + self.min_jitter_sigma
                jitter_e[i_epoch] = \
                    uniform_deviate() * (self.max_jitter_e - \
                                             self.min_jitter_e) + self.min_jitter_e
                jitter_beta[i_epoch] = uniform_deviate()*180.0
                charge_sigma[i_epoch] = \
                    uniform_deviate() * (self.max_charge_sigma - \
                                             self.min_charge_sigma) + self.min_charge_sigma
                charge_e1[i_epoch] = \
                    uniform_deviate() * (self.max_charge_e1 - \
                                             self.min_charge_e1) + self.min_charge_e1

        # Make schema for catalog.
        schema = [("opt_psf_lam_over_diam", float),
                  ("opt_psf_obscuration", float),
                  ("opt_psf_n_struts", int),
                  ("opt_psf_strut_angle", float),
                  ("opt_psf_pad_factor", float)]
        for aber in self.use_aber:
            schema.append((self.opt_schema_pref+aber, float))
        if self.obs_type == "ground":
            schema.append(("atmos_psf_fwhm", float))
            schema.append(("atmos_psf_e", float))
            schema.append(("atmos_psf_beta", float))
        else:
            schema.append(("opt_psf_jitter_sigma", float))
            schema.append(("opt_psf_jitter_e", float))
            schema.append(("opt_psf_jitter_beta", float))
            schema.append(("opt_psf_charge_sigma", float))
            schema.append(("opt_psf_charge_e1", float))

        # Prepare the dict that this function must return.
        psf_dict = dict(schema=schema, lam_over_diam=lam_over_diam, obscuration=obscuration,
                        n_struts=n_struts, strut_angle=strut_angle, pad_factor=pad_factor)
        for aber in self.use_aber:
            psf_dict[aber] = aber_dict[aber]
        if self.obs_type == "ground":
            psf_dict["atmos_psf_fwhm"]=atmos_psf_fwhm
            psf_dict["atmos_psf_e"]=atmos_psf_e
            psf_dict["atmos_psf_beta"]=atmos_psf_beta
        else:
            psf_dict["jitter_sigma"]=jitter_sigma
            psf_dict["jitter_e"]=jitter_e
            psf_dict["jitter_beta"]=jitter_beta
            psf_dict["charge_sigma"]=charge_sigma
            psf_dict["charge_e1"]=charge_e1

        return psf_dict

    def generateEpochParameters(self, rng, subfield_index, epoch_index, field_parameters):
        # At the field level we already pre-determined the results for the different epochs, so just
        # return the result for this epoch.

        psf_dict = \
            dict(lam_over_diam=field_parameters["lam_over_diam"][epoch_index],
                 obscuration=field_parameters["obscuration"][epoch_index],
                 n_struts=field_parameters["n_struts"][epoch_index],
                 strut_angle=field_parameters["strut_angle"][epoch_index],
                 pad_factor=field_parameters["pad_factor"][epoch_index],
                 schema=field_parameters["schema"])
        for aber in self.use_aber:
            psf_dict[aber]=field_parameters[aber][epoch_index]
        if self.obs_type == "ground":
            psf_dict["atmos_psf_fwhm"] = field_parameters["atmos_psf_fwhm"][epoch_index]
            psf_dict["atmos_psf_e"] = field_parameters["atmos_psf_e"][epoch_index]
            psf_dict["atmos_psf_beta"] = field_parameters["atmos_psf_beta"][epoch_index]
        else:
            psf_dict["jitter_sigma"] = field_parameters["jitter_sigma"][epoch_index]
            psf_dict["jitter_e"] = field_parameters["jitter_e"][epoch_index]
            psf_dict["jitter_beta"] = field_parameters["jitter_beta"][epoch_index]
            psf_dict["charge_sigma"] = field_parameters["charge_sigma"][epoch_index]
            psf_dict["charge_e1"] = field_parameters["charge_e1"][epoch_index]

        return psf_dict

    def generateCatalog(self, rng, catalog, parameters, offsets, normalized):
        psf_parameters = parameters["psf"]
        for record in catalog:
            record["opt_psf_lam_over_diam"] = psf_parameters["lam_over_diam"]
            record["opt_psf_obscuration"] = psf_parameters["obscuration"]
            record["opt_psf_n_struts"] = psf_parameters["n_struts"]
            record["opt_psf_strut_angle"] = psf_parameters["strut_angle"]
            record["opt_psf_pad_factor"] = psf_parameters["pad_factor"]
            for aber in self.use_aber:
                record[self.opt_schema_pref+aber] = psf_parameters[aber]
            if self.obs_type == "ground":
                record["atmos_psf_fwhm"] = psf_parameters["atmos_psf_fwhm"]
                record["atmos_psf_e"] = psf_parameters["atmos_psf_e"]
                record["atmos_psf_beta"] = psf_parameters["atmos_psf_beta"]
            else:
                record["opt_psf_jitter_sigma"] = psf_parameters["jitter_sigma"]
                record["opt_psf_jitter_e"] = psf_parameters["jitter_e"]
                record["opt_psf_jitter_beta"] = psf_parameters["jitter_beta"]
                record["opt_psf_charge_sigma"] = psf_parameters["charge_sigma"]
                record["opt_psf_charge_e1"] = psf_parameters["charge_e1"]

    def makeConfigDict(self):
        d = {
            'type' : 'OpticalPSF',
            'lam_over_diam' : { 'type' : 'Catalog', 'col' : 'opt_psf_lam_over_diam' },
            'obscuration' : { 'type' : 'Catalog', 'col' : 'opt_psf_obscuration' },
            'nstruts' : { 'type' : 'Catalog', 'col' : 'opt_psf_n_struts' },
            'strut_angle' : { 'type' : 'Deg',
                              'theta' : { 'type' : 'Catalog', 'col' : 'opt_psf_strut_angle' }
                            },
            'defocus' : { 'type' : 'Catalog', 'col' : 'opt_psf_defocus' },
            'astig1' : { 'type' : 'Catalog', 'col' : 'opt_psf_astig1' },
            'astig2' : { 'type' : 'Catalog', 'col' : 'opt_psf_astig2' },
            'coma1' : { 'type' : 'Catalog', 'col' : 'opt_psf_coma1' },
            'coma2' : { 'type' : 'Catalog', 'col' : 'opt_psf_coma2' },
            'trefoil1' : { 'type' : 'Catalog', 'col' : 'opt_psf_trefoil1' },
            'trefoil2' : { 'type' : 'Catalog', 'col' : 'opt_psf_trefoil2' },
            'spher' : { 'type' : 'Catalog', 'col' : 'opt_psf_spher' },
            'pad_factor' : { 'type' : 'Catalog', 'col' : 'opt_psf_pad_factor' },
            'suppress_warning' : True
        }
        if self.obs_type == 'ground':
            d = {
                'type' : 'Convolve',
                'items' : [
                    {
                        'type' : 'Kolmogorov',
                        'fwhm' : { 'type' : 'Catalog', 'col' : 'atmos_psf_fwhm' },
                        'ellip' : { 
                            'type' : 'EBeta', 
                            'e' : { 'type' : 'Catalog', 'col' : 'atmos_psf_e' },
                            'beta' : { 'type' : 'Deg',
                                       'theta' : { 'type' : 'Catalog', 'col' : 'atmos_psf_beta' } }
                        }
                    },
                    d  # The above constructed OpticalPSF
                ]
            } 
        else:
            d = {
                'type' : 'Convolve',
                'items' : [
                    {
                        'type' : 'Gaussian',
                        'sigma' : { 'type' : 'Catalog', 'col' : 'opt_psf_jitter_sigma' },
                        'ellip' : { 
                            'type' : 'EBeta', 
                            'e' : { 'type' : 'Catalog', 'col' : 'opt_psf_jitter_e' },
                            'beta' : { 'type' : 'Deg',
                                       'theta' : { 'type' : 'Catalog', 
                                                   'col' : 'opt_psf_jitter_beta' } }
                        }
                    },
                    {
                        'type' : 'Gaussian',
                        'sigma' : { 'type' : 'Catalog', 'col' : 'opt_psf_charge_sigma' },
                        'ellip' : { 
                            'type' : 'E1E2',
                            'e1' : { 'type' : 'Catalog', 'col' : 'opt_psf_charge_e1' },
                            'e2' : 0
                        }
                    },
                    d  # The above constructed OpticalPSF
                ]
            } 

        # The parameters here are constant, so set index=0 for all Catalog entries:
        UseZeroIndex(d)

        return d

    def makeGalSimObject(self, record, parameters):
        # Size is specified in arcsec.

        # Choose which aberrations to use based on what should be in catalog
        # We have to start by making a dict with all aberrations that are nonzero.
        aber_dict = dict()
        # Then populate those that we use.
        for aber in self.use_aber:
            aber_dict[aber] = record[self.opt_schema_pref+aber]

        optical_psf = galsim.OpticalPSF(record["opt_psf_lam_over_diam"],
                                        obscuration=record["opt_psf_obscuration"],
                                        nstruts=record["opt_psf_n_struts"],
                                        strut_angle=record["opt_psf_strut_angle"]*galsim.degrees,
                                        pad_factor=record["opt_psf_pad_factor"],
                                        suppress_warning=True,
                                        **aber_dict)
        if self.obs_type == "space":
            jitter_psf = galsim.Gaussian(sigma=record["opt_psf_jitter_sigma"])
            e = record["opt_psf_jitter_e"]
            if e > 0.:
                jitter_psf.applyShear(e=e,
                                      beta = record["opt_psf_jitter_beta"]*galsim.degrees)
            charge_psf = galsim.Gaussian(sigma=record["opt_psf_charge_sigma"])
            e1 = record["opt_psf_charge_e1"]
            if e1 > 0.:
                charge_psf.applyShear(e1=e1)
            return galsim.Convolve(jitter_psf, charge_psf, optical_psf)
        else:
            atmos_psf = galsim.Kolmogorov(fwhm = record["atmos_psf_fwhm"])
            atmos_psf.applyShear(e = record["atmos_psf_e"], beta =
                                 record["atmos_psf_beta"]*galsim.degrees)
            return galsim.Convolve(atmos_psf, optical_psf)

class VariablePSFBuilder(PSFBuilder):
    """PSFBuilder for experiments with variable PSF across the images.

    This class basically functions as follows:
    (1) When initialized, it sets up some basic parameters like the number of tiles across the
        field-of-view.
    (2) generateFieldParameters draws all random numbers that are needed for generation of the PSF
        models.  These will be 2d arrays (indexed by tile and epoch).
    (3) generateEpochParameters just takes the generateFieldParameters outputs and extracts the ones
        for the epoch we care about.
    (4) generateCatalog does the heavy lifting: make (or take from the cache) the array of
        OpticalPSFModels and atmospheric PSFs to get the aberrations etc. interpolated to various
        positions.
    """
    # Set some basic parameters here.  We are sticking a bit more closely to the optical designs
    # that give us our variable PSFs compared to how we did things in ConstPSFBuilder.

    # For lam/diam, WFIRST has a range for the different bands, but with 2.4m mirror and 1 micron,
    # it is 0.086".  This is the right sampling for our single-epoch images, and appropriately
    # undersampled for our multi-epoch case.
    #
    # For the ground, let's consider typical wavelengths between 500-900 nm, and diam
    # of 4m.  This gives a total range for lam/diam = 0.026 -> 0.046.  Let's use 0.035".
    lam_over_diam = {
        "space" : 0.086, # arcsec
        "ground" : 0.035, # arcsec
        }
    diam = {
        "space" : 2.4, # m
        "ground" : 4., #m
        }
    lam = {
        "space" : 1.e9*lam_over_diam["space"]*diam["space"]*np.pi/(3600.*180.), # nm
        "ground" : 1.e9*lam_over_diam["ground"]*diam["ground"]*np.pi/(3600.*180.), # nm
        }

    # Use the value for obscuration by the secondary for our realistic PSF models.
    obscuration = {
        "space" : 0.28,
        "ground" : 0.35,
        }

    # Define the size of the tiles to be used to make the PSF models.  The considerations are that
    # we want to fit an integer number of tiles in a 10x10 degree field, but we also don't want to
    # stretch the scale of our pre-existing models excessively.  The choices given below lead to
    # stretching that can be as much as 15%, which is not excessive.
    n_tile_linear = {
        "space" : 20, # 20 x 20 tiles, each 0.5 x 0.5 degrees
        "ground" : 5, # 5 x 5 tiles, each 2 x 2 degrees
        }

    # We have to define the set of aberrations that we'll use, and some names that will go into the
    # schema to specify them.  We know we need all of these for the realistic PSF branch.
    use_aber = ["defocus", "astig1", "astig2", "coma1", "coma2", "trefoil1", "trefoil2", "spher"]
    n_aber = len(use_aber)
    opt_schema_pref = "opt_psf_"

    # For space, additional aberrations will be drawn randomly
    # according to some RMS that we define.
    space_rms_additional_aber = 0.075
    # For ground, our additional aberrations come from defocus, tilt, and decenter.
    # According to Hironao's notes, some expected values that Aaron Roodman gave him are
    # - decenter(dx, dy): 0.15 mm
    # - tilt(tx, ty): 10 arcsec
    # We found typical value of defocus is likely
    # - defocus(dz): 0.06 mm
    # We will draw randomly from Gaussian distributions with these values of sigma.
    ground_tilt_rms = 10. # arcsec, per component (i.e., x and y)
    ground_decenter_rms = 0.15 # mm, per commponent (k.e., dx and dy)
    ground_defocus_rms = 0.06 # mm along z

    # Struts: we choose an integer number of struts corresponding to our models.  Based on the
    # discussion about accuracy of rendering PSFs with struts for ground-based sims on
    # https://github.com/rmandelb/great3-private/pull/53, we use n_struts=0 for ground-based sims,
    # unlike our DECam optical model which has 4 struts.
    n_struts = {
        "space" : 6,
        "ground" : 0,
        }

    # This is the default.  But option here to go larger if desired.
    pad_factor = 1.5

    # For space-based PSFs only: include jitter and charge diffusion.  The plan as stated on #42,
    # based on discussion with Lance Miller and several WFIRST people, is to make jitter a Gaussian
    # with RMS of 0.005-0.015" per axis, e=0-0.3 (uniform distribution), and a random direction.  So
    # let's make a total sigma in the range sqrt(2)*(0.005-0.015).
    min_jitter_sigma = np.sqrt(2.)*0.005 # arcsec
    max_jitter_sigma = np.sqrt(2.)*0.005
    min_jitter_e = 0.
    max_jitter_e = 0.3
    # And for charge diffusion, the plan is a Gaussian with sigma=0.1-0.2pix (uniform
    # distribution), with e=0-0.2 (uniform distribution), always the same direction.  Let's use
    # always the same pixel scale for determining that sigma, 0.05".  For simplicity since we always
    # want the same direction of ellipicity, let's say that it's always in +e1.
    min_charge_sigma = 0.1*0.05 # arcsec
    max_charge_sigma = 0.2*0.05 # arcsec
    min_charge_e1 = 0.
    max_charge_e1 = 0.2

    # Atmospheric PSF info:
    atmos_scheme_pref = "atmos_psf_"
    # Draw seeing from a distribution:
    fwhm_arcsec = 0.05 + 0.10*np.arange(17)
    #Old version is commented out.  We use a new version that has only 0.5-0.85
    #freq = (0., 0., 0., 7.5, 19., 20., 17., 13., 9., 5., 3.5, 2., 1., 1., 0.5, 0.0, 0.0)
    freq = (0., 0., 0., 0.0, 0.0, 20., 17., 13., 9., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    fwhm_scatter = 0.1 # fractional scatter in FWHM for tiles within a single field
    # Later on we will draw from this distribution using:
    # dd = galsim.DistDeviate(uniform_deviate, 
    #                        function=galsim.LookupTable(fwhm_arcsec, freq, interpolant='linear'))

    # For atmospheric PSF anisotropy, we have some basic numbers that we will use to make PSF
    # anisotropy fields.
    # For A, take a range based on imSim for 20s exposures: 1e-4 to 8e-4 (could debate about what
    # distribution to use, but for now, flat isn't too bad).  This has to be rescaled according to
    # the exposure time and diameter, the latter of which we fix now to 4m.
    min_A = 1.e-4 * (6.7/4.) # will eventually be x 20 / t_exp
    max_A = 8.e-4 * (6.7/4.)
    min_t_exp = 30. # s
    max_t_exp = 180. # s
    # For theta_0, take a range based on imSim, from 0.1-2 degree (again, flat distribution).
    min_theta_0 = 360. # arcsec
    max_theta_0 = 7200. # arcsec
    # Set up empty cache for optical PSF and atmospheric PSF information, to be used for the entire
    # field.
    cached_atmos = None
    cached_optics = None
    cached_field = None
    # set up grid needed to make filenames and interpolate
    log_min_theta_0 = np.log10(min_theta_0)
    log_max_theta_0 = np.log10(max_theta_0)
    n_theta_0 = 11
    dlog_theta_0 = (log_max_theta_0 - log_min_theta_0) / (n_theta_0 - 1)
    theta_0_grid = np.logspace(np.log10(min_theta_0), np.log10(max_theta_0), n_theta_0)

    def __init__(self, obs_type, multiepoch, shear_type, atmos_ps_dir):
        # First, a sanity check: we require nrows == ncols for the GalSim lensing engine (in order
        # to have a square grid), so we check for that before attempting any calculations for
        # variable PSF, ground-based sims.
        if obs_type == "ground":
            if constants.nrows != constants.ncols:
                raise NotImplementedError("Atmospheric PSF calculations require nrows=ncols")

        self.obs_type = obs_type
        self.multiepoch = multiepoch
        self.variable_psf = True # by definition, so it doesn't get passed in as an arg
        self.shear_type = shear_type
        self.atmos_ps_dir = atmos_ps_dir

        # Now that we've specified an obs_type, we can figure out the tile structure:
        self.n_tiles = self.n_tile_linear[self.obs_type]**2
        # Define coordinates for tiles, in degrees, without any subfield offset
        self.tile_size_deg = constants.image_size_deg / self.n_tile_linear[self.obs_type]
        self.tile_x_min = []
        self.tile_y_min = []
        for i_tile in range(self.n_tiles):
            self.tile_x_min.append((i_tile % self.n_tile_linear[self.obs_type]) *
                                   self.tile_size_deg)
            self.tile_y_min.append((i_tile / self.n_tile_linear[self.obs_type]) *
                                   self.tile_size_deg)

    def generateFieldParameters(self, rng, field_index):
        # Define some constants that we will need throughout.
        if self.multiepoch:
            n_epochs = constants.n_epochs
        else:
            n_epochs = 1

        # Make RNG
        uniform_deviate = galsim.UniformDeviate(rng)
        gaussian_deviate = galsim.GaussianDeviate(rng)
        if self.obs_type == "ground":
            dist_deviate = galsim.DistDeviate(
                rng,
                function = galsim.LookupTable(self.fwhm_arcsec, self.freq, interpolant='linear')
                )

        # Get the stellar density to be used for the field, uniformly distributed between the values
        # in constants.py.  Use this to also calculate the number of stars to be used in a subfield
        # image.
        star_density = constants.min_star_density + \
            uniform_deviate()*(constants.max_star_density - constants.min_star_density)
        n_star_ideal = ( 3600. * star_density * constants.image_size_deg**2 /
                         constants.n_subfields_per_field[self.shear_type][self.variable_psf] )
        n_star_linear = int(np.ceil(np.sqrt(n_star_ideal)))


        # Get the basic parameters per epoch.  The per-observation parameters have two dimensions:
        # n_tiles x n_epochs.  The ones that are overall the same for a given observation type are
        # just numbers that get tacked onto the dict at the very end.
        if self.obs_type == "ground":
            atmos_psf_fwhm = np.zeros((self.n_tiles, n_epochs))
            atmos_psf_pk_amp = np.zeros((self.n_tiles, n_epochs))
            atmos_psf_pk_theta0 = np.zeros((self.n_tiles, n_epochs))
            x_tilt = np.zeros((self.n_tiles, n_epochs))
            y_tilt = np.zeros((self.n_tiles, n_epochs))
            x_decenter = np.zeros((self.n_tiles, n_epochs))
            y_decenter = np.zeros((self.n_tiles, n_epochs))
            dz = np.zeros((self.n_tiles, n_epochs))
        else:
            jitter_sigma = np.zeros((self.n_tiles, n_epochs))
            jitter_e = np.zeros((self.n_tiles, n_epochs))
            jitter_beta = np.zeros((self.n_tiles, n_epochs))
            charge_sigma = np.zeros((self.n_tiles, n_epochs))
            charge_e1 = np.zeros((self.n_tiles, n_epochs))
            additional_aber_seed = np.zeros((self.n_tiles, n_epochs))

        if self.obs_type == "ground":
            # Define exposure time for those exposures in this entire field, not per tile
            t_exp = uniform_deviate() * (self.max_t_exp - self.min_t_exp) + self.min_t_exp

        # Get numbers for various parameters that will be used to make a variable OpticalPSFModel
        # that is different per tile and epoch:
        for i_epoch in range(n_epochs):
            for i_tile in range(self.n_tiles):

                if self.obs_type == "ground":
                    # Some atmosphere parameters: PSF FWHM, P(k) numbers
                    if i_tile == 0:
                        atmos_psf_fwhm[i_tile, i_epoch] = dist_deviate()
                    else:
                        atmos_psf_fwhm[i_tile, i_epoch] = atmos_psf_fwhm[0, i_epoch] * \
                            (1.0 + self.fwhm_scatter*gaussian_deviate())
                    atmos_psf_pk_amp[i_tile, i_epoch] = \
                        (uniform_deviate()*(self.max_A-self.min_A) + self.min_A) * (20./t_exp)
                    atmos_psf_pk_theta0[i_tile, i_epoch] = \
                        uniform_deviate()*(self.max_theta_0-self.min_theta_0) + self.min_theta_0
                    # Stochastic parameters for optical PSF model
                    x_tilt[i_tile, i_epoch] = self.ground_tilt_rms*gaussian_deviate()
                    y_tilt[i_tile, i_epoch] = self.ground_tilt_rms*gaussian_deviate()
                    x_decenter[i_tile, i_epoch] = self.ground_decenter_rms*gaussian_deviate()
                    y_decenter[i_tile, i_epoch] = self.ground_decenter_rms*gaussian_deviate()
                    dz[i_tile, i_epoch] = self.ground_defocus_rms*gaussian_deviate()
                else:
                    # For space, we just need the stochastic parameters of the optical PSF model,
                    # the jitter, and the charge diffusion.  We also have to define a positive
                    # random seed to pass in for the additional aberrations.  We'll just make this
                    # some random number between 1 and 1000000
                    jitter_sigma[i_tile, i_epoch] = \
                        uniform_deviate() * (self.max_jitter_sigma - \
                                                 self.min_jitter_sigma) + self.min_jitter_sigma
                    jitter_e[i_tile, i_epoch] = \
                        uniform_deviate() * (self.max_jitter_e - \
                                                 self.min_jitter_e) + self.min_jitter_e
                    jitter_beta[i_tile, i_epoch] = uniform_deviate()*180.0
                    charge_sigma[i_tile, i_epoch] = \
                        uniform_deviate() * (self.max_charge_sigma - \
                                                 self.min_charge_sigma) + self.min_charge_sigma
                    charge_e1[i_tile, i_epoch] = \
                        uniform_deviate() * (self.max_charge_e1 - \
                                                 self.min_charge_e1) + self.min_charge_e1
                    additional_aber_seed[i_tile, i_epoch] = int(uniform_deviate()*1000000)+1

        # Make schema for catalog.
        schema = [("opt_psf_lam_over_diam", float),
                  ("opt_psf_obscuration", float),
                  ("opt_psf_n_struts", int),
                  ("opt_psf_strut_angle", float),
                  ("opt_psf_pad_factor", float)]
        for aber in self.use_aber:
            schema.append((self.opt_schema_pref+aber, float))
        if self.obs_type == "ground":
            schema.append(("atmos_psf_fwhm", float))
            schema.append(("atmos_psf_e1", float))
            schema.append(("atmos_psf_e2", float))
        else:
            schema.append(("opt_psf_jitter_sigma", float))
            schema.append(("opt_psf_jitter_e", float))
            schema.append(("opt_psf_jitter_beta", float))
            schema.append(("opt_psf_charge_sigma", float))
            schema.append(("opt_psf_charge_e1", float))
        # Save information about tile and position on tile.
        schema.append(("tile_index", int)) # in 1d tile listing
        schema.append(("x_tile_index", int)) # in x direction
        schema.append(("y_tile_index", int)) # in y direction
        schema.append(("tile_x_pos_deg", float)) # x offset from bottom left corner of tile, deg
        schema.append(("tile_y_pos_deg", float)) # y offset from bottom left corner of tile, deg
        # These are the position within the field, in degrees; for stars, these are the only way
        # to get the true position (since the gridded image positions mean nothing) whereas for
        # galaxies these should be computable by participants based on position on the grid +
        # subfield offsets.
        schema.append(("x_field_true_deg", float))
        schema.append(("y_field_true_deg", float))
        # Finally, our schema must include a star S/N
        schema.append(("star_snr", float))

        # Prepare the dict that this function must return.
        psf_dict = dict(schema=schema, lam_over_diam=self.lam_over_diam[self.obs_type],
                        obscuration=self.obscuration[self.obs_type],
                        n_struts=self.n_struts[self.obs_type], strut_angle=0.,
                        pad_factor=self.pad_factor,
                        star_density=star_density, n_star_linear=n_star_linear)
        if self.obs_type == "ground":
            psf_dict["atmos_psf_fwhm"]=atmos_psf_fwhm
            psf_dict["atmos_psf_pk_amp"]=atmos_psf_pk_amp
            psf_dict["atmos_psf_pk_theta0"]=atmos_psf_pk_theta0
            psf_dict["x_tilt"]=x_tilt
            psf_dict["y_tilt"]=y_tilt
            psf_dict["x_decenter"]=x_decenter
            psf_dict["y_decenter"]=y_decenter
            psf_dict["dz"]=dz
        else:
            psf_dict["jitter_sigma"]=jitter_sigma
            psf_dict["jitter_e"]=jitter_e
            psf_dict["jitter_beta"]=jitter_beta
            psf_dict["charge_sigma"]=charge_sigma
            psf_dict["charge_e1"]=charge_e1
            psf_dict["additional_aber_seed"]=additional_aber_seed

        return psf_dict

    def generateEpochParameters(self, rng, subfield_index, epoch_index, field_parameters):
        # At the field level we already pre-determined the results for the different epochs, so just
        # return the result for this epoch.

        psf_dict = \
            dict(lam_over_diam=field_parameters["lam_over_diam"],
                 obscuration=field_parameters["obscuration"],
                 n_struts=field_parameters["n_struts"],
                 strut_angle=field_parameters["strut_angle"],
                 pad_factor=field_parameters["pad_factor"],
                 star_density=field_parameters["star_density"],
                 n_star_linear=field_parameters["n_star_linear"],
                 schema=field_parameters["schema"])
        if self.obs_type == "ground":
            psf_dict["atmos_psf_fwhm"] = field_parameters["atmos_psf_fwhm"][:,epoch_index]
            psf_dict["atmos_psf_pk_amp"] = field_parameters["atmos_psf_pk_amp"][:,epoch_index]
            psf_dict["atmos_psf_pk_theta0"] = field_parameters["atmos_psf_pk_theta0"][:,epoch_index]
            psf_dict["x_tilt"] = field_parameters["x_tilt"][:,epoch_index]
            psf_dict["y_tilt"] = field_parameters["y_tilt"][:,epoch_index]
            psf_dict["x_decenter"] = field_parameters["x_decenter"][:,epoch_index]
            psf_dict["y_decenter"] = field_parameters["y_decenter"][:,epoch_index]
            psf_dict["dz"] = field_parameters["dz"][:,epoch_index]
        else:
            psf_dict["jitter_sigma"] = field_parameters["jitter_sigma"][:,epoch_index]
            psf_dict["jitter_e"] = field_parameters["jitter_e"][:,epoch_index]
            psf_dict["jitter_beta"] = field_parameters["jitter_beta"][:,epoch_index]
            psf_dict["charge_sigma"] = field_parameters["charge_sigma"][:,epoch_index]
            psf_dict["charge_e1"] = field_parameters["charge_e1"][:,epoch_index]
            psf_dict["additional_aber_seed"] = \
                field_parameters["additional_aber_seed"][:,epoch_index]

        return psf_dict

    def generateCatalog(self, rng, catalog, parameters, offsets, normalized):
        # Here's where most of the work happens.
        #
        # Need to import the optical PSF modules for space or ground depending on obs_type.
        import sys
        sys.path.append('../psfs')
        if self.obs_type == "ground":
            import ground_optical_psf
        else:
            import space_optical_psf

        # If this is the first subfield in a field, then loop over tiles, making an OpticalPSFModel
        # for each one.  Save it so we have a list of OpticalPSFModels for which we can get PSF info
        # as a function of position.  We should save one for each epoch as well, so generateCatalog
        # has a lot of caching to do for the first subfield in a field (or, the first subfield that
        # we choose to run), but then for the other subfields in the field, it can just use the
        # cached values.
        epoch_index = parameters["epoch_index"]
        field_index = parameters["field_index"]
        psf_parameters = parameters["psf"]

        # Check if there is already a cache in place.  If not, prepare to set one up.  If yes, then
        # just construct the index we'll use to access it.
        # This is the index in the list we'll use to save the lists of OpticalPSFModel or
        # galsim.PowerSpectrum (for atmosphere) objects.  The list index for all subfields/fields
        # within the field comes from the epoch index.  Once we access based on the list_index, then
        # we get a list of the results for each tile.
        list_index = epoch_index
        if self.cached_field != field_index:
            self.cached_field = field_index
            # First, do optical PSF stuff, which has to happen regardless of whether it's ground or
            # space.  Try to do this as much as possible in a way that is independent of observation
            # type, though at some point, different functions must be called.
            self.cached_optics = []
            tmp_list = []
            for i_tile in range(self.n_tiles):
                if self.obs_type == "ground":
                    # Note: strut information is not needed, given that we don't actually make
                    # the PSFs, we just get the level of aberrations.
                    new_model = \
                        ground_optical_psf.OpticalPSFModel(
                            position_list_filename = \
                            '../psfs/ground_optical_psf_zernike_coefficients_41x41/ZEMAXInput.dat',
                            lam = self.lam[self.obs_type],
                            diameter = self.diam[self.obs_type],
                            obscuration = self.obscuration[self.obs_type],
                            pad_factor = self.pad_factor,
                            dz = psf_parameters["dz"][i_tile],
                            dx = psf_parameters["x_decenter"][i_tile],
                            dy = psf_parameters["y_decenter"][i_tile],
                            tx = psf_parameters["x_tilt"][i_tile],
                            ty = psf_parameters["y_tilt"][i_tile])
                else:
                    # Note: strut information is not needed, given that we don't actually make
                    # the PSFs, we just get the level of aberrations.
                    new_model = \
                        space_optical_psf.OpticalPSFModel(
                            filename = \
                            '../psfs/afta_wfirst_example_psf_exaggerated.fields_and_coefs.fits',
                            lam = self.lam[self.obs_type],
                            diameter = self.diam[self.obs_type],
                            obscuration = self.obscuration[self.obs_type],
                            pad_factor = self.pad_factor,
                            rms = self.space_rms_additional_aber,
                            seed = int(psf_parameters["additional_aber_seed"][i_tile]))
                tmp_list.append(new_model)
            self.cached_optics.append(tmp_list)

            # And now, for ground-based observations, we must build up the cache for the atmospheric
            # PSF galsim.PowerSpectrum objects.
            if self.obs_type == "ground":
                self.cached_atmos = []
                tmp_list = [] # list for galsim.PowerSpectrum instances

                import os
                # this is actually a list (all tiles)
                atmos_psf_pk_amp = psf_parameters["atmos_psf_pk_amp"] 
                atmos_psf_pk_theta0 = psf_parameters["atmos_psf_pk_theta0"]

                # To begin, we have to read in some tabulated P(k), one per tile.
                # Could cache these to avoid rereading, but it's really not a lot of overhead, so
                # going the lazy route for now.
                for i_tile in range(self.n_tiles):
                    # Find nearest theta_0 value for which we have tabulated results.
                    log_theta_0 = np.log10(atmos_psf_pk_theta0[i_tile])
                    theta_0_index = \
                        int(round((log_theta_0 - self.log_min_theta_0)/self.dlog_theta_0))
                    theta_0_str = str(int(round(self.theta_0_grid[theta_0_index])))
                    infile = os.path.join(self.atmos_ps_dir, 'Pk'+theta_0_str+'.dat')
                    dat = np.loadtxt(infile).transpose()
                    # get k, P(k) in 1/arcsec, arcsec^2.  The latter should be used for both the E
                    # and B mode power.
                    k = dat[0]
                    pk = atmos_psf_pk_amp[i_tile]*np.pi*dat[1]
                    tab_pk = galsim.LookupTable(k, pk, x_log=True, f_log=True)
                    ps = galsim.PowerSpectrum(tab_pk, tab_pk, units = galsim.arcsec)
                    # Define the grid on which we want to get the PSF anisotropies.
                    # See comments in shear.py, VariableShearBuilder for more explanation.
                    # We just have an additional parameter, tile_fac, that accounts for the fact
                    # that each tile is some fraction of the overall FOV.
                    tile_fac = 1. / self.n_tile_linear[self.obs_type]
                    # We've assumed a square grid.  In the __init__ routine for this builder, we
                    # have already checked that constants.nrows == constants.ncols, so that's okay.
                    n_grid = int(np.ceil(constants.subfield_grid_subsampling * constants.nrows * tile_fac))
                    grid_spacing = self.tile_size_deg / n_grid
                    grid_center_zerod = 0.5 * self.tile_size_deg
                    ps.buildGrid(grid_spacing = grid_spacing,
                                 ngrid = n_grid,
                                 units = galsim.degrees,
                                 rng = rng,
                                 center = (grid_center_zerod, grid_center_zerod),
                                 kmin_factor = 3)
                    tmp_list.append(ps)
                self.cached_atmos.append(tmp_list)

        # Now we finished building up the cache (if it was necessary; otherwise we can just access
        # the right entry for this subfield / epoch using list_index which was calculated above).
        # To do anything with this, like get the optical PSF and atmospheric PSF (for ground) model
        # parameters at the position of each object, we first need to do some calculations to get
        # positions for the galaxies within the field in degrees.  For stars, this was already done
        # for us, so we should check whether those fields in the catalog are populated, and only do
        # the calculation if they were not.
        if not np.any(catalog["x_field_true_deg"] != 0):
            # We must take into account the position of the galaxy within the subfield, as well as
            # subfield offsets within the field.  We are not taking into account sub-pixel xdither
            # and ydither between epochs, i.e., the PSF model is assumed to be constant on pixel
            # scales.  These positions within the field are used to choose (a) which tile we are on
            # and (b) location within the tile.  Then, for optical PSFs, we apply a stretching
            # factor to account for the fact that we're kind of faking the size of the optical PSF
            # model to fit conveniently onto these tiles.  Code is based on that in shear.py,
            # VariableShearBuilder.
            xsize = constants.xsize[self.obs_type][self.multiepoch]
            ysize = constants.ysize[self.obs_type][self.multiepoch]
            # Below, we define object indices that range from 0 to constants.nrows-1
            x_ind = (catalog["x"]+1+0.5*xsize)/xsize-1
            y_ind = (catalog["y"]+1+0.5*ysize)/ysize-1
            # Turn this into (x, y) positions within the subfield, in degrees.
            x_pos = x_ind * constants.image_size_deg / constants.nrows
            y_pos = y_ind * constants.image_size_deg / constants.ncols
            # But now we have to add the subfield offset.  These are calculated as a fraction of the
            # separation between galaxies, so we have to convert to degrees.
            x_pos += offsets[0] * constants.image_size_deg / constants.nrows
            y_pos += offsets[1] * constants.image_size_deg / constants.ncols
            catalog["x_field_true_deg"] = x_pos
            catalog["y_field_true_deg"] = y_pos
        else:
            # If this is a star catalog that already has true field positions, simply grab from the
            # catalog for later use.
            x_pos = catalog["x_field_true_deg"]
            y_pos = catalog["y_field_true_deg"]

        # Get x, y tile indices
        x_tile_ind = (x_pos / self.tile_size_deg).astype(int)
        y_tile_ind = (y_pos / self.tile_size_deg).astype(int)
        catalog["x_tile_index"] = x_tile_ind
        catalog["y_tile_index"] = y_tile_ind
        # Convert to our single tile index (as defined for each object)
        tile_ind = x_tile_ind.astype(int) + self.n_tile_linear[self.obs_type]*y_tile_ind.astype(int)
        catalog["tile_index"] = tile_ind
        cat_indices = np.arange(len(tile_ind))
        # Do the calculations on a tile-by-tile basis
        optical_psf_model_list = self.cached_optics[list_index]
        for iter_tile_ind in range(max(tile_ind)+1):
            if len(tile_ind[tile_ind == iter_tile_ind]) > 0:
                # take the OpticalPSFModel object
                optical_psf_model = optical_psf_model_list[iter_tile_ind]

                # First, if we are making space data, we should assign jitter and charge diffusion
                # parameters based on which tile we live on.  These numbers are the same for all
                # objects on the tile.
                if self.obs_type == 'space':
                    catalog["opt_psf_jitter_sigma"][tile_ind == iter_tile_ind] = \
                        psf_parameters["jitter_sigma"][iter_tile_ind]
                    catalog["opt_psf_jitter_e"][tile_ind == iter_tile_ind] = \
                        psf_parameters["jitter_e"][iter_tile_ind]
                    catalog["opt_psf_jitter_beta"][tile_ind == iter_tile_ind] = \
                        psf_parameters["jitter_beta"][iter_tile_ind]
                    catalog["opt_psf_charge_sigma"][tile_ind == iter_tile_ind] = \
                        psf_parameters["charge_sigma"][iter_tile_ind]
                    catalog["opt_psf_charge_e1"][tile_ind == iter_tile_ind] = \
                        psf_parameters["charge_e1"][iter_tile_ind]

                # Then we get the Zernike coefficients which vary for all objects within the tile.
                # For the objects on this tile, find the position with respect to the lower left
                # corner of the tile.
                tile_x_pos = x_pos[tile_ind == iter_tile_ind] - self.tile_x_min[iter_tile_ind]
                tile_y_pos = y_pos[tile_ind == iter_tile_ind] - self.tile_y_min[iter_tile_ind]
                catalog["tile_x_pos_deg"][tile_ind == iter_tile_ind] = tile_x_pos
                catalog["tile_y_pos_deg"][tile_ind == iter_tile_ind] = tile_y_pos

                # Apply stretching factor so that all the data fit into the optical PSF model.  (We
                # pretend the tiles are bigger than we are, so to interpolate on the tile, we must
                # decrease the distance with respect to the lower left corner of the tile.)
                tile_x_pos *= (optical_psf_model.xmax - optical_psf_model.xmin) / self.tile_size_deg
                tile_y_pos *= (optical_psf_model.ymax - optical_psf_model.ymin) / self.tile_size_deg
                # Add tile (x_min, ymin) to get position on the tile.
                tile_x_pos += optical_psf_model.xmin
                tile_y_pos += optical_psf_model.ymin
                # Interpolate to get the OpticalPSF parameters.  Have to do this one by one, so
                # track indices carefully.
                tile_cat_indices = cat_indices[tile_ind == iter_tile_ind]
                for ind_on_tile in range(len(tile_cat_indices)):
                    coefs = optical_psf_model.get_zernike_coefficients(tile_x_pos[ind_on_tile],
                                                                       tile_y_pos[ind_on_tile])
                    this_cat_ind = tile_cat_indices[ind_on_tile]
                    catalog[self.opt_schema_pref+"defocus"][this_cat_ind] = coefs[0]
                    catalog[self.opt_schema_pref+"astig1"][this_cat_ind] = coefs[1]
                    catalog[self.opt_schema_pref+"astig2"][this_cat_ind] = coefs[2]
                    catalog[self.opt_schema_pref+"coma1"][this_cat_ind] = coefs[3]
                    catalog[self.opt_schema_pref+"coma2"][this_cat_ind] = coefs[4]
                    catalog[self.opt_schema_pref+"trefoil1"][this_cat_ind] = coefs[5]
                    catalog[self.opt_schema_pref+"trefoil2"][this_cat_ind] = coefs[6]
                    catalog[self.opt_schema_pref+"spher"][this_cat_ind] = coefs[7]

        # For ground-based PSFs, we use the cached galsim.PowerSpectrum objects, one per
        # tile.  They each have a saved grid of shears that we can use to make PSF
        # ellipticities and modified sizes as a function of position.  This has to include all
        # possible positions in all subfields.
        if self.obs_type == "ground":
            for iter_tile_ind in range(max(tile_ind)+1):
                if len(tile_ind[tile_ind == iter_tile_ind]) > 0:
                    tile_ps = self.cached_atmos[list_index][iter_tile_ind]
                    tile_x_pos = catalog["tile_x_pos_deg"][tile_ind == iter_tile_ind]
                    tile_y_pos = catalog["tile_y_pos_deg"][tile_ind == iter_tile_ind]
                    g1, g2 = tile_ps.getShear(pos = (tile_x_pos, tile_y_pos),
                                              units = galsim.degrees)
                    catalog["atmos_psf_e1"][tile_ind == iter_tile_ind] = 2*g1
                    catalog["atmos_psf_e2"][tile_ind == iter_tile_ind] = 2*g2
                    kappa = tile_ps.getConvergence(pos = (tile_x_pos, tile_y_pos),
                                                   units = galsim.degrees)
                    kappa /= 2.
                    mu = 1./((1.-kappa)**2 - (g1**2 + g2**2))
                    catalog["atmos_psf_fwhm"][tile_ind == iter_tile_ind] = \
                        psf_parameters["atmos_psf_fwhm"][iter_tile_ind] * mu
        
        # Finally, we can save the information that goes in the catalog that is the same for all
        # objects.  We need this to make the optical PSF either for space or for ground.
        for record in catalog:
            record["opt_psf_lam_over_diam"] = psf_parameters["lam_over_diam"]
            record["opt_psf_obscuration"] = psf_parameters["obscuration"]
            record["opt_psf_n_struts"] = psf_parameters["n_struts"]
            record["opt_psf_strut_angle"] = psf_parameters["strut_angle"]
            record["opt_psf_pad_factor"] = psf_parameters["pad_factor"]

    def makeConfigDict(self):
        d = {
            'type' : 'OpticalPSF',
            'lam_over_diam' : { 'type' : 'Catalog', 'col' : 'opt_psf_lam_over_diam' },
            'obscuration' : { 'type' : 'Catalog', 'col' : 'opt_psf_obscuration' },
            'nstruts' : { 'type' : 'Catalog', 'col' : 'opt_psf_n_struts' },
            'strut_angle' : { 'type' : 'Deg',
                              'theta' : { 'type' : 'Catalog', 'col' : 'opt_psf_strut_angle' }
                            },
            'defocus' : { 'type' : 'Catalog', 'col' : 'opt_psf_defocus' },
            'astig1' : { 'type' : 'Catalog', 'col' : 'opt_psf_astig1' },
            'astig2' : { 'type' : 'Catalog', 'col' : 'opt_psf_astig2' },
            'coma1' : { 'type' : 'Catalog', 'col' : 'opt_psf_coma1' },
            'coma2' : { 'type' : 'Catalog', 'col' : 'opt_psf_coma2' },
            'trefoil1' : { 'type' : 'Catalog', 'col' : 'opt_psf_trefoil1' },
            'trefoil2' : { 'type' : 'Catalog', 'col' : 'opt_psf_trefoil2' },
            'spher' : { 'type' : 'Catalog', 'col' : 'opt_psf_spher' },
            'pad_factor' : { 'type' : 'Catalog', 'col' : 'opt_psf_pad_factor' },
            'suppress_warning' : True
        }
        if self.obs_type == 'ground':
            d = {
                'type' : 'Convolve',
                'items' : [
                    {
                        'type' : 'Kolmogorov',
                        'fwhm' : { 'type' : 'Catalog', 'col' : 'atmos_psf_fwhm' },
                        'ellip' : { 
                            'type' : 'E1E2', 
                            'e1' : { 'type' : 'Catalog', 'col' : 'atmos_psf_e1' },
                            'e2' : { 'type' : 'Catalog', 'col' : 'atmos_psf_e2' }
                        }
                    },
                    d  # The above constructed OpticalPSF
                ]
            } 
        else:
            d = {
                'type' : 'Convolve',
                'items' : [
                    {
                        'type' : 'Gaussian',
                        'sigma' : { 'type' : 'Catalog', 'col' : 'opt_psf_jitter_sigma' },
                        'ellip' : { 
                            'type' : 'EBeta', 
                            'e' : { 'type' : 'Catalog', 'col' : 'opt_psf_jitter_e' },
                            'beta' : { 'type' : 'Deg',
                                       'theta' : { 'type' : 'Catalog', 
                                                   'col' : 'opt_psf_jitter_beta' } }
                        }
                    },
                    {
                        'type' : 'Gaussian',
                        'sigma' : { 'type' : 'Catalog', 'col' : 'opt_psf_charge_sigma' },
                        'ellip' : { 
                            'type' : 'E1E2',
                            'e1' : { 'type' : 'Catalog', 'col' : 'opt_psf_charge_e1' },
                            'e2' : 0
                        }
                    },
                    d  # The above constructed OpticalPSF
                ]
            } 

        return d

    def makeGalSimObject(self, record, parameters):
        # Size is specified in arcsec.

        # Choose which aberrations to use based on what should be in catalog
        # We have to start by making a dict with all aberrations that are nonzero.
        aber_dict = dict()
        # Then populate those that we use.
        for aber in self.use_aber:
            aber_dict[aber] = record[self.opt_schema_pref+aber]

        optical_psf = galsim.OpticalPSF(record["opt_psf_lam_over_diam"],
                                        obscuration=record["opt_psf_obscuration"],
                                        nstruts=record["opt_psf_n_struts"],
                                        strut_angle=record["opt_psf_strut_angle"]*galsim.degrees,
                                        pad_factor=record["opt_psf_pad_factor"],
                                        suppress_warning=True,
                                        **aber_dict)
        if self.obs_type == "space":
            jitter_psf = galsim.Gaussian(sigma=record["opt_psf_jitter_sigma"])
            e = record["opt_psf_jitter_e"]
            if e > 0.:
                jitter_psf.applyShear(e=e,
                                      beta = record["opt_psf_jitter_beta"]*galsim.degrees)
            charge_psf = galsim.Gaussian(sigma=record["opt_psf_charge_sigma"])
            e1 = record["opt_psf_charge_e1"]
            if e1 > 0.:
                charge_psf.applyShear(e1=e1)
            return galsim.Convolve(jitter_psf, charge_psf, optical_psf)
        else:
            atmos_psf = galsim.Kolmogorov(fwhm = record["atmos_psf_fwhm"])
            atmos_psf.applyShear(e1 = record["atmos_psf_e1"],
                                 e2 = record["atmos_psf_e2"])
            return galsim.Convolve(atmos_psf, optical_psf)
