import os
import datetime
import getpass
import numpy
import galsim

import great3sims.mapper
import great3sims.psf
import great3sims.noise
import great3sims.shear
import great3sims.galaxies
from . import constants

class SimBuilder(object):

    # Per-experiment parameters; should be overridden by per-experiment subclasses (see customize())
    experiment = None
    real_galaxy = False
    variable_psf = False
    multiepoch = False

    @staticmethod
    def customize(experiment, real_galaxy=None, variable_psf=None, multiepoch=None):
        """Create a custom subclass of SimBuilder with class variables overridden from
        their default (control experiment) values.
        """
        cls = type(experiment, (SimBuilder,), dict())
        cls.experiment = experiment
        if real_galaxy is not None:
            cls.real_galaxy = real_galaxy
        if variable_psf is not None:
            cls.variable_psf = variable_psf
        if multiepoch is not None:
            cls.multiepoch = multiepoch
        return cls

    def __init__(self, root, obs_type, shear_type, gal_dir, ps_dir, atmos_ps_dir, public_dir,
                 truth_dir, preload=False, nproc=-1):
        """Initialize a builder for the given obs_type and shear_type.

        @param[in] root         Root directory for generated files
        @param[in] obs_type     Type of observation to simulate: either "ground" or "space"
        @param[in] shear_type   Type of shear field to generate: either "constant" or "variable"
        @param[in] gal_dir      Directory with real galaxy catalog information
        @param[in] ps_dir       Directory with tabulated iCosmo shear power spectra.
        @param[in] atmos_ps_dir Directory with tabulated atmospheric PSF anisotropy power spectra.
        @param[in] public_dir   Directory for placing files to be distributed publicly.
        @param[in] truth_dir    Directory containing files used for metric evaluation.
        @param[in] preload      Preload the RealGalaxyCatalog images to speed up generation of large
                                numbers of real galaxies?  Note that for parametric galaxy branches,
                                the catalog is never preloaded. (default = False)
        @param[in] nproc        How many processes to use in the config file.  (default = -1)
        """
        self.obs_type = obs_type
        self.shear_type = shear_type
        self.public_dir = public_dir
        self.truth_dir = truth_dir
        self.preload = preload
        self.psf_builder = great3sims.psf.makeBuilder(obs_type=obs_type,
                                                      variable_psf=self.variable_psf,
                                                      multiepoch=self.multiepoch,
                                                      shear_type=self.shear_type,
                                                      atmos_ps_dir=atmos_ps_dir)
        self.shear_builder = great3sims.shear.makeBuilder(shear_type=shear_type, obs_type=obs_type,
                                                          multiepoch=self.multiepoch, ps_dir=ps_dir)
        self.galaxy_builder = great3sims.galaxies.makeBuilder(real_galaxy=self.real_galaxy,
                                                              obs_type=obs_type,
                                                              shear_type=shear_type,
                                                              multiepoch=self.multiepoch,
                                                              gal_dir=gal_dir,
                                                              preload=preload)
        self.noise_builder = great3sims.noise.makeBuilder(obs_type=obs_type,
                                                          multiepoch=self.multiepoch,
                                                          variable_psf = self.variable_psf)
        self.mapper = great3sims.mapper.Mapper(root, self.experiment, obs_type, shear_type)
        self.n_epochs = constants.n_epochs if self.multiepoch else 1
        self.nproc = nproc

    def writeParameters(self, seed):
        """Generate and write the metaparameters of the builder.

        This creates four types of files:
          - a single dict of parameters for the branch (=experiment+obs_type+shear_type combination)
          - a dict of parameters for each field
          - a dict of parameters for each subfield
          - a dict of parameters for each subfield+epoch combination

        The given seed (an integer) is the first of a sequence of sequential integers used to seed
        random number generators for different steps; these seeds will be saved in the metaparameter
        dictionaries to make simulation generation deterministic (and suitable for parallelization)
        after this step.  The last seed number used is returned.
        """
        # First get some basic quantities for this branch: metadata, pixel scale, and the grouping
        # of subfields into fields with identical PSF and shear fields
        metadata = {"timestamp": str(datetime.datetime.now()), "user": getpass.getuser()}
        pixel_scale = constants.pixel_scale[self.obs_type][self.multiepoch]
        n_subfields = constants.n_subfields
        n_subfields_per_field = constants.n_subfields_per_field[self.shear_type][self.variable_psf]
        if n_subfields % n_subfields_per_field != 0:
            raise RuntimeError("%d subfields does not divide evenly into %d fields!" % \
                                   (n_subfields, n_subfields_per_field) )
        n_fields = n_subfields / n_subfields_per_field

        # Also check some things about the deep fields.
        n_deep_subfields = constants.n_deep_subfields
        if n_deep_subfields % n_subfields_per_field != 0:
            raise RuntimeError("%d deep subfields does not divide evenly into %d fields!" %\
                                   (n_deep_subfields, n_subfields_per_field) )
        n_reg_subfields = n_subfields - n_deep_subfields
        if constants.deep_frac > float(n_deep_subfields) / n_reg_subfields:
            raise RuntimeError("%d subfields is insufficient for %f deep fraction!" % \
                                   (n_deep_subfields, constants.deep_frac) )
        n_deep_fields = n_deep_subfields / n_subfields_per_field
        n_reg_fields = n_fields - n_deep_fields
        # A note regarding these deep fields: in most respects they are to be like the regular
        # fields.  For example, galaxy selection is the same: we only make images of galaxies that
        # would be observed at S/N>=20 in the regular fields.  PSF and shear determination is
        # carried out in the same way.  The only difference comes in the noise builder, when
        # choosing the level of noise to add (all the way down at the level of epochs).  Even then,
        # we use the same typical noise variance for all stored quantities - so that galaxy
        # selection can use them - and use a multiplicative factor to tell it to add less noise.

        # put together the basic parameters to be stored in metaparameters file
        self.parameters = {"metadata": metadata, "seed": seed, "pixel_scale": pixel_scale,
                           "n_fields": n_fields}
        self.mapper.write(self.parameters, 'parameters', self.parameters)
        rng = galsim.UniformDeviate(seed)
        base_schema = [("index", int), ("x", int), ("y", int),
                       ("xshift", float), ("yshift", float),
                       ("xmin", int), ("xmax", int), ("ymin", int), ("ymax", int)]
        shear_schema = [("g1", float), ("g2", float), ("mu", float),
                       ("x_field_pos", int), ("y_field_pos", int), ("ID", int)]
        seed += 1  # we could also draw random integers to set seeds, but based on a discussion, it
                   # does not seem necessary

        # Start a separate sequence for noise_seed.  The second number is just a random
        # number that is much larger than the plausible number of items in the original
        # seed sequence.  (It is phi * 10^6)
        noise_seed = seed + 1618033

        # Now we begin to loop over successively smaller units - starting with field, then subfield,
        # then epoch.
        for field_index in xrange(n_fields):
            # A given field has the same shear and per-epoch PSF.
            # The builders can decide what format to use for the results of generateFieldParameters.
            # We don't actually care what that format is out here - we will just pass it along to
            # generateSubfieldParameters or generateEpochParameters.  So we require internal
            # consistency between the various generate*Parameters, but we should be able to switch
            # parameter selection between the field/subfield/epoch layer without modifying this
            # code, just modifying the builders themselves (in psf.py, shear.py, noise.py, or
            # galaxies.py).
            field_shear_parameters = self.shear_builder.generateFieldParameters(rng, field_index)
            field_psf_parameters = self.psf_builder.generateFieldParameters(rng, field_index)
            field_parameters = {
                "shear": field_shear_parameters,
                "psf": field_psf_parameters,
                "field_index": field_index,
                "subfield_offsets": self.generateSubfieldOffsets(
                                             rng, n_subfields_per_field,
                                             constants.subfield_grid_subsampling),
                "metadata": metadata,
                "seed": seed,
                }
            seed += 1
            self.mapper.write(field_parameters, 'field_parameters', field_parameters)

            field_min_subfield = field_index * n_subfields_per_field
            field_max_subfield = field_min_subfield + n_subfields_per_field - 1
            for subfield_index in xrange(field_min_subfield, field_max_subfield+1):
                # A given subfield has the same shear (already determined at field level) and
                # galaxies.  But to allow for flexibility later on, we'll have a
                # generateSubfieldParameters for the shear builder anyway; it has to take the output
                # of generateFieldParameters as an input.  For now, it's a no-op, but that might not
                # be the case later on.
                subfield_parameters = {
                    "shear": self.shear_builder.generateSubfieldParameters(rng, subfield_index,
                                                                           field_shear_parameters),
                    "galaxy": self.galaxy_builder.generateSubfieldParameters(rng, subfield_index),
                    "subfield_index": subfield_index,
                    "subfield_offset": \
                        field_parameters["subfield_offsets"][subfield_index-field_min_subfield],
                    "field_index": field_index,
                    "metadata": metadata,
                    "seed": seed,
                    }
                seed += 1
                subfield_parameters["subfield_schema"] = \
                    (base_schema + shear_schema + subfield_parameters["galaxy"]["schema"])
                self.mapper.write(subfield_parameters, 'subfield_parameters', subfield_parameters)
                for epoch_index in xrange(self.n_epochs):
                    # Each epoch has its own PSF (already determined at field level) and noise.
                    # But the galaxies and shears were determined at the subfield level, so nothing
                    # new is needed here.
                    epoch_parameters = dict(subfield_parameters)
                    epoch_parameters["psf"] = \
                        self.psf_builder.generateEpochParameters(rng, subfield_index, epoch_index,
                                                                 field_psf_parameters)
                    # Determine multiplying factor for the sky variance based on whether we are in a
                    # deep field or not.  Note that sky variance changes due to single
                    # vs. multiepoch images are handled directly in the noise_builder, so they do
                    # not need to be included here.
                    if subfield_index < n_reg_subfields:
                        noise_mult = 1.
                    else:
                        noise_mult = constants.deep_variance_mult
                    if self.obs_type == "ground":
                        epoch_parameters["noise"] = \
                            self.noise_builder.generateEpochParameters(
                                    rng, subfield_index, epoch_index,
                                    field_parameters["psf"]["atmos_psf_fwhm"], noise_mult)
                    else:
                        epoch_parameters["noise"] = \
                            self.noise_builder.generateEpochParameters(rng, subfield_index,
                                                                       epoch_index, None,
                                                                       noise_mult)
                    epoch_parameters["epoch_index"] = epoch_index
                    epoch_parameters["subfield_index"] = subfield_index
                    epoch_parameters["field_index"] = field_index
                    epoch_parameters["seed"] = seed
                    seed += 1
                    epoch_parameters["noise_seed"] = noise_seed
                    noise_seed += constants.nrows * constants.ncols
                    epoch_parameters["epoch_schema"] = (epoch_parameters["subfield_schema"]
                                                        + epoch_parameters["psf"]["schema"])
                    epoch_parameters["star_schema"] = \
                        base_schema + epoch_parameters["psf"]["schema"]
                    # xdither, ydither say the amount of dithering between epochs of this subfield.
                    # In contrast, epoch_offset (a tuple) is the amount of offsetting between this
                    # subfield and the first one in the field.  This is not truly a per-epoch
                    # parameter, however, it is included here so that the per-epoch catalog maker
                    # (specifically, for PSFs) will be able to do its job.
                    if self.multiepoch:
                        epoch_parameters["xdither"] = (2.0 * rng() - 1.0) * \
                                                      constants.epoch_shift_max
                        epoch_parameters["ydither"] = (2.0 * rng() - 1.0) * \
                                                      constants.epoch_shift_max
                    else:
                        epoch_parameters["xdither"] = 0.0
                        epoch_parameters["ydither"] = 0.0
                    epoch_parameters["epoch_offset"] = subfield_parameters["subfield_offset"]
                    self.mapper.write(epoch_parameters, 'epoch_parameters', epoch_parameters)
        return seed

    def writeSubfieldCatalog(self, subfield_index):
        """Given the subfield index, load the corresponding metaparameters and generate a catalog of
        galaxies and shear values.
        """
        subfield_parameters = self.mapper.read("subfield_parameters", subfield_index=subfield_index)
        field_parameters = self.mapper.read(
            "field_parameters",
            field_index = ( subfield_index / 
                            constants.n_subfields_per_field[self.shear_type][self.variable_psf] ) 
        )
        catalog = numpy.zeros(constants.nrows * constants.ncols,
                              dtype=numpy.dtype(subfield_parameters["subfield_schema"]))
        index = 0
        rng = galsim.UniformDeviate(subfield_parameters["seed"])
        for row in xrange(constants.nrows):
            for col in xrange(constants.ncols):
                # The numbers below are all in pixels.
                sx = (2.0*rng() - 1.0) * constants.centroid_shift_max
                sy = (2.0*rng() - 1.0) * constants.centroid_shift_max
                record = catalog[index]
                record['index'] = index
                record['xmin'] = col * constants.xsize[self.obs_type][self.multiepoch]
                record['ymin'] = row * constants.ysize[self.obs_type][self.multiepoch]
                record['xmax'] = (col + 1) * constants.xsize[self.obs_type][self.multiepoch] - 1
                record['ymax'] = (row + 1) * constants.ysize[self.obs_type][self.multiepoch] - 1
                record['x'] = (record['xmin'] + record['xmax']) / 2
                record['y'] = (record['ymin'] + record['ymax']) / 2
                record['xshift'] = sx
                record['yshift'] = sy
                index += 1

        # Determine multiplying factor for the sky variance based on whether we are in a deep field
        # or not.  Note that sky variance changes due to single vs. multiepoch images do not need to
        # be included at all at this stage, because we want to impose our cuts based on whether the
        # S/N would be 20 in a single combined image, not based on its value in the individual epoch
        # images.  Also note that while the noise variance parameter output into the
        # epoch_parameters files by the noise builder already includes this noise_mult for the deep
        # fields (and any factors due to single vs. multiepoch which we do *not* want here), we have
        # to recalculate the deep field noise multiplying factor because we're just going to use the
        # noise_builder.typical_variance which does not include any of those factors.
        if subfield_index < constants.n_subfields - constants.n_deep_subfields:
            noise_mult = 1.
        else:
            noise_mult = constants.deep_variance_mult

        # We give it a value for seeing to use when selecting galaxies.  This calculation becomes
        # more complex in the multi-epoch case since we have to decide on a relevant effective FWHM.
        if self.obs_type == "space":
            effective_seeing = None
        else:
            # Note: really we care about the full PSF FWHM, not just the atmospheric part.  However,
            # we use the seeing as a proxy for it, so we don't have to start generating images.  If
            # this seems really worrisome, we could make some simple sims, derive approximate rules
            # for total PSF size including optics as well (which will mainly affect really
            # good-seeing images), and use those instead of just the seeing.
            if not self.multiepoch and not self.variable_psf:
                # This is just a scalar value.
                effective_seeing = field_parameters["psf"]["atmos_psf_fwhm"]
            else:
                # This is a 1d numpy array of length n_epochs.
                seeing = field_parameters["psf"]["atmos_psf_fwhm"]
                effective_seeing = 1. / numpy.mean(1./seeing)
        self.galaxy_builder.generateCatalog(rng, catalog, subfield_parameters,
                                            self.noise_builder.typical_variance, noise_mult,
                                            effective_seeing)
        self.shear_builder.generateCatalog(rng, catalog, subfield_parameters["shear"],
                                           subfield_parameters["subfield_offset"], subfield_index)
        self.mapper.write(catalog, "subfield_catalog", subfield_parameters)

    def writeEpochCatalog(self, subfield_index, epoch_index):
        """Given the subfield and epoch indices, load the epoch metaparameters and a
        previously-generated subfield catalog, add per-object PSF parameter information (and
        possibly shift the centroids) to create and save an epoch catalog.
        """
        epoch_parameters = self.mapper.read('epoch_parameters', subfield_index=subfield_index,
                                            epoch_index=epoch_index)
        subfield_catalog = self.mapper.read('subfield_catalog', epoch_parameters)
        epoch_catalog = numpy.zeros(constants.nrows * constants.ncols,
                                    dtype=numpy.dtype(epoch_parameters["epoch_schema"]))
        rng = galsim.UniformDeviate(epoch_parameters["seed"])
        for name, _ in epoch_parameters["subfield_schema"]:
            epoch_catalog[name] = subfield_catalog[name]
        self.psf_builder.generateCatalog(rng, epoch_catalog, epoch_parameters,
                                         epoch_parameters["epoch_offset"], normalized=True)
        self.mapper.write(epoch_catalog, "epoch_catalog", epoch_parameters)
        xsize = constants.xsize[self.obs_type][self.multiepoch]
        ysize = constants.ysize[self.obs_type][self.multiepoch]
        if self.variable_psf:
            # We have to write catalog entries indicating the "true" position within the **field**
            # (not subfield).  We also have to write x, y entries that are used to place the objects
            # on an image grid, the parameters of which are to be defined here.  Finally, let's
            # define a DistDeviate to use to draw random values of S/N.  I decided on the function
            # to use via the super-duper-scientific method of reading points off of Chihway's second
            # dN/dmag plot in
            # https://github.com/rmandelb/great3-private/issues/2#issuecomment-20482730 for the
            # range of magnitudes where that function is roughly linear (i.e., excluding the bits
            # where it's not linear since stars are often saturated and not used for PSF
            # estimation).  On this plot, dN/dmag ~ 1000 (mag - 18).  To get dN/d(S/N), I used
            # dN/d(S/N) ~ dN/dmag dmag/d(S/N),
            # plus the fact that at fixed sky level, we can say
            # mag = -2.5 log10(S/N) + A.
            # I assumed S/N(mag=25)=25, so A=28.5.  Putting this together, we have
            # dN/d(S/N) ~ (mag - 18) / S/N
            #           ~ (-2.5 log10(S/N) + 10.5) / S/N.
            dist_deviate = galsim.DistDeviate(
                rng,
                function = lambda x : (-2.5*numpy.log10(x)+10.5)/x, x_min=25., x_max=400.)

            # First, determine how many stars, which will determine the size of the catalog to
            # write.  Note that this is for a single epoch/subfield, which means the number of stars
            # should be 1/20 of those for the entire field.
            n_star_linear = epoch_parameters["psf"]["n_star_linear"]
            star_catalog = numpy.zeros(n_star_linear * n_star_linear,
                                       dtype=numpy.dtype(epoch_parameters["star_schema"]))
            index = 0
            for row in xrange(n_star_linear):
                for col in xrange(n_star_linear):
                    # The numbers below are in pixels, and are used to define image positions.
                    sx = (2.0*rng() - 1.0) * constants.centroid_shift_max
                    sy = (2.0*rng() - 1.0) * constants.centroid_shift_max
                    record = star_catalog[index]
                    record['index'] = index
                    record['xmin'] = col * constants.xsize[self.obs_type][self.multiepoch]
                    record['ymin'] = row * constants.ysize[self.obs_type][self.multiepoch]
                    record['xmax'] = (col + 1) * constants.xsize[self.obs_type][self.multiepoch] - 1
                    record['ymax'] = (row + 1) * constants.ysize[self.obs_type][self.multiepoch] - 1
                    record['x'] = (record['xmin'] + record['xmax']) / 2
                    record['y'] = (record['ymin'] + record['ymax']) / 2
                    # Here are some numbers that we'll only compute if it's the first epoch,
                    # otherwise pull from cache:
                    if epoch_index == 0:
                        record['xshift'] = sx
                        record['yshift'] = sy
                        # But these numbers are the true positions within the field.
                        record['x_field_true_deg'] = (constants.image_size_deg-1.) * rng() + \
                            epoch_parameters["epoch_offset"][0] * constants.image_size_deg / \
                            constants.nrows
                        record['y_field_true_deg'] = (constants.image_size_deg-1.) * rng() + \
                            epoch_parameters["epoch_offset"][1] * constants.image_size_deg / \
                            constants.ncols
                        # Finally, let's make a S/N value for this star (per epoch).  In other
                        # words, the star would have S/N=record['star_snr'] in a single-epoch
                        # branch, or in each of the images in a multi-epoch branch.  This differs
                        # from our handling of galaxies, where the flux is split up among the
                        # different epochs.  The rationale (possibly weak?) is this: most people
                        # will choose some range of stars to use for PSF estimation, and put that
                        # selection into their star-finder before estimating PSF.  That means that
                        # if they have one long exposure vs. a few shorter ones, and they have to
                        # estimate the PSF per exposure, they would legitimately select different
                        # sets of stars in the two cases.  (e.g., in the one long exposure, some
                        # stars might be saturated and unusable, whereas for fewer short exposures,
                        # they could be used) We are glossing over the slight difference in star
                        # number density that should also occur, and just using the fact that the
                        # S/N distribution should be similar.
                        record['star_snr'] = dist_deviate()
                    index += 1
            if epoch_index == 0.:
                self.cached_xshift = star_catalog['xshift']
                self.cached_yshift = star_catalog['yshift']
                self.cached_x_field_true_deg = star_catalog['x_field_true_deg']
                self.cached_y_field_true_deg = star_catalog['y_field_true_deg']
                self.cached_star_snr = star_catalog['star_snr']
            else:
                star_catalog['xshift'] = self.cached_xshift
                star_catalog['yshift'] = self.cached_yshift
                star_catalog['x_field_true_deg'] = self.cached_x_field_true_deg
                star_catalog['y_field_true_deg'] = self.cached_y_field_true_deg
                star_catalog['star_snr'] = self.cached_star_snr

        else:
            star_catalog = numpy.zeros(constants.nx_constpsf * constants.ny_constpsf,
                                       dtype=numpy.dtype(epoch_parameters["star_schema"]))
            index = 0
            for yc in xrange(constants.ny_constpsf):
                for xc in xrange(constants.nx_constpsf):
                    # The numbers below are in pixel units.
                    record = star_catalog[index]
                    record["xmin"] = xc * xsize
                    record["ymin"] = yc * ysize
                    record["xmax"] = (xc + 1) * xsize - 1
                    record["ymax"] = (yc + 1) * ysize - 1
                    record["index"] = index
                    record["x"] = (record["xmin"] + record["xmax"]) / 2
                    record["y"] = (record["ymin"] + record["ymax"]) / 2
                    if index > 0:
                        sx = (2.0*rng() - 1.0) * constants.centroid_shift_max
                        sy = (2.0*rng() - 1.0) * constants.centroid_shift_max
                        record["xshift"] = sx
                        record["yshift"] = sy
                    else:
                        record["xshift"] = 0
                        record["yshift"] = 0
                    index += 1
        self.psf_builder.generateCatalog(rng, star_catalog, epoch_parameters,
                                         epoch_parameters["epoch_offset"], normalized=False)
        self.mapper.write(star_catalog, "star_catalog", epoch_parameters)

    def writeStarTestCatalog(self, subfield_min, subfield_max):
        """Given a range of subfield and epoch indices, write a test catalog for generating star
        images to check for oddities."""
        n_subfields = subfield_max + 1 - subfield_min

        # Define a final catalog with some additional information about subfield / epoch.  That
        # way if we find a problematic star image in the data cube, we know where it came from.
        epoch_parameters = self.mapper.read('epoch_parameters', subfield_index=subfield_min,
                                            epoch_index=0)
        test_schema = epoch_parameters["star_schema"]
        test_schema.append(("subfield", int))
        test_schema.append(("epoch", int))
        test_schema.append(("star_catalog_entry", int))
        test_catalog = numpy.zeros(n_subfields * self.n_epochs, dtype=numpy.dtype(test_schema))

        # Now loop over subfields and epochs
        test_ind = 0
        for subfield_index in xrange(subfield_min, subfield_max+1):
            for epoch_index in xrange(self.n_epochs):
                # Read in star catalog
                epoch_parameters = self.mapper.read('epoch_parameters',
                                                    subfield_index=subfield_index,
                                                    epoch_index=epoch_index)
                star_catalog = self.mapper.read("star_catalog", epoch_parameters)
                # What we do with the catalog depends on whether it's constant or variable PSF.  We
                # have to choose which star in the catalog we want to use differently in these cases.
                if not self.variable_psf:
                    # Transfer first entry (non-offset one) to test_catalog.
                    star_cat_ind = 0
                else:
                    # Find the index of the most aberrated PSF.
                    tot_aber = numpy.zeros(len(star_catalog))
                    for aber in self.psf_builder.use_aber:
                        tot_aber += star_catalog[self.psf_builder.opt_schema_pref+aber]**2
                    star_cat_ind = numpy.argmax(tot_aber)

                # We cannot just do
                #   test_catalog[test_ind] = star_catalog[0]
                # because the schema are not identical.  Instead we loop over the schema entries for
                # star_catalog, and transfer the information for each one over to test_catalog.
                for schema_entry in epoch_parameters["star_schema"]:
                    test_catalog[schema_entry[0]][test_ind] = \
                        star_catalog[schema_entry[0]][star_cat_ind]
                test_catalog["subfield"][test_ind] = subfield_index
                test_catalog["epoch"][test_ind] = epoch_index
                test_catalog["star_catalog_entry"][test_ind] = star_cat_ind
                test_ind += 1
        # Write to file.
        self.mapper.write(test_catalog, "star_test_catalog", epoch_parameters)

    def writeConfig(self, experiment, obs_type, shear_type, subfield_min, subfield_max):

        # Build the dictionary, which we'll output with yaml.dump()
        # We start with the PSF dict, which has much in common with the gal dict.
        # After we write that out, we'll change what has to change and add the gal field.
        d = {}

        # We'll use this same format_items several times below:
        format_items = [
            { 'type' : 'Sequence', 
              'first' : subfield_min, 'last' : subfield_max,
              'repeat' : self.n_epochs 
            },
            { 'type' : 'Sequence', 'nitems' : self.n_epochs } 
        ]

        #
        # (1) First make the config files that will create the PSF images.
        #

        # The input field:
        d['input'] = {
            'catalog' : {
                'file_name' : { 
                    'type' : 'FormattedStr',
                    'format' : 'star_catalog-%03d-%1d.fits',
                    'items' : format_items
                 },
                'dir' : self.mapper.dir
            },

            'dict' : { 
                'file_name' : { 
                    'type' : 'FormattedStr',
                    'format' : 'epoch_parameters-%03d-%1d.yaml',
                    'items' : format_items
                },
                'dir' : self.mapper.dir
            }
        }

        # The output field:
        d['output'] = {
            'type' : 'Fits',

            'file_name' : { 
                'type' : 'FormattedStr',
                'format' : 'starfield_image-%03d-%1d.fits',
                'items' : format_items
            },
            'dir' : self.mapper.dir,

            'nfiles' : self.n_epochs*(subfield_max - subfield_min + 1),
        }

        # The image field:
        if self.variable_psf:
            nx = {'type' : 'Dict', 'key' : 'psf.n_star_linear'}
            ny = {'type' : 'Dict', 'key' : 'psf.n_star_linear'}
        else:
            nx = constants.nx_constpsf
            ny = constants.ny_constpsf
        d['image'] = {
            'type' : 'Tiled',

            'nx_tiles' : nx,
            'ny_tiles' : ny,
            'stamp_xsize' : constants.xsize[self.obs_type][self.multiepoch],
            'stamp_ysize' : constants.ysize[self.obs_type][self.multiepoch],

            'offset' : {
                'type' : 'XY',
                'x' : {'type' : 'Catalog', 'col' : 'xshift' },
                'y' : {'type' : 'Catalog', 'col' : 'yshift' },
            },

            'pixel_scale' : constants.pixel_scale[self.obs_type][self.multiepoch],

            'random_seed' : {
                'type' : 'Sequence',
                'first' : { 'type' : 'Dict', 'key' : 'noise_seed' },
                'nitems' : {
                    'type' : 'Eval',
                    'str' : 'nx*ny',
                    'inx' : { 'type' : 'Current', 'key' : 'image.nx_tiles' },
                    'iny' : { 'type' : 'Current', 'key' : 'image.ny_tiles' }
                }
            },

            'index_convention' : 'python',
        }

        if self.variable_psf:
            # The variable psf images are rather large, so parallelize at the image level.
            d['image']['nproc'] = self.nproc
        else:
            # The constant psf images are small enough that it is probably better to
            # parallelize at the file level.
            d['output']['nproc'] = self.nproc

        # Delegate the basic 'psf' dict to psf_builder
        d['psf'] = self.psf_builder.makeConfigDict()

        if self.variable_psf:
            d['psf']['signal_to_noise'] = { 'type' : 'Catalog', 'col' : 'star_snr' }
            d['image']['noise'] = {
                'type' : 'Gaussian',
                'variance' : { 'type' : 'Dict', 'key' : 'noise.variance' }
            }

        # Set up the file name for the yaml config file:
        experiment_letter = experiment[0]
        obs_letter = obs_type[0]
        shear_letter = shear_type[0]

        file_name = os.path.join(self.mapper.root,
                                 experiment_letter + obs_letter + shear_letter + '_psf.yaml')
        print 'Write PSF config dict to ',file_name

        import yaml
        # This workaround makes sure we avoid anchors and aliases, since they make it a bit 
        # less readable.  (The default likes to alias the format_items variable I used above.)
        # c.f. http://pyyaml.org/ticket/91
        Dumper = yaml.SafeDumper
        Dumper.ignore_aliases = lambda self, data: True

        with open(file_name,'w') as f:
            yaml.dump(d, f, indent=4, Dumper=Dumper)

        # 
        # (2) Make config files for the galaxy images
        #

        # Start with the above dict and make appropriate changes as necessary.
        if self.variable_psf:
            del d['psf']['signal_to_noise']
        else:
            d['image']['noise'] = {
                'type' : 'Gaussian',
                'variance' : { 'type' : 'Dict', 'key' : 'noise.variance' }
            }
        d['input']['catalog']['file_name']['format'] = 'epoch_catalog-%03d-%1d.fits'
        d['output']['file_name']['format'] = 'image-%03d-%1d.fits'
        nx = constants.ncols
        ny = constants.nrows
        d['image']['nx_tiles'] = nx
        d['image']['ny_tiles'] = ny
        d['image']['random_seed']['nitems'] = nx*ny
        d['image']['offset'] = {
            'type' : 'XY',
            'x' : { 'type' : 'Eval',
                    'str' : 'xdither + xshift',
                    'fxdither' : { 'type' : 'Dict', 'key' : 'xdither' },
                    'fxshift' : { 'type' : 'Catalog', 'col' : 'xshift' }
                  },
            'y' : { 'type' : 'Eval',
                    'str' : 'ydither + yshift',
                    'fydither' : { 'type' : 'Dict', 'key' : 'ydither' },
                    'fyshift' : { 'type' : 'Catalog', 'col' : 'yshift' }
                  }
        }
        d['image']['gsparams'] = { 'maximum_fft_size' : 10240 }

        # Delegate the basic 'gal' dict to galaxy_builder
        d['gal'] = self.galaxy_builder.makeConfigDict()

        # Modifications to gal dict:
        d['gal']['shear'] = {
            'type' : 'G1G2',
            'g1' : { 'type' : 'Catalog', 'col' : 'g1' },
            'g2' : { 'type' : 'Catalog', 'col' : 'g2' }
        }
        d['gal']['magnification'] = { 'type' : 'Catalog', 'col' : 'mu' }

        if not self.variable_psf:
            # The galaxy images are large, so parallelize at the image level (unlike for the star
            # fields, for which we already had set up to parallelize at the file level).
            d['image']['nproc'] = self.nproc
            del d['output']['nproc']

        if self.real_galaxy:
            # Need to add the RealGalaxyCatalog to input
            d['input']['real_catalog'] = {
                'dir' : os.path.abspath(self.galaxy_builder.gal_dir),
                'file_name' : self.galaxy_builder.rgc_file,
                'preload' : self.preload
            }


        file_name = os.path.join(self.mapper.root,
                                 experiment_letter + obs_letter + shear_letter + '.yaml')
        print 'Write gal config dict to ',file_name

        with open(file_name,'w') as f:
            yaml.dump(d, f, indent=4, Dumper=Dumper)

        #
        # (3) Finally, make a config file for the "StarTest" images
        #

        # Easiest to just start over here.
        d = {}

        # The input field:
        d['input'] = {
            'catalog' : {
                'file_name' : 'star_test_catalog.fits',
                'dir' : self.mapper.dir
            },
        }

        # The output field:
        d['output'] = {
            'type' : 'DataCube',
            'file_name' : 'star_test_images.fits',
            'dir' : self.mapper.dir,
            'nimages' : self.n_epochs*(subfield_max - subfield_min + 1),
            # The star_test images are all small, so parallelize at the file level.
            'nproc' : self.nproc
        }

        # The image field:
        d['image'] = {
            # Make them 2x the normal size.
            'xsize' : 2 * constants.xsize[self.obs_type][self.multiepoch],
            'ysize' : 2 * constants.ysize[self.obs_type][self.multiepoch],
            'pixel_scale' : constants.pixel_scale[self.obs_type][self.multiepoch],
        }

        # Delegate the basic 'psf' dict to psf_builder
        d['psf'] = self.psf_builder.makeConfigDict(use_zero_index=False)

        # Set up the file name for the yaml config file:
        experiment_letter = experiment[0]
        obs_letter = obs_type[0]
        shear_letter = shear_type[0]

        file_name = os.path.join(self.mapper.root,
                                 experiment_letter + obs_letter + shear_letter + '_star_test.yaml')
        print 'Write star test config dict to ',file_name

        with open(file_name,'w') as f:
            yaml.dump(d, f, indent=4, Dumper=Dumper)


    def writeGalImage(self, subfield_index, epoch_index):
        """This method builds and writes the galaxy images to disk.
        """
        epoch_parameters = self.mapper.read("epoch_parameters", subfield_index=subfield_index,
                                            epoch_index=epoch_index)
        epoch_catalog = self.mapper.read("epoch_catalog", epoch_parameters)
        seed = epoch_parameters["noise_seed"]

        pixel_scale = constants.pixel_scale[self.obs_type][self.multiepoch]
        xsize = constants.xsize[self.obs_type][self.multiepoch]
        ysize = constants.ysize[self.obs_type][self.multiepoch]

        # Setup full image on which to place the postage stamps
        galaxy_image = galsim.ImageF(constants.ncols * xsize, constants.nrows * ysize,
                                     scale=pixel_scale)
        galaxy_image.setOrigin(0,0)

        # maximum sizes for padding RealGalaxy objects with noise
        max_xsize = xsize + 2*(constants.centroid_shift_max + constants.epoch_shift_max)
        max_ysize = ysize + 2*(constants.centroid_shift_max + constants.epoch_shift_max)

        # The GSObjects that are returned by the builders are in arcsec, so our galsim.Pixel needs
        # to use arcsec as well.  However, some of the later manipulations (shifts of the centroid
        # within the image carried out by the draw() function) will be in pixels.
        pixel = galsim.Pixel(pixel_scale)

        # We sometimes need a bit larger FFT than usual.
        params = galsim.GSParams(maximum_fft_size=10240)

        cached_psf_obj = None
        for record in epoch_catalog:
            rng = galsim.UniformDeviate(seed)
            seed = seed + 1

            # Build PSF (or take cached value if possible)
            if not self.variable_psf:
                if cached_psf_obj is None:
                    psf = self.psf_builder.makeGalSimObject(record, epoch_parameters["psf"])
                    cached_psf_obj = psf
                else:
                    psf = cached_psf_obj
            else:
                psf = self.psf_builder.makeGalSimObject(record, epoch_parameters["psf"])

            # Build galaxy
            galaxy = self.galaxy_builder.makeGalSimObject(
                record, epoch_parameters["galaxy"], xsize=max_xsize, ysize=max_ysize, rng=rng)
            galaxy.applyLensing(g1=record['g1'], g2=record['g2'], mu=record['mu'])
            final = galsim.Convolve([psf, pixel, galaxy], gsparams=params)

            # Apply both offsets
            offset = galsim.PositionD(epoch_parameters['xdither'] + record['xshift'],
                                      epoch_parameters['ydither'] + record['yshift'])

            # Draw postage stamp
            bbox = galsim.BoundsI(
                xmin=int(record['xmin']), ymin=int(record['ymin']),
                xmax=int(record['xmax']), ymax=int(record['ymax']),
            )
            stamp = galaxy_image.subImage(bbox)
            # Draw into the postage stamp.
            final.draw(stamp, normalization='f', dx=pixel_scale, offset=offset)

            # Apply whitening if necessary:
            if hasattr(final, 'noise'):
                current_var = final.noise.applyWhiteningTo(stamp)
            else:
                current_var = 0.

            # The lines below are diagnostics that can be used to check that the actual S/N is
            # fairly consistent with the estimated one.
            if False:
                # G08 is the best possible S/N estimate:
                #   S = sum W(x,y) I(x,y) / sum W(x,y)
                #   N^2 = Var(S) = sum W(x,y)^2 Var(I(x,y)) / (sum W(x,y))^2
                # with W(x,y) = I(x,y), so
                #   S = sum I^2(x,y) / sum I(x,y)
                #   N^2 = noise variance * sum I^2(x,y) / (sum I(x,y))^2
                #   S/N = sqrt(sum I^2(x,y)) / sqrt(noise variance)
                actual_sn_g08 = \
                    numpy.sqrt((stamp.array**2).sum() / float(epoch_parameters['noise']['variance']))
                try:
                    res = stamp.FindAdaptiveMom()
                    aperture_noise = numpy.sqrt(float(epoch_parameters['noise']['variance']) * \
                                                    2.*numpy.pi*(res.moments_sigma**2))
                    # The number below is the flux S/N within an elliptical Gaussian filter.  My
                    # guess is that it will be somewhere below the optimal actual_sn_g08 but not too
                    # horrible.
                    sn_ellip_gauss = res.moments_amp / aperture_noise
                    # We also want to estimate the S/N on the size, using an unweighted estimator
                    #   S = Sum I(x,y) [(x-x_c)^2 + (y-y_c)^2]
                    #   N^2 = (noise variance) * Sum [(x-x_c)^2 + (y-y_c)^2]^2
                    # For this, we use the centroid estimate from the adaptive moments.  But we also
                    # have to set up the grid of x, y values for the postage stamp, according to the
                    # same exact convention as used for adaptive moments, which is that the center
                    # of the first pixel is 1.
                    if stamp.array.shape[0] != stamp.array.shape[1]:
                        raise RuntimeError
                    min = 1.
                    max = float(stamp.array.shape[0]+1)
                    x_pix, y_pix = numpy.meshgrid(numpy.arange(min, max, 1.),
                                                  numpy.arange(min, max, 1.))
                    dx_pix = x_pix - res.moments_centroid.x
                    dy_pix = y_pix - res.moments_centroid.y
                    sn_size = numpy.sum(stamp.array * (dx_pix**2 + dy_pix**2)) / \
                        numpy.sqrt(float(epoch_parameters['noise']['variance']) * \
                                       numpy.sum((dx_pix**2 + dy_pix**2)**2))
                except:
                    sn_ellip_gauss = -10.
                    sn_size = -10.
                print 'SN: ', record['gal_sn'], actual_sn_g08, sn_ellip_gauss, sn_size
            self.noise_builder.addNoise(rng, epoch_parameters['noise'], stamp, current_var)

        self.mapper.write(galaxy_image, "image", epoch_parameters)


    def writePSFImage(self, subfield_index, epoch_index):
        """This method builds and writes the star field images to disk.
        """
        epoch_parameters = self.mapper.read("epoch_parameters", subfield_index=subfield_index,
                                            epoch_index=epoch_index)
        star_catalog = self.mapper.read("star_catalog", epoch_parameters)
        seed = epoch_parameters["noise_seed"]

        pixel_scale = constants.pixel_scale[self.obs_type][self.multiepoch]
        xsize = constants.xsize[self.obs_type][self.multiepoch]
        ysize = constants.ysize[self.obs_type][self.multiepoch]
        pixel = galsim.Pixel(pixel_scale)

        # Make PSF image
        if self.variable_psf:
            n_star_linear = epoch_parameters["psf"]["n_star_linear"]
            star_image = galsim.ImageF(n_star_linear * xsize,
                                       n_star_linear * ysize,
                                       scale=pixel_scale)
        else:
            star_image = galsim.ImageF(constants.nx_constpsf * xsize,
                                       constants.ny_constpsf * ysize,
                                       scale=pixel_scale)
        star_image.setOrigin(0, 0)

        cached_psf_obj = None
        for record in star_catalog:
            rng = galsim.UniformDeviate(seed)
            seed = seed + 1

            bbox = galsim.BoundsI(
                xmin=int(record['xmin']), ymin=int(record['ymin']),
                xmax=int(record['xmax']), ymax=int(record['ymax']),
            )
            stamp = star_image.subImage(bbox)

            # Build PSF (or take cached value if possible)
            if not self.variable_psf:
                if cached_psf_obj is None:
                    psf = self.psf_builder.makeGalSimObject(record, epoch_parameters["psf"])
                    cached_psf_obj = psf
                else:
                    psf = cached_psf_obj
            else:
                psf = self.psf_builder.makeGalSimObject(record, epoch_parameters["psf"])

            final = galsim.Convolve([psf, pixel])
            offset = galsim.PositionD(record['xshift'],record['yshift'])
            final.draw(stamp, normalization='f', offset=offset)
            if self.variable_psf:
                self.noise_builder.addStarImageNoise(
                    rng, epoch_parameters['noise'], record['star_snr'], stamp)

        self.mapper.write(star_image, "starfield_image", epoch_parameters)

    def writeStarParameters(self, subfield_index, epoch_index):
        """This method writes out a dict for the PSF shapes needed for metric calculation.
        """
        # Only do this for constant shear, not variable shear!  The star shapes are needed to create
        # the metric for the constant shear branch fits to (m, c) values.
        if self.shear_type == "variable":
            return

        epoch_parameters = self.mapper.read("epoch_parameters", subfield_index=subfield_index,
                                            epoch_index=epoch_index)
        star_catalog = self.mapper.read("star_catalog", epoch_parameters)
        seed = epoch_parameters["noise_seed"]

        star_image = self.mapper.read("starfield_image", epoch_parameters)

        starshape_parameters = None
        if not self.variable_psf:
            for record in star_catalog:
                bbox = galsim.BoundsI(
                    xmin=int(record['xmin']), ymin=int(record['ymin']),
                    xmax=int(record['xmax']), ymax=int(record['ymax']),
                    )
                stamp = star_image.subImage(bbox)
                if starshape_parameters is None:
                    try:
                        shape_res = galsim.hsm.FindAdaptiveMom(stamp)
                        star_g1 = shape_res.observed_shape.g1
                        star_g2 = shape_res.observed_shape.g2
                    except:
                        star_g1 = -10.
                        star_g2 = -10.
                    starshape_parameters = {"psf_g1": star_g1,
                                            "psf_g2": star_g2,
                                            "subfield_index": subfield_index,
                                            "epoch_index": epoch_index}
        # For variable PSF, choose a random subset of the stars to measure.  We can just use the
        # first N in the catalog, since they correspond to completely random positions in the field.
        if self.variable_psf:
            # Barney suggested using 1% of the stars.  Use ceil() to make sure that we don't end
            # up with zero for test runs with few stars.
            n_star_use = int(numpy.ceil(0.01*len(star_catalog)))
            sub_catalog = star_catalog[0:n_star_use]
            for record in sub_catalog:
                bbox = galsim.BoundsI(
                    xmin=int(record['xmin']), ymin=int(record['ymin']),
                    xmax=int(record['xmax']), ymax=int(record['ymax']),
                    )
                stamp = star_image.subImage(bbox)
                try:
                    shape_res = galsim.hsm.FindAdaptiveMom(stamp)
                    star_g1 = shape_res.observed_shape.g1
                    star_g2 = shape_res.observed_shape.g2
                    if starshape_parameters is None:
                        starshape_parameters = {"psf_g1": star_g1,
                                                "psf_g2": star_g2,
                                                "subfield_index": subfield_index,
                                                "epoch_index": epoch_index,
                                                "n_star_actual": 1}
                    else:
                        starshape_parameters["psf_g1"] += star_g1
                        starshape_parameters["psf_g2"] += star_g2
                        starshape_parameters["n_star_actual"] += 1
                except:
                    pass
            # Now normalize average shapes.
            if starshape_parameters != None:
                if starshape_parameters["n_star_actual"] > 1:
                    starshape_parameters["psf_g1"] /= starshape_parameters["n_star_actual"]
                    starshape_parameters["psf_g2"] /= starshape_parameters["n_star_actual"]
            else:
                starshape_parameters = {"psf_g1": -10.,
                                        "psf_g2": -10.,
                                        "subfield_index": subfield_index,
                                        "epoch_index": epoch_index,
                                        "n_star_actual": 0}

        self.mapper.write(starshape_parameters, "starshape_parameters", starshape_parameters)

    def packagePublic(self, subfield_min, subfield_max):
        """This method packages up the public outputs."""
        import shutil
        import tarfile

        # First, do some basic calculations related to the deep fields.  If they are included in the
        # range of requested subfields, i.e., (subfield_min, ..., subfield_max), then things are
        # slightly more complicated: we have to choose which of those subfields to actually package
        # up, since we only want the deep fields to be a certain fraction of the field overall.
        # Also, we change the name on output, because we want to be completely clear that these are
        # not typical images.
        n_deep_subfields = constants.n_deep_subfields
        n_reg_subfields = constants.n_subfields - n_deep_subfields
        n_deep_to_output = int(round(constants.deep_frac * n_reg_subfields))
        max_reg_subfield = n_reg_subfields - 1
        min_deep_subfield = n_reg_subfields
        max_deep_subfield = n_reg_subfields + n_deep_to_output - 1
        if subfield_max > max_deep_subfield:
            import warnings
            deep_warning = "Requested range of subfields includes extra deep fields. Adjusting."
            warnings.warn(deep_warning)
            subfield_max = max_deep_subfield

        # Define the output directory (and create if necessary).  Use a mapper for this.
        public_mapper = great3sims.mapper.Mapper(self.public_dir, self.experiment, self.obs_type,
                                                 self.shear_type)
        # Also define an absolute path to the root directory structure, because we're going to be
        # moving around, so self.mapper (which uses relative path) won't work.
        root_rel_mapper = great3sims.mapper.Mapper(os.path.abspath(self.mapper.root),
                                                   self.experiment, self.obs_type, self.shear_type)

        # Zipping / tarring.  Open tarfile at the start, then add the files as they are created.
        tarfile_name = os.path.join(self.public_dir,
                                    self.experiment+'-'+self.obs_type+'-'+self.shear_type+'.tar.gz')
        tar = tarfile.open(tarfile_name, "w:gz")

        # Now we want to move into public_dir.  The reason for this is that we don't want the files
        # in public_dir that go into the tarfile to have public_dir/ at the start of their paths,
        # otherwise when untarred they end up in public_dir/public_dir/... which is kind of silly.
        saved_path = os.getcwd()
        os.chdir(self.public_dir)
        sub_mapper = great3sims.mapper.Mapper('.', self.experiment, self.obs_type, self.shear_type)
        for subfield_index in xrange(subfield_min, subfield_max+1):
            # Loop over galaxy and star field images for the defined set of subfields and epochs,
            # and copy over without modification of any sort.
            tmp_dict = {"subfield_index" : subfield_index}
            if subfield_index > max_reg_subfield:
                    tmp_dict["deep_subfield_index"] = subfield_index - n_reg_subfields

            for epoch_index in xrange(self.n_epochs):
                tmp_dict["epoch_index"] = epoch_index
                if subfield_index <= max_reg_subfield:
                    outfile = root_rel_mapper.copyTo(sub_mapper, 'image', tmp_dict)
                    tar.add(outfile)
                    outfile = root_rel_mapper.copyTo(sub_mapper, 'starfield_image', tmp_dict)
                    tar.add(outfile)
                else:
                    outfile = root_rel_mapper.copyTo(sub_mapper, 'image', tmp_dict,
                        new_template = "deep_image-%(deep_subfield_index)03d-%(epoch_index)1d.fits")
                    tar.add(outfile)
                    outfile = root_rel_mapper.copyTo(sub_mapper, 'starfield_image', tmp_dict,
                        new_template = \
                             "deep_starfield_image-%(deep_subfield_index)03d-%(epoch_index)1d.fits")
                    tar.add(outfile)

            # Loop over galaxy catalogs for each subfield, and copy only the information we want to
            # be public.  For now, let's stick with 'x', 'y', 'ID'.  We could consider giving the
            # information on position within the field, but for now, that's encoded in ID rather
            # than being given explicitly in the catalog. Also, the subfield offset will be output
            # separately, so that can be used with x / y to get positions in degrees, assuming the
            # user understands how to go from x / y to degrees in a single subfield.
            # However... if this is a variable_psf branch, then we want to explicitly output
            # information about the location within the field, the tile, and location within the
            # tile.  
            # Unfortunately this is in the epoch_catalog files rather than subfield_catalog
            # ("unfortunately" because it's the same for all epochs and means we have to read in
            # multiple files) so in that case we have to read in that info as well, and merge the
            # bits of info together.  Bah.
            gal_use_cols = [('x', int), ('y', int), ('ID', int)]
            if self.variable_psf:
                gal_epoch_use_cols = [('x_tile_index', int), ('y_tile_index', int),
                                      ('tile_x_pos_deg', float), ('tile_y_pos_deg', float),
                                      ('x_field_true_deg', float), ('y_field_true_deg', float)]
            if subfield_index <= max_reg_subfield:
                tmp_dict = {"subfield_index" : subfield_index}
                if self.variable_psf:
                    tmp_dict["epoch_index"] = 0
                    outfile = root_rel_mapper.mergeSub(sub_mapper, 'subfield_catalog',
                                                      'epoch_catalog', tmp_dict, gal_use_cols,
                                                      gal_epoch_use_cols,
                                                      new_template =
                                                      "galaxy_catalog-%(subfield_index)03d")
                else:
                    outfile = root_rel_mapper.copySub(sub_mapper, 'subfield_catalog', tmp_dict,
                                                      gal_use_cols,
                                                      new_template =
                                                      "galaxy_catalog-%(subfield_index)03d")
            else:
                tmp_dict["deep_subfield_index"] = subfield_index - n_reg_subfields
                if self.variable_psf:
                    tmp_dict["epoch_index"] = 0
                    outfile = root_rel_mapper.mergeSub(sub_mapper, 'subfield_catalog',
                                                      'epoch_catalog', tmp_dict, gal_use_cols,
                                                      gal_epoch_use_cols,
                                                      new_template =
                                                      "deep_galaxy_catalog-%(subfield_index)03d")
                outfile = root_rel_mapper.copySub(sub_mapper, 'subfield_catalog', tmp_dict,
                                                  gal_use_cols,
                                                  new_template =
                                                  "deep_galaxy_catalog-%(deep_subfield_index)03d")
            tar.add(outfile)
            # ... and also copy to text file that gets added
            outfile_no_ext = os.path.splitext(outfile)[0]
            great3sims.mapper.fitsToTextCatalog(outfile_no_ext)
            tar.add(outfile_no_ext + '.txt')

            # Loop over star catalogs, and copy only the information that we want to be public.  For
            # constant PSF branches, this is just 'x' and 'y'.  For variable PSF branches, we want
            # the same additional info as was given for the galaxies re: position within the field
            # and tiles.
            if self.variable_psf:
                star_use_cols = [('x', int), ('y', int), ('x_tile_index', int), 
                                       ('y_tile_index', int), ('tile_x_pos_deg', float),
                                       ('tile_y_pos_deg', float), ('x_field_true_deg', float),
                                       ('y_field_true_deg', float)]
            else:
                star_use_cols = [('x', int), ('y', int)]
            if subfield_index <= max_reg_subfield:
                tmp_dict = {"subfield_index" : subfield_index, "epoch_index" : 0}
                outfile = root_rel_mapper.copySub(
                    sub_mapper, 'star_catalog', tmp_dict, star_use_cols,
                    new_template="star_catalog-%(subfield_index)03d")
            else:
                tmp_dict["deep_subfield_index"] = subfield_index - n_reg_subfields
                tmp_dict["epoch_index"] = 0
                outfile = root_rel_mapper.copySub(
                    sub_mapper, 'star_catalog', tmp_dict, star_use_cols,
                    new_template="deep_star_catalog-%(deep_subfield_index)03d")
            # ... and also copy to text file that gets added 
            outfile_no_ext = os.path.splitext(outfile)[0]
            great3sims.mapper.fitsToTextCatalog(outfile_no_ext)
            tar.add(outfile)
            tar.add(outfile_no_ext+'.txt')

            # We can also give some overall information about subfield offsets.  For now, we use the
            # subfield_parameters file, extract subfield_offset [which currently is in units of
            # separation between galaxies in the grid], convert that to degrees, and output that as
            # a yaml file and a text file.  This is redundant, but since the files are tiny it
            # doesn't seem like a big issue to support both formats.
            template, reader, writer = root_rel_mapper.mappings['subfield_parameters']
            in_path = os.path.join(root_rel_mapper.full_dir, template % tmp_dict)
            subfield_params = great3sims.mapper.readDict(in_path)
            mult_val = constants.image_size_deg / constants.nrows
            offset_parameters = {
                "offset_deg_x": mult_val * subfield_params['subfield_offset'][0],
                "offset_deg_y": mult_val * subfield_params['subfield_offset'][1]
                }
            if subfield_index <= max_reg_subfield:
                template = "subfield_offset-%(subfield_index)03d"
            else:
                template = "deep_subfield_offset-%(deep_subfield_index)03d"
            outfile = os.path.join(sub_mapper.full_dir, template % tmp_dict)
            great3sims.mapper.writeDict(offset_parameters, outfile)
            great3sims.mapper.writeDict(offset_parameters, outfile, type='txt')
            tar.add(outfile + '.yaml')
            tar.add(outfile + '.txt')

            # Finally, for each epoch, we need to give information about the size of the dithering
            # with respect to the first epoch.  For the sake of having a consistent data format for
            # all branches, we include this even for single epoch branches (for which the dithers
            # are all 0 since each image IS the first and only epoch).  We extract the xdither and
            # ydither from the epoch_parameters, and write a file for each subfield and epoch.
            # These files are
            # tiny, so let's write one for each subfield and epoch, in both yaml and txt format.
            for epoch_index in xrange(self.n_epochs):
                tmp_dict["epoch_index"] = epoch_index
                template, reader, writer = root_rel_mapper.mappings['epoch_parameters']
                in_path = os.path.join(root_rel_mapper.full_dir, template % tmp_dict)
                epoch_params = great3sims.mapper.readDict(in_path)
                dither_parameters = {
                    "xdither_pixels": epoch_params['xdither'],
                    "ydither_pixels": epoch_params['ydither']
                    }
                if subfield_index <= max_reg_subfield:
                    template = "epoch_dither-%(subfield_index)03d-%(epoch_index)1d"
                else:
                    template = "deep_epoch_dither-%(deep_subfield_index)03d-%(epoch_index)1d"
                outfile = os.path.join(sub_mapper.full_dir, template % tmp_dict)
                great3sims.mapper.writeDict(dither_parameters, outfile)
                great3sims.mapper.writeDict(dither_parameters, outfile, type='txt')
                tar.add(outfile + '.yaml')
                tar.add(outfile + '.txt')

        # Close the tarfile.  Now that we're done with it, we can go back to our original working
        # directory.
        tar.close()
        os.chdir(saved_path)

        # Delete the files / dirs, just keep the tarfiles.
        shutil.rmtree(public_mapper.full_dir)

    def packageTruth(self, subfield_min, subfield_max):
        """This method packages up the true shear values and PSF ellipticities for metric
        calculations.
        """
        import shutil
        import tarfile

        # First, do some basic calculations related to the deep fields.  If they are included in the
        # range of requested subfields, i.e., (subfield_min, ..., subfield_max), then emit a warning
        # that we are excluding them because they are not used for metric evaluation.
        n_deep_subfields = constants.n_deep_subfields
        n_reg_subfields = constants.n_subfields - n_deep_subfields
        if subfield_max > n_reg_subfields - 1:
            import warnings
            deep_warning = "Requested range of subfields includes deep fields. Adjusting."
            warnings.warn(deep_warning)
            subfield_max = n_reg_subfields - 1

        # Define the output directory (and create if necessary).  Use a mapper for this.
        truth_mapper = great3sims.mapper.Mapper(self.truth_dir, self.experiment, self.obs_type,
                                                self.shear_type)
        # Also define an absolute path to the root directory structure, because we're going to be
        # moving around, so self.mapper (which uses relative path) won't work.
        root_rel_mapper = great3sims.mapper.Mapper(os.path.abspath(self.mapper.root),
                                                   self.experiment, self.obs_type, self.shear_type)

        # Zipping / tarring.  Open tarfile at the start, then add the files as they are created.
        tarfile_name = os.path.join(self.truth_dir,
                                    self.experiment+'-'+self.obs_type+'-'+self.shear_type+'.tar.gz')
        tar = tarfile.open(tarfile_name, "w:gz")

        # Now we want to move into truth_dir.  The reason for this is that we don't want the files
        # in truth_dir that go into the tarfile to have truth_dir/ at the start of their paths,
        # otherwise when untarred they end up in truth_dir/truth_dir/... which is kind of silly.
        saved_path = os.getcwd()
        os.chdir(self.truth_dir)
        sub_mapper = great3sims.mapper.Mapper('.', self.experiment, self.obs_type, self.shear_type)

        # First, we copy over the star test catalog and images.
        # Make the old, new target filenames for the star test catalog:
        template, reader, writer = root_rel_mapper.mappings['star_test_catalog']
        infile = os.path.join(root_rel_mapper.full_dir, template % {}) + '.fits'
        outfile = os.path.join(sub_mapper.full_dir, template % {}) + '.fits'
        shutil.copy2(infile, outfile)
        tar.add(outfile)
        template, reader, writer = root_rel_mapper.mappings['star_test_images']
        infile = os.path.join(root_rel_mapper.full_dir, template % {}) + '.fits'
        outfile = os.path.join(sub_mapper.full_dir, template % {}) + '.fits'
        shutil.copy2(infile, outfile)
        tar.add(outfile)

        # Now do all the per-subfield stuff.
        for subfield_index in xrange(subfield_min, subfield_max+1):
            tmp_dict = {"subfield_index" : subfield_index}

            # If variable shear, then loop over subfield catalogs and copy over just the ID and the
            # per-galaxy reduced shear.
            if self.shear_type == 'variable':
                use_cols = [('ID', int), ('g1', float), ('g2', float),
                            ('g1_intrinsic', float), ('g2_intrinsic', float)]
                outfile = root_rel_mapper.copySub(sub_mapper, 'subfield_catalog', tmp_dict,
                                                  use_cols,
                                                  new_template =
                                                  "galaxy_catalog-%(subfield_index)03d")
                tar.add(outfile)

                # We can also give some overall information about subfield offsets.  This is
                # necessary for variable shear sims in order to convert the positions in each
                # subfield to a position in the field.  For now, we use the subfield_parameters
                # file, extract subfield_offset [which currently is in units of separation between
                # galaxies in the grid], convert that to degrees, and output that as a yaml file and
                # a text file.  This is redundant, but since the files are tiny it doesn't seem like
                # a big issue to support both formats.
                template, reader, writer = root_rel_mapper.mappings['subfield_parameters']
                in_path = os.path.join(root_rel_mapper.full_dir, template % tmp_dict)
                subfield_params = great3sims.mapper.readDict(in_path)
                mult_val = constants.image_size_deg / constants.nrows
                offset_parameters = {
                    "offset_deg_x": mult_val * subfield_params['subfield_offset'][0],
                    "offset_deg_y": mult_val * subfield_params['subfield_offset'][1]
                    }
                template = "subfield_offset-%(subfield_index)03d"
                outfile = os.path.join(sub_mapper.full_dir, template % tmp_dict)
                great3sims.mapper.writeDict(offset_parameters, outfile)
                great3sims.mapper.writeDict(offset_parameters, outfile, type='txt')
                tar.add(outfile + '.yaml')
                tar.add(outfile + '.txt')

            else:
                # If constant shear, then take the subfield g1, g2 and write as yaml/text.
                template, reader, writer = root_rel_mapper.mappings['subfield_parameters']
                in_path = os.path.join(root_rel_mapper.full_dir, template % tmp_dict)
                subfield_params = great3sims.mapper.readDict(in_path)
                shear_params = {
                    "g1": subfield_params['shear']['g1'],
                    "g2": subfield_params['shear']['g2']
                    }
                template = "shear_params-%(subfield_index)03d"
                outfile = os.path.join(sub_mapper.full_dir, template % tmp_dict)
                great3sims.mapper.writeDict(shear_params, outfile)
                great3sims.mapper.writeDict(shear_params, outfile, type='txt')
                tar.add(outfile + '.yaml')
                tar.add(outfile + '.txt')

                # If this branch has constant shear, then we need the PSF (star) shape parameters.
                if self.shear_type == "constant":
                    # And loop over epochs, copying over the PSF shape parameter files.  Can't use
                    # copyTo() because we want to write as .txt in addition to .yaml
                    for epoch_index in xrange(self.n_epochs):
                        tmp_dict["epoch_index"] = epoch_index
                        template, reader, writer = root_rel_mapper.mappings['starshape_parameters']
                        in_path = os.path.join(root_rel_mapper.full_dir, template % tmp_dict)
                        starshape_params = great3sims.mapper.readDict(in_path)
                        outfile = os.path.join(sub_mapper.full_dir, template % tmp_dict)
                        great3sims.mapper.writeDict(starshape_params, outfile)
                        great3sims.mapper.writeDict(starshape_params, outfile, type='txt')
                        tar.add(outfile + '.yaml')
                        tar.add(outfile + '.txt')

        # Close the tarfile.  Now that we're done with it, we can go back to our original working
        # directory.
        tar.close()
        os.chdir(saved_path)

        # Delete the files / dirs, just keep the tarfiles.
        shutil.rmtree(truth_mapper.full_dir)

    def generateSubfieldOffsets(self, rng, n_subfields_per_field, subfield_grid_subsampling):
        """A utility to decide about the offsets between subfields in a field.

        Currently, the offsets are required to be regular, so that the shear and PSF builders can
        simply make an overly dense grid compared to what's needed for a subfield, and use a subset
        of its grid points.  We assume that the options for offsetting are on an
        subfield_grid_subsampling x subfield_grid_subsampling grid.  This then gives
        subfield_grid_subsampling^2 possible locations.  We then choose a random
        n_subfields_per_field-1 of those options for the subfields that are not the first in the
        field.

        Offsets between subfields in a field are first specified as integers and are defined as the
        offset with respect to the first subfield in the field.  Then, we divide by the amount of
        grid subsampling for subfields in a field.  For example, if subfield_grid_subsampling = 7,
        then the minimum and maximum values for x/y offsets are 0/7 and 6/7.  Because of our
        definition with respect to the first subfield, that subfield must get offset (0, 0).

        Results for offsets for each subfield are returned as a list of tuples, where the length of
        the list is determined by n_subfields_per_field.
        """
        # Check for the simplest case, and do it.
        if n_subfields_per_field == 1:
            offsets = []
            offsets.append((0., 0.))
            return offsets
        else:
            # first do everything in terms of ints, for easier comparison
            int_offsets = []
            int_offsets.append((0, 0))
            for i_subfield in range(n_subfields_per_field-1):
                test_tuple = (0, 0)
                # Make sure we end up with a unique one that does not exist in int_offsets
                while test_tuple in int_offsets:
                    test_tuple = (numpy.int(numpy.floor(rng()*subfield_grid_subsampling)),
                                  numpy.int(numpy.floor(rng()*subfield_grid_subsampling)))
                int_offsets.append(test_tuple)
            offsets = [(float(int_offset[0])/subfield_grid_subsampling,
                        float(int_offset[1])/subfield_grid_subsampling)
                       for int_offset in int_offsets]
            return offsets
