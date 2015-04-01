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
"""
This file contains the SimBuilder class that does the heavy lifting; it coordinates all the steps
that are needed to generate parameters, catalogs, config files, images, and image packages for a
given GREAT3 branch (experiment, data type, shear type).
"""
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
        their default (control experiment) values depending on what keyword arguments are set.
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

    def __init__(self, root, obs_type, shear_type, gal_dir, ps_dir, opt_psf_dir, atmos_ps_dir,
                 public_dir, truth_dir, preload=False, nproc=-1, gal_pairs=True):
        """Initialize a builder for the given `obs_type` and `shear_type`.

        @param[in] root         Root directory for generated files.
        @param[in] obs_type     Type of observation to simulate: either "ground" or "space".
        @param[in] shear_type   Type of shear field to generate: either "constant" or "variable".
        @param[in] gal_dir      Directory with real galaxy catalog information.
        @param[in] ps_dir       Directory with tabulated iCosmo shear power spectra.
        @param[in] opt_psf_dir  Directory with the optical PSF models for ground and space
                                variable PSF simulations.
        @param[in] atmos_ps_dir Directory with tabulated atmospheric PSF anisotropy power spectra.
        @param[in] public_dir   Directory for placing files to be distributed publicly.
        @param[in] truth_dir    Directory for placing files to be used for metric evaluation.
        @param[in] preload      Preload the RealGalaxyCatalog images to speed up generation of large
                                numbers of real galaxies?  Note that for parametric galaxy branches,
                                the catalog is never preloaded. [default = False]
        @param[in] nproc        How many processes to use in the config file.  [default = -1]
        @param[in] gal_pairs    For constant shear branches, should it use 90 degree rotated pairs
                                to cancel out shape noise, or not?  This option is ignored for
                                variable shear branches. [default: True]
        """
        self.obs_type = obs_type
        self.shear_type = shear_type
        self.public_dir = public_dir
        self.truth_dir = truth_dir
        self.preload = preload
        # Below we initialize the builders for the PSF, shear, galaxy population, and noise field.
        # They each require various bits of information as appropriate (e.g., only the PSF builder
        # needs to know where information about atmospheric PSFs lives).
        self.psf_builder = great3sims.psf.makeBuilder(obs_type=obs_type,
                                                      variable_psf=self.variable_psf,
                                                      multiepoch=self.multiepoch,
                                                      shear_type=self.shear_type,
                                                      opt_psf_dir=opt_psf_dir,
                                                      atmos_ps_dir=atmos_ps_dir)
        self.shear_builder = great3sims.shear.makeBuilder(shear_type=shear_type, obs_type=obs_type,
                                                          multiepoch=self.multiepoch, ps_dir=ps_dir)
        self.galaxy_builder = great3sims.galaxies.makeBuilder(real_galaxy=self.real_galaxy,
                                                              obs_type=obs_type,
                                                              shear_type=shear_type,
                                                              multiepoch=self.multiepoch,
                                                              gal_dir=gal_dir,
                                                              preload=preload,
                                                              gal_pairs=gal_pairs)
        self.noise_builder = great3sims.noise.makeBuilder(obs_type=obs_type,
                                                          multiepoch=self.multiepoch,
                                                          variable_psf = self.variable_psf)
        # We also initialize a mapper, which assists with i/o for this branch.  It knows how to make
        # directory and file names depending on the branch, and what types of files need to be
        # output for that branch.
        self.mapper = great3sims.mapper.Mapper(root, self.experiment, obs_type, shear_type)
        # And store some additional necessary information.
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

        # Also check some things about the deep fields.  We need to figure out how many to make.
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
        # selection can use them - and use a multiplicative factor to tell it to add less noise for
        # the deep field.

        # Put together the basic parameters to be stored in the metaparameters file.
        self.parameters = {"metadata": metadata, "seed": seed, "pixel_scale": pixel_scale,
                           "n_fields": n_fields}
        # And use the mapper to write them out.
        self.mapper.write(self.parameters, 'parameters', self.parameters)
        # Now, initialize a RNG with the required seed.  We'll use that to set up seeds for
        # everything else that follows.
        rng = galsim.UniformDeviate(seed)
        # We also have to set up schema for catalogs and such.  We'll set up a basic schema here, as
        # well as schema to be used for shear-related parameters.
        base_schema = [("index", int), ("x", int), ("y", int),
                       ("xshift", float), ("yshift", float),
                       ("xmin", int), ("xmax", int), ("ymin", int), ("ymax", int)]
        shear_schema = [("g1", float), ("g2", float), ("mu", float),
                       ("x_field_pos", int), ("y_field_pos", int), ("ID", int)]
        seed += 1  # We could also draw random integers to set seeds, but based on a discussion, it
                   # does not seem necessary

        # Start a separate sequence for noise_seed, which will be used to make the noise fields in
        # the images.  The second number is just a random number that is much larger than the
        # plausible number of items in the original seed sequence.  (It is phi * 10^6.)
        noise_seed = seed + 1618033

        # Now we begin to loop over successively smaller units - starting with field, then subfield,
        # then epoch.  For each of these, we will generate the necessary parameters.
        for field_index in xrange(n_fields):
            # A given field has the same shear and per-epoch PSF.  Thus, we set up the shear and PSF
            # parameters at the field level.
            #
            # The builders can decide what format to use for the results of
            # generateFieldParameters().  We don't actually care what that format is out here - we
            # will just pass it along to generateSubfieldParameters() or generateEpochParameters(),
            # which are other methods of the builders.  So we require internal consistency between
            # the various generate*Parameters() methods, but we should be able to switch parameter
            # selection between the field/subfield/epoch layer without modifying this code, just
            # modifying the builders themselves (in psf.py, shear.py, noise.py, or galaxies.py).
            field_shear_parameters = self.shear_builder.generateFieldParameters(rng, field_index)
            field_psf_parameters = self.psf_builder.generateFieldParameters(rng, field_index)
            # Now we set up a `field_parameters` dict which includes the field parameters we just
            # generated for the shear and PSF, as well as other basic info like "which field is
            # this", the offsets of subfields within the field (see call to
            # generateSubfieldOffsets() below), the random seed, and some metadata.
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
            # Use the mapper to write the field parameters to file.
            self.mapper.write(field_parameters, 'field_parameters', field_parameters)

            # Now loop over subfields within this field.
            field_min_subfield = field_index * n_subfields_per_field
            field_max_subfield = field_min_subfield + n_subfields_per_field - 1
            for subfield_index in xrange(field_min_subfield, field_max_subfield+1):
                # A given subfield has the same shear (already determined at field level) and
                # galaxies.  But to allow for flexibility later on, we'll have a
                # generateSubfieldParameters() for the shear builder anyway; it has to take the
                # output of generateFieldParameters() as an input.  For now, it's a no-op, returning
                # the input as output - but that might not be the case later on.
                #
                # As at the field level, we then collect the subfield parameters (which are the
                # field-level shear parameters and the galaxy parameters determined at the subfield
                # level), plus some other necessary bits of data.
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
                # Include schema for the subfield-level parameters.
                subfield_parameters["subfield_schema"] = \
                    (base_schema + shear_schema + subfield_parameters["galaxy"]["schema"])
                # Use the mapper to write the subfield parameters to file.
                self.mapper.write(subfield_parameters, 'subfield_parameters', subfield_parameters)
                # Finally, loop over the epoch.
                for epoch_index in xrange(self.n_epochs):
                    # Each epoch has its own PSF (already determined at field level) and noise.
                    # But the galaxies and shears were determined at the subfield level, so nothing
                    # new is needed here.
                    epoch_parameters = dict(subfield_parameters)
                    epoch_parameters["psf"] = \
                        self.psf_builder.generateEpochParameters(rng, subfield_index, epoch_index,
                                                                 field_psf_parameters)
                    # We also determine noise-related parameters at the epoch level.  So, here we
                    # determine a multiplying factor for the sky variance based on whether we are in
                    # a deep field or not.  Note that changes in the sky variance due to single
                    # vs. multiepoch images are handled directly in the noise_builder (which knows
                    # what branch we're in), so they do not need to be included here.
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
                    # `xdither`, `ydither` are the amount of dithering between epochs of this
                    # subfield.  In contrast, `epoch_offset` (a tuple) is the amount of offsetting
                    # between this subfield and the first one in the field, i.e., it's the same as
                    # `subfield_offset`.  This is not truly a per-epoch parameter, however, it is
                    # included here so that the per-epoch catalog maker (specifically, for PSFs)
                    # will be able to do its job.
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
        """Given the subfield index, load the corresponding metaparameters (previously generated by
        writeParameters()), and use them to generate a catalog of galaxies and shear values.
        """
        # Read the subfield and field parameters.
        subfield_parameters = self.mapper.read("subfield_parameters", subfield_index=subfield_index)
        field_parameters = self.mapper.read(
            "field_parameters",
            field_index = ( subfield_index / 
                            constants.n_subfields_per_field[self.shear_type][self.variable_psf] ) 
        )
        # Set up a catalog with the appropriate schema.
        catalog = numpy.zeros(constants.nrows * constants.ncols,
                              dtype=numpy.dtype(subfield_parameters["subfield_schema"]))
        index = 0
        # Use the given seed to initialize a RNG.
        rng = galsim.UniformDeviate(subfield_parameters["seed"])
        # Loop over the galaxies and generate some basic numbers like where they belong in the image
        # corresponding to this subfield, what is their centroid shift compared to the nominal
        # position, etc.
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
        # or not.  Note that sky variance changes due to single vs. multiepoch images are not
        # included at this stage, because we want to impose our cuts based on whether the S/N would
        # be 20 in a single combined image, not based on its value in the individual epoch images.
        # Also note that while the noise variance parameter output into the epoch_parameters files
        # by the noise builder already includes this noise_mult for the deep fields (and any factors
        # due to single vs. multiepoch which we do *not* want here), we have to recalculate the deep
        # field noise multiplying factor because we're just going to use the
        # noise_builder.typical_variance which does not include any of those factors.
        if subfield_index < constants.n_subfields - constants.n_deep_subfields:
            noise_mult = 1.
        else:
            noise_mult = constants.deep_variance_mult

        # We give the galaxy catalog generation routine a value for seeing to use when selecting
        # galaxies.  This calculation becomes more complex in the multi-epoch case since we have to
        # decide on a relevant effective FWHM.
        if self.obs_type == "space":
            effective_seeing = None
        else:
            # Note: really we care about the full PSF FWHM, not just the atmospheric part.  However,
            # we use the seeing as a proxy for it, so we don't have to start generating images.  If
            # this seems really worrisome, we could make some simple sims, derive approximate rules
            # for total PSF size including optics as well (which will mainly affect really
            # good-seeing images), and use those instead of just the atmospheric seeing.
            if not self.multiepoch and not self.variable_psf:
                # For single epoch images with a constant PSF, the FWHM is just a single (scalar)
                # value for the entire subfield.
                effective_seeing = field_parameters["psf"]["atmos_psf_fwhm"]
            else:
                # This is a 1d numpy array of FWHM values.  We determine a single effective seeing
                # value from it.
                seeing = field_parameters["psf"]["atmos_psf_fwhm"]
                effective_seeing = 1. / numpy.mean(1./seeing)
        # The galaxy builder generates a catalog given some noise variance information (which
        # determines galaxy S/N), and the effective seeing (for selecting galaxies that are
        # resolved).
        self.galaxy_builder.generateCatalog(rng, catalog, subfield_parameters,
                                            self.noise_builder.typical_variance, noise_mult,
                                            effective_seeing)
        # The shear builder generates shear values for the catalog.
        self.shear_builder.generateCatalog(rng, catalog, subfield_parameters["shear"],
                                           subfield_parameters["subfield_offset"], subfield_index)
        # The mapper writes out the galaxy catalog for this subfield in the appropriate directory
        # and file.
        self.mapper.write(catalog, "subfield_catalog", subfield_parameters)

    def writeEpochCatalog(self, subfield_index, epoch_index):
        """Given the subfield and epoch indices, load the epoch metaparameters and a
        previously-generated subfield catalog.  Then, add per-object PSF parameter information (and
        possibly shift the centroids) to create and save an epoch catalog.
        """
        # First, read in the stored parameter information and subfield catalog.
        epoch_parameters = self.mapper.read('epoch_parameters', subfield_index=subfield_index,
                                            epoch_index=epoch_index)
        subfield_catalog = self.mapper.read('subfield_catalog', epoch_parameters)
        # Make a catalog for this epoch, according to the stored schema.
        epoch_catalog = numpy.zeros(constants.nrows * constants.ncols,
                                    dtype=numpy.dtype(epoch_parameters["epoch_schema"]))
        # Initialize a RNG with the given seed.
        rng = galsim.UniformDeviate(epoch_parameters["seed"])
        # First, just carry over values from the subfield catalog into the epoch catalog.
        for name, _ in epoch_parameters["subfield_schema"]:
            epoch_catalog[name] = subfield_catalog[name]
        # Then, generate the PSF information for the catalog, given the parameters for this epoch.
        self.psf_builder.generateCatalog(rng, epoch_catalog, epoch_parameters,
                                         epoch_parameters["epoch_offset"], normalized=True)
        # Write out the epoch-level catalog for the galaxies.
        self.mapper.write(epoch_catalog, "epoch_catalog", epoch_parameters)

        # We also need a star catalog for this epoch.  To set up the positions in the catalog,
        # figure out the size that each star postage stamp will cover.  (These should be the same as
        # the sizes covered by galaxy postage stamps.)
        xsize = constants.xsize[self.obs_type][self.multiepoch]
        ysize = constants.ysize[self.obs_type][self.multiepoch]
        # Then write star catalog entries, which depends on whether this is a variable PSF branch or
        # not.
        if self.variable_psf:
            # We have to write catalog entries indicating the "true" star position within the
            # **field** (not subfield).  We also have to write x, y entries that are used to place
            # the objects on an image grid, the parameters of which are to be defined here.
            #
            # Finally, let's define a galsim.DistDeviate to use to draw random values of S/N for
            # each star.  The distribution of S/N values comes from a catalog of main sequence stars
            # for galactic longitude=180 degrees, in the I band, from the LSST ImSim.  Reading
            # numbers off a plot from Chihway Chang and assuming that the S/N of an I=25 star is 25,
            # the function we get is
            # dN/d(S/N) ~ (mag - 18) / S/N
            #           ~ (-2.5 log10(S/N) + 10.5) / S/N.
            # We assume that the distribution of star S/N values can go from 25 to 400, with higher
            # S/N stars being excluded due to saturation.
            dist_deviate = galsim.DistDeviate(
                rng,
                function = lambda x : (-2.5*numpy.log10(x)+10.5)/x, x_min=25., x_max=400.)

            # First, determine how many stars should be in this subfield, which will determine the
            # size of the catalog to write.  Note that this is for a single subfield, which means
            # the number of stars should be 1/20 of those for the entire field.
            n_star_linear = epoch_parameters["psf"]["n_star_linear"]
            star_catalog = numpy.zeros(n_star_linear * n_star_linear,
                                       dtype=numpy.dtype(epoch_parameters["star_schema"]))
            index = 0
            # Loop over the entries in the star catalog.
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
                    # Here are some numbers that we'll only compute if it's the first epoch.  Since
                    # we want to represent the same star population in each epoch, we store them in
                    # the cache after generating the first epoch, and read from the cache if it's
                    # not the first epoch.
                    if epoch_index == 0:
                        record['xshift'] = sx
                        record['yshift'] = sy
                        # But these numbers are the true positions within the field.
                        record['x_field_true_deg'] = constants.image_size_deg * rng()
                        record['y_field_true_deg'] = constants.image_size_deg * rng()
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
                        # they could be used.) We are glossing over the slight difference in star
                        # number density that should also occur, and just using the fact that the
                        # S/N distribution should be similar.
                        record['star_snr'] = dist_deviate()
                    index += 1
            if epoch_index == 0:
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
            # For constant PSF branches, the star catalog generation is much simpler.
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
                    # Here are some numbers that we will save in a cache, only computing for the
                    # first epoch in a field, so as to preserve the random subpixel shifts between
                    # stars in the starfield.  This way simple coaddition will work for all stars in
                    # one of the constant PSF starfields.
                    if epoch_index == 0:
                        if index > 0:
                            sx = (2.0*rng() - 1.0) * constants.centroid_shift_max
                            sy = (2.0*rng() - 1.0) * constants.centroid_shift_max
                            record["xshift"] = sx
                            record["yshift"] = sy
                        else:
                            record["xshift"] = 0
                            record["yshift"] = 0
                    index += 1
            if epoch_index == 0:
                self.cached_xshift = star_catalog['xshift']
                self.cached_yshift = star_catalog['yshift']
            else:
                star_catalog['xshift'] = self.cached_xshift
                star_catalog['yshift'] = self.cached_yshift
        # Given the basic catalog information generated above, make the catalog of PSF parameters
        # (e.g., optical PSF and atmospheric PSF parameters) for each star.
        self.psf_builder.generateCatalog(rng, star_catalog, epoch_parameters,
                                         epoch_parameters["epoch_offset"], normalized=False)
        # Write the star catalog to file in the appropriate location and format.
        self.mapper.write(star_catalog, "star_catalog", epoch_parameters)

    def writeStarTestCatalog(self, subfield_min, subfield_max):
        """Given a range of subfield and epoch indices, write a test catalog for generating star
        images to check for oddities.  We want a small set of objects (suitable for eyeballing)
        chosen in a representative way."""
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

        # Now loop over subfields and epochs.
        test_ind = 0
        for subfield_index in xrange(subfield_min, subfield_max+1):
            for epoch_index in xrange(self.n_epochs):
                # Read in the epoch parameters and star catalog.
                epoch_parameters = self.mapper.read('epoch_parameters',
                                                    subfield_index=subfield_index,
                                                    epoch_index=epoch_index)
                star_catalog = self.mapper.read("star_catalog", epoch_parameters)
                # What we do with the catalog depends on whether it's a constant or variable PSF
                # branch.  We have to choose which star in the catalog we want to use differently in
                # these cases.
                if not self.variable_psf:
                    # For constant PSF, just transfer the first entry (i.e., the non-offset one) to
                    # the test catalog.
                    star_cat_ind = 0
                else:
                    # For variable PSF, find the index of the most aberrated PSF in the catalog.  In
                    # general, more aberrated PSFs are the ones that might have the most numerical
                    # difficulties in the rendering process.
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
        # Write to the appropriate file and directory as specified by the mapper.
        self.mapper.write(test_catalog, "star_test_catalog", epoch_parameters)

    def writeConfig(self, experiment, obs_type, shear_type, subfield_min, subfield_max):
        """This function writes yaml-style config files that can be used by GalSim to automatically
        generate the galaxy, star, and test images for this branch and range of subfields."""

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
            'noclobber' : True,
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
            # The galaxy images are large, so parallelize at the image level (unlike for the small
            # star fields for constant PSF, for which we already had set up to parallelize at the
            # file level).
            d['image']['nproc'] = self.nproc
            del d['output']['nproc']

        if self.real_galaxy:
            # Need to add the RealGalaxyCatalog to input.
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
        # (3) Finally, make a config file for the "StarTest" images.
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
            'nproc' : self.nproc,
            'noclobber' : True
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
        """This method builds and writes the galaxy images for a given subfield and epoch to disk.
        It was not used for generation of the GREAT3 simulations, since the GalSim config interface
        allows for the work done here to be parallelized much more.  However, for testing and
        generation of small images this method can be useful.
        """
        # Read in the epoch parameters and catalog.
        epoch_parameters = self.mapper.read("epoch_parameters", subfield_index=subfield_index,
                                            epoch_index=epoch_index)
        epoch_catalog = self.mapper.read("epoch_catalog", epoch_parameters)
        seed = epoch_parameters["noise_seed"]

        # Define basic numbers like pixel scale and size of postage stamps.
        pixel_scale = constants.pixel_scale[self.obs_type][self.multiepoch]
        xsize = constants.xsize[self.obs_type][self.multiepoch]
        ysize = constants.ysize[self.obs_type][self.multiepoch]

        # Set up the full image on which to place the postage stamps.
        galaxy_image = galsim.ImageF(constants.ncols * xsize, constants.nrows * ysize,
                                     scale=pixel_scale)
        galaxy_image.setOrigin(0,0)

        # Define the maximum sizes for padding RealGalaxy objects with noise.  This is necessary for
        # 'real_galaxy' and 'full' branches.
        max_xsize = xsize + 2*(constants.centroid_shift_max + constants.epoch_shift_max)
        max_ysize = ysize + 2*(constants.centroid_shift_max + constants.epoch_shift_max)

        # The GSObjects that are returned by the builders have scale sizes in arcsec, so our
        # galsim.Pixel needs to use arcsec as well.  However, some of the later manipulations
        # (shifts of the centroid within the image carried out by the draw() function) will be in
        # pixels.
        pixel = galsim.Pixel(pixel_scale)

        # We sometimes need a larger FFT than allowed by the default GSParams, so define a new
        # GSParams that will allow the larger FFT.
        params = galsim.GSParams(maximum_fft_size=10240)

        # We'll make a cache for the PSF object, since in constant PSF branches the PSF is the same
        # for all galaxies.  This way, we only build the PSF object once.
        cached_psf_obj = None
        # Loop over the objects in the galaxy catalog for this epoch.
        for record in epoch_catalog:
            # Make the RNG.
            rng = galsim.UniformDeviate(seed)
            seed = seed + 1

            # Build PSF (or take cached value if possible).  If we have to build the PSF for the
            # constant PSF branch, then save it to the cache.
            if not self.variable_psf:
                if cached_psf_obj is None:
                    psf = self.psf_builder.makeGalSimObject(record, epoch_parameters["psf"])
                    cached_psf_obj = psf
                else:
                    psf = cached_psf_obj
            else:
                psf = self.psf_builder.makeGalSimObject(record, epoch_parameters["psf"])

            # Build galaxy, apply the lensing shear and magnification, and convolve with the PSF.
            galaxy = self.galaxy_builder.makeGalSimObject(
                record, epoch_parameters["galaxy"], xsize=max_xsize, ysize=max_ysize, rng=rng)
            galaxy.applyLensing(g1=record['g1'], g2=record['g2'], mu=record['mu'])
            final = galsim.Convolve([psf, pixel, galaxy], gsparams=params)

            # Apply both offsets: the one related to dithering between epochs (same for all objects
            # in a given epoch), and the random shift of this particular object from the center of
            # the postage stamp.
            offset = galsim.PositionD(epoch_parameters['xdither'] + record['xshift'],
                                      epoch_parameters['ydither'] + record['yshift'])

            # Define postage stamp.
            bbox = galsim.BoundsI(
                xmin=int(record['xmin']), ymin=int(record['ymin']),
                xmax=int(record['xmax']), ymax=int(record['ymax']),
            )
            stamp = galaxy_image.subImage(bbox)
            # Draw into the postage stamp.
            final.draw(stamp, normalization='f', dx=pixel_scale, offset=offset)

            # Apply whitening if necessary (i.e., for 'real_galaxy' and 'full' branches, which use
            # real HST images).
            if hasattr(final, 'noise'):
                current_var = final.noise.applyWhiteningTo(stamp)
            else:
                current_var = 0.

            # The lines below are diagnostics that can be used to check that the actual S/N is
            # fairly consistent with the estimated one.  Turn it to True if you want to run this
            # code.
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
                    # of the first pixel is 1.  I do not like this estimator because if we make the
                    # postage stamp larger (with white space) then S doesn't change but N^2
                    # changes.  So let's instead use a weighted version:
                    #   S = Sum W(x,y) I(x,y) [(x-x_c)^2 + (y-y_c)^2] / Sum W(x,y)
                    #   N^2 = (noise variance) * Sum W^2(x,y) [(x-x_c)^2 + (y-y_c)^2]^2 /
                    #                                      (Sum W(x,y))^2
                    # Use W(x,y) = I(x,y),
                    #   S = Sum I(x,y)^2 [(x-x_c)^2 + (y-y_c)^2] / Sum I(x,y)
                    #   N^2 = (noise variance) * Sum I^2(x,y) [(x-x_c)^2 + (y-y_c)^2]^2 /
                    #                                      (Sum I(x,y))^2
                    #   S/N = Sum I(x,y)^2 [(x-x_c)^2 + (y-y_c)^2] /
                    #         sqrt[(noise variance) * Sum I^2(x,y) [(x-x_c)^2 + (y-y_c)^2]^2]
                    if stamp.array.shape[0] != stamp.array.shape[1]:
                        raise RuntimeError
                    min = 1.
                    max = float(stamp.array.shape[0]+1)
                    x_pix, y_pix = numpy.meshgrid(numpy.arange(min, max, 1.),
                                                  numpy.arange(min, max, 1.))
                    dx_pix = x_pix - (res.moments_centroid.x - (res.image_bounds.xmin-1))
                    dy_pix = y_pix - (res.moments_centroid.y - (res.image_bounds.ymin-1))
                    sn_size = numpy.sum(stamp.array**2 * (dx_pix**2 + dy_pix**2)) / \
                        numpy.sqrt(float(epoch_parameters['noise']['variance']) * \
                                   numpy.sum(stamp.array**2 * (dx_pix**2 + dy_pix**2)**2))
                except:
                    sn_ellip_gauss = -10.
                    sn_size = -10.
                print 'SN: ', record['gal_sn'], actual_sn_g08, sn_ellip_gauss, sn_size, \
                    record['bulge_n'], record['bulge_hlr'], record['bulge_flux']

            # Now, actually add the noise to this postage stamp.
            self.noise_builder.addNoise(rng, epoch_parameters['noise'], stamp, current_var)

        # Write the entire big image to the appropriate file, as determined by the mapper.
        self.mapper.write(galaxy_image, "image", epoch_parameters)


    def writePSFImage(self, subfield_index, epoch_index):
        """This method builds and writes the star field images for a particular subfield and epoch
        to disk.  It was not used for generation of the GREAT3 simulations, since the GalSim config
        interface allows for the work done here to be parallelized much more.  However, for testing
        and generation of small images this method can be useful.
        """
        # Read in the epoch parameters and star catalog.
        epoch_parameters = self.mapper.read("epoch_parameters", subfield_index=subfield_index,
                                            epoch_index=epoch_index)
        star_catalog = self.mapper.read("star_catalog", epoch_parameters)
        seed = epoch_parameters["noise_seed"]

        # Define basic numbers like pixel scale and size of postage stamps.
        pixel_scale = constants.pixel_scale[self.obs_type][self.multiepoch]
        xsize = constants.xsize[self.obs_type][self.multiepoch]
        ysize = constants.ysize[self.obs_type][self.multiepoch]
        # Make a galsim.Pixel representing the top-hat pixel.
        pixel = galsim.Pixel(pixel_scale)

        # Set up the image for the star field.  Its size depends on whether this is a constant PSF
        # or variable PSF branch.
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

        # Set up a cache for the galsim.GSObject corresponding to this star.  This is useful for
        # constant PSF branches, for which the star is the same in each star image (just shifted).
        cached_psf_obj = None
        for record in star_catalog:
            rng = galsim.UniformDeviate(seed)
            seed = seed + 1

            # Define the bounds of this postage stamp.
            bbox = galsim.BoundsI(
                xmin=int(record['xmin']), ymin=int(record['ymin']),
                xmax=int(record['xmax']), ymax=int(record['ymax']),
            )
            stamp = star_image.subImage(bbox)

            # Build PSF (or take cached value if possible).
            if not self.variable_psf:
                if cached_psf_obj is None:
                    psf = self.psf_builder.makeGalSimObject(record, epoch_parameters["psf"])
                    cached_psf_obj = psf
                else:
                    psf = cached_psf_obj
            else:
                psf = self.psf_builder.makeGalSimObject(record, epoch_parameters["psf"])

            # Convolve with the pixel response.
            final = galsim.Convolve([psf, pixel])
            offset = galsim.PositionD(record['xshift'],record['yshift'])
            # Draw into the postage stamp, including the centroid shift with the draw() method
            # (rather than actually shifting the GSObject).  The draw() `offset` keyword takes pixel
            # units, rather than arcsec.
            final.draw(stamp, normalization='f', offset=offset)
            # Only the variable PSF branches have noisy star fields, so add noise in that case.
            if self.variable_psf:
                self.noise_builder.addStarImageNoise(
                    rng, epoch_parameters['noise'], record['star_snr'], stamp)

        # Write the entire star field image to file.
        self.mapper.write(star_image, "starfield_image", epoch_parameters)

    def writeStarParameters(self, subfield_index, epoch_index):
        """This method writes out a dict for the PSF shapes needed for metric calculation.  The
        metric calculation for constant shear requires us to know the direction of PSF anisotropy,
        so we measure some of the star shapes from the star images.  (We cannot do this based on the
        catalogs since the stars are the composition of an aberrated optical PSF and an atmospheric
        PSF, for which the composite shape is not obvious.)
        """
        # Only do this for constant shear, not variable shear!  The star shapes are needed to create
        # the metric for the constant shear branch fits to (m, c) values.
        if self.shear_type == "variable":
            return

        # Read in epoch parameters and a star catalog.
        epoch_parameters = self.mapper.read("epoch_parameters", subfield_index=subfield_index,
                                            epoch_index=epoch_index)
        star_catalog = self.mapper.read("star_catalog", epoch_parameters)
        seed = epoch_parameters["noise_seed"]

        # Read in the star image.
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
            # We will use 1% of the stars.  Use ceil() to make sure that we don't end up with zero
            # for test runs with few stars.
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

        # Write the results to the appropriate file using the mapper.
        self.mapper.write(starshape_parameters, "starshape_parameters", starshape_parameters)

    def packagePublic(self, subfield_min, subfield_max):
        """This method packages up the public outputs (no truth values) into a single big tarfile
        for this branch.  We can choose to use a subset of the subfields if we wish."""
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
                else:
                    # If this line is not in an 'else' statement, then it will simply overwrite the
                    # correct FITS catalog with one that lacks some important columns that are
                    # necessary for variable PSF experiments!  Oops.  And then those incorrect
                    # catalogs will be converted to the text versions, which will also be wrong.
                    outfile = root_rel_mapper.copySub(sub_mapper, 'subfield_catalog', tmp_dict,
                                                      gal_use_cols,
                                                      new_template =
                                                      "deep_galaxy_catalog-%(deep_subfield_index)03d")
            tar.add(outfile)
            # ... and also copy to text file that gets added to the tarball.
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
            # ... and also copy to text file that gets added to the tarball.
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
            # These files are tiny, so let's write one for each subfield and epoch, in both yaml and
            # txt format.
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
        of its grid points.  (This eliminates the need to interpolate the shears, which is useful
        since interpolation of shear fields was not tested in GalSim until after the GREAT3 sims
        were made.)  We assume that the options for offsetting are on an subfield_grid_subsampling x
        subfield_grid_subsampling grid.  This then gives subfield_grid_subsampling^2 possible
        locations.  We then choose a random n_subfields_per_field-1 of those options for the
        subfields that are not the first in the field.

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
