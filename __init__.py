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
from .builder import SimBuilder

builders = {
    "control": SimBuilder.customize("control"),
    "real_galaxy": SimBuilder.customize("real_galaxy", real_galaxy=True),
    "variable_psf": SimBuilder.customize("variable_psf", variable_psf=True),
    "multiepoch": SimBuilder.customize("multiepoch", multiepoch=True),
    "full": SimBuilder.customize("full", real_galaxy=True, variable_psf=True, multiepoch=True),
    }

def run(root, experiments=None, obs_type=None, shear_type=None, seed=10, steps=None,
        subfield_min=0, subfield_max=(constants.n_subfields-1),
        gal_dir='/Users/rmandelb/great3/data-23.5', ps_dir='inputs/ps/tables',
        opt_psf_dir = '../inputs/optical-psfs',
        atmos_ps_dir = '../inputs/atmospsf/pk_math',
        public_dir='public', truth_dir='truth', preload=False, nproc=-1):
    """Top-level driver for GREAT3 simulation code.

    This driver parses the input parameters to decide what work must be done.  Here are the
    definitions for some of the terms used throughout this driver and the rest of the code:

    experiment - one of the 5 GREAT3 experiments ("control", "real_galaxy", "variable_psf",
    "multiepoch", "full")

    obs_type - type of observation ("ground", "space")

    shear_type - type of input shear field ("variable", "constant")

    branch - a set of simulations corresponding to a particular combination of (experiment,
    obs_type, shear_type).  There are 5 x 2 x 2 = 20 branches.

    steps - the distinct stages of building the simulations:
        'metaparameters' = building the global parameters of the catalogs for a given branch
        'catalogs' = writing the per-subfield and per-epoch catalogs listing the characteristics 
                       of each object
        'config' = writing the config files for galsim_yaml to building the images
        'gal_images' = making the galaxy images within this script
        'psf_images' = making the PSF images within this script
        'star_params' = measure the HSM shapes for the PSF images
        'packages' = packing up the tarballs to be distributed publicly and for metric evaluation
    In order to carry out the later stages without the former, the script must have
    already produced output files for the first stages in the specified directory.

    field - a slightly more than 10x10 degree area of simulated sky that samples the same underlying
    shear field and PSF pattern (in the case of varying PSF).

    subfield - a region corresponding to a particular field, containing just a subset of the
    galaxies in that simulated area of sky purely for convenience sake (so the galaxies do not
    overlap, etc.).  It is 10x10 degrees in extent, and compared to the other subfields in its
    field, it is slightly offset.

    epoch - a single observation of a given subfield at some epoch (i.e., in the single epoch
    experiments, there is one observation per subfield, whereas in the multi-epoch case, there are
    >1, as determined in constants.py)

    tile - a 2x2 degree subset of a subfield on which PSF properties are determined for the variable
    PSF case.  This terminology is only needed by the PSF builders, not __init__.py or builder.py

    mapper - a class that manages I/O within a branch (specifying file formats and so on)

    Within a given branch, the hierarchy is branch -> field -> subfield -> epoch.  The numbers that
    determine how many of each is in a branch are in constants.py (number of subfields per branch,
    number of subfields per field, number of epochs per subfield).

    @param[in] root          root directory for generated files
    @param[in] experiments   sequence of simulation experiments to build
                             (None == all 5 experiments)
    @param[in] obs_type      observation type for simulations, one of 'ground' or 'space'
                             (None == both)
    @param[in] shear_type    shear field type for simulations, one of 'constant' or 'variable'
                             (None == both)
    @param[in] seed          first random number seed (ignored if the steps parameter does not
                             include 'metaparameters'); random number generators will be seeded with
                             this number *and* the integers immediately after it.  Note that seed
                             should not be 0, since this does not result in a repeatable random
                             number sequence.
    @param[in] steps         a sequence of steps to simulate; valid entries are 'metaparameters',
                             'catalogs', 'config', 'gal_images', 'psf_images', 'star_params',
                             'packages' (None == all steps). 
                             Previous steps must have already been run with the same root.
    @param[in] subfield_min  Minimum range "subfield index" to limit the number of catalogs/images
                             produced e.g. when testing.
    @param[in] subfield_max  Maximum range "subfield index" to limit the number of catalogs/images
                             produced e.g. when testing.
    @param[in] gal_dir       Directory containing real galaxy catalogs.
    @param[in] ps_dir        Directory containing cosmological power spectrum inputs from iCosmo.
    @param[in] opt_psf_dir   Directory containing the optical PSF models for ground and space
                             variable PSF simulations.
    @param[in] atmos_ps_dir  Directory containing atmospheric power spectrum inputs based on
                             Mathematica numerical integration.
    @param[in] public_dir    Directory containing files to be distributed publicly.
    @param[in] truth_dir     Directory containing files used for metric evaluation.
    @param[in] preload       Preload the RealGalaxyCatalog images to speed up generation of large
                             numbers of real galaxies? (default=False)  Note that for parametric
                             galaxy branches, the catalog is never preloaded.
    @param[in] nproc         How many processes to use in the config file.  (default = -1)
    """
    import sys

    # Select experiments based on keywords, or do all of them if no experiment was specified.
    if experiments is None:
        experiments = builders.keys()
    elif isinstance(experiments,basestring):
        experiments = [experiment]

    # Specify ground or space observations according to keywords, or do both if none was specified.
    if obs_type is None:
        obs_types = ["ground", "space"]
    elif isinstance(obs_type,basestring):
        obs_types = [obs_type]
    else:
        obs_types = obs_type

    # Specify shear type according to keywords, or do both if none was specified.
    if shear_type is None:
        shear_types = ["constant", "variable"]
    elif isinstance(shear_type,basestring):
        shear_types = [shear_type]
    else:
        shear_types = shear_type

    # Choose which step of the simulation process to do, or do all if no steps were specified.
    if steps is None:
        steps = ["metaparameters", "catalogs", "config", "gal_images", "psf_images", "star_params",
                 "packages"]
    elif isinstance(steps,basestring):
        steps = [steps]

    # Based on all the above choices, decide which branch to simulate, where a branch is a
    # combination of experiment [5] x observation type [2] x shear type [2] = a maximum of 20
    # possible branches.
    branches = [ (experiment, obs_type, shear_type,
                  builders[experiment](root, obs_type, shear_type,
                                       gal_dir, ps_dir, opt_psf_dir, atmos_ps_dir, public_dir,
                                       truth_dir, preload, nproc))
                        for experiment in experiments
                        for obs_type in obs_types
                        for shear_type in shear_types ]

    # Now actually do the requested work for each step of the process.
    if 'metaparameters' in steps:
        for experiment, obs_type, shear_type, builder in branches:
            print "Generating metaparameters for %s / %s / %s" % (experiment, obs_type, shear_type)
            builder.writeParameters(seed)
            if subfield_min>0 or subfield_max<constants.n_subfields-1:
                import warnings
                warnings.warn('Regenerating metaparameters for all subfields,' +
                              ' not just the chosen subset')

    if 'catalogs' in steps:
        for experiment, obs_type, shear_type, builder in branches:
            print "Generating catalogs for %s / %s / %s" % (experiment, obs_type, shear_type)
            for subfield_index in xrange(subfield_min, subfield_max+1):
                sys.stdout.write('\r')
                sys.stdout.write("  subfield %d" % subfield_index)
                sys.stdout.flush()
                builder.writeSubfieldCatalog(subfield_index)
                for epoch_index in xrange(builder.n_epochs):
                    builder.writeEpochCatalog(subfield_index, epoch_index)
            builder.writeStarTestCatalog(subfield_min, subfield_max)
            sys.stderr.write("\n")

    if 'config' in steps:
        for experiment, obs_type, shear_type, builder in branches:
            print "Generating config files for %s / %s / %s" % (experiment, obs_type, shear_type)
            builder.writeConfig(experiment, obs_type, shear_type, subfield_min, subfield_max)

    if 'gal_images' in steps:
        for experiment, obs_type, shear_type, builder in branches:
            print "Generating galaxy images for %s / %s / %s" % (experiment, obs_type, shear_type)
            for subfield_index in xrange(subfield_min, subfield_max+1):
                for epoch_index in xrange(builder.n_epochs):
                    sys.stdout.write('\r')
                    sys.stdout.write("  subfield %d / epoch %d" % (subfield_index, epoch_index))
                    sys.stdout.flush()
                    builder.writeGalImage(subfield_index, epoch_index)
            sys.stderr.write("\n")
            if hasattr(builder.galaxy_builder,'rgc'):
                builder.galaxy_builder.rgc.close()

    if 'psf_images' in steps:
        for experiment, obs_type, shear_type, builder in branches:
            print "Generating PSF images for %s / %s / %s" % (experiment, obs_type, shear_type)
            for subfield_index in xrange(subfield_min, subfield_max+1):
                for epoch_index in xrange(builder.n_epochs):
                    sys.stdout.write('\r')
                    sys.stdout.write("  subfield %d / epoch %d" % (subfield_index, epoch_index))
                    sys.stdout.flush()
                    builder.writePSFImage(subfield_index, epoch_index)
            sys.stderr.write("\n")

    if 'star_params' in steps:
        for experiment, obs_type, shear_type, builder in branches:
            print "Measuring star parameters for %s / %s / %s" % (experiment, obs_type, shear_type)
            for subfield_index in xrange(subfield_min, subfield_max+1):
                for epoch_index in xrange(builder.n_epochs):
                    sys.stdout.write('\r')
                    sys.stdout.write("  subfield %d / epoch %d" % (subfield_index, epoch_index))
                    sys.stdout.flush()
                    builder.writeStarParameters(subfield_index, epoch_index)
            sys.stderr.write("\n")

    if 'packages' in steps:
        for experiment, obs_type, shear_type, builder in branches:
            print "Packaging data for %s / %s / %s" % (experiment, obs_type, shear_type)
            builder.packagePublic(subfield_min, subfield_max)
            builder.packageTruth(subfield_min, subfield_max)
            sys.stderr.write("\n")
