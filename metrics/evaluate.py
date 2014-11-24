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
"""@file evaluate.py

Module containing functions for the evaluation of the constant and variable shear-type branches
of the GREAT3 Challenge.

For the functions used to directly evaluate the GREAT3 metrics, see q_constant() and q_variable().

This code requires the GREAT3 metric evaluation truth data to have been unpacked into the directory
specified by the path in TRUTH_DIR, set below.  Please update TRUTH_DIR to match the location of the
truth data on your system.

For information about getting this data, please see the following page on the great3-public wiki:
https://github.com/barnabytprowe/great3-public/wiki/Metric-evaluation-scripts-and-truth-data

Constant shear branches
-----------------------
Each submission file (one per branch) should take the format of a 3-column ASCII catalog, e.g.:

    #  SUBFIELD_INDEX  G1  G2
    0   -.26664  0.11230
    1   -.13004  0.48103
    ...

or similar.  The hashed header/comment can be omitted, and almost all formats for the numbers are
fine.  The main criterion to be satisfied is that after

    >>> data = np.loadtxt(submission_file)

the `data` object must be a NumPy array with shape `(NSUBFIELDS, 3)`, where `NSUBFIELDS` is the
total number of subfields in the branch (currently fixed at 200).

In addition, the array slice `data[:, 0]` must be the subfield number in ascending order from `0` to
`NSUBFIELDS - 1`.  The array slices `data[:, 1]` and `data[:, 2]` must be the corresponding
estimates of mean shear g1 and g2 in each subfield.

In these details the submission should match the output of the helper script `presubmission.py`
available at https://github.com/barnabytprowe/great3-public .
 
Variable shear branches
-----------------------
Each submission file (one per branch) should take the format of an ASCII catalog with a minimum of
3 columns as in the following example:

    # FIELD_INDEX  THETA [degrees]  MAP_E
    0  0.0246  2.5650e-06
    0  0.0372  4.1300e-06
    ...

The `FIELD_INDEX` will be an integer between 0 and 9.  `THETA` should be a sequence of floats giving
the annular bin centres in degrees; these are logarithmically spaced between the minimum separation
considered (0.01 degrees) and the maximum (10.0 degrees).  `MAP_E` is the E-mode aperture mass
dispersion.  `FIELD`, `THETA` (and the thus corresponding `MAP_E` entries) must be ordered as in the
output of `presubmission.py`.

The hashed header/comment can be omitted.  Additional columns can be present provided that the
location and order of the three described above  Additional columns can be present provided that the
location and order of the three described above.  An example of this is the output of 
`presubmission.py` for variable shear branches, which also append columns for the B-mode aperture
mass dispersion and a (shot noise only) error estimate.

After

    >>> data = np.loadtxt(submission_file)

the `data` object must be a NumPy array with shape `(NFIELDS * NBINS_THETA, n)`, where `NFIELDS` is
the total number of fields in the branch (currently fixed at 10), `NBINS_THETA` is the number of
annular bins in angle used to estimate Map_E in each field (currently fixed at 15), and `n >= 3`.

As mentioned, in these details the submission should match the output of the helper script
`presubmission.py` available at https://github.com/barnabytprowe/great3-public .
"""

import os
import sys
import logging
import numpy as np

try:
    import great3sims
    import great3sims.mapper
except ImportError:
    path, module = os.path.split(__file__)
    sys.path.append(os.path.join(path, "..")) # Appends the folder great3-public/ to sys.path
    import great3sims
    import great3sims.mapper

try:
    import g3metrics
except ImportError:
    path, module = os.path.split(__file__)
    sys.path.append(os.path.join(path, "..", "metrics")) # Appends the great3-private/metrics
                                                               # folder to path
    import g3metrics

TRUTH_DIR = "/great3/beta/truth"  # Root folder in which the truth values are upacked

NFIELDS = 10     # Total number of fields
NSUBFIELDS = 200 # Total number of subfields, not necessarily equal to the number of subfields made
                 # in mass_produce as that script also generates the deep fields
NSUBFIELDS_PER_FIELD = NSUBFIELDS / NFIELDS
NGALS_PER_SUBFIELD = 10000 # 100x100 galaxies per subfield

CFID = 2.e-4
MFID = 2.e-3

XMAX_GRID_DEG = 10.0 # Maximum image spatial extent in degrees
DX_GRID_DEG = 0.1    # Grid spacing in degrees

THETA_MIN_DEG = 0.02 # Minimum and maximum angular scales for logarithmic bins used to calculate the
THETA_MAX_DEG = 10.0 # aperture mass disp. - MUST match specs given to participants - in degrees
NBINS_THETA = 15     # Number of logarithmic bins theta for the aperture mass dispersion

EXPECTED_THETA = np.array([ # Array of theta values expected in submissions, good to 3 d.p.
    0.0246,  0.0372,  0.0563,  0.0853,  0.1290,  0.1953,  0.2955, 0.4472,  0.6768,  1.0242,  1.5499,
    2.3455,  3.5495,  5.3716,  8.1289] * NFIELDS)

USEBINS = np.array([ # Which of the theta bins above to actually use in calculating the metric?
    False,   False,   False,   True,    True,    True,    True,   True,    True,    True,    True,
    True,    True,    False,   False] * NFIELDS) # Note the *NFIELDS to match per-field theta layout

STORAGE_DIR = "./metric_calculation_products" # Folder into which to store useful intermediate
                                              # outputs of metric calculations (e.g. rotation files,
                                              # dicts, mapE tables) which need be calculated only
                                              # once

SUBFIELD_DICT_FILE_PREFIX = "subfield_dict_"
GTRUTH_FILE_PREFIX = "gtruth_"
ROTATIONS_FILE_PREFIX = "rotations_"
OFFSETS_FILE_PREFIX = "offsets_"
MAPESHEAR_FILE_PREFIX = "mapEshear_"
MAPEINT_FILE_PREFIX = "mapEint_"
MAPEOBS_FILE_PREFIX = "mapEobs_"

# These constant normalization factors come from a run of ~1000 sims done on 6 Jan 2014, modified on
# 30 Jan 2014 to bring space and ground into agreement at high bias
NORMALIZATION_CONSTANT_SPACE = 1.232
NORMALIZATION_CONSTANT_GROUND = NORMALIZATION_CONSTANT_SPACE

NORMALIZATION_VARIABLE_SPACE = 0.0001837  # Factor comes from tests with
                                          # tabulate_variable_shear_metric_rev1.py on 1000 runs and
                                          # NOISE_SIGMA = 0.10, 6 Jan 2015, with sigma2_min = 4.e-8
NORMALIZATION_VARIABLE_GROUND = NORMALIZATION_VARIABLE_SPACE  # Bring space=ground at high bias

# Values of sigma2_min to adopt as the defaults for the Q_c and Q_v metrics, as of 30 Dec 2013.
# These parameters add a damping
SIGMA2_MIN_CONSTANT_GROUND = 4.     # 2**2
SIGMA2_MIN_CONSTANT_SPACE = 1.      # 1**2
SIGMA2_MIN_VARIABLE_GROUND = 9.e-8  # [2 * 1.e-3]**2
SIGMA2_MIN_VARIABLE_SPACE = 4.e-8   # [3 * 1.e-3]**2


def get_generate_const_truth(experiment, obs_type, truth_dir=TRUTH_DIR, storage_dir=STORAGE_DIR,
                             logger=None):
    """Get or generate arrays of subfield_index, g1true, g2true, each of length `NSUBFIELDS`.

    If the gtruth file has already been built for this constant shear branch, loads and returns the
    saved copies.

    If the array of truth values has not been built, or is older than the first entry in the set of
    shear_params files, the arrays are built first, saved to file, then returned.

    @param experiment     Experiment for this branch, one of 'control', 'real_galaxy',
                          'variable_psf', 'multiepoch', 'full'
    @param obs_type       Observation type for this branch, one of 'ground' or 'space'
    @param storage_dir    Directory from/into which to load/store rotation files
    @param truth_dir      Root directory in which the truth information for the challenge is stored
    @param logger         Python logging.Logger instance, for message logging
    @return subfield_index, g1true, g2true
    """
    gtruefile = os.path.join(storage_dir, GTRUTH_FILE_PREFIX+experiment[0]+obs_type[0]+"c.asc")
    mapper = great3sims.mapper.Mapper(truth_dir, experiment, obs_type, "constant")
    use_stored = True
    if not os.path.isfile(gtruefile):
        use_stored = False
        if logger is not None:
            logger.info(
                "First build of shear truth tables using values from "+
                os.path.join(mapper.full_dir, "shear_params-*.yaml"))
    else:
        # Compare timestamps for the gtruefile and the first shear_params file
        # (subfield = 000) for this branch.  If the former is older than the latter, or this file,
        # force rebuild...
        gtruemtime = os.path.getmtime(gtruefile)
        shear_params_file = os.path.join(mapper.full_dir, "shear_params-000.yaml")
        shear_params_mtime = os.path.getmtime(shear_params_file)
        if gtruemtime < shear_params_mtime or gtruemtime < os.path.getmtime(__file__):
            use_stored = False
            if logger is not None:
                logger.info(
                    "Updating out-of-date shear truth tables using newer values from "+
                    os.path.join(mapper.full_dir, "shear_params-*.yaml"))
    # Then load or build (and save) the array of truth values per subfield
    if use_stored:
        if logger is not None:
            logger.info("Loading shear truth tables from "+gtruefile) 
        gtruedata = np.loadtxt(gtruefile)
    else: 
        params_prefix = os.path.join(mapper.full_dir, "shear_params-")
        import yaml
        # Check to see if this is a variable_psf or full branch, in which case we only need the
        # first entry from each set of subfields
        if experiment in ("variable_psf", "full"):
            gtruedata = np.empty((NFIELDS, 3))
            gtruedata[:, 0] = np.arange(NFIELDS)
            subfield_index_targets = range(0, NSUBFIELDS, NSUBFIELDS_PER_FIELD)
        else:
            gtruedata = np.empty((NSUBFIELDS, 3))
            gtruedata[:, 0] = np.arange(NSUBFIELDS)
            subfield_index_targets = range(NSUBFIELDS)
        # Then loop over the required subfields reading in the shears
        for i, subfield_index in enumerate(subfield_index_targets):

            params_file = params_prefix+("%03d" % subfield_index)+".yaml"
            with open(params_file, "rb") as funit:
                gdict = yaml.load(funit)
                gtruedata[i, 1] = gdict["g1"]
                gtruedata[i, 2] = gdict["g2"]

        if logger is not None:
            logger.info("Saving shear truth table to "+gtruefile)
        if not os.path.isdir(storage_dir):
            os.mkdir(storage_dir)
        with open(gtruefile, "wb") as fout:
            fout.write("# True shears for "+experiment+"-"+obs_type+"-constant\n")
            fout.write("# subfield_index  g1true  g2true\n")
            np.savetxt(fout, gtruedata, fmt=" %4d %+.18e %+.18e")
    return (gtruedata[:, 0]).astype(int), gtruedata[:, 1], gtruedata[:, 2]

def get_generate_const_rotations(experiment, obs_type, storage_dir=STORAGE_DIR, truth_dir=TRUTH_DIR,
                                 logger=None):
    """Get or generate an array of rotation angles for Q_const calculation.

    If the rotation file has already been built for this constant shear branch, loads and returns an
    array of rotation angles to align with the PSF.  This array is of shape `(NSUBFIELDS,)`, having
    averaged over the `n_epochs` epochs in the case of multi-epoch branches.

    If the rotation file has not been built, or is older than the first entry in the set of
    starshape_parameters files, the array of rotations is built, saved to file, then returned.

    @param experiment     Experiment for this branch, one of 'control', 'real_galaxy',
                          'variable_psf', 'multiepoch', 'full'
    @param obs_type       Observation type for this branch, one of 'ground' or 'space'
    @param storage_dir    Directory from/into which to load/store rotation files
    @param truth_dir      Root directory in which the truth information for the challenge is stored
    @param logger         Python logging.Logger instance, for message logging
    @return rotations     An array containing all the rotation angles, in radians
    """
    import great3sims
    rotfile = os.path.join(storage_dir, ROTATIONS_FILE_PREFIX+experiment[0]+obs_type[0]+"c.asc")
    mapper = great3sims.mapper.Mapper(truth_dir, experiment, obs_type, "constant")
    use_stored = True
    if not os.path.isfile(rotfile):
        use_stored = False
        if logger is not None:
            logger.info(
                "First build of rotations file using starshape_parameters from "+
                mapper.full_dir)
    else:
        # Then compare timestamps for the rotation file and the first starshape_parameters file
        # (subfield = 000, epoch =0) for this branch.  If the former is older than the latter, or
        # this file, force rebuild...
        rotmtime = os.path.getmtime(rotfile)
        starshape_file_template, _ ,_ = mapper.mappings['starshape_parameters']  
        starshape_file = os.path.join(
            mapper.full_dir, starshape_file_template % {"epoch_index": 0, "subfield_index": 0})
        starshapemtime = os.path.getmtime(starshape_file+".yaml")
        if rotmtime < starshapemtime or rotmtime < os.path.getmtime(__file__):
            use_stored = False
            if logger is not None:
                logger.info(
                    "Updating out-of-date rotations file using newer starshape_parameters from "+
                    mapper.full_dir)
    # Then load / build as required 
    if use_stored is True:
        if logger is not None:
            logger.info("Loading rotations from "+rotfile) 
        rotations = np.loadtxt(rotfile)[:, 1] # First column is just subfield indices
    else:
        # To build we must loop over all the subfields and epochs
        # First work out if the experiment is multi-exposure and has multiple epochs
        if experiment in ("multiepoch", "full"):
            import great3sims.constants
            n_epochs = great3sims.constants.n_epochs
        else:
            n_epochs = 1
        # Setup the array for storing the PSF values from which rotations are calculated
        psf_g1 = np.empty((NSUBFIELDS, n_epochs))
        psf_g2 = np.empty((NSUBFIELDS, n_epochs))
        mean_psf_g1 = np.empty(NSUBFIELDS)
        mean_psf_g2 = np.empty(NSUBFIELDS)
        for subfield_index in range(NSUBFIELDS):

            n_ignore = 0 # Counter for how many epochs had flagged, bad PSF g1/g2 values
            for epoch_index in range(n_epochs):

                starshape_parameters = mapper.read(
                    "starshape_parameters",
                    data_id={"epoch_index": epoch_index, "subfield_index": subfield_index})
                star_g1 = starshape_parameters["psf_g1"]
                star_g2 = starshape_parameters["psf_g2"]
                # Test for flagged failures (these do happen rarely and are given the value
                # psf_g1=psf_g2=-10.0, see writeStarParameters in great3sims/builders.py)
                # If the psf ellipticities are failed, we just ignore these for the (m, c) calcs
                if star_g1 > -9.9 and star_g2 > -9.9:
                    psf_g1[subfield_index, epoch_index] = star_g1
                    psf_g2[subfield_index, epoch_index] = star_g2
                else:
                    n_ignore += 1
                    psf_g1[subfield_index, epoch_index] = 0.
                    psf_g2[subfield_index, epoch_index] = 0.

            # Calculate the mean across the epochs in this subfield taking any flagged values into
            # account
            n_eff = n_epochs - n_ignore
            if n_eff > 0:
                mean_psf_g1[subfield_index] = (psf_g1[subfield_index, :]).sum() / float(n_eff) 
                mean_psf_g2[subfield_index] = (psf_g2[subfield_index, :]).sum() / float(n_eff)
            else:
                mean_psf_g1[subfield_index] = 0. # This is safe in np.arctan2() -> 0. 
                mean_psf_g2[subfield_index] = 0.

        if experiment in ("variable_psf", "full"):
            # Average over all subfields per field
            final_psf_g1 = np.empty(NFIELDS)
            final_psf_g2 = np.empty(NFIELDS)
            for i in range(NFIELDS):

                final_psf_g1[i] = np.mean(
                    mean_psf_g1[i * NSUBFIELDS_PER_FIELD: (i + 1) * NSUBFIELDS_PER_FIELD])
                final_psf_g2[i] = np.mean(
                    mean_psf_g2[i * NSUBFIELDS_PER_FIELD: (i + 1) * NSUBFIELDS_PER_FIELD])

        else:
            final_psf_g1 = mean_psf_g1
            final_psf_g2 = mean_psf_g2
 
        rotations = .5 * np.arctan2(final_psf_g2, final_psf_g1)
        # We have built rotations, but then save this file as ascii for use next time
        if logger is not None:
            logger.info("Saving rotations to "+rotfile)    
        if not os.path.isdir(storage_dir):
            os.mkdir(storage_dir)
        with open(rotfile, "wb") as fout:
            fout.write("# Rotations for "+experiment+"-"+obs_type+"-constant\n")
            fout.write("# subfield_index  rotation [radians]\n")
            np.savetxt(fout, np.array((np.arange(len(rotations)), rotations)).T, fmt=" %4d %+.18f")
    return rotations

def get_generate_variable_offsets(experiment, obs_type, storage_dir=STORAGE_DIR,
                                  truth_dir=TRUTH_DIR, logger=None):
    """Get or generate arrays of subfield_index, offset_deg_x, offset_deg_y, each of length
    `NSUBFIELDS`.

    If the offsets file has already been built for this variable shear branch, loads and returns the
    saved arrays.

    If the arrays of offset values have not been built, or are older than the first entry in the set
    of subfield_offset files, the arrays are built first, saved to file, then returned.

    @param experiment     Experiment for this branch, one of 'control', 'real_galaxy',
                          'variable_psf', 'multiepoch', 'full'
    @param obs_type       Observation type for this branch, one of 'ground' or 'space'
    @param storage_dir    Directory from/into which to load/store rotation files
    @param truth_dir      Root directory in which the truth information for the challenge is stored
    @param logger         Python logging.Logger instance, for message logging
    @return subfield_index, offset_deg_x, offset_deg_y
    """
    offsetfile = os.path.join(storage_dir, OFFSETS_FILE_PREFIX+experiment[0]+obs_type[0]+"v.asc") 
    mapper = great3sims.mapper.Mapper(truth_dir, experiment, obs_type, "variable")
    use_stored = True
    if not os.path.isfile(offsetfile):
        use_stored = False
        if logger is not None:
            logger.info(
                "First build of offsets file using subfield_offset files from "+
                mapper.full_dir)
    else:
        # Then compare timestamps for the offsets file and the first file
        # (subfield = 000) for this branch.  If the former is older than the latter, or
        # this file, force rebuild...
        offsetmtime = os.path.getmtime(offsetfile)
        subfield_offset_file = os.path.join(mapper.full_dir, "subfield_offset-000.yaml")
        subfield_offset_mtime = os.path.getmtime(subfield_offset_file)
        if offsetmtime < subfield_offset_mtime or offsetmtime < os.path.getmtime(__file__):
            use_stored = False 
            if logger is not None:
                logger.info(
                    "Updating out-of-date offset file using newer values from "+
                    os.path.join(mapper.full_dir, "subfield_offset-*.yaml"))
    # Then load / build as required 
    if use_stored is True:
        if logger is not None:
            logger.info("Loading offsets from "+offsetfile) 
        offsets = np.loadtxt(offsetfile)
    else:
        offsets_prefix = os.path.join(mapper.full_dir, "subfield_offset-") 
        offsets = np.empty((NSUBFIELDS, 3))
        import yaml
        offsets[:, 0] = np.arange(NSUBFIELDS)
        for i in range(NSUBFIELDS):

            offsets_file = offsets_prefix+("%03d" % i)+".yaml"
            with open(offsets_file, "rb") as funit:
                offsetdict = yaml.load(funit)
                offsets[i, 1] = offsetdict["offset_deg_x"]
                offsets[i, 2] = offsetdict["offset_deg_y"]

        if logger is not None:
            logger.info("Saving offset file to "+offsetfile)
        if not os.path.isdir(storage_dir):
            os.mkdir(storage_dir)
        with open(offsetfile, "wb") as fout:
            fout.write("# Subfield offsets for "+experiment+"-"+obs_type+"-variable\n")
            fout.write("# subfield_index  offset_deg_x  offset_deg_y\n")
            np.savetxt(fout, offsets, fmt=" %4d %.18e %.18e")
    return (offsets[:, 0]).astype(int), offsets[:, 1], offsets[:, 2]

def run_corr2(x, y, e1, e2, w, min_sep=THETA_MIN_DEG, max_sep=THETA_MAX_DEG, nbins=NBINS_THETA,
              cat_file_suffix='_temp.fits', params_file_suffix='_corr2.params',
              m2_file_suffix='_temp.m2', xy_units='degrees', sep_units='degrees',
              corr2_executable='corr2'):
    """Copied from presubmission.py
    """
    import pyfits
    import subprocess
    import tempfile
    # Create temporary, unique files for I/O
    catfile = tempfile.mktemp(suffix=cat_file_suffix)
    paramsfile = tempfile.mktemp(suffix=params_file_suffix)
    m2file = tempfile.mktemp(suffix=m2_file_suffix)
    # Write the basic corr2.params to temp location
    print_basic_corr2_params(paramsfile, min_sep=min_sep, max_sep=max_sep, nbins=nbins,
                             xy_units=xy_units, sep_units=sep_units, fits_columns=True)
    # Use fits binary table for faster I/O. (Converting to/from strings is slow.)
    # First, make the data into np arrays
    x_array = np.asarray(x).flatten()
    y_array = np.asarray(y).flatten()
    g1_array = np.asarray(e1).flatten()
    g2_array = np.asarray(e2).flatten()
    w_array = np.asarray(w).flatten()
    # Then, mask out the >= 10 values
    use_mask = np.logical_and.reduce([g1_array<10.,g2_array<10.])
    # And finally make the FITS file
    x_col = pyfits.Column(name='x', format='1D', array=x_array[use_mask])
    y_col = pyfits.Column(name='y', format='1D', array=y_array[use_mask])
    g1_col = pyfits.Column(name='g1', format='1D', array=g1_array[use_mask])
    g2_col = pyfits.Column(name='g2', format='1D', array=g2_array[use_mask])
    w_col = pyfits.Column(name='w', format='1D', array=w_array[use_mask])
    cols = pyfits.ColDefs([x_col, y_col, g1_col, g2_col, w_col])
    table = pyfits.new_table(cols)
    phdu = pyfits.PrimaryHDU()
    hdus = pyfits.HDUList([phdu, table])
    hdus.writeto(catfile, clobber=True)
    subprocess.Popen([
        corr2_executable, str(paramsfile), 'file_name='+str(catfile), 'm2_file_name='+str(m2file)
        ]).wait()
    results = np.loadtxt(m2file)
    os.remove(paramsfile)
    os.remove(catfile)
    os.remove(m2file)
    return results

def print_basic_corr2_params(outfile, min_sep=THETA_MIN_DEG, max_sep=THETA_MAX_DEG,
                             nbins=NBINS_THETA, xy_units='degrees', sep_units='degrees',
                             fits_columns=False):
    """Write a bare-bones corr2.params file (used by corr2) to the file named outfile.
    """
    with open(outfile, 'wb') as fout:
        if fits_columns:
            fout.write("# Column description\n")
            fout.write("x_col = x\n")
            fout.write("y_col = y\n")
            fout.write("g1_col = g1\n")
            fout.write("g2_col = g2\n")
            fout.write("w_col = w\n")
            fout.write("\n")
            fout.write("# File info\n")
            fout.write("file_type=FITS")
        else:
            fout.write("# Column description\n")
            fout.write("x_col = 1\n")
            fout.write("y_col = 2\n")
            fout.write("g1_col = 3\n")
            fout.write("g2_col = 4\n")
            fout.write("w_col = 5\n")
        fout.write("\n")
        fout.write(
            "# Assume sign conventions for gamma were correct in the catalog passed to "+
            "presubmission.py\n")
        fout.write("flip_g1 = false\n")
        fout.write("flip_g2 = false\n")
        fout.write("\n")
        fout.write("# Describe the parameters of the requested correlation function\n")
        fout.write('min_sep=%f\n'%min_sep)
        fout.write('max_sep=%f\n'%max_sep)
        fout.write('nbins=%f\n'%nbins)
        fout.write('x_units='+str(xy_units)+'\n')
        fout.write('y_units='+str(xy_units)+'\n')
        fout.write('sep_units='+str(sep_units)+'\n')
        fout.write('\n')
        fout.write("# verbose specifies how much progress output the code should emit.\n")
        fout.write("verbose = 0\n")
        fout.write("\n")

def get_generate_variable_truth(experiment, obs_type, storage_dir=STORAGE_DIR, truth_dir=TRUTH_DIR,
                                logger=None, corr2_exec="corr2", make_plots=False,
                                file_prefixes=("galaxy_catalog",), suffixes=("",),
                                mape_file_prefix=MAPESHEAR_FILE_PREFIX, output_xy_prefix=None):
    """Get or generate an array of truth map_E vectors for all the fields in this branch.

    If the map_E truth file has already been built for this variable shear branch, loads and returns
    the saved copies.

    If the array of truth values has not been built, or is older than the first entry in the set of
    galaxy_catalog files, the arrays are built first, saved to file, then returned.

    @param experiment        Experiment for this branch, one of 'control', 'real_galaxy',
                             'variable_psf', 'multiepoch', 'full'
    @param obs_type          Observation type for this branch, one of 'ground' or 'space'
    @param storage_dir       Directory from/into which to load/store rotation files
    @param truth_dir         Root directory in which the truth info for the challenge is stored
    @param logger            Python logging.Logger instance, for message logging
    @param corr2_exec        Path to Mike Jarvis' corr2 exectuable
    @param make_plots        Generate plotting output
    @param file_prefixes     Tuple containing one or more prefixes for file type in which to load
                             up shears, summing shears when `len(file_prefixes) >= 2`
                             [default = `("galaxy_catalog",)`]
    @param suffixes          Load up shear from entries "g1"+suffixes[0] and "g2"+suffixes[0] in the
                             `file_prefixes[0]`-type files, then add "g1"+suffixes[1] from
                             `file_prefixes[1]`-type files, etc.  Must be same length as 
                             `file_prefixes` tuple [default = `("",)`]
    @param mape_file_prefix  Prefix for output filename
    @param output_xy_prefix  Filename prefix (and switch if not None) for x-y position debug output
    @return field, theta, map_E, map_B, maperr
    """
    # Sanity check on suffixes & prefixes
    if len(suffixes) != len(file_prefixes):
        raise ValueError("Input file_prefixes and suffixes kwargs must be same length.")
    # Build basic x and y grids to use for coord positions: note we do this here rather than as
    # needed later so as to check the dimensions (meshgrid is very quick anyway)
    xgrid_deg, ygrid_deg = np.meshgrid(
        np.arange(0., XMAX_GRID_DEG, DX_GRID_DEG), np.arange(0., XMAX_GRID_DEG, DX_GRID_DEG))
    xgrid_deg = xgrid_deg.flatten() # Flatten these - the default C ordering corresponds to the way
    ygrid_deg = ygrid_deg.flatten() # the true shears are ordered too, which is handy
    if len(xgrid_deg) != NGALS_PER_SUBFIELD:
        raise ValueError(
            "Dimensions of xgrid_deg and ygrid_deg do not match NGALS_PER_SUBFIELD.  Please check "+
            "the values of XMAX_GRID_DEG and DX_GRID_DEG in evaluate.py.")
    # Define storage file and check for its existence and/or age
    mapEtruefile = os.path.join(
        storage_dir, mape_file_prefix+experiment[0]+obs_type[0]+"v.asc")
    mapper = great3sims.mapper.Mapper(truth_dir, experiment, obs_type, "variable") 
    use_stored = True
    if not os.path.isfile(mapEtruefile):
        use_stored = False
        if logger is not None:
            logger.info(
                "First build of map_E truth file using "+str(file_prefixes)+" files from "+
                mapper.full_dir)
    else:
        # Then compare timestamps for the mapE file and the newest file_prefixes[:]-000.fits file
        # (subfield = 000) for this branch.  If the former is older than the latter, or
        # this file, force rebuild...
        mapEmtime = os.path.getmtime(mapEtruefile)
        catalogmtime = 0 # Set earliest possible T
        for prefix in file_prefixes:

            catalog_file = os.path.join(mapper.full_dir, prefix+"-000.fits")
            tmpmtime = os.path.getmtime(catalog_file)
            if tmpmtime > catalogmtime: catalogmtime = tmpmtime

        if mapEmtime < catalogmtime or mapEmtime < os.path.getmtime(__file__):
            use_stored = False
            if logger is not None:
                logger.info(
                    "Updating out-of-date map_E file using newer "+str(file_prefixes)+" files "+
                    "from "+mapper.full_dir)
    # Then load / build as required 
    if use_stored is True:
        if logger is not None:
            logger.info("Loading truth map_E from "+mapEtruefile)
        data = np.loadtxt(mapEtruefile)
        field, theta, map_E, map_B, maperr = (
            data[:, 0].astype(int), data[:, 1], data[:, 2], data[:, 3], data[:, 4])
    else:
        # Define the field array, then theta and map arrays in which we'll store the results
        field = np.arange(NBINS_THETA * NFIELDS) / NBINS_THETA
        theta = np.empty(NBINS_THETA * NFIELDS)
        map_E = np.empty(NBINS_THETA * NFIELDS)
        map_B = np.empty(NBINS_THETA * NFIELDS)
        maperr = np.empty(NBINS_THETA * NFIELDS)
        # Load the offsets
        subfield_indices, offset_deg_x, offset_deg_y = get_generate_variable_offsets(
            experiment, obs_type, storage_dir=storage_dir, truth_dir=truth_dir, logger=logger)
        # Setup some storage arrays into which we'll write
        xfield = np.empty((NGALS_PER_SUBFIELD, NSUBFIELDS_PER_FIELD)) 
        yfield = np.empty((NGALS_PER_SUBFIELD, NSUBFIELDS_PER_FIELD)) 
        # Loop over fields
        import pyfits
        for ifield in range(NFIELDS):

            # Read in all the shears in this field and store
            g1 = np.zeros((NGALS_PER_SUBFIELD, NSUBFIELDS_PER_FIELD))
            g2 = np.zeros((NGALS_PER_SUBFIELD, NSUBFIELDS_PER_FIELD)) 
            for jsub in range(NSUBFIELDS_PER_FIELD):

                # Build the x,y grid using the subfield offsets
                isubfield_index = jsub + ifield * NSUBFIELDS_PER_FIELD
                xfield[:, jsub] = xgrid_deg + offset_deg_x[isubfield_index]
                yfield[:, jsub] = ygrid_deg + offset_deg_y[isubfield_index]
                # If requested (by setting output_xy_prefix) then write these xy out for diagnostic
                if output_xy_prefix is not None:
                    output_xy_filename = output_xy_prefix+("-sub%03d" % isubfield_index)+".asc"
                    print "Writing "+output_xy_filename+" as requested..."
                    with open(output_xy_filename, 'wb') as fout:
                        fout.write("# x  y\n")
                        np.savetxt(fout, np.array((xfield[:, jsub], yfield[:, jsub])).T)
                # Then loop over the supplied file_prefixes and g1/g2 suffixes, summing shears
                for prefix, suffix in zip(file_prefixes, suffixes):

                    galcatfile = os.path.join(
                        mapper.full_dir, (prefix+"-%03d.fits" % isubfield_index))
                    truedata = pyfits.getdata(galcatfile)
                    if len(truedata) != NGALS_PER_SUBFIELD:
                        raise ValueError(
                            "Number of records in "+galcatfile+" (="+str(len(truedata))+") is not "+
                            "equal to NGALS_PER_SUBFIELD (="+str(NGALS_PER_SUBFIELD)+")")
                    # Use the correct rule for shear addition, best (most safely) evaluated using
                    # arrays of complex numbers, see Schneider 2006 eq 12
                    gtoaddc = truedata["g1"+suffix] + truedata["g2"+suffix]*1j 
                    gpriorc = g1[:, jsub] + g2[:, jsub]*1j
                    gfinalc = (gpriorc + gtoaddc) / (1. + gtoaddc.conj() * gpriorc)
                    g1[:, jsub] = gfinalc.real
                    g2[:, jsub] = gfinalc.imag 

            # If requested (by setting output_xy_prefix) then write these xy out for diagnostic
            if output_xy_prefix is not None:
                output_xy_filename = output_xy_prefix+("-%03d" % ifield)+".asc"
                with open(output_xy_filename, 'wb') as fout:
                    fout.write("# x  y\n")
                    np.savetxt(fout, np.array((xfield.flatten(), yfield.flatten())).T)

            # Having got the x,y and g1, g2 for all the subfields in this field, flatten and use
            # to calculate the map_E
            map_results = run_corr2(
                xfield.flatten(), yfield.flatten(), g1.flatten(), g2.flatten(),
                np.ones(NGALS_PER_SUBFIELD * NSUBFIELDS_PER_FIELD), min_sep=THETA_MIN_DEG,
                max_sep=THETA_MAX_DEG, nbins=NBINS_THETA, corr2_executable=corr2_exec,
                xy_units="degrees", sep_units="degrees")
            theta[ifield * NBINS_THETA: (ifield + 1) * NBINS_THETA] = map_results[:, 0] 
            map_E[ifield * NBINS_THETA: (ifield + 1) * NBINS_THETA] = map_results[:, 1]     
            map_B[ifield * NBINS_THETA: (ifield + 1) * NBINS_THETA] = map_results[:, 2]
            maperr[ifield * NBINS_THETA: (ifield + 1) * NBINS_THETA] = map_results[:, 5]

        # Save these in ASCII format
        if logger is not None:
            logger.info("Saving truth map_E file to "+mapEtruefile)
        with open(mapEtruefile, "wb") as fout:
            fout.write("# True aperture mass statistics for "+experiment+"-"+obs_type+"-variable\n")
            fout.write("# field_index  theta [deg]  map_E  map_B  maperr\n")
            np.savetxt(
                fout, np.array((field, theta, map_E, map_B, maperr)).T,
                fmt=" %2d %.18e %.18e %.18e %.18e")
    if make_plots and not use_stored: # No point plotting if already built!
        import matplotlib.pyplot as plt
        plt.figure(figsize=(10, 8))
        plt.subplot(211)
        for ifield in range(NFIELDS):
            plt.semilogx(
                theta[ifield * NBINS_THETA: (ifield + 1) * NBINS_THETA],
                map_E[ifield * NBINS_THETA: (ifield + 1) * NBINS_THETA], label="Field "+str(ifield))
        plt.ylim(-2.e-5, 2.e-5)
        plt.title(mapEtruefile+" E-mode")
        plt.ylabel("Ap. Mass Dispersion")
        plt.axhline(ls="--", color="k")
        plt.legend()
        plt.subplot(212)
        for ifield in range(NFIELDS):
            plt.semilogx(
                theta[ifield * NBINS_THETA: (ifield + 1) * NBINS_THETA],
                map_B[ifield * NBINS_THETA: (ifield + 1) * NBINS_THETA], label="Field "+str(ifield))
        plt.ylim(-2.e-5, 2.e-5)
        plt.title(mapEtruefile+" B-mode")
        plt.xlabel("Theta [degrees]")
        plt.ylabel("Ap. Mass Dispersion")
        plt.axhline(ls="--", color="k")
        plt.legend()
        plotfile = mapEtruefile.rstrip("asc")+"png"
        if logger is not None:
            logger.info("Saving plot output to "+plotfile)
        plt.savefig(plotfile)
    # Then return
    return field, theta, map_E, map_B, maperr

def q_constant(submission_file, experiment, obs_type, storage_dir=STORAGE_DIR, truth_dir=TRUTH_DIR,
               logger=None, normalization=None, sigma2_min=None, just_q=False, cfid=CFID, mfid=MFID,
               pretty_print=False, flip_g1=False, flip_g2=False, plot=False, ignore_fields=None):
    """Calculate the Q_c for a constant shear branch submission.

    @param submission_file  File containing the user submission.
    @param experiment       Experiment for this branch, one of 'control', 'real_galaxy',
                            'variable_psf', 'multiepoch', 'full'
    @param obs_type         Observation type for this branch, one of 'ground' or 'space'
    @param storage_dir      Directory from/into which to load/store rotation files
    @param truth_dir        Root directory in which the truth information for the challenge is
                            stored
    @param logger           Python logging.Logger instance, for message logging
    @param normalization    Normalization factor for the metric (default `None` uses either
                            `NORMALIZATION_CONSTANT_GROUND` or `NORMALIZATION_CONSTANT_SPACE`
                            depending on `obs_type`)
    @param sigma2_min       Damping term to put into the denominator of QZ1 metric (default `None`
                            uses either `SIGMA2_MIN_CONSTANT_GROUND` or `SIGMA2_MIN_CONSTANT_SPACE`
                            depending on `obs_type`)
    @param just_q           Set `just_q = True` (default is `False`) to only return Q_c rather than
                            the default behaviour of returning a tuple including best fitting c+,
                            m+, cx, mx, etc.
    @param cfid             Fiducial, target c value
    @param mfid             Fiducial, target m value
    @param ignore_fields    List or tuple of fields to ignore.  If None, use all fields.
    @return The metric Q_c, & optionally best fitting c+, m+, cx, mx, sigc+, sigcm+, sigcx, sigmx.
    """
    if not os.path.isfile(submission_file):
        raise ValueError("Supplied submission_file '"+submission_file+"' does not exist.")
    # If the sigma2_min is not changed from None, set using defaults based on obs_type
    if sigma2_min is None:
        if obs_type == "ground":
            sigma2_min = SIGMA2_MIN_CONSTANT_GROUND
        elif obs_type == "space":
            sigma2_min = SIGMA2_MIN_CONSTANT_SPACE
        else:
            raise ValueError("Default sigma2_min cannot be set as obs_type not recognised")
    # If the normalization is not changed from None, set using defaults based on obs_type
    if normalization is None:
        if obs_type == "ground":
            normalization = NORMALIZATION_CONSTANT_GROUND
        elif obs_type == "space":
            normalization = NORMALIZATION_CONSTANT_SPACE
        else:
            raise ValueError("Default sigma2_min cannot be set as obs_type not recognised")
    # Load the submission and label the slices we're interested in
    if logger is not None:
        logger.info("Calculating Q_c metric for "+submission_file)
    data = np.loadtxt(submission_file)
    subfield = data[:, 0]  
    g1sub = data[:, 1]
    g2sub = data[:, 2]
    if flip_g1: g1sub = -g1sub
    if flip_g2: g2sub = -g2sub
    # Load up the rotations, then rotate g1 & g2 in the correct sense.
    # NOTE THE MINUS SIGNS!  This is because we need to rotate the coordinates *back* into a frame
    # in which the primary direction of the PSF is g1, and the orthogonal is g2
    try: # Put this in a try except block to handle funky submissions better
        rotations = get_generate_const_rotations(
            experiment, obs_type, truth_dir=truth_dir, storage_dir=storage_dir, logger=logger)
        g1srot = g1sub * np.cos(-2. * rotations) - g2sub * np.sin(-2. * rotations)
        g2srot = g1sub * np.sin(-2. * rotations) + g2sub * np.cos(-2. * rotations)
        # Load the truth
        _, g1truth, g2truth = get_generate_const_truth(
             experiment, obs_type, truth_dir=truth_dir, storage_dir=storage_dir, logger=logger)
        # Rotate the truth in the same sense, then use the g3metrics.fitline routine to
        # perform simple linear regression
        g1trot = g1truth * np.cos(-2. * rotations) - g2truth * np.sin(-2. * rotations)
        g2trot = g1truth * np.sin(-2. * rotations) + g2truth * np.cos(-2. * rotations)
        # Decide which subfields to use / ignore
        use = np.ones_like(g1srot, dtype=bool)
        if ignore_fields is not None:
            for field_index in ignore_fields:
                use[field_index] = False
        Q_c, c1, m1, c2, m2, sigc1, sigm1, sigc2, sigm2 = g3metrics.metricQZ1_const_shear(
            g1srot[use], g2srot[use], g1trot[use], g2trot[use],
            cfid=cfid, mfid=mfid, sigma2_min=sigma2_min)
        Q_c *= normalization
        if plot:
            import matplotlib.pyplot as plt
            plt.subplot(211)
            plt.plot(g1trot, g1srot - g1trot, 'k+')
            plt.plot(
                [min(g1trot), max(g1trot)], [m1 * min(g1trot) + c1, m1 * max(g1trot) + c1],
                'b-', label="c+  = %+.5f +/- %.5f \nm+  = %+.5f +/- %.5f" % (c1, sigc1, m1, sigm1))
            plt.xlim()
            plt.xlabel("gtrue_+")
            plt.ylabel("(gsub - gtrue)_+")
            plt.ylim(-0.015, 0.015)
            plt.title(os.path.split(submission_file)[-1])
            plt.axhline(ls='--', color='k')
            plt.legend()
            plt.subplot(212)
            plt.plot(g2trot, g2srot - g2trot, 'k+')
            plt.plot(
                [min(g2trot), max(g2trot)], [m2 * min(g2trot) + c2, m2 * max(g2trot) + c2],
                'r-', label="cx  = %+.5f +/- %.5f \nmx  = %+.5f +/- %.5f" % (c2, sigc2, m2, sigm2))
            plt.xlabel("gtrue_x")
            plt.ylabel("(gsub - gtrue)_x")
            plt.ylim(-0.015, 0.015)
            plt.axhline(ls='--', color='k')
            plt.legend()
            if type(plot) == str:
                print "Saving plot to "+plot
                plt.savefig(plot)
            plt.show()

    except Exception as err:
        # Something went wrong... We'll handle this silently setting all outputs to zero but warn
        # the user via any supplied logger; else raise
        Q_c, c1, m1, c2, m2, sigc1, sigm1, sigc2, sigm2 = 0, 0, 0, 0, 0, 0, 0, 0, 0
        print err
        if logger is not None:
            logger.warn(err.message)
        else:
            raise err # ...Raise exception if there is no logger
    # Then return
    if just_q:
        ret = Q_c
    else:
        if pretty_print:
            print
            print "Evaluated results for submission "+str(submission_file)
            print "Using sigma2_min = "+str(sigma2_min)
            print
            print "Q_c =  %.4f" % Q_c
            print "c+  = %+.5f +/- %.5f" % (c1, sigc1)
            print "cx  = %+.5f +/- %.5f" % (c2, sigc2)
            print "m+  = %+.5f +/- %.5f" % (m1, sigm1)
            print "mx  = %+.5f +/- %.5f" % (m2, sigm2)
            print
        ret = (Q_c, c1, m1, c2, m2, sigc1, sigm1, sigc2, sigm2)
    return ret

def q_variable(submission_file, experiment, obs_type, normalization=None, truth_dir=TRUTH_DIR,
               storage_dir=STORAGE_DIR, logger=None, corr2_exec="corr2", poisson_weight=False,
               usebins=USEBINS, fractional_diff=False, squared_diff=False, sigma2_min=None):
    """Calculate the Q_v for a variable shear branch submission.

    @param submission_file  File containing the user submission.
    @param experiment       Experiment for this branch, one of 'control', 'real_galaxy',
                            'variable_psf', 'multiepoch', 'full'
    @param obs_type         Observation type for this branch, one of 'ground' or 'space'
    @param normalization    Normalization factor for the metric, default will be set differently for
                            obs_type='space' and obs_type='ground'
    @param storage_dir      Directory from/into which to load/store rotation files
    @param truth_dir        Root directory in which the truth information for the challenge is
                            stored
    @param logger           Python logging.Logger instance, for message logging
    @param corr2_exec       Path to Mike Jarvis' corr2 exectuable
    @param poisson_weight   If `True`, use the relative Poisson errors in each bin of map_E
                            to form an inverse variance weight for the difference metric
                            [default = `False`]
    @param usebins          An array the same shape as EXPECTED_THETA specifying which bins to
                            use in the calculation of Q_v [default = `USEBINS`].  If set to `None`,
                            uses all bins
    @param fractional_diff  Use a |fractional|, rather than absolute difference in metric
    @param squared_diff     Use the squared, rather than the absolute difference in metric
    @param sigma2_min       Damping term to put into the denominator of metric (default `None`
                            uses either `SIGMA2_MIN_VARIABLE_GROUND` or `SIGMA2_MIN_VARIABLE_SPACE`
                            depending on `obs_type`)
    @return The metric Q_v
    """
    if not os.path.isfile(submission_file):
        raise ValueError("Supplied submission_file '"+submission_file+"' does not exist.")
    # Set the default normalization based on whether ground or space data
    if normalization is None:
        if obs_type == "ground":
            normalization = NORMALIZATION_VARIABLE_GROUND
        elif obs_type == "space":
            normalization = NORMALIZATION_VARIABLE_SPACE
        else:
            raise ValueError("Default normalization cannot be set as obs_type not recognised")
    # If the sigma2_min is not changed from `None`, set using defaults based on `obs_type`
    if sigma2_min is None:
        if obs_type == "ground":
            sigma2_min = SIGMA2_MIN_VARIABLE_GROUND
        elif obs_type == "space":
            sigma2_min = SIGMA2_MIN_VARIABLE_SPACE
        else:
            raise ValueError("Default sigma2_min cannot be set as obs_type not recognised")
    # Load the submission and label the slices we're interested in
    if logger is not None:
        logger.info("Calculating Q_v metric for "+submission_file)
    data = np.loadtxt(submission_file)
    # We are stating that we want at least 4 and up to 5 columns, so check for this
    if data.shape not in ((NBINS_THETA * NFIELDS, 4), (NBINS_THETA * NFIELDS, 5)):
        raise ValueError("Submission "+str(submission_file)+" is not the correct shape!")
    # Extract the salient parts of the submission from data
    field_sub = data[:, 0].astype(int)
    theta_sub = data[:, 1]
    map_E_sub = data[:, 2]
    # Load/generate the truth shear signal
    field_shear, theta_shear, map_E_shear, _, maperr_shear = get_generate_variable_truth(
        experiment, obs_type, truth_dir=truth_dir, storage_dir=storage_dir, logger=logger,
        corr2_exec=corr2_exec, mape_file_prefix=MAPESHEAR_FILE_PREFIX, suffixes=("",),
        make_plots=False)
    # Then generate the intrinsic only map_E, useful for examinging plots, including the maperr
    # (a good estimate of the relative Poisson errors per bin) which we will use to provide a weight
    field_int, theta_int, map_E_int, _, maperr_int = get_generate_variable_truth(
        experiment, obs_type, truth_dir=truth_dir, storage_dir=storage_dir, logger=logger,
        corr2_exec=corr2_exec, mape_file_prefix=MAPEINT_FILE_PREFIX, suffixes=("_intrinsic",),
        make_plots=False)
    # Then generate the theory observed = int + shear combined map signals - these are our reference
    # Note this uses the new functionality of get_generate_variable_truth for adding shears
    field_ref, theta_ref, map_E_ref, _, maperr_ref = get_generate_variable_truth(
        experiment, obs_type, truth_dir=truth_dir, storage_dir=storage_dir, logger=logger,
        corr2_exec=corr2_exec, mape_file_prefix=MAPEOBS_FILE_PREFIX,
        file_prefixes=("galaxy_catalog", "galaxy_catalog"), suffixes=("_intrinsic", ""),
        make_plots=False)
    # Set up the weight
    if poisson_weight:
        weight = max(maperr_int**2) / maperr_int**2 # Inverse variance weight
    else:
        weight = np.ones_like(map_E_ref)
    # Set up the usebins to use if `usebins == None` (use all bins)
    if usebins is None:
        usebins = np.repeat(True, NBINS_THETA * NFIELDS)
    # Get the total number of active bins per field
    nactive = sum(usebins) / NFIELDS
    try: # Put this in a try except block to handle funky submissions better
        np.testing.assert_array_almost_equal( # Sanity check out truth / expected theta bins
            theta_shear, EXPECTED_THETA, decimal=3,
            err_msg="BIG SNAFU! Truth theta does not match the EXPECTED_THETA, failing...")
        np.testing.assert_array_equal(
            field_sub, field_ref, err_msg="User field array does not match truth.")
        np.testing.assert_array_almost_equal(
            theta_sub, theta_ref, decimal=3, err_msg="User theta array does not match truth.")
        # The definition of Q_v is so simple there is no need to use the g3metrics version
        # NOTE WE ARE TRYING A NEW DEFINITION OF Q_v THAT IS NOT SO SIMPLE
        Q_v_fields = np.zeros(nactive) # To store diffs averaged over fields, per bin
        if not fractional_diff:
            for i in range(nactive): # Sum over all fields for each bin, nactive being the stride

                Q_v_fields[i] = np.sum(((weight * (map_E_sub - map_E_ref))[usebins])[i::nactive])

        else:
            for i in range(nactive): # Sum over all fields for each bin, nactive being the stride

                Q_v_fields[i] = np.sum(
                    ((weight * (map_E_sub - map_E_ref) / map_E_ref)[usebins])[i::nactive])

        # Then take the weighted average abs(Q_v_fields)
        if not squared_diff:
            Q_v = normalization / (
                sigma2_min + (np.sum(np.abs(Q_v_fields)) / np.sum(weight[usebins])))
        else:
            Q_v = normalization / (
                sigma2_min + (np.sum(Q_v_fields**2) / np.sum(weight[usebins])))
    except Exception as err:
        Q_v = 0. # If the theta or field do not match, let's be strict and force Q_v...
        if logger is not None:
            logger.warn(err.message) # ...But let's warn if there's a logger!
        else:                        # ...And raise the exception if not
            raise err
    # Then return Q_v
    return Q_v

def map_diff_func(cm_array, mapEsub, maperrsub, mapEref, mapEunitc):
        """Difference between an m,c model of a biased aperture mass statistic submission and the 
        submission itself, as a vector corresponding to the theta vector.

        The model of the biased submission is simply:

            mapEmodel = mapEunitc * c^2 + mapEref * (1 + 2 * m + m^2)

        where c, m = cm_array[0], cm_array[1].  This code returns

            (mapEmodel - mapEsub) / maperrsub

        for the use of scipy.optimize.leastsq() within q_variable_by_mc().
        """
        ret = (
            mapEunitc * cm_array[0]**2 + mapEref * (1. + 2. * cm_array[1] + cm_array[1]**2)
            - mapEsub) / maperrsub
        return ret

def q_variable_by_mc(submission_file, experiment, obs_type, map_E_unitc, normalization=None,
                     truth_dir=TRUTH_DIR, storage_dir=STORAGE_DIR, logger=None, usebins=None,
                     corr2_exec="corr2", sigma2_min=None, cfid=CFID, mfid=MFID, just_q=False,
                     pretty_print=False):
    """Calculate the Q_v for a variable shear branch submission, using a best-fitting m and c model
    of submission biases to evaluate the score.  Experimental metric, not used in the GREAT3
    challenge due to the difficulty of reliably modelling m & c in simulation tests.

    @param submission_file  File containing the user submission.
    @param experiment       Experiment for this branch, one of 'control', 'real_galaxy',
                            'variable_psf', 'multiepoch', 'full'
    @param obs_type         Observation type for this branch, one of 'ground' or 'space'
    @param normalization    Normalization factor for the metric, default will be set differently for
                            obs_type='space' and obs_type='ground'
    @param storage_dir      Directory from/into which to load/store rotation files
    @param truth_dir        Root directory in which the truth information for the challenge is
                            stored
    @param logger           Python logging.Logger instance, for message logging
    @param usebins          An array the same shape as EXPECTED_THETA specifying which bins to
                            use in the calculation of Q_v [default = `USEBINS`].  If set to `None`,
                            uses all bins
    @param corr2_exec       Path to Mike Jarvis' corr2 exectuable
    @param sigma2_min       Damping term to put into the denominator of metric (default `None`
                            uses either `SIGMA2_MIN_VARIABLE_GROUND` or `SIGMA2_MIN_VARIABLE_SPACE`
                            depending on `obs_type`)
    @param cfid             Fiducial, target c value
    @param mfid             Fiducial, target m value
    @param just_q           Set `just_q = True` (default is `False`) to only return Q_v rather than
                            the default behaviour of returning a tuple including best fitting |c|,
                            m, uncertainties etc.
    @return The metric Q_v
    """
    if not os.path.isfile(submission_file):
        raise ValueError("Supplied submission_file '"+submission_file+"' does not exist.")
    # Set the default normalization based on whether ground or space data
    if normalization is None:
        if obs_type == "ground":
            normalization = NORMALIZATION_CONSTANT_GROUND
        elif obs_type == "space":
            normalization = NORMALIZATION_CONSTANT_SPACE
        else:
            raise ValueError("Default normalization cannot be set as obs_type not recognised")
    # If the sigma2_min is not changed from `None`, set using defaults based on `obs_type`
    if sigma2_min is None:
        if obs_type == "ground":
            sigma2_min = SIGMA2_MIN_VARIABLE_GROUND
        elif obs_type == "space":
            sigma2_min = SIGMA2_MIN_VARIABLE_SPACE
        else:
            raise ValueError("Default sigma2_min cannot be set as obs_type not recognised")
    # Load the submission and label the slices we're interested in
    if logger is not None:
        logger.info("Calculating Q_v metric (by m & c) for "+submission_file)
    data = np.loadtxt(submission_file)
    # We are stating that we want at least 4 and up to 5 columns, so check for this
    if data.shape not in ((NBINS_THETA * NFIELDS, 4), (NBINS_THETA * NFIELDS, 5)):
        raise ValueError("Submission "+str(submission_file)+" is not the correct shape!")
    # Extract the salient parts of the submission from data
    field_sub = data[:, 0].astype(int)
    theta_sub = data[:, 1]
    map_E_sub = data[:, 2]
    # Load/generate the truth shear signal
    field_shear, theta_shear, map_E_shear, _, maperr_shear = get_generate_variable_truth(
        experiment, obs_type, truth_dir=truth_dir, storage_dir=storage_dir, logger=logger,
        corr2_exec=corr2_exec, mape_file_prefix=MAPESHEAR_FILE_PREFIX, suffixes=("",),
        make_plots=False)
    # Then generate the intrinsic only map_E, useful for examinging plots, including the maperr
    # (a good estimate of the relative Poisson errors per bin) which we will use to provide a weight
    field_int, theta_int, map_E_int, _, maperr_int = get_generate_variable_truth(
        experiment, obs_type, truth_dir=truth_dir, storage_dir=storage_dir, logger=logger,
        corr2_exec=corr2_exec, mape_file_prefix=MAPEINT_FILE_PREFIX, suffixes=("_intrinsic",),
        make_plots=False)
    # Then generate the theory observed = int + shear combined map signals - these are our reference
    # Note this uses the new functionality of get_generate_variable_truth for adding shears
    field_ref, theta_ref, map_E_ref, _, maperr_ref = get_generate_variable_truth(
        experiment, obs_type, truth_dir=truth_dir, storage_dir=storage_dir, logger=logger,
        corr2_exec=corr2_exec, mape_file_prefix=MAPEOBS_FILE_PREFIX,
        file_prefixes=("galaxy_catalog", "galaxy_catalog"), suffixes=("_intrinsic", ""),
        make_plots=False)
    # Set up the usebins to use if `usebins == None` (use all bins)
    if usebins is None:
        usebins = np.repeat(True, NBINS_THETA * NFIELDS)
    # Get the total number of active bins per field
    nactive = sum(usebins) / NFIELDS
    if True: # Put this in a try except block to handle funky submissions better
        np.testing.assert_array_almost_equal( # Sanity check out truth / expected theta bins
            theta_shear, EXPECTED_THETA, decimal=3,
            err_msg="BIG SNAFU! Truth theta does not match the EXPECTED_THETA, failing...")
        np.testing.assert_array_equal(
            field_sub, field_ref, err_msg="User field array does not match truth.")
        np.testing.assert_array_almost_equal(
            theta_sub, theta_ref, decimal=3, err_msg="User theta array does not match truth.")
        # Use optimize.leastsq to find the best fitting linear bias model params, and covariances
        import scipy.optimize
        optimize_results = scipy.optimize.leastsq(
            map_diff_func, np.array([0., 0.]),
            args=(
                map_E_sub[usebins],
                maperr_ref[usebins],
                map_E_ref[usebins],  # Note use of ref errors: this will appropriately
                                              # weight different bins and is not itself noisy
                map_E_unitc[usebins]), full_output=True)
        csub = optimize_results[0][0]
        msub = optimize_results[0][1]
        map_E_model = map_E_unitc * csub**2 + map_E_ref * (1. + 2. * msub + msub**2)
        residual_variance = np.var(
            ((map_E_sub - map_E_model) / maperr_ref)[usebins], ddof=1)
        if optimize_results[1] is not None:
            covcm = optimize_results[1] * residual_variance
            sigcsub = np.sqrt(covcm[0, 0])
            sigmsub = np.sqrt(covcm[1, 1])
            covcm = covcm[0, 1]
        else:
            sigcsub = 0.
            sigmsub = 0.
            covcm = 0.
        # Then we define the Q_v
        Q_v = 2449. * normalization / np.sqrt(
            (csub / cfid)**2 + (msub / mfid)**2 + sigma2_min)
    try:
        pass
    except Exception as err:
        Q_v = 0. # If the theta or field do not match, let's be strict and force Q_v...
        if logger is not None:
            logger.warn(err.message) # ...But let's warn if there's a logger!
        else:                        # ...And raise the exception if not
            raise err
    # Then return
    # Then return
    if just_q:
        ret = Q_v
    else:
        if pretty_print:
            print "Evaluated results for submission "+str(submission_file)
            print "Using sigma2_min = "+str(sigma2_min)
            print "Q_v =  %.4f" % Q_v
            print "|c| = %+.5f +/- %.5f" % (csub, sigcsub)
            print " m  = %+.5f +/- %.5f" % (msub, sigmsub)
            print "Cov(0, 0) = %+.4e" % sigcsub**2
            print "Cov(0, 1) = %+.4e" % covcm
            print "Cov(1, 1) = %+.4e" % sigmsub**2
        ret = (Q_v, csub, msub, sigcsub, sigmsub, covcm)
    return ret
