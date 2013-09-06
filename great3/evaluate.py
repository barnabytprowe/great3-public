"""@file evaluate.py

Module containing functions for the evaluation of the constant and variable shear-type branches
of the GREAT3 Challenge.

Constant shear branches
-----------------------
Each submission file (one per branch) should take the format of a 3-column ASCII catalog, e.g.:

    # SUBFIELD  G1  G2
    0   -.25424  0.23111
    1   -.05111  0.37123
    ...

or similar.  The hashed header/comment can be omitted, and almost all formats for the numbers are
fine.  The main criterion to be satisfied is that after

    >>> data = np.loadtxt(submission_file)

the `data` object must be a NumPy array with shape `(NSUBFIELDS, 3)`, where `NSUBFIELDS` is the
total number of subfields in the branch (currently fixed at 200).

In addition, the array slice `data[:, 0]` must be the subfield number in ascending order from `0` to
`NSUBFIELDS - 1`.  The array slices `data[:, 1]` and `data[:, 2]` must be the corresponding
estimates of mean shear g1 and g2 in each subfield.
 
Variable shear branches
-----------------------
Each submission file (one per branch?) CHECK THESE DETAILS WITH MELANIE'S CODE OUTPUT, TODO.
"""

import os
import sys
import logging
import numpy as np
import great3sims
import great3sims.mapper
sys.path.append(sys.path.pop(0)) # Return the path back to normal (i.e. with ./ as the first entry)
try:
    import g3metrics
except ImportError:
    path, module = os.path.split(__file__)
    sys.path.append(os.path.join(path, "..", "..", "metrics")) # Appends the great3-private/metrics
                                                               # folder to path
    import g3metrics

NFIELDS = 10     # Total number of fields
NSUBFIELDS = 200 # Total number of subfields, not necessarily equal to the number of subfields made
                 # in mass_produce as that script also generates the deep fields
NSUBFIELDSPERFIELD = NSUBFIELDS / NFIELDS

CFID = 2.e-4
MFID = 2.e-3

TRUTH_SUBFIELD_DICT = {} # A dictionary containing the mapping between subfields containing the
                         # same applied shear [one of NFIELDS pairs of independent (g1, g2) values]

STORAGE_DIR = "./metric_calculation_products" # Folder into which to store useful intermediate
                                              # outputs of metric calculations (e.g. rotation files,
                                              # dicts, mapE tables) which need be calculated only
                                              # once
TRUTH_DIR = "/Users/browe/great3/truth"       # Root folder in which the truth values are unpacked
                                              # (admin only)

SUBFIELD_DICT_FILE_PREFIX = "subfield_dict_"
GTRUTH_FILE_PREFIX = "gtruth_"
ROTATIONS_FILE_PREFIX = "rotations_"



def get_generate_const_truth(experiment, obs_type, truth_dir=TRUTH_DIR, storage_dir=STORAGE_DIR,
                             logger=None):
    """Get or generate arrays of subfield_index, g1true, g2true, each of length `NSUBFIELDS`.

    If the gtruth file has already been built for this constant shear branch, loads and
    returns the saved copies.

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
    gtruefile = os.path.join(storage_dir, GTRUTH_FILE_PREFIX+experiment[0]+obs_type[0]+".asc")
    mapper = great3sims.mapper.Mapper(truth_dir, experiment, obs_type, 'constant')
    use_stored = True
    if not os.path.isfile(gtruefile):
        use_stored = False
        if logger:
            logger.info(
                "First build of shear truth tables using values from "+
                os.path.join(mapper.full_dir, "shear_params-*.yaml"))
    else:
        # Compare timestamps for the gtruefile and the first shear_params file
        # (subfield = 000) for this branch.  If the former is older than the latter, or this file,
        # force rebuild...
        gtruemtime = os.path.getmtime(gtruefile)
        shear_params_file0 = os.path.join(mapper.full_dir, "shear_params-000.yaml")
        shear_params_mtime = os.path.getmtime(shear_params_file0)
        if gtruemtime < shear_params_mtime or gtruemtime < os.path.getmtime(__file__):
            use_stored = False
            if logger:
                logger.info(
                    "Updating out-of-date shear truth tables using newer values from "+
                    os.path.join(mapper.full_dir, "shear_params-*.yaml"))
    # Then load or build (and save) the array of truth values per subfield
    if use_stored:
        if logger:
            logger.info("Loading shear truth tables from "+gtruefile) 
        gtruedata = np.loadtxt(gtruefile)
    else: 
        params_prefix = os.path.join(mapper.full_dir, "shear_params-") 
        gtruedata = np.empty((NSUBFIELDS, 3))
        import yaml
        gtruedata[:, 0] = np.arange(NSUBFIELDS)
        for i in range(NSUBFIELDS):

            params_file = params_prefix+("%03d" % i)+".yaml"
            with open(params_file, 'rb') as funit:
                gdict = yaml.load(funit)
                gtruedata[i, 1] = gdict['g1']
                gtruedata[i, 2] = gdict['g2']

        if logger:
            logger.info("Saving shear truth table to "+gtruefile)
        if not os.path.isdir(storage_dir):
            os.mkdir(storage_dir)
        with open(gtruefile, 'wb') as fout:
            fout.write("#  True shears for "+experiment+"-"+obs_type+"-constant\n")
            fout.write("#  subfield_index  g1true  g2true\n")
            np.savetxt(fout, gtruedata, fmt=" %4d %+.18e %+.18e")
    return (gtruedata[:, 0]).astype(int), gtruedata[:, 1], gtruedata[:, 2]

def get_generate_const_subfield_dict(experiment, obs_type, storage_dir=STORAGE_DIR,
                                     truth_dir=TRUTH_DIR, logger=None):
    """Get or generate a dict mapping which subfields contain which of NFIELDS independent truth
    shear values.

    If the subfield_dict file has already been built for this constant shear branch, loads and
    returns the saved copy.

    If the subfield_dict has not been built, or is older than the first entry in the set of
    shear_params files, the subfield_dict is built first, saved to file, then returned.

    @param experiment     Experiment for this branch, one of 'control', 'real_galaxy',
                          'variable_psf', 'multiepoch', 'full'
    @param obs_type       Observation type for this branch, one of 'ground' or 'space'
    @param storage_dir    Directory from/into which to load/store rotation files
    @param truth_dir      Root directory in which the truth information for the challenge is stored
    @param logger         Python logging.Logger instance, for message logging
    @return               The subfield_dict (see code below for details)
    """
    import cPickle

    subfield_dict_file = os.path.join(
        storage_dir, SUBFIELD_DICT_FILE_PREFIX+experiment[0]+obs_type[0]+".pkl")
    mapper = great3sims.mapper.Mapper(truth_dir, experiment, obs_type, 'constant')
    use_stored = True
    if not os.path.isfile(subfield_dict_file):
        use_stored = False
        if logger:
            logger.info(
                "First build of shear-subfield dictionary using values from "+
                os.path.join(mapper.full_dir, "shear_params-*.yaml"))
    else:
        # Compare timestamps for the subfield_dict file and the first shear_params file
        # (subfield = 000) for this branch.  If the former is older than the latter, or this file, 
        # force rebuild...
        dictmtime = os.path.getmtime(subfield_dict_file)
        shear_params_file0 = os.path.join(mapper.full_dir, "shear_params-000.yaml")
        shear_params_mtime = os.path.getmtime(shear_params_file0)
        if dictmtime < shear_params_mtime or dictmtime < os.path.getmtime(__file__):
            use_stored = False 
            if logger:
                logger.info(
                    "Updating out-of-date shear-subfield dictionary using newer values from "+
                    os.path.join(mapper.full_dir, "shear_params-*.yaml"))
    # Then load or build (and save) the subfield_dict
    if use_stored is True:
        if logger:
            logger.info("Loading shear-subfield dictionary from "+subfield_dict_file)  
        with open(subfield_dict_file, 'rb') as funit:
            subfield_dict = cPickle.load(funit)
    else:
        _, g1true, g2true = get_generate_const_truth(experiment, obs_type, logger=logger)
        # First of all get all the unique (g1, g2) values, in the order in which they first
        # appear in the arrays (so as not to mess up pairs) as suggested in
        # http://stackoverflow.com/questions/12926898/numpy-unique-without-sort
        g1unique, g1unique_indices = np.unique(g1true, return_index=True) 
        g2unique, g2unique_indices = np.unique(g2true, return_index=True)
        # Put back into first-found order
        g1unique_unsorted = [g1true[index] for index in sorted(g1unique_indices)]
        g2unique_unsorted = [g2true[index] for index in sorted(g2unique_indices)]
        # Sanity check
        if len(g1unique_unsorted) != NFIELDS or len(g2unique_unsorted) != NFIELDS:
            raise ValueError(
                "Number of unique shear values != NFIELDS!\n"+
                "NFIELDS = "+str(NFIELDS)+"; n_unique = "+str(len(g1unique_indices)))
        # Then construct the subfield dict by looping over these unique values and saving the
        # g1 and g2 values to an ordered list, along with a corresponding ordered list of
        # NumPy arrays containing indices of the subfields this shear was applied to 
        g1dict = {"value": g1unique_unsorted, "subfield_indices":[]}
        g2dict = {"value": g2unique_unsorted, "subfield_indices":[]}
        ifullset = np.arange(NSUBFIELDS, dtype=int)
        for g1trueval, g2trueval in zip(g1dict["value"], g2dict["value"]):

            ig1subset = ifullset[g1true == g1trueval]
            ig2subset = ifullset[g2true == g2trueval]
            if tuple(ig1subset) != tuple(ig2subset): # Tuple comparison
                raise ValueError(
                    "Unique values of truth g1 and g2 do not correspond pairwise!")
            if len(ig1subset) != NSUBFIELDSPERFIELD:
                raise ValueError(
                    "Fields of truth shear values do not contain equal numbers of subfields!")
            g1dict["subfield_indices"].append(ig1subset)
            g2dict["subfield_indices"].append(ig2subset)

        # Then put both g1 and g2 into the subfield_dict
        subfield_dict = {"g1": g1dict, "g2": g2dict}
        # Save the resulting subfield_dict
        if logger:
            logger.info("Saving subfield_dict to "+subfield_dict_file)
        if not os.path.isdir(storage_dir):
            os.mkdir(storage_dir) 
        with open(subfield_dict_file, 'wb') as fout:
            cPickle.dump(subfield_dict, fout)
    # Then return
    return subfield_dict

def get_generate_const_rotations(experiment, obs_type, storage_dir=STORAGE_DIR,
                                 truth_dir=TRUTH_DIR, logger=None):
    """Get or generate an array of rotation angles for Q_const calculation.

    If the rotation file has already been built for this constant shear branch, loads and returns an
    array of rotation angles to align with the PSF.  This array is of shape `(NSUBFIELDS, n_epochs)`
    where the number of epochs `n_epochs` is determined from the experiment name using the mapper.

    If the rotation file has not been built, or is older than the first entry in the set of
    starshape_parameters files, the array of rotations is built, saved to file, then returned.

    @param experiment     Experiment for this branch, one of 'control', 'real_galaxy',
                          'variable_psf', 'multiepoch', 'full'
    @param obs_type       Observation type for this branch, one of 'ground' or 'space'
    @param storage_dir    Directory from/into which to load/store rotation files
    @param truth_dir      Root directory in which the truth information for the challenge is stored
    @return               An array containing all the rotation angles, in radians
    """
    rotfile = os.path.join(storage_dir, ROTATIONS_FILE_PREFIX+experiment[0]+obs_type[0]+".asc")
    mapper = great3sims.mapper.Mapper(truth_dir, experiment, obs_type, 'constant')
    use_stored = True
    if not os.path.isfile(rotfile):
        use_stored = False
        if logger:
            logger.info(
                "First build of rotations file using starshape_parameters from "+
                mapper.full_dir)
    else:
        # Then compare timestamps for the rotation file and the first starshape_parameters file
        # (subfield = 000, epoch =0) for this branch.  If the former is older than the latter, or
        # this file, force rebuild...
        rotmtime = os.path.getmtime(rotfile)
        starshape_file_template, _ ,_ = mapper.mappings['starshape_parameters']  
        starshape_file00 = os.path.join(
            mapper.full_dir, starshape_file_template % {"epoch_index": 0, "subfield_index": 0})
        starshapemtime = os.path.getmtime(starshape_file00+".yaml")
        if rotmtime < starshapemtime or rotmtime < os.path.getmtime(__file__):
            use_stored = False
            logger.info(
                    "Updating out-of-date rotations file using newer starshape_parameters from "+
                    mapper.full_dir)
    # Then load / build as required 
    if use_stored is True:
        rotations = np.loadtxt(rotfile)
        if logger:
            logger.info("Loading rotations from "+rotfile) 
    else:
        # To build we must loop over all the subfields and epochs
        # First work out if the experiment is multi-exposure and has multiple epochs
        if experiment in ("multiepoch", "full"):
            import great3sims.constants
            n_epochs = great3sims.constants.n_epochs
        else:
            n_epochs = 1
        # Setup the array for storing the rotation values
        rotations = np.empty((NSUBFIELDS, n_epochs))
        # TODO: Change how this is handled for MULTIEPOCH experiments!
        output_header = "#  Rotations for "+experiment+"-"+obs_type+"-constant\n"+"#  epoch"
        for epoch_index in range(n_epochs):
            output_header+="  "+str(epoch_index)
            for subfield_index in range(NSUBFIELDS):
                starshape_parameters = mapper.read(
                    "starshape_parameters", data_id={
                        "epoch_index": epoch_index, "subfield_index": subfield_index})
                rotations[subfield_index, epoch_index] = .5 * np.arctan2(
                    starshape_parameters['psf_g2'], starshape_parameters['psf_g1'])
        # We have built rotations, but then save this file as ascii for use next time
        if logger:
            logger.info("Saving rotations to "+rotfile)    
        if not os.path.isdir(storage_dir):
            os.mkdir(storage_dir)
        with open(rotfile, 'wb') as fout:
            fout.write(output_header+"\n")
            np.savetxt(fout, rotations, fmt=" %+.18f" * n_epochs)
    if len(rotations.shape) > 1:
        if rotations.shape[1] == 1:
            rotations = rotations.flatten()
    return rotations

def Q_const(submission_file, experiment, obs_type, truth_dir=TRUTH_DIR, storage_dir=STORAGE_DIR,
            logger=None):
    """Calculate the Q_c for a constant shear branch submission.

    @param submission_file  File containing the user submission.
    @param experiment       Experiment for this branch, one of 'control', 'real_galaxy',
                            'variable_psf', 'multiepoch', 'full'
    @param obs_type         Observation type for this branch, one of 'ground' or 'space'
    @param storage_dir      Directory from/into which to load/store rotation files
    @param truth_dir        Root directory in which the truth information for the challenge is
                            stored
    @return The metric Q_const, & best fitting c1, m1, c2, m2, sigc1, sigcm1, sigc2, sigm2
    """
    if not os.path.isfile(submission_file):
        raise ValueError("Supplied submission_file '"+submission_file+"' does not exist.")
    
    # Load the submission and label the slices we're interested in
    data = np.loadtxt(submission_file)
    subfield = data[:, 0]  
    g1sub = data[:, 1]
    g2sub = data[:, 2]
    # Load up the rotations, then rotate g1 & g2 in the correct sense
    # NOTE THE MINUS SIGN!  This is because we need to rotated the coordinates back into a frame
    # in which the primary direction of the PSF is g1, and the orthogonal is g2
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
    Q_c = g3metrics.metricQZ1_const_shear(g1srot, g2srot, g1trot, g2trot, cfid=CFID, mfid=MFID)
    return Q_c 
