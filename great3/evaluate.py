"""@file evaluate.py

Module containing functions for the evaluation of the constant and variable shear-type branches
of the GREAT3 Challenge.

Constant shear branches
-----------------------
Each submission file (one per branch) should take the format of a 3-column ASCII catalog, e.g.:

    # SUBFIELD  G1  G2
    1   -.25424  0.23111
    2   -.05111  0.37123
    ...

or similar.  The hashed header/comment can be omitted, and almost all formats for the numbers are
fine.  The main criterion to be satisfied is that after

    >>> data = np.loadtxt(submission_file)

the `data` object must be a NumPy array with shape `(NSUBFIELDS, 3)`, where `NSUBFIELDS` is the
total number of subfields in the branch (currently fixed at 200).

In addition, the array slice `data[:, 0]` must be the subfield number in ascending order from 1 to 
`NSUBFIELDS`.  The array slices `data[:, 1]` and `data[:, 2]` must be the corresponding estimates
of mean shear g1 and g2 in each subfield.
 
Variable shear branches
-----------------------
Each submission file (one per branch?) CHECK THESE DETAILS WITH MELANIE'S CODE OUTPUT.
"""

import os
import sys
import numpy as np
sys.path.append(os.path.join("..", ".."))
import great3.mapper
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

TRUTH_SUBFIELD_DICT = {} # A dictionary containing the mapping between subfields containing the
                         # same applied shear [one of NFIELDS pairs of independent (g1, g2) values]

STORAGE_DIR = "./metric_calculation_products" # Folder into which to store useful intermediate
                                              # outputs of metric calculations (e.g. rotation files,
                                              # dicts, mapE tables) which need be calculated only
                                              # once
TRUTH_DIR = "/Users/browe/great3/truth"

def get_generate_const_truth_subfield_dict(experiment, obs_type, naive=False, diff_tol=1.e-14,
                                           storage_dir=STORAGE_DIR, truth_dir=TRUTH_DIR):
    """Get or generate a dict mapping which subfields contain which of NFIELDS independent truth
    shear values.
    """
    subfield_dict_file = os.path.join(
        STORAGE_DIR, "subfield_dict_"+experiment[0]+obs_type[0]+".yaml") 
    use_stored = True
    if not os.path.isfile(subfield_dict_file):
        use_stored = False
    else:
        # Compare timestamps for the subfield_dict file and the first shear_params file
        # (subfield = 000) for this branch.  If the former is older than the latter, force
        # rebuild...
        dictmtime = os.path.getmtime(subfield_dict_file)
        mapper = great3.mapper.Mapper(truth_dir, experiment, obs_type, 'constant')
        shear_params_file0 = os.path.join(mapper.full_dir, "shear_params-000.yaml")
        shear_params_mtime = os.path.getmtime(shear_params_file0)
        if dictmtime < shear_params_mtime:
            use_stored = False 
    # Then load or build (and save) the subfield_dict
    if use_stored is True:
        import yaml
        with open(subfield_dict_file, 'rb') as funit:
            subfield_dict = yaml.load(funit)
    else:
        if not naive:
            raise NotImplementedError(
                "Sorry, only the naive true shear <-> subfield mapping is coded!")
        else: # Hackily make this dict by hand for now...
            for i in range(NFIELDS):
                subfield_dict.update( # In the current dataset this is just sequential
                    {i: range(i * NSUBFIELDSPERFIELD, (i + 1) * NSUBFIELDSPERFIELD)})
            import yaml
            with open(subfield_dict_file, 'wb') as fout:
                yaml.dump(subfield_dict, fout)
    # Then return
    return subfield_dict

def get_generate_const_rotations(experiment, obs_type, storage_dir=STORAGE_DIR,
                                 truth_dir=TRUTH_DIR):
    """Get or generate an array of rotation angles for Q_const calculation.

    If the rotation file has already been built for this constant shear branch, simply returns an
    array of rotation angles to align with the PSF.  This array is of shape `(NSUBFIELDS, n_epochs)`
    where the number of epochs `n_epochs` is determined from the experiment name using the mapper.

    If the rotation file has not been built, or is older than the first entry in the set of
    starshape_parameters files, the rotation file is built first.

    @param experiment     Experiment for this branch, one of 'control', 'real_galaxy', 'real_psf',
                          'multiepoch', 'full'
    @param obs_type       Observation type for this branch, one of 'ground' or 'space'
    @param storage_dir    Directory from/into which to load/store rotation files
    @param truth_dir      Root directory in which the truth information for the challenge is stored
    @return               An array containing all the rotation angles, in radians
    """
    rotfile = os.path.join(storage_dir, "rotations_"+experiment[0]+obs_type[0]+".asc")
    use_stored = True
    if not os.path.isfile(rotfile):
        use_stored = False
    else:
        # Then compare timestamps for the rotation file and the first starshape_parameters file
        # (subfield = 000, epoch =0) for this branch.  If the former is older than the latter, force
        # rebuild...
        rotmtime = os.path.getmtime(rotfile)
        mapper = great3.mapper.Mapper(truth_dir, experiment, obs_type, 'constant')
        starshape_file_template, _ ,_ = mapper.mappings['starshape_parmeters']  
        starshape_file00 = os.path.join(
            mapper.full_dir, starshape_file_template % {"epoch_index": 0, "subfield_index": 0})
        starshapemtime = os.path.getmtime(starshape_file_00)
        if rotmtime < starshapemtime:
            use_stored = False
    # Then load / build as required 
    if use_stored is True:
        rotations = np.loadtxt(rotfile)
    else:
        # To build we must loop over all the subfields and epochs
        # First work out if the experiment is multi-exposure and has multiple epochs
        if experiment in ("multiepoch", "full"):
            import great3.constants
            n_epochs = great3.constants.n_epochs
        else:
            n_epochs = 1
        # Setup the array for storing the rotation values
        rotations = np.empty((NSUBFIELDS, n_epochs))
        # Then build a mapper for this branch to find the stored StarParameters as written by the
        # builder defined in great3.builder... Note shear_type = constant always!
        mapper = great3.mapper.Mapper(sim_root_dir, experiment, obs_type, 'constant')
        output_header = "# Rotations for "+experiment+"-"+obs_type+"-constant\n"+"# epoch"
        for epoch_index in range(n_epochs):
            output_header+=" "+str(epoch_index)
            for subfield_index in range(NSUBFIELDS):
                starshape_parameters = mapper.read(
                    "starshape_parameters", data_id={
                        "epoch_index": epoch_index, "subfield_index": subfield_index})
                rotations[subfield_index, epoch_index] = .5 * np.arctan2(
                    starshape_parameters['psf_g2'], starshape_parameters['psf_g1'])
        # We have built rotations, but then save this file as ascii for use next time
        if not os.path.isdir(storage_dir):
            os.mkdir(storage_dir)
        with open(rotfile, 'wb') as fout:
            fout.write(output_header+"\n")
            np.savetxt(fout, rotations, fmt="%e20.16 " * n_epochs)
    return rotations

def get_generate_const_truth(experiment, obs_type):
    """Get or generate an array containing the true g1, g2 per subfield

    (NSUBFIELDS, 2)
    """


def Q_const(submission_file, experiment, obs_type):
    """Calculate the Q_c for a constant shear branch submission.

    @param submission_file  File containing the user submission.
    @param experiment       Experiment for this branch, one of 'control', 'real_galaxy', 'real_psf',
                            'multiepoch', 'full'
    @param obs_type         Observation type for this branch, one of 'ground' or 'space'
    @return                 The metric Q_const
    """
    if not os.path.isfile(submission_file):
        raise ValueError("Supplied submission_file '"+submission_file+"' does not exist.")
    
    # Load the submission and label the slices we're interested in
    data = np.loadtxt(submission_file)
    subfield = data[:, 0]  # Or should this be field?  Was confused slightly...
    g1sub = data[:, 1]
    g2sub = data[:, 2]
    # Load up the rotations, then rotate g1 & g2 in the correct sense
    g1rot = +g1sub * np.cos(rotation) + g2sub * np.sin(rotation)
    g2rot = -g1sub * np.sin(rotation) + g2sub * np.cos(rotation)
    # Rotate the truth in the same sense, then use the g3metrics.fitline routine to
    # perform simple linear regression
    g1true = +truth[:, 0] * np.cos(rotation) + truth[:, 1] * np.sin(rotation)
    g2true = -truth[:, 0] * np.sin(rotation) + truth[:, 1] * np.cos(rotation)
