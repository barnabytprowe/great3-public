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
import great3.mass_produce as mass_produce
try:
    import g3metrics
except ImportError:
    path, module = os.path.split(__file__)
    sys.path.append(os.path.join(path, "..", "..", "metrics")) # Appends the great3-private/metrics
                                                               # folder to path
    import g3metrics

NSUBFIELDS = 200 # Total number of subfields, not necessarily equal to the number of subfields made
                 # in mass_produce as that script also generates the deep fields

def get_generate_rotations(experiment, obs_type, rot_dir="./rotations",
                           sim_root_dir="/great3/public/truth"):
    """If the rotation file has already been built for this constant shear branch, simply returns an
    array of rotation angles to align with the PSF.  This array is of shape `(NSUBFIELDS, n_epochs)`
    where the number of epochs `n_epochs` is determined from the experiment name using the mapper.

    If the rotation file has not been built, does this first.

    @param[in] experiment    experiment for this branch, one of 'control', 'real_galaxy',
                             'real_psf', 'multiepoch', 'full'
    @param[in] obs_type      observation type for this branch, one of 'ground' or 'space'
    @return An array containing all the rotation angles, in radians.

    TODO: Would like to work out a way to compare timestamps or something so that we can ensure that
    the rotation files get rebuilt if necessary after a regeneration of the sims.
    """
    rotfile = os.path.join(rot_dir, "g3rot_"+experiment[0]+obs_type[0]+".asc")
    if os.path.isfile(os.path.join(rotfile)):
        rotations = np.loadtxt(rotfile)
    else:
        # If we haven't already built the rotation file, loop over all the subfields and epochs
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
        import great3.mapper
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
        with open(rotfile, 'wb') as fout:
            fout.write(output_header+"\n")
            np.savetxt(fout, rotations, fmt="%e20.16 " * n_epochs)
    return rotations

def const_branch(submission, rotation, truth):
    """Calculate the Q_c for a constant shear branch submission.

    @param submission_file  File containing the user submission.
    @param rotation         NumPy array containing the rotation angles that need to be applied for
                            each subfield before fitting `m_i` and `c_i`.  Supplied as radians.
    @param truth            NumPy array containing the true shears in each field.
    @return                 The metric Q_c.
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
