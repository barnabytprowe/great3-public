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
Each submission file (one per branch?) CHECK THESE DETAILS WITH MELANIE'S CODE OUTPUT
"""

import os
import numpy as np
try:
    import g3metrics
except ImportError:
    import sys
    path, module = os.path.split(__file__)
    sys.path.append(os.path.join(path, "..", "..", "metrics")) # Appends the great3-private/metrics
                                                               # folder to path
    import g3metrics

NFIELDS = 10      # This is the total number of fields in a branch
NSUBFIELDS = 200  # This is the total number of subfields in a branch (so that each field consists
                  # of NSUBFIELDS/NFIELDS subfields)

def get_generate_rotations(branchname, path="."):
    """If the rotation file has already been built for this `branchname`, simply returns an array
    of rotation angles to align with the PSF.  If the rotation file has not been built, does this
    first.
    """
    if not os.path.isfile(os.path.join(branchname)



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
