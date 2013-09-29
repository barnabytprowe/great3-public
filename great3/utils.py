import os
import sys
import numpy as np
sys.path.append(os.path.join("..", "..", "validation"))
import constants
import evaluate

def load_constant_submission(filename):
    """Load a constant shear branch submission from `filename`.

    Does no checking.

    @return subfield_index, g1, g2
    """
    data = np.loadtxt(filename)
    subfield_index, g1, g2 = data[:, 0].astype(int), data[:, 1], data[:, 2]
    return subfield_index, g1, g2

def load_variable_submission(filename):
    """Load a variable shear branch submission from `filename`.
    @return field_index, theta, map_E, map_B, maperr

    Does no real checking.

    If map_B and/or map_err are not present in the data, these are returned as zeroed float arrays
    of length equal to that of field_index, theta and map_E.
    """
    data = np.loadtxt(filename)
    field_index, theta, map_E = data[:, 0].astype(int), data[:, 1], data[:, 2]
    if data.shape[1] >= 4: 
        map_B = data[:, 3]
        if data.shape[1] == 5:
            maperr = data[:, 4]
        else:
            maperr = np.zeros_like(data[:, 1])
    else:
        map_B = np.zeros_like(data[:, 1])
        maperr = np.zeros_like(data[:, 1])
    return field_index, theta, map_E, map_B, maperr

def verify_subfield_indices(subfield_index):
    """Returns `True` if input subfield indices match expectation.
    """
    reference = range(evaluate.NSUBFIELDS)
    checklist = [subitem == refitem for subitem, refitem in zip(subfield_index, reference)]
    return all(checklist)

def verify_field_indices(field_index):
    """Returns `True` if input field indices match expectation.
    """
    reference = np.arange(evaluate.NBINS_THETA * evaluate.NFIELDS, dtype=int) / evaluate.NBINS_THETA
    checklist = [fielditem == refitem for fielditem, refitem in zip(field_index, reference)]
    return all(checklist)

def verify_theta(theta, decimal=3):
    """Returns `True` if input theta array matches expectation to within `decimal` decimal places.
    """
    try:
        np.testing.assert_array_almost_equal(
            np.asarray(theta), evaluate.EXPECTED_THETA, decimal=decimal)
    except AssertionError as err:
        retval = False
    else:
        retval = True 
    return retval

def get_verify_submission_basics(filename):
    """Check the basics of a submission: i) file exists; ii) numpy.loadtxt(filename) is a 2D array.

    If successful, returns numpy.loadtxt(filename).

    Otherwise raises an IOError if i) fails or ValueError if ii) fails.
    """
    # Try loading
    try: 
        data = np.loadtxt(filename)
    except IOError as err:
        err_msg = "Cannot open submission file: "+str(err.message)   
        raise IOError(err_msg)
    except Exception as err:
        err.message = "Trouble reading submission file ("+str(filename)+"): "+err.message
        raise err
    # Then take the length of the data.shape (if loadtxt worked this will not raise an exception)
    ndim = len(data.shape)
    if ndim != 2:
        raise ValueError(
            "Expecting submission to be a 2D table, but numpy.loadtxt() returned an array of "+
            "dimension = "+str(ndim))
    return data 

def verify_constant_submission(filename, verbose=False, raise_exception_on_fail=True,
                               check_variable=False):
    """Check whether a given submission (stored in `filename`) conforms to constant shear type
    submission expectations.

    Returns `True` if the submission looks like a valid constant shear submission, `False`
    otherwise.
    """
    # Start off assuming that this is not a valid constant submission, set True only after all hoops
    # leapt through
    constant = False 
    error = False  # Gets set if this looks like neither a variable or constant submission

    try: # Try just loading the submission
        data = get_verify_submission_basics(filename)
    except Exception as err:
        if raise_exception_on_fail:
            raise err
        elif verbose:
            print err.message
    else: # If successfully loaded, process the resultant 2D array and check its dimensions and
          # content 
        nrows = data.shape[0]
        ncols = data.shape[1]
        # Check the number of rows first, as this *must* be NSUBFIELDS
        if nrows == evaluate.NSUBFIELDS:
            if ncols == 3:
                subfield_index = data[:, 0].astype(int)
                if verify_subfield_indices(subfield_index):
                    constant = True # If ncols, nrows, & subfield_index are right, we're good
                else:
                    error = True
                    err_msg = "Submission "+str(filename)+" has the correct number of rows and "+\
                        "columns for a constant shear branch submission but the subfield "+\
                        "indices (first column) are wrong."
            else:
                error = True
                err_msg = "Submission "+str(filename)+" has correct number of rows but the wrong "+\
                    "number of columns (ncols="+str(ncols)+") to be a valid constant shear "+\
                    "branch submission!"
        else:
            if check_variable:
                if verbose:
                    print "Submission ("+str(filename)+") is not constant, checking whether "+\
                        "variable instead."
                variable = verify_variable_submission(
                    filename, verbose=verbose, raise_exception_on_fail=False, check_constant=False)
                if variable:
                    if verbose:
                        print "Submission ("+str(filename)+") looks like a valid variable shear "+\
                            "branch submission."
                else:
                    error = True
                    err_msg = "Submission ("+str(filename)+") looks like neither a constant nor "+\
                        "a variable shear branch submission!"
            else:
                error = True
                err_msg = "Submission ("+str(filename)+") does not look like a "+\
                    "a constant shear branch submission!"
    if error:
        if raise_exception_on_fail:
            raise TypeError(err_msg)
        elif verbose:
            print err_msg
    elif constant:
        if verbose:
            print "Submission ("+str(filename)+") looks like a valid constant shear branch "+\
                "submission."
    # Return 
    return constant

def verify_variable_submission(filename, verbose=False, raise_exception_on_fail=True,
                               check_constant=False):
    """Check whether a given submission (stored in `filename`) conforms to variable shear type
    submission expectations.

    Returns `True` if the submission looks like a valid variable shear submission, `False`
    otherwise.
    """
    # Start off assuming that this is not a valid variable submission, set True only after all hoops
    # leapt through
    variable = False 
    error = False  # Gets set if this looks like neither a variable or constant submission

    try: # Try just loading the submission
        data = get_verify_submission_basics(filename)
    except Exception as err:
        if raise_exception_on_fail:
            raise err
        elif verbose:
            print err.message
    else: # If successfully loaded, process the resultant 2D array and check its dimensions and
          # content 
        nrows = data.shape[0]
        ncols = data.shape[1]
        # Check the number of rows first, as this *must* be NBINS_THETA * NFIELDS
        if nrows == evaluate.NBINS_THETA * evaluate.NFIELDS:
            if (ncols >= 3) and (ncols <= 5):
                field_index = data[:, 0].astype(int)
                if verify_field_indices(field_index):
                    if verify_theta(data[:, 1]):
                        variable = True # If ncols, nrows, field_index, theta are right, we're good
                    else:
                        error = True
                        err_msg = "Submission "+str(filename)+" has the correct number of rows "+\
                            "and columns, and the correct field indices, for a variable shear "+\
                            "branch submission but the theta values (second column) are wrong at "+\
                            "at 3 decimal places."
                else:
                    error = True
                    err_msg = "Submission "+str(filename)+" has the correct number of rows and "+\
                        "columns a for a variable shear branch submission but the field "+\
                        "indices (first column) are wrong." 
            else:
                error = True
                err_msg = "Submission "+str(filename)+" has correct number of rows but the wrong "+\
                    "number of columns (ncols="+str(ncols)+") to be a valid variable shear "+\
                    "branch submission!"
        else:
            if check_constant:
                if verbose:
                    print "Submission ("+str(filename)+") is not variable, checking whether "+\
                        "constant instead."
                constant = verify_constant_submission(
                    filename, verbose=verbose, raise_exception_on_fail=False, check_variable=False)
                if constant:
                    if verbose:
                        print "Submission ("+str(filename)+") looks like a valid constant shear "+\
                            "branch submission."
                else:
                    error = True
                    err_msg = "Submission ("+str(filename)+") looks like neither a constant nor "+\
                        "a variable shear branch submission!"
            else:
                error = True
                err_msg = "Submission ("+str(filename)+") does not look like a "+\
                    "a variable shear branch submission!"
    if error:
        if raise_exception_on_fail:
            raise TypeError(err_msg)
        elif verbose:
            print err_msg
    elif variable:
        if verbose:
            print "Submission ("+str(filename)+") looks like a valid variable shear branch "+\
                "submission."
    # Return 
    return variable


if __name__ == "__main__":

    from sys import argv
    result = verify_constant_submission(argv[1], verbose=True, raise_exception_on_fail=True,
        check_variable=True)
    print result
    result = verify_variable_submission(argv[1], verbose=True, raise_exception_on_fail=True,
        check_constant=True)
    print result
