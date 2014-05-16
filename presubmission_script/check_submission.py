#!/usr/bin/env python
#
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
"""Check a submission to see whether it will be correctly read by GREAT3, and print summary
information if so.

Comments in the input submission files will be ignored if prefixed with the '#' symbol.

Usage depends on whether your Python installation is earlier than version 2.3.  (Note: the code has
not been tested with Python 3, or extensively checked with earlier Python versions, please raise an
Issue on the GitHub Issues pages if you encounter problems and we will try to fix them as soon as
possible.) 

On Python 2.3 - 2.7
-------------------
usage: check_submission.py [-h] [-c] [-v] submission

positional arguments:
  submission  Filename of GREAT3 submission to be checked

optional arguments:
  -h, --help  show this help message and exit
  -c          Force checking as a constant branch submission only
  -v          Force checking as a variable branch submission only

On earlier Python versions
--------------------------
usage: check_submission.py submission

positional arguments:
  submission  Filename of GREAT3 submission to be checked

Dependencies
------------
This script requires NumPy in the Python install used, and will fail immediately if `import numpy`
fails.
"""

import os
import sys
try:
    import numpy
except ImportError as err:
    print str(err)
    print "Very sorry, check_submission.py requires NumPy!  Please install NumPy and try again"
    exit(1) 

# Some constants describing the GREAT3 simulation datasets, used to check the submission
NSUBFIELDS = 200 # Total number of subfields per branch
NFIELDS = 10     # Number of fields (each a set of subfields) in variable / variable_psf branches
NSUBFIELDS_PER_FIELD = NSUBFIELDS / NFIELDS
NBINS_THETA = 15 # Number of angular bins of aperture mass dispersion in variable submissions
EXPECTED_THETA = numpy.array([ # Array of theta values expected in variable subs, good to 3 d.p.
    0.0246, 0.0372,  0.0563,  0.0853,  0.129 ,  0.1953,  0.2955, 0.4472,  0.6768,  1.0242,  1.5499,
    2.3455,  3.5495,  5.3716, 8.1289] * NFIELDS)


def verify_subfield_indices(subfield_index):
    """Returns `True` if input subfield indices match expectation.
    """
    reference = range(NSUBFIELDS)
    checklist = [subitem == refitem for subitem, refitem in zip(subfield_index, reference)]
    return all(checklist)

def verify_field_indices(field_index):
    """Returns `True` if input field indices match expectation.
    """
    reference = numpy.arange(NBINS_THETA * NFIELDS, dtype=int) / NBINS_THETA
    checklist = [fielditem == refitem for fielditem, refitem in zip(field_index, reference)]
    return all(checklist)

def verify_theta(theta, decimal=3):
    """Returns `True` if input theta array matches expectation to within `decimal` decimal places.
    """
    try:
        numpy.testing.assert_array_almost_equal(
            numpy.asarray(theta), EXPECTED_THETA, decimal=decimal)
    except AssertionError as err:
        retval = False
    else:
        retval = True 
    return retval

def get_verify_submission_basics(filename):
    """Check the basics of a submission: i) file exists; ii) numpy.loadtxt(filename) is a 2D array.

    If successful, returns data = numpy.loadtxt(filename).

    Otherwise raises an IOError if i) fails or ValueError if ii) fails.
    """
    # Try loading
    try: 
        data = numpy.loadtxt(filename)
    except IOError as err:
        err_msg = "Cannot open submission file ("+str(filename)+"): "+str(err)
        raise IOError(err_msg)
    except Exception as err:
        err_msg = "Trouble reading submission file ("+str(filename)+"): "+str(err)+"\n"+\
            "Does "+str(filename)+" contain a valid data table?"
        raise err.__class__(err_msg)
    # Then take the length of the data.shape (if loadtxt worked this will not raise an exception)
    ndim = len(data.shape)
    if ndim != 2:
        raise ValueError(
            "Expecting submission ("+str(filename)+") to be a 2D table, but numpy.loadtxt() "+\
            "returned an array of dimension = "+str(ndim)+"!")
    return data 

def verify_constant_submission(data, filename, verbose=False, raise_exception_on_fail=True,
                               check_variable=False):
    """Check whether a given submission (supplied as the array data) conforms to constant shear type
    submission expectations.

    Assumes that the data array has been returned by get_verify_submission_basics().
 
    Returns True if the submission looks like a valid constant shear submission, False otherwise.
    """
    # Start off assuming that this is not a valid constant submission, set True only after all hoops
    # leapt through
    constant = False 
    error = False  # Gets set if this looks like neither a variable or constant submission
    # Get the shape parameters
    nrows = data.shape[0]
    ncols = data.shape[1]
    # Check the number of rows first, as this *must* be NSUBFIELDS
    if nrows == NSUBFIELDS:
        if ncols == 3:
            subfield_index = data[:, 0].astype(int)
            if verify_subfield_indices(subfield_index):
                constant = True # If ncols, nrows, & subfield_index are right, we're good
            else:
                error = True
                err_msg = "Submission "+str(filename)+" has the correct number of rows and "+\
                    "columns for a constant shear branch submission but the subfield "+\
                    "indices (first column) are wrong!"
        else:
            error = True
            err_msg = "Submission "+str(filename)+" has correct number of rows but the wrong "+\
                "number of columns (ncols="+str(ncols)+") to be a valid constant shear "+\
                "branch submission!"
    else:
        if check_variable:
            if verbose:
                print "Submission ("+str(filename)+") does not look like a valid constant "+\
                    "submission, checking whether variable instead..."
            variable = verify_variable_submission(
                data, filename, verbose=False, raise_exception_on_fail=False,
                check_constant=False)
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
            err_msg = "Submission ("+str(filename)+") does not look like a valid constant shear "+\
                "branch submission!"
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

def verify_variable_submission(data, filename, verbose=False, raise_exception_on_fail=True,
                               check_constant=False):
    """Check whether a given submission (supplied as the array data) conforms to variable shear type
    submission expectations.

    Assumes that the data array has been returned by get_verify_submission_basics().

    Returns True if the submission looks like a valid variable shear submission, False otherwise.
    """
    # Start off assuming that this is not a valid variable submission, set True only after all hoops
    # leapt through
    variable = False 
    error = False  # Gets set if this looks like neither a variable or constant submission
    # Get the shape parameters
    nrows = data.shape[0]
    ncols = data.shape[1]
    # Check the number of rows first, as this *must* be NBINS_THETA * NFIELDS
    if nrows == NBINS_THETA * NFIELDS:
        if (ncols >= 4) and (ncols <= 5):
            field_index = data[:, 0].astype(int)
            if verify_field_indices(field_index):
                if verify_theta(data[:, 1]):
                    variable = True # If ncols, nrows, field_index, theta are right, we're good
                else:
                    error = True
                    err_msg = "Submission "+str(filename)+" has the correct number of rows "+\
                        "and columns, and the correct field indices, for a variable shear "+\
                        "branch submission but the theta values (second column) are wrong at "+\
                        "3 decimal places!"
            else:
                error = True
                err_msg = "Submission "+str(filename)+" has the correct number of rows and "+\
                    "columns for a variable shear branch submission but the field "+\
                    "indices (first column) are wrong!" 
        else:
            error = True
            err_msg = "Submission "+str(filename)+" has correct number of rows but the wrong "+\
                "number of columns (ncols="+str(ncols)+") to be a valid variable shear branch "+\
                "submission!"
    else:
        if check_constant:
            if verbose:
                print "Submission ("+str(filename)+") does not look like a valid variable "+\
                    "submission, checking whether constant instead..."
            constant = verify_constant_submission(
                data, filename, verbose=False, raise_exception_on_fail=False, 
                check_variable=False)
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
            err_msg = "Submission ("+str(filename)+") does not look like a variable shear branch "+\
                "submission!"
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

def print_constant_summary(data, filename):
    """Print a summary of data for a verified constant branch submission.
    """
    print "Summary of "+str(filename)+":"
    print "g1 = "+str(data[:, 1])
    print "g2 = "+str(data[:, 2])
    print "End of summary for successful constant shear submission "+str(filename)+", returning..."
    return

def print_variable_summary(data, filename):
    """Print a summary of data for a verified variable branch submission.
    """
    print "Summary of "+str(filename)+":"
    print "field = "+str(data[:, 0])
    print "theta = "+str(data[:, 1])
    print "map_E = "+str(data[:, 2])
    print "map_B = "+str(data[:, 3]) 
    print "End of summary for successful variable shear submission "+str(filename)+", returning..."
    return


# Main body of program
if __name__ == "__main__":

    description = \
    """Checks a submission to see whether it will be correctly read by GREAT3, and
    prints summary information if so
    """
    # Test for & set up the argument parser if present, otherwise do basic I/O 
    try:
        import argparse
    except ImportError as err:
        try:
            import optparse
        except:
            has_optparse = False
        else:
            has_optparse = True
        has_argparse = False
    else:
        has_argparse = True
        has_optparse = False # Doesn't matter if it has optparse--argparse takes precedence
    # Get submission filename 
    if has_argparse:
        parser = argparse.ArgumentParser(description=description)
        parser.add_argument(
            "submission", type=str, help="Filename of GREAT3 submission to be checked")
        parser.add_argument(
            "-c", action="store_true", help="Force checking as a constant branch submission only")
        parser.add_argument(
            "-v", action="store_true", help="Force checking as a variable branch submission only")
        args = parser.parse_args()
        if args.c and args.v:
            parser.print_help()
            print "ArgumentError: You cannot specify both -c and -v on the command line"
            exit(1)        
        submission = args.submission
    elif has_optparse:
        parser = optparse.OptionParser()
        parser.add_option(
            "-c", action="store_true", help="Force checking as a constant branch submission only")
        parser.add_option(
            "-v", action="store_true", help="Force checking as a variable branch submission only")
        args, submission = parser.parse_args()
        if len(submission)>1:
            parser.print_help()
            print "ArgumentError: You cannot analyze more than 1 submission file"
            exit(1)
        elif len(submission)==0:
            parser.print_help()
            print "ArgumentError: You must specify a submission file to check"
            exit(1)
        submission = submission[0]
        if args.c and args.v:
            parser.print_help()
            print "ArgumentError: You cannot specify both -c and -v on the command line"
            exit(1)        
    else:
        from sys import argv
        if len(argv) != 2:
            print "usage: check_submission.py submission"
            print
            print description
            print
            print "positional arguments:"
            print "  submission  Filename of GREAT3 submission to be checked" 
            exit(1)
        submission = argv[1]
    # OK, first of all we verify the basics and load the data
    data = get_verify_submission_basics(submission)
    # Then we do a few checks... First of all if we have been told to only inspect one
    # type of submission, check only that type and be verbose about it.
    if args.c:
        constant = verify_constant_submission(data, submission, check_variable=False, verbose=True)
        if constant: # If successful print summary
            print_constant_summary(data, submission)
    elif args.v:
        variable = verify_variable_submission(data, submission, check_constant=False, verbose=True)
        if variable: # If successful print summary
            print_variable_summary(data, submission)
    else:
        # First silently check (quick) and then work out what to call with what verbosity to get
        # most useful info
        constant = verify_constant_submission(
            data, submission, check_variable=False, verbose=False, raise_exception_on_fail=False)
        variable = verify_variable_submission(
            data, submission, check_constant=False, verbose=False, raise_exception_on_fail=False)
        if constant and not variable: # This looks constant, print success message and summary
            constant = verify_constant_submission(
                data, submission, check_variable=False, verbose=True)
            print_constant_summary(data, submission)
        elif variable and not constant: # This looks variable, print success message and summary
            variable = verify_variable_submission(
                data, submission, check_constant=False, verbose=True)
            print_variable_summary(data, submission)
        elif variable and constant: # Uh-oh something went very wrong, abort with prejudice!
            raise RuntimeError(
                "Submission seems to be of both constant and variable shear type, total fail!\n"+
                "Please report on the GitHub great3-public repository Issues pages!")
        else:
            constant = verify_constant_submission(
                data, submission, check_variable=False, verbose=True, raise_exception_on_fail=False)
            variable = verify_variable_submission(
                data, submission, check_constant=False, verbose=True, raise_exception_on_fail=False)
            raise RuntimeError(
                "Checked for both variable and constant shear type characteristics without "+
                "success.  See info above.")
