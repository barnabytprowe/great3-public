#!/usr/bin/env python

#    presubmission.py: a script to generate summary statistics for submission
#    to the GREAT3 Challenge (http://great3challenge.info)
#    Copyright (C) 2013  Melanie Simet & the GREAT3 Team:
#    https://github.com/barnabytprowe/great3-public
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
presubmission.py: a script to generate summary statistics for submission to the GREAT3 challenge
(http://great3challenge.info).

General information for the use of this script can be found in the README, which should have been
distributed with this file.  If you are missing the README, it can be downloaded from the
great3-public repository (http://github.com/barnabytprowe/great3-public).
"""
import math
import sys
import branch_id
import os
from collections import defaultdict
has_fits_handler = False
try:
    import pyfits as fits_handler
    # If this works, import its dependency NumPy
    import numpy
    has_fits_handler = True
except ImportError as err:
    # System does not have PyFITS installed; try AstroPy
    try:
        import astropy.io.fits as fits_handler
        # If AstroPy is found, also import NumPy
        import numpy
        has_fits_handler = True
    except ImportError as err:
        # Neither pyfits nor astropy found
        has_fits_handler = False

# Various parameters relating to the size of the input images and structure of the subfields
image_size_deg = 10. # Image size in degrees
nrows = 100 # Number of galaxies in one row of the image
subfield_grid_subsampling = 7 # To help convert galaxy ID to field position
ngalaxies_per_subfield = nrows**2
nsubfields_per_field = {}
nfields_per_branch = {}
for branch in branch_id.branch_names: # Set up the field/subfield structure per branch
    if 'constant' in branch and not ('variable_psf' in branch or 'full' in branch):
        nsubfields_per_field[branch] = 1
        nfields_per_branch[branch] = 200
    else:
        nsubfields_per_field[branch] = 20
        nfields_per_branch[branch] = 10

# Parameters for the correlation function
min_sep = 0.02 # in degrees
max_sep = 10.0 # in degrees
nbins = 15

def get_field_id(subfield_list,branch = 'control_ground_constant'):
    """
    Return a list of field IDs corresponding to the subfield IDs in the input list.
    """
    if 'constant' in branch and not ('variable_psf' in branch or 'full' in branch):
        return subfield_list # In this case, "subfield" list is actually field list.
    else:
        if hasattr(subfield_list,'__getitem__'): # If a list, return a list
            return [s/nsubfields_per_field[branch] for s in subfield_list]
        return s/nsubfields_per_field[branch] # Otherwise, return a single data point

def get_branch(file_list):
    """
    Search the current directory and the file_list for branch names; raise a warning for conflicts
    and, if a single branch can be determined, return its name.
    """
    path = os.getcwd() # Current directory
    path_branch = branch_id.branch_regex.findall(path) # Find unique character strings that match
    path_branch = set(path_branch)                     # branch names in the current directory name
    if len(path_branch)>1:
        path_branch = 'multi' # Flag to indicate that more than one branch name was matched
    elif len(path_branch)==0:
        path_branch = ''
    else:
        path_branch = list(path_branch)[0]
    file_branch = branch_id.branch_regex.findall(''.join(file_list)) # Repeat for filenames
    file_branch = set(file_branch)
    if len(file_branch)>1:
        file_branch = 'multi'
    elif len(file_branch)==0:
        file_branch = ''
    else:
        file_branch = list(file_branch)[0] 
    # Now, run through the possible branches that were found.  If only a single branch was found,
    # return that name; otherwise, if the path indicated multiple branches but the files only one,
    # or vice versa, print a warning and continue with the single indicated branch from the files
    # (or path); if multiple branches were found in every test, raise an error suggesting the 
    # command-line option -b. (Null strings are handled in generate_submission_files().)
    error_message = ('Attempted to determine branch name from current directory and filenames and '
                     'found conflicting information.  Please make sure you are analyzing only one '
                     'branch; if you are, please pass -b branchname as a command-line option and '
                     'run this script again.')
    if path_branch == 'multi': 
        if file_branch == 'multi' or not file_branch:
            raise RuntimeError(error_message)
        else:
            print ('WARNING: found multiple branch names in the current directory path, '
                   'but only one in given filenames; proceeding using the branch name '
                   'determined from the filenames: %s\n'%file_branch)
            # This function will match control/ground/constant as well as control-ground-constant,
            # so subdirectories can be matched, but all the dictionaries are indexed with 
            # strings using hyphens to delineate the parts of the branch name.  This function
            # reformats the branch name with hyphens, just in case.
            return branch_id.formatted_branch_name(file_branch)
    elif file_branch=='multi':
        if not path_branch:
            raise RuntimeError(error_message)
        else:
            print ('WARNING: found multiple branch names in the given filenames, but '
                   'only one in the current directory path; proceeding using the branch '
                   'name determined from the directory: %s\n'%path_branch)[0] 
            return branch_id.formatted_branch_name(path_branch)
    elif path_branch and file_branch:
        if path_branch==file_branch:
            return branch_id.formatted_branch_name(file_branch)
        else:
            raise RuntimeError(error_message)
    # We know there is at most one, so this "or" will return the non-empty string if it exists.
    return branch_id.formatted_branch_name(file_branch or path_branch) 
                           
def get_gridpoint(gid):
    """
    Given a 3-digit galaxy ID subnumber, return the angular position
    """
    return int(gid)*image_size_deg / (nrows * subfield_grid_subsampling)
    
def get_column_number(number=-1,name='',default_index=-1):
    """
    Return the desired column number, taking defaults into account, and print a diagnostic.
    """
    try:
        number = int(number)
    except:
        pass
    if number==-1:
        print '  Column containing', name, 'is', default_index, '(default)'
        return default_index
    print '  Column containing', name, 'is', number
    return number
 
def get_list_from_fitsfile(arr,column,filename):
    try:
        return list(arr.field(column))
    except KeyError:
        raise KeyError('Column %s not found in file %s'%(column,filename))
    except IndexError:
        raise IndexError('Column %i not found in file %s'%(column,filename))
 
def readfitsfile(filename, hdu_number=1):
    """
    Read in a FITS file and return its data.  This requires either the package pyfits or the
    package astropy.io.fits--otherwise it will raise an error.
    """
    try:
        obj = fits_handler.open(filename)
    except:
        raise IOError('FITS file %s could not be opened.  Please check for existence of file '
                      '& for proper extension and format.'%filename)
    try:
        data = obj[hdu_number].data
    except IndexError:
        raise IndexError('Requested HDU %i could not be found in file %s.  Please check for '
                         'correct HDU number.'%(hdu_number,filename))
    obj.close()
    return data
                              
def readfile(filename, comment_identifier='#'):
    """
    Read in an ASCII file containing a table, remove comments, and return its data as a list, each
    element of which is a list of the fields of a single line of the table.  Columns should be 
    separated with whitespace (no comma-separated values).
    """
    lenc = len(comment_identifier)
    with open(filename) as f:
        return [line.split() for line in f if line[0:lenc]!=comment_identifier]
    
def count_items(items):
    """
    Given a list of items, return a dict whose keys are the elements of the list and whose values
    are the number of times those elements appear.
    """
    # A defaultdict is essentially a Python dictionary with automatic initialization of previously
    # unknown keys.  Rough testing on my laptop says this is ~25% faster than running a set()
    # operation to generate all possible keys and then initializing the dict by hand.  A 
    # defaultdict of type int initializes all dictionary entries to 0.  This could be dangerous if 
    # we were accessing specific keys later, but since we just loop over the included keys, it 
    # shouldn't cause problems.
    nitems = defaultdict(int)
    for item in items:
        nitems[item]+=1
    return nitems

def run_corr2_ascii(x, y, e1, e2, w, min_sep=min_sep, max_sep=max_sep, nbins=nbins,
                    cat_file_suffix='_temp.cat', params_file_suffix='_corr2.params',
                    m2_file_suffix='_temp.m2', xy_units='degrees', sep_units='degrees',
                    corr2_executable='corr2'):
    """ 
    Stolen from some scripts by Barney Rowe.  This runs the corr2 code on a temporary
    catalog generated from lists of positions, shears, and weights, and returns the contents
    of the output file from corr2.
    """
    import subprocess
    import tempfile
    # Create temporary, unique files for I/O
    catfile = tempfile.mktemp(suffix=cat_file_suffix)
    paramsfile = tempfile.mktemp(suffix=params_file_suffix)
    m2file = tempfile.mktemp(suffix=m2_file_suffix)
    # Write the basic corr2.params to temp location
    print_basic_corr2_params(paramsfile, min_sep=min_sep, max_sep=max_sep, nbins=nbins,
                             xy_units=xy_units, sep_units=sep_units)
    # Write catfile and run corr2 
    f = open(catfile, 'wb')
    for (xi, yi, e1i, e2i, wi) in zip(x, y, e1, e2, w):
        if e1i < 10. and e2i < 10.:
            f.write('%e  %e  %e  %e %e\n' % (xi, yi, e1i, e2i, wi))
    f.close()
    subprocess.Popen([
        corr2_executable, str(paramsfile), 'file_name='+str(catfile), 'm2_file_name='+str(m2file),
        ]).wait()
    os.remove(paramsfile)
    os.remove(catfile)
    if not os.path.isfile(m2file):
        raise RuntimeError('Corr2 output file does not exist--this usually indicates an error '
                           'within corr2 itself.  Please check output stream for error '
                           'messages.')
    results = readfile(m2file)
    os.remove(m2file)
    return results

def run_corr2(x, y, e1, e2, w, min_sep=min_sep, max_sep=max_sep, nbins=nbins,
              cat_file_suffix='_temp.cat', params_file_suffix='_corr2.params',
              m2_file_suffix='_temp.m2', xy_units='degrees', sep_units='degrees', 
              corr2_executable='corr2'):
    import subprocess
    import tempfile
    # Create temporary, unique files for I/O
    catfile = tempfile.mktemp(suffix=cat_file_suffix)
    paramsfile = tempfile.mktemp(suffix=params_file_suffix)
    m2file = tempfile.mktemp(suffix=m2_file_suffix)
    # Write the basic corr2.params to temp location
    print_basic_corr2_params(paramsfile, min_sep=min_sep, max_sep=max_sep, nbins=nbins,
                             xy_units=xy_units, sep_units=sep_units,fits_columns=True)
    # Use fits binary table for faster I/O. (Converting to/from strings is slow.)
    # First, make the data into numpy arrays
    x_array = numpy.asarray(x).flatten()
    y_array = numpy.asarray(y).flatten()
    g1_array = numpy.asarray(e1).flatten()
    g2_array = numpy.asarray(e2).flatten()
    w_array = numpy.asarray(w).flatten()
    # Then, mask out the >= 10 values
    use_mask = numpy.logical_and.reduce([g1_array<10.,g2_array<10.])
    # And finally make the FITS file
    x_col = fits_handler.Column(name='x', format='1D', array=x_array[use_mask])
    y_col = fits_handler.Column(name='y', format='1D', array=y_array[use_mask])
    g1_col = fits_handler.Column(name='g1', format='1D', array=g1_array[use_mask])
    g2_col = fits_handler.Column(name='g2', format='1D', array=g2_array[use_mask])
    w_col = fits_handler.Column(name='w', format='1D', array=w_array[use_mask])
    cols = fits_handler.ColDefs([x_col, y_col, g1_col, g2_col, w_col])
    table = fits_handler.new_table(cols)
    phdu = fits_handler.PrimaryHDU()
    hdus = fits_handler.HDUList([phdu, table])
    hdus.writeto(catfile, clobber=True)
    subprocess.Popen([
        corr2_executable, str(paramsfile), 'file_name='+str(catfile), 'm2_file_name='+str(m2file)
        ]).wait()
    os.remove(paramsfile)
    os.remove(catfile)
    if not os.path.isfile(m2file):
        raise RuntimeError('Corr2 output file does not exist--this usually indicates an error '
                           'within corr2 itself.  Please check output stream for error '
                           'messages.')
    results = readfile(m2file)
    os.remove(m2file)
    return results

def print_basic_corr2_params(outfile, min_sep=min_sep, max_sep=max_sep, nbins=nbins,
              xy_units='degrees', sep_units='degrees',fits_columns=False):
    """
    Write a bare-bones corr2.params file (used by corr2) to the file named outfile.
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
 
def print_variable_summary(file_object, field, x_coordinates, y_coordinates, g1, g2, weight, corr2):
    """
    Print the correlation-function output of the corr2 code for the given field to the given
    file_object.  x_coordinates and y_coordinates should be grid positions IN DEGREES, g1, g2
    and weight the shears and weights for each galaxy.
    """
    if has_fits_handler:
        results = run_corr2(
            x_coordinates, y_coordinates, g1, g2, weight, min_sep=min_sep, max_sep=max_sep, 
            nbins=nbins, xy_units='degrees', sep_units='degrees', corr2_executable=corr2)
    else:
        results = run_corr2_ascii(
            x_coordinates, y_coordinates, g1, g2, weight, min_sep=min_sep, max_sep=max_sep, 
            nbins=nbins, xy_units='degrees', sep_units='degrees', corr2_executable=corr2)
    for line in results:
        file_object.write(str(field)+' '+line[0]+' '+line[1]+' '+line[2]+' '+line[5]+'\n')
    
def print_constant_summary(file_object,field,g1,g2,weight):    
    """
    Print the average g1 and g2 shears for the field to the given file_object.
    """
    g1_weight = [a[0]*a[1] for a in zip(g1,weight) if a[0]<9.9]
    sumg1weight = sum([a[1] for a in zip(g1,weight) if a[0]<9.9])
    g2_weight = [a[0]*a[1] for a in zip(g2,weight) if a[0]<9.9]
    sumg2weight = sum([a[1] for a in zip(g2,weight) if a[0]<9.9])
    file_object.write(str(field)+' '+
                      str(sum(g1_weight)/sumg1weight)+' '+
                      str(sum(g2_weight)/sumg2weight)+'\n')    
    
def generate_submission_files(args):
    """
    Take command-line input parsed into the object 'args' and generate a submission file for the
    branch covered by the input shear file(s).  This function performs the following tasks:
        - Read in the given file(s)
        - Determine which field each galaxy belongs to based on the subfield it comes from
        - Determine the nearest gridpoint/fiducial galaxy position based on the given galaxy ID,
          if the branch type is variable shear
        - Check that the right number of galaxies are included in each subfield and field
        - Call the appropriate functions to write the summary statistics for either constant or 
          variable shear branches
    """
    # Check for existence of catalogs
    nfiles = len(args.file_list)
    if nfiles==0:
        raise RuntimeError('No shear catalogs passed to presubmission.py.  Please include '
                           'at least one shear catalog.')
    # Check branch type
    if not args.branch:
        args.branch = get_branch(args.file_list)
        if not args.branch:
            raise RuntimeError('No branch name was passed to presubmission.py, and the name could '
                               'not be determined from current path or shear catalog filenames.  '
                               'Please rerun this script with the option -b branchname set.')
    else:
        if not args.branch in branch_id.branch_names:
            raise ValueError('Branch name %s given with command-line option -b does not appear to '
                             'be a valid branch name.  Please pass one of the following with -b: '
                             '%s.'%(args.branch, ' '.join(branch_id.branch_names)))
    print 'Running presubmission.py on',
    if nfiles>1:
        print 'files',
    else:
        print 'file',
    for file in args.file_list[:3]:
        print file,
    if nfiles>3:
        print '...',
    print ''
    # Pull column numbers out of the argument list and print them to the terminal
    if args.use_fortran_convention:
        fmod=1
    else:
        fmod=0
    if 'constant' in args.branch:
        print '  Using constant-shear analysis'
    elif 'variable' in args.branch:
        print '  Using variable-shear analysis'
    id_column = get_column_number(args.id_column,name='galaxy ids', default_index=0+fmod)
    g1_column = get_column_number(args.g1_column,name='g1', default_index=1+fmod)
    g2_column = get_column_number(args.g2_column,name='g2', default_index=2+fmod)
    if args.weight_column:
        weight_column = get_column_number(args.weight_column, name='weight', 
                               default_index=None)
    else:
        print '  Using uniform weighting'
        weight_column = None
    if args.use_fortran_convention:
        # If we're using FORTRAN conventions, the column numbers are too large by 1 right now;
        # subtract 1 from the non-string values.  (Strings may be used for indexing FITS files.)
        print '  Using Fortran-style (starting-from-1) column numbering'
        if not isinstance(id_column,str): 
            id_column-=1
        if not isinstance(g1_column,str):
            g1_column-=1
        if not isinstance(g2_column,str):
            g2_column-=1
        if weight_column and not isinstance(weight_column,str):
            weight_column-=1
    # Determine the largest column we'll be requesting, to check the files for proper size later.
    max_column = 0
    for col in (id_column,g1_column,g1_column,weight_column):
        if col and not isinstance(col,str):
            max_column = max(max_column,col)

    # Check if all of the files are fits files--due to differences in column-index handling,
    # it's not possible to mix the two.
    if args.use_fits:
        is_fitsfile = True
    else:
        # This makes a list of True/False values, True if the extension looks like FITS
        file_type = [
            (os.path.splitext(f)[1].lower()=='.fits' or os.path.splitext(f)[1].lower()=='.fit') 
            for f in args.file_list]
        # This makes a list of only the True values, so we can count the number of FITS files
        is_fitsfile = [f for f in file_type if f]
        if len(is_fitsfile)==len(file_type):
            is_fitsfile=True
        elif len(is_fitsfile)==0:
            is_fitsfile=False
        else:
            raise RuntimeError('Some, but not all, requested shear catalogs are FITS files.  '
                               'Please call this script with only ASCII or only FITS files.')
    # Read in all the shear information
    shear_info = []
    if is_fitsfile:
        if nfiles==1:
            print 'Catalog file determined to be FITS file...'
        else:
            print 'Catalog files determined to be FITS files...'
        if not has_fits_handler:
            raise IOError('No FITS handler found.  Please install pyfits or astropy and rerun, or '
                            'convert catalogs to ASCII format.')
        for ifn, filename in enumerate(args.file_list):
            # So that we don't have to deal with processing issues later with NumPy arrays vs
            # lists of lists, read the FITS file(s) and then convert the data we want into
            # a list of lists.  Once this is done, rewrite the column numbers so the rest of
            # the program works properly.  (This is why the script can't handle both FITS
            # and ASCII files at once.)
            if ifn%20==0 or ifn==nfiles-1:
                print 'Reading file', str(ifn+1)+'/'+str(nfiles), '...'
            newdata = readfitsfile(filename)
            gid = get_list_from_fitsfile(newdata, id_column, filename)
            g1 = get_list_from_fitsfile(newdata, g1_column, filename)
            g2 = get_list_from_fitsfile(newdata, g2_column, filename)
            if weight_column:
                w = get_list_from_fitsfile(newdata, weight_column, filename)
            if weight_column:
                shear_info+=zip(gid,g1,g2,w)
            else:
                shear_info+=zip(gid,g1,g2)
        id_column = 0
        g1_column = 1
        g2_column = 2
        if weight_column:
            weight_column = 3
    else:
        if (isinstance(id_column,str) or isinstance(g1_column,str) or 
            isinstance(g2_column,str) or isinstance(weight_column,str)):
            raise RuntimeError('Some column IDs given as strings, which is not allowed for ASCII '
                               'catalog files.  Please use integer column numbers.')
        for ifn,filename in enumerate(args.file_list):
            if ifn%20==0 or ifn==nfiles-1:
                print 'Reading file', str(ifn+1)+'/'+str(nfiles), '...'
            newdata = readfile(filename, comment_identifier=args.comment_identifier)
            shear_info+=newdata
            min_column = min([len(si) for si in shear_info]) # Returns shortest row in the table
            if max_column >= min_column:
                raise RuntimeError('Requested column %i is greater than the number of columns in '
                                   'file %s'%(max_column+fmod,filename))
    print 'All files read.'
    print 'Checking for correct galaxy IDs... (this may take a while)'
    # Check that each submitted branch has all its galaxies submitted.
    # Start by checking that the correct total numbre of galaxies appears.
    ngals = nfields_per_branch[args.branch]*nsubfields_per_field[args.branch]*ngalaxies_per_subfield
    if len(shear_info)!=ngals:
        raise RuntimeError('Incorrect number of galaxies: %i instead of the expected %i.  Please '
                           'check your shear catalogs for completeness and make sure you have '
                           "only included one branch's worth of galaxies, then run this script "
                           'again.'%(len(shear_info),ngals))
    # To parse the next line: take the id_column of shear_info, turn it into an int (via float and
    # round in case it was read in as a float from a FITS file), zero-pad it, take the first 3 
    # characters of the resulting string, and turn those into an int.  (Just in case the catalog 
    # formatting has stripped the leading 0s!)
    galaxy_id = ['%09i'%int(round(float(si[id_column]))) for si in shear_info]
    subfield_id = [int(gid[:3]) for gid in galaxy_id]
    field_id = get_field_id(subfield_id,branch=args.branch)
    nitems = count_items(field_id)       # count_items returns a defaultdict of the count of each
    fields = [field for field in nitems] # item in the list that gets passed to it
    # Check that the correct number of fields was found in count_items.
    # (These intermediate steps are to help pin down the problem for people who have the correct
    # number of galaxies, but they're distributed in some odd way.)
    if len(fields)!=nfields_per_branch[args.branch]:
        raise RuntimeError('Incorrect number of fields--please check your shear catalogs and '
                           'try again.')
    for field in nitems: # remember, nitems is currently the dict of ngalaxies per field      
        # Check that each field has the correct number of galaxies
        if nitems[field]!=nsubfields_per_field[args.branch]*ngalaxies_per_subfield:
            raise RuntimeError('Incorrect number of subfields in field %s--please check your '
                               'shear catalogs and try again.'%field)
    nitems = count_items(subfield_id) 
    # Check that each subfield has the correct number of galaxies
    for subfield in nitems:           
        if nitems[subfield]!=ngalaxies_per_subfield:
            raise RuntimeError('Incorrect number of galaxies in subfield %s--please check your '
                               'shear catalogs and try again.'%subfield)
    # Now, check uniqueness of galaxy IDs                               
    if len(set(galaxy_id))!=len(galaxy_id):
        raise RuntimeError('Duplicate galaxy IDs found.  Please check your catalog for uniqueness '
                           'of galaxy IDs.')
    # and, finally, check that the set of galaxy IDs is the expected set of galaxy IDs.
    if set(galaxy_id)!=set(branch_id.expected_galaxy_ids(args.branch)):
        raise RuntimeError('Unexpected galaxy IDs found.  Please check that your galaxy IDS are '
                           'the same as the catalogs distributed with the Great3 simulations, '
                           'and that you have passed the correct branch name.')

    # See the comments in count_items for a description of what defaultdict does
    # This makes a set of dicts x, y, etc.; the keys are the field IDs, and the values are a list of
    # the quantity x, y, etc for the galaxies in that field. 
    x = defaultdict(list)
    y = defaultdict(list)
    g1 = defaultdict(list)
    g2 = defaultdict(list)
    weight = defaultdict(list)

    for i,si in enumerate(shear_info):
        g1[field_id[i]].append(float(si[g1_column]))
        g2[field_id[i]].append(float(si[g2_column]))
    if 'constant' not in args.branch: # Only need x and y positions for variable-shear branches
        for i,si in enumerate(shear_info):
            idnum = galaxy_id[i]
            x[field_id[i]].append(get_gridpoint(idnum[3:6]))
            y[field_id[i]].append(get_gridpoint(idnum[-3:]))
    if weight_column:
        for i,si in enumerate(shear_info):
            weight[field_id[i]].append(float(si[weight_column]))
    else:
        for field in g1:
            weight[field] = [1.]*len(g1[field])
    # Finally, compute the summary statistics.
    with open(args.outfile,'w') as file_object:
        if 'constant' in args.branch:
            for field in fields:
                print 'Computing summary statistics for field', field, '...'
                print_constant_summary(file_object,field,g1[field],g2[field],weight[field])
        # Check if "variable" is in the branch name, apart from where it appears in "variable_psf"
        # if that is part of the branch name
        elif 'variable' in args.branch.replace('variable_psf',''):
            for field in fields:
                print 'Computing summary statistics for field', field, '...'
                print_variable_summary(file_object,field,x[field],y[field],
                                       g1[field],g2[field],weight[field],
                                       corr2=args.corr2)
        else:
            # This should never happen, but just in case...
            raise RuntimeError('Could not determine analysis type from branch name %s!  Please '
                               'make sure that the branch name contains either the word '
                               '"constant" or the word "variable" (not part of the string '
                               '"variable_psf").'%args.branch)
    print 'All summary statistics calculated successfully!  Output written to', args.outfile
        
if __name__=='__main__':
    try: 
        # This is our preferred option for analyzing command-line arguments, but it doesn't
        # work for Python 2.6 and earlier
        import argparse
        parser = argparse.ArgumentParser(description='Create a set of submission files for the '
                                                     'Great3 challenge from one or more shear '
                                                     'catalogs')
        parser.add_argument('--version', action='version', 
                            version='great3-public: presubmission.py 0.1')    
        # The next few lines read in the column numbers describing the tables, if given.
        # Note that the add_argument default and the help default look different, but the correct
        # defaults are set in generate_submission_files() -- the default here is just so the main
        # program diagnostics know that the column was not passed as a command-line argument.
        parser.add_argument('-i','--id-column', 
                            help='Number of the column (or FITS column name) containing galaxy '
                                 'IDs (default: 0)',
                            default='-1', type=str, dest='id_column')
        parser.add_argument('-g1','--g1-column', 
                            help='Number of the column (or FITS column name) containing g1 shape '
                                 'of galaxies (default: 1)',
                            default='-1', type=str, dest='g1_column')
        parser.add_argument('-g2','--g2-column', 
                            help='Number of the column (or FITS column name) containing g1 shape '
                                 'of galaxies (default: 2)',
                            default='-1', type=str, dest='g2_column')
        parser.add_argument('-w','--weight-column', 
                            help='Number of the column (or FITS column name) containing weights '
                                 'for galaxies (default: not used)',
                            type=str, dest='weight_column')
        parser.add_argument('-b','--branch', 
                        help='Which branch the catalog(s) apply to (default: determined from path '
                            'or filename)',
                        default=None, dest='branch')
        # This syntax creates a variable use_fortran_convention which is True if -f is passed or
        # False if it is not
        parser.add_argument('-f','--use-fortran-convention', 
                          action='store_const', const=True, default=False,
                          help='Whether to use Fortran-numbered columns (ie first column is 1, '
                               'not 0; default is False)',
                          dest='use_fortran_convention')
        parser.add_argument('-c','--comment-identifier', 
                          help='What symbol to use to indicate comment lines in the shear '
                               'catalogs (default: #)',
                          default='#', dest='comment_identifier')
        parser.add_argument('-hdu','--hdu-number', 
                          help='Which extension of the FITS file contains the shear-catalog table '
                               '(default: 1)',
                          default='1', type=int, dest='hdu_number')
        parser.add_argument('-o','--output-file', 
                          help='Name of file to write the summary statistics to',
                          default='great3_submission.txt', type=str, dest='outfile')
        parser.add_argument('--use-fits', 
                          help='Force the script to process all catalogs as FITS files '
                               '(default: False)',
                          action='store_const', const=True, default=False, dest='use_fits')
        parser.add_argument('-c2','--corr2', 
                          help="The location of the executable for Mike Jarvis's corr2 code "
                               '(default: corr2)',
                          default='corr2', type=str, dest='corr2')
        # Finally, collect any filenames given
        parser.add_argument('file_list',nargs='*',help='One or more names of files containing '
                                                       'shear catalogs')
        # Parse the command-line options and pass everything to generate_submission_scripts
        args = parser.parse_args()
    except ImportError:
        # Use optparse.  It's less optimal in that we can't use single-hyphen flags with
        # more than one character, so only the long version of flags -g1, -g2, -hdu, and -c2
        # are accepted.
        import optparse
        parser = optparse.OptionParser()
        parser.add_option('-i','--id-column', 
                          help='Number of the column (or FITS column name) containing galaxy '
                               'IDs (default: 0)',
                          default='-1', type=str, dest='id_column')
        parser.add_option('--g1-column', 
                          help='Number of the column (or FITS column name) containing g1 shape '
                               'of galaxies (default: 1)',
                          default='-1', type=str, dest='g1_column')
        parser.add_option('--g2-column', 
                          help='Number of the column (or FITS column name) containing g1 shape '
                               'of galaxies (default: 2)',
                          default='-1', type=str, dest='g2_column')
        parser.add_option('-w','--weight-column', 
                          help='Number of the column (or FITS column name) containing weights '
                               'for galaxies (default: not used)',
                          type=str, dest='weight_column')
        parser.add_option('-b','--branch', 
                          help='Which branch the catalog(s) apply to (default: determined from '+
                               'path or filename)',
                          default=None, dest='branch')
        parser.add_option('-f','--use-fortran-convention', 
                          action='store_const', const=True, default=False,
                          help='Whether to use Fortran-numbered columns (ie first column is 1, '
                               'not 0; default is False)',
                          dest='use_fortran_convention')
        parser.add_option('-c','--comment-identifier', 
                          help='What symbol to use to indicate comment lines in the shear '
                               'catalogs (default: #)',
                          default='#', dest='comment_identifier')
        parser.add_option('--hdu-number', 
                          help='Which extension of the FITS file contains the shear-catalog table '
                               '(default: 1)',
                          default='1', type=int, dest='hdu_number')
        parser.add_option('-o','--output-file', 
                          help='Name of file to write the summary statistics to',
                          default='great3_submission.txt', type=str, dest='outfile')
        parser.add_option('--use-fits', 
                          help='Force the script to process all catalogs as FITS files '
                               '(default: False)',
                          action='store_const', const=True, default=False, dest='use_fits')
        parser.add_option('--corr2', 
                          help="The location of the executable for Mike Jarvis's corr2 code "
                               '(default: corr2)',
                          default='corr2', type=str, dest='corr2')
        parser.add_option('--filexyz123',default=[],dest='file_list') # to create args.file_list
        args = parser.parse_args()
        # optparse gives a tuple (named_stuff, positional_stuff), so hack it to include
        # the positional arguments in the first part of the tuple
        args[0].file_list = args[1] 
        args = args[0]
    generate_submission_files(args)
    
