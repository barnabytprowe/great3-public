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
"""Helper classes and encapsulation for I/O
"""

import os
import galsim
import shutil

# We'll read and write a lot of dictionaries below, using yaml to save them right now.  But all the
# I/O goes through these functions, so we can switch to some other format easily.

def readDict(path, yaml_dict = True):
    """Read a dict from disk; 'path' should not include file extension.
    """
    if yaml_dict:
        import yaml
        with open(path + ".yaml") as stream:
            return yaml.load(stream)
    else:
        import cPickle
        with open(path + ".p") as stream:
            return cPickle.load(stream)

def writeDict(d, path, type = 'yaml', comment_pref = '#'):
    """Write a dict to disk; 'path' should not include file extension.

    Type is either yaml dict ('yaml'), pickle ('p'), or text ('txt').
    """
    if type == 'yaml':
        import yaml
        with open(path + ".yaml", 'w') as stream:
            yaml.dump(d, stream)
    elif type == 'p':
        import cPickle
        with open(path + ".p", 'w') as stream:
            cPickle.dump(d, stream, protocol=2)
    elif type == 'txt':
        with open(path + ".txt", 'w') as f:
            # First print one line with keys (as a comment).
            # Note that we want to do this alphabetically so e.g. x comes before y.
            list_keys = sorted(d)
            f.write(comment_pref + " ")
            for key in list_keys:
                f.write(key + " ")
            f.write('\n')
            # Then print one line with values
            for key in list_keys:
                f.write(str(d[key]) + " ")
            f.write('\n')
    else:
        raise NotImplementedError('%s type of file is not supported by writeDict'%type)

# We'll also be writing structured NumPy arrays, as that's what we use for catalogs.
# Similarly, these functions encapsulate how we want to save them.

def readCatalog(path, fits_catalog = True):
    """Read a catalog (structured NumPy array) from disk.  'path' should not include file
    extensions.

    Technically, when reading in a fits catalog, the result is not a NumPy array; it is a
    pyfits.FITSrec, which has all the methods we need to use the data in the same way as we would
    use a NumPy array.

    @param[in] fits_catalog        Read a FITS catalog?  If True, file extension is assumed to be
                                   .fits.  If False, it is assumed to be a pickle dump with
                                   extention .p.  [Default: fits_catalog = True.]
    """
    if fits_catalog:
        import pyfits
        return pyfits.getdata(path + ".fits")
    else:
        import cPickle
        with open(path + ".p") as stream:
            return cPickle.load(stream)

def writeCatalog(catalog, path, type = 'fits', comment_pref = '#', format = None):
    """Write a catalog (structured NumPy array) to disk; 'path' should not include file extension.

    @param[in] type        Type of file to write.  options are 'fits', 'p' (pickle dump), 'txt' 
                           [Default: 'fits']
    """
    import numpy as np
    for record in catalog:
        for val in record:
            if np.isnan(val) or np.isinf(val):
                raise RuntimeError("NaN/Inf values found in catalog!")
    if type == 'fits':
        import pyfits
        pyfits.writeto(path + ".fits", catalog, clobber = True)
    elif type == 'p':
        import cPickle
        with open(path + ".p", 'w') as stream:
            cPickle.dump(catalog, stream, protocol=2)
    elif type == 'txt':
        import tempfile
        # First print lines with column names.  This goes into a separate file for now.
        f1, tmp_1 = tempfile.mkstemp(dir='.')
        with open(tmp_1, 'w') as f:
            name_list = catalog.names
            for name in name_list:
                f.write(comment_pref + " " + name + "\n")

        # Then save catalog itself.
        f2, tmp_2 = tempfile.mkstemp(dir='.')
        if format is not None:
            np.savetxt(tmp_2, catalog, fmt=format)
        else:
            np.savetxt(tmp_2, catalog)

        # Finally, concatenate, and remove tmp files
        with open(path + '.txt', 'wb') as destination:
            shutil.copyfileobj(open(tmp_1, 'rb'), destination)
            shutil.copyfileobj(open(tmp_2, 'rb'), destination)
        # Need to close the tempfiles opened by mkstemp.  cf:
        # http://stackoverflow.com/questions/9944135/how-do-i-close-the-files-from-tempfile-mkstemp
        os.close(f1)
        os.close(f2)
        os.remove(tmp_1)
        os.remove(tmp_2)
    else:
        raise ValueError("Invalid catalog type requested!  Options are fits, p, txt.")

def fitsToTextCatalog(path, comment_pref = '#', format = None):
    """Function to copy a catalog that already exists as FITS to a text file with the same name and
    a .txt extension."""
    import pyfits
    cat = pyfits.getdata(path + ".fits")
    writeCatalog(cat, path, type = 'txt', comment_pref = comment_pref, format = format)

class Mapper(object):
    """Class that manages I/O within an experiment, setting field names and encapsulating
    which file formats to use for different datasets.

    Loosely based on the mapper/butler classes in the LSST/HSC software pipeline,
    but waaaaay simplified.
    """

    # A dictionary of {dataset-name: (path-template, reader, writer)}
    mappings = {
        "parameters": ("parameters", readDict, writeDict),
        "field_parameters": ("field_parameters-%(field_index)03d", readDict, writeDict),
        "subfield_parameters": ("subfield_parameters-%(subfield_index)03d", readDict, writeDict),
        "epoch_parameters": \
            ("epoch_parameters-%(subfield_index)03d-%(epoch_index)1d", readDict, writeDict),
        "subfield_catalog": ("subfield_catalog-%(subfield_index)03d", readCatalog, writeCatalog),
        "epoch_catalog": \
            ("epoch_catalog-%(subfield_index)03d-%(epoch_index)1d", readCatalog, writeCatalog),
        "star_catalog": \
            ("star_catalog-%(subfield_index)03d-%(epoch_index)1d", readCatalog, writeCatalog),
        "star_test_catalog": \
            ("star_test_catalog", readCatalog, writeCatalog),
        "image": ("image-%(subfield_index)03d-%(epoch_index)1d.fits",
                             galsim.fits.read, galsim.fits.write),
        "starfield_image" : ("starfield_image-%(subfield_index)03d-%(epoch_index)1d.fits",
                             galsim.fits.read, galsim.fits.write),
        "star_test_images" : ("star_test_images", galsim.fits.readCube, galsim.fits.writeCube),
        "starshape_parameters" : ("starshape_parameters-%(subfield_index)03d-%(epoch_index)1d",
                                  readDict, writeDict)
    }

    def __init__(self, root, experiment, obs_type, shear_type):
        """Initialize a Mapper with the given root and experiment parameters.

        @param[in] root        Root for the entire simulation set.
        @param[in] experiment  Experiment parameter: "control", "real_galaxy", "variable_psf",
                               "multiepoch", or "full"
        @param[in] obs_type    Type of observation to simulate: either "ground" or "space"
        @param[in] shear_type  Type of shear field: "constant" or "variable"
        """
        self.root = root
        self.dir = os.path.join(experiment, obs_type, shear_type)
        self.full_dir = os.path.join(root, experiment, obs_type, shear_type)
        if not os.path.exists(self.full_dir):
            os.makedirs(os.path.abspath(self.full_dir))

    def read(self, dataset, data_id=None, **kwds):
        """Read a dataset from disk.

        @param[in] dataset    dataset name to get; one of the keys in self.mappings
        @param[in] data_id    a dict of values with which to expand the path template
                              (the first value in self.mappings).
        Additional keyword arguments are included in the data_id dict.

        @return the loaded object.
        """
        if data_id is None: data_id = dict()
        data_id.update(kwds)
        template, reader, writer = self.mappings[dataset]
        return reader(os.path.join(self.full_dir, template % data_id))

    def write(self, obj, dataset, data_id=None, **kwds):
        """Write a dataset to disk.

        @param[in] obj        object to write
        @param[in] dataset    dataset name to get; one of the keys in self.mappings
        @param[in] data_id    a dict of values with which to expand the path template
                              (the first value in self.mappings).
        Additional keyword arguments are included in the data_id dict.
        """
        if data_id is None: data_id = dict()
        data_id.update(kwds)
        template, reader, writer = self.mappings[dataset]
        return writer(obj, os.path.join(self.full_dir, template % data_id))

    def copyTo(self, other_mapper, dataset, data_id, new_template = None):
        """Copy files in the directory structure defined by this mapper to the same location in a
        directory structure defined by some other mapper.

        @param[in] other_mapper    The mapper defining the directory structure to which the file
                                   should be copied.
        @param[in] dataset         dataset name to get; one of the keys in self.mappings
        @param[in] data_id         a dict of values with which to expand the path template
                                   (the first value in self.mappings).
        @param[out] outfile        The new file name.
        """
        template, reader, writer = self.mappings[dataset]
        infile = os.path.join(self.full_dir, template % data_id)
        if new_template is None:
            outfile = os.path.join(other_mapper.full_dir, template % data_id)
        else:
            outfile = os.path.join(other_mapper.full_dir, new_template % data_id)
        # Use of shutil.copy2 preserves some of the file metadata.
        shutil.copy2(infile, outfile)
        return outfile

    def copySub(self, other_mapper, dataset, data_id, use_cols, new_template=None):
        """Copy subsets of files in the directory structure defined by this mapper to the same
        location in a directory structure defined by some other mapper.

        @param[in] other_mapper    The mapper defining the directory structure to which the file
                                   should be copied.
        @param[in] dataset         dataset name to get; one of the keys in self.mappings
        @param[in] data_id         a dict of values with which to expand the path template
                                   (the first value in self.mappings).
        @param[in] use_cols        a list of columns to copy over (i.e., neglect the others)
        @param[in] new_template    naming template to use for output catalog, if different from
                                   previous.
        @param[out] outfile        The new file name.
        """
        import numpy
        import pyfits
        # read in the catalog
        template, reader, writer = self.mappings[dataset]
        infile = os.path.join(self.full_dir, template % data_id)
        incat = readCatalog(infile)

        # choose the subset of data to save
        outcat = numpy.zeros(len(incat),
                             dtype=numpy.dtype(use_cols))
        for col in use_cols:
            outcat[col[0]] = incat[col[0]]

        # write to output file
        if new_template is None:
            new_template = template
        outfile = os.path.join(other_mapper.full_dir, new_template % data_id)
        pyfits.writeto(outfile + ".fits", outcat, clobber = True)
        return outfile+'.fits'

    def mergeSub(self, other_mapper, dataset, dataset_2, data_id, use_cols, use_cols_2,
                 new_template=None):
        """Copy subsets of files in the directory structure defined by this mapper to the same
        location in a directory structure defined by some other mapper.  This routine merges
        different columns from different catalogs instead of just taking a subset of one
        catalog.

        @param[in] other_mapper    The mapper defining the directory structure to which the file
                                   should be copied.
        @param[in] dataset         dataset name to get; one of the keys in self.mappings
        @param[in] dataset_2       2nd dataset name to get; one of the keys in self.mappings
        @param[in] data_id         a dict of values with which to expand the path template
                                   (the first value in self.mappings).
        @param[in] use_cols        a list of columns to copy over (i.e., neglect the others)
        @param[in] use_cols_2      a list of columns to copy over (i.e., neglect the others) for
                                   second dataset
        @param[in] new_template    naming template to use for output catalog, if different from
                                   previous.
        @param[out] outfile        The new file name.
        """
        import numpy
        import pyfits
        # read in the catalog
        template, reader, writer = self.mappings[dataset]
        infile = os.path.join(self.full_dir, template % data_id)
        incat = readCatalog(infile)

        # choose the subset of data to save, combined across both datasets. Note: this list
        # manipulation is necessary to avoid overwriting use_cols.
        all_cols = []
        all_cols.extend(use_cols)
        all_cols.extend(use_cols_2)

        # read in the second catalog
        template, reader, writer = self.mappings[dataset_2]
        infile = os.path.join(self.full_dir, template % data_id)
        incat_2 = readCatalog(infile)

        # Now make a catalog for everything
        outcat = numpy.zeros(len(incat),
                             dtype=numpy.dtype(all_cols))
        for col in use_cols:
            outcat[col[0]] = incat[col[0]]
        for col in use_cols_2:
            outcat[col[0]] = incat_2[col[0]]

        # write to output file
        if new_template is None:
            new_template = template
        outfile = os.path.join(other_mapper.full_dir, new_template % data_id)
        pyfits.writeto(outfile + ".fits", outcat, clobber = True)
        return outfile+'.fits'
