#!/usr/bin/env python

import os
import sys
import glob
import pyfits
import constants
sys.path.append("..")
import great3sims.mapper

def file_count(filename, sex_exec="/usr/local/bin/sex", silent=False,
	           config_filename="sex.config"):
    """Run SExtractor on the given `filename`.

    @param sex_exec        Executable file for SExtractor (default=`/usr/local/bin/sex`)
    @param silent          Run silently (stdout & stderr pipe to temporary file, default=`True`)
    @param config_filename Filename for the SExtractor configuation file (default=`sex.config`)

    @return Number of objects detected in the file
    """
    import subprocess
    import tempfile
    # Make a unique (safe) tempfile to which to write the FITS-format catalog
    _, catfile = tempfile.mkstemp(suffix=".fits")
    if silent:
        stdoutfile = tempfile.mkstemp()
        with open(stdoutfile, "wb") as ftmp:
            subprocess.check_call(
                [sex_exec, filename, "-c", config_filename, "-CATALOG_NAME", catfile],
                stdout=ftmp, stderr=ftmp)
        os.remove(stdoutfile)
    else:
    	print filename
    	print sex_exec
    	print config_filename
    	print catfile
        subprocess.check_call([sex_exec, filename, "-c", config_filename, "-CATALOG_NAME", catfile])
    nobs = len(pyfits.getdata(catfile))
    os.remove(catfile)
    return nobs