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
        subprocess.check_call([sex_exec, filename, "-c", config_filename, "-CATALOG_NAME", catfile])
    nobs = len(pyfits.getdata(catfile))
    os.remove(catfile)
    return nobs

def count_all(root_dir, experiments=constants.experiments, obs_types=constants.obs_types,
              shear_types=constants.shear_types):
    """Check (including detection) objects in all FITS files in branches for the experiments,
    obs_types and shear_types specified,  within the specified root directory.

    Assumes all FITS files end in suffix .fits.

    @return good  If all checks pass, returns a list of all the object totals that were counted by
                  this tool, alongside the filename in a (nojects, filename) tuple.  Any files
                  found to contain neither 9 or 10 0000 objects are judged as failed.  All failed
                  filenames are printed and the function raises a ValueError exception with all
                  the failed filenames listed in the error message.
    """
    # Set storage lists
    good = []
    bad = []
    found = {}
    found_obs = []
    found_shears = []
    for experiment in experiments:

        for obs_type in obs_types:

            for shear_type in shear_types:

                # Check for all files matching *.fits
                mapper = great3sims.mapper.Mapper(root_dir, experiment, obs_type, shear_type)
                fitsfiles = glob.glob(os.path.join(mapper.full_dir, "*.fits"))
                # If any fits files found, check
                if len(fitsfiles) > 0:
                    if experiment not in found_exps:
                        found[experiment] = {}
                        if obs_type not in found[obs_type]:
                            found[experiment][obs_type] = {}
                            if shear_type not in found[experiment][obs_type]:
                                found[experiment][obs_type][shear_type] = True
                                print "Counting objects in FITS files in "+str(mapper.full_dir)
                                for fitsfile in fitsfiles:

                                    nobs = count_file good.append(fitsfile)


    print "Ran fitsverify on all FITS files in "+root_dir
    print "Found FITS files for experiments "+str(found_exps)
    print "Found FITS files for obs_types "+str(found_obs)
    print "Found FITS files for shear_types "+str(found_shears)
    print "Found "+str(len(good))+" FITS files total"

if __name__ == "__main__"

    print "blah"

