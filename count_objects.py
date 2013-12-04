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

    @return good  If all checks pass, returns a dictionay containing the object totals that were
                  counted by this tool, itemized by the filename.  Any files
                  found to contain neither 9 or 10 0000 objects are judged as failed.  All failed
                  filenames are printed and the function raises a ValueError exception with all
                  the failed filenames listed in the error message.
    """
    # Set storage lists
    good = {}
    bad = {}
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

                                    nobs = count_file(fitsfile)
                                    if nobs == 3 and "starfield_image" in fitsfile:
                                        good[fitsfile] = nobs
                                    elif nobs == 10000 and "image-" in fitsfile:
                                        good[fitsfile] = nobs
                                    else:
                                        bad[fitsfile] = nobs
 
    print "Ran SExtractor on all FITS files in "+root_dir
    print "Verified object totals in FITS files for the following branches: \n"+str(found)
    print "Verified "+str(len(good))+" FITS files total"
    if len(bad) > 0:
        ret = bad
        message = "The following files failed FITS object counting:\n"
        for filename, nobs in bad.iteritems():

            message += filename+" with "+str(nobs)+" objects\n" 

        raise ValueError(message)
    else:
        print "All files verified successfully!"
        retl = good
    return ret


if __name__ == "__main__"

    good_public = count_all(constants.public_dir)
