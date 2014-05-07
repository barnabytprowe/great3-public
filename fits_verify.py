#!/usr/bin/env python

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
"""@file fitsverify.py 

Use the fitsverify exectutable to check basic FITS properties for all the files in a GREAT3 release 
specified by the `constants.py` module in this same directory.
"""

import os
import sys
import glob
import constants
sys.path.append("..")
import great3sims.mapper


def file_good(filename, fitsverify_exec="/usr/local/bin/fitsverify", silent=True):
    """Run fitsverify on the given `filename`.

    @param fitsverify_exec Executable file for fitsverify (default=`fitsverify`).
    @param silent          Run silently (stdout & stderr pipe to temporary file, default=`True`).

    @return True if fits file is verified, false if there was an error.
    """
    import subprocess
    if silent:
        import tempfile
        tmpfile = tempfile.mktemp()
        with open(tmpfile, "wb") as ftmp:
            retcode = subprocess.call(
                [fitsverify_exec, "-q", "-e", filename], stdout=ftmp, stderr=ftmp)
        os.remove(tmpfile)
    else:
        retcode = subprocess.call([fitsverify_exec, "-q", "-e", filename])
    if retcode == 0:
        ret = True
    else:
        ret = False
    return ret

def check_all(root_dir, experiments=constants.experiments, obs_types=constants.obs_types,
              shear_types=constants.shear_types):
    """Check all FITS files are good in branches for the experiments, obs_types and shear_types
    specified,  within the specified root directory.

    Assumes all FITS files end in suffix .fits.

    @return good  If all checks pass, returns a list of all the filenames that were found by
                  this check and passed fitsverify.  If any checks failed, prints all failed
                  filenames and raises a TypeError exception with all the failed filenames listed
                  in the error message.
    """
    # Set storage lists
    good = []
    bad = []
    found_exps = []
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
                    if experiment not in found_exps: found_exps.append(experiment)
                    if obs_type not in found_obs: found_obs.append(obs_type)
                    if shear_type not in found_shears: found_shears.append(shear_type)
                    print "Verifying FITS files in "+str(mapper.full_dir)
                    for fitsfile in fitsfiles:

                        if file_good(fitsfile):
                            good.append(fitsfile)
                        else:
                            bad.append(fitsfile)

    print "Ran fitsverify on all FITS files in "+root_dir
    print "Found FITS files for experiments "+str(found_exps)
    print "Found FITS files for obs_types "+str(found_obs)
    print "Found FITS files for shear_types "+str(found_shears)
    print "Found "+str(len(good))+" FITS files total"
    if len(bad) > 0:
        retlist = bad
        message = "The following files failed FITS verify:\n"+str(bad)
        raise TypeError(message)
    else:
        print "All files verified successfully!"
        retlist = good
    print
    return retlist


if __name__ == "__main__":

    good_truth = check_all(constants.truth_dir)
    good_public = check_all(constants.public_dir)
