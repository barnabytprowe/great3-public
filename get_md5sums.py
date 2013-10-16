#!/usr/bin/env python

import os
import constants

def get_md5sum(filename, md5sum_exec="md5sum", silent=True):
    """Run md5sum on the given filename.

    @return A string containing "<MD5 checksum> experiment/obs_type/shear_type/"
    """
    import subprocess
    import tempfile
    tmpfile = tempfile.mktemp()
    with open(tmpfile, "wb") as ftmp:
        retcode = subprocess.check_call(
            [md5sum_exec, filename], stdout=ftmp)
    with open(tmpfile, "rb") as ftmp:
        retstring = ftmp.readline()
    os.remove(tmpfile) # Tidy up
    if not silent:
        print retstring
    return retstring

def collate_all(root_dir, outfile, experiments=constants.experiments, obs_types=constants.obs_types,
                shear_types=constants.shear_types, md5sum_exec="md5sum"):
    """Put together a file containing an <MD5 checksum> <filename> entry for every file in the
    all the experiment/obs_type/shear_type folders specified.
    """
    import sys
    import glob
    sys.path.append("..")
    import great3sims.mapper
    # Save the folders in which files were found
    found_exps = []
    found_obs = []
    found_shears = []
    with open(outfile, "wb") as fout:
        for experiment in experiments:

            for obs_type in obs_types:

                for shear_type in shear_types:

                    # Get *all* the files in this folder
                    mapper = great3sims.mapper.Mapper(root_dir, experiment, obs_type, shear_type) 
                    allfiles = glob.glob(os.path.join(mapper.full_dir, "*"))
                    if len(allfiles) > 0:
                        if experiment not in found_exps: found_exps.append(experiment)
                        if obs_type not in found_obs: found_obs.append(obs_type)
                        if shear_type not in found_shears: found_shears.append(shear_type)
                        print "Getting md5sum for all files in "+str(mapper.full_dir)
                        for checkfile in allfiles:

                            fout.write(get_md5sum(checkfile, silent=True))

    print "Ran md5sum on all files in "+root_dir
    print "Found files for experiments "+str(found_exps)
    print "Found files for obs_types "+str(found_obs)
    print "Found files for shear_types "+str(found_shears)
    print "Full md5sum list written to "+str(outfile)
    return


if __name__ == "__main__":

    import argparse

    description = \
    """Get all the md5sums for the all the branches in the experiments, obs_types, shear_types
    lists in validation/constants.py.

    Writes a full md5sum output file with each entry in the format <MD5 checksum> <filename> to
    the specified outfile.
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "outfile", type=str,
        help="Filename for ASCII file containing all md5sums calculated")
    parser.add_argument(
        '-root_dir', default=str(constants.public_dir),
        help="Root directory for the GREAT3 release for which you want to calculate md5sums "+
        "[default = "+str(constants.public_dir)+"]")
    parser.add_argument(
        "-md5sum", default="md5sum", help="Path to md5sum executable [default = md5sum]") 
    args = parser.parse_args()
    collate_all(args.root_dir, args.outfile, md5sum_exec=args.md5sum)

