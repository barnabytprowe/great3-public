#!/usr/bin/env python
#
# Copyright 2013 Rachel Mandelbaum & the GREAT3 Team:
# https://github.com/barnabytprowe/great3-public
#
# This file is example code released as part of the GREAT3 lensing competition.
#
# The GREAT3 example code is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This example code is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#

"""
Example module that can drive the shear estimation script simple_shear.py for an entire branch.

This Python module provides the machinery to run the simple shear estimation routine,
simple_shear.py, on all subfields in a branch.  Currently it allows two options: farming out all the
processes to run directly, or setting up scripts to run them all on a cluster using the PBS system.


DEPENDENCIES

Python 2.6.x-2.7.x
NumPy 1.6+
PBS [optional]


USAGE

Executing the script looks something like this:

./mass_produce_shear.py --output_prefix my_prefix /sims/control/ground/constant ./output 10

This will process all subfields for the branch located in '/sims/control/ground/constant', placing
all outputs in the './output' directory.  Use of the output_prefix option means the output files
will be named 'my_prefix-*.fits' (where * runs from 000 through 199); the default prefix is
'output_catalog'.  The last required argument (which here is '10') indicates the maximum number of
processes to be launched simultaneously.

For more help with the command-line interface, run the following command:

./mass_produce_shear.py --help


KNOWN ISSUES / LIMITATIONS

This code currently only permits launching jobs from the command line or via PBS.  If you have some
other batch processing system like condor, you will have to modify the function that farms out the
jobs for use on your own system.  Even for use of PBS, it is likely that the header information in
the PBS scripts (used to specify things like queues) will have to be edited in order for this script
to be usable.

"""

import sys
import time
import os
import optparse
import numpy
import subprocess

verbose = False

def log(msg):
    if verbose:
        print msg

def check_done(process_str, sleep_time=60):
    """Check the PBS queue for processes with names that contain process_str.

    Sleep for `sleep_time' seconds between checks, and return only when no jobs are found.
    """
    ind_found = 1000
    while ind_found > 0:
        if ind_found < 1000:
            # only sleep the full time before checking for processes if it's not the first time
            # through this loop
            time.sleep(sleep_time)
        else:
            # otherwise, sleep a short time just to let the queue register any jobs that were
            # submitted / ended very recently
            time.sleep(1)
        res = subprocess.check_output('qstat')
        # Find index of first occurrence of process_str in output of qstat.  We don't have to find
        # all occurrences, since we only care to know that we have reached a point where all jobs
        # are done.
        ind_found = res.find(process_str)

def MassProduce(sim_dir, output_dir, n_process, output_prefix=None, output_type=None,
                clobber=None, sn_weight=None, calib_factor=None, subfield_min=0, subfield_max=199,
                verbose_shear=False):
    """Main driver for launching jobs directly through the shell.

    Required arguments:

      sim_dir ------ simulation directory, containing the GREAT3 images and catalogs for a single
                     branch

      output_dir --- directory in which outputs should be placed

      n_process ---- maximum number of processes to launch simultaneously

    Optional arguments:

      output_prefix --- Prefix for output catalog; the subfield (as a 3-digit number) will be
                        appended

      output_type ----- Type for output catalogs: fits (default) or ascii.

      clobber --------- Overwrite pre-existing output files?  Default: true.

      sn_weight ------- Apply S/N-dependent weighting to each object?  Default: true.  If false, all
                        objects get an equal weight.

      calib_factor ---- Multiplicative calibration factor to apply to all shears; otherwise, use
                        default in simple_shear.py.

      subfield_min ---- Minimum subfield index to process; default = 0 (first subfield in the
                        branch).

      subfield_max ---- Maximum subfield index to process; default = 199 (last subfield in the
                        branch).

      verbose_shear --- Print progress statements while estimating the shear for each subfield.
    """
    t1 = time.time()

    n_subfields = subfield_max + 1 - subfield_min
    log("Preparing to analyze %d subfields."%n_subfields)

    process_list = []
    for subfield_index in range(subfield_min, subfield_max+1):
        # Check number of jobs running, and wait until it's below the limit.
        n_jobs = len(process_list)
        while n_jobs > n_process - 1:
            # Too many jobs, so wait 3 seconds, and then check whether all jobs that were launched
            # still exist.
            time.sleep(3)
            n_popped = 0
            for process_index in range(len(process_list)):
                # Make sure we don't fall off the end of the list due to some of the jobs having
                # been popped off after completion.
                if process_index >= n_jobs - n_popped:
                    n_jobs = len(process_list)
                    break
                # Check for return code from this process
                process_list[process_index].poll()
                # If there is a return code, that means the process ended, so remove it from the
                # list of active processes.
                if process_list[process_index].returncode is not None:
                    process_list.pop(process_index)
                    n_popped += 1
            n_jobs = len(process_list)

        # Launch job for this subfield, and add it to the list of active processes.
        log("...Launching job for subfield %d."%subfield_index)
        command_str = "./simple_shear.py %d %s %s "%(subfield_index, sim_dir, output_dir)
        if verbose_shear is False: command_str += "--quiet "
        if output_prefix is not None: command_str += "--output_prefix %s "%output_prefix
        if output_type is not None: command_str += "--output_type %s "%output_type
        if clobber is False: command_str += "--no_clobber "
        if sn_weight is False: command_str += "--no_weight "
        if calib_factor is not None: command_str += "--calib_factor %f "%calib_factor
        p = subprocess.Popen(command_str, shell=True, close_fds=True)
        process_list.append(p)

    # Now wait until all jobs are done.
    n_jobs = len(process_list)
    while n_jobs > 0:
        time.sleep(3)
        n_popped = 0
        for process_index in range(len(process_list)):
            # Make sure we don't fall off the end of the list.
            if process_index >= n_jobs - n_popped:
                n_jobs = len(process_list)
                break
            # Check for return code from this process
            process_list[process_index].poll()
            # If there is a return code, that means the process ended, so remove it from the list of
            # active processes.
            if process_list[process_index].returncode is not None:
                process_list.pop(process_index)
                n_popped += 1
        n_jobs = len(process_list)
        
    log("Total time: %f s."%(time.time()-t1))

def MassProducePBS(sim_dir, output_dir, n_process, output_prefix=None, output_type=None,
                   clobber=None, sn_weight=None, calib_factor=None, subfield_min=0,
                   subfield_max=199, verbose_shear=False, coadd=False, variable_psf_dir=""):
    """Main driver for launching jobs via the PBS batch system.

    Required arguments:

      sim_dir ------ simulation directory, containing the GREAT3 images and catalogs for a single
                     branch

      output_dir --- directory in which outputs should be placed

      n_process ---- maximum number of processes to launch simultaneously

    Optional arguments:

      output_prefix --- Prefix for output catalog; the subfield (as a 3-digit number) will be
                        appended

      output_type ----- Type for output catalogs: fits (default) or ascii.

      clobber --------- Overwrite pre-existing output files?  Default: true.

      sn_weight ------- Apply S/N-dependent weighting to each object?  Default: true.  If false, all
                        objects get an equal weight.

      calib_factor ---- Multiplicative calibration factor to apply to all shears; otherwise, use
                        default in simple_shear.py.

      subfield_min ---- Minimum subfield index to process; default = 0 (first subfield in the
                        branch).

      subfield_max ---- Maximum subfield index to process; default = 199 (last subfield in the
                        branch).

      verbose_shear --- Print progress statements while estimating the shear for each subfield.

      coadd ----------- Are we processing the outputs of the coaddition script, coadd_multiepoch.py?
                        If so, set this to True.  Default: false.

      variable_psf_dir- Directory in which to find PSF model outputs from psf_models.py, if this is
                        a variable PSF branch.  Default value of "" indicates that this is not a
                        variable_psf branch.

    """
    t1 = time.time()

    n_subfields = subfield_max + 1 - subfield_min
    log("Preparing to analyze %d subfields, launching jobs via PBS."%n_subfields)

    # Figure out how to divide up the subfields into the requested number of processes.
    x = numpy.arange(n_process).astype(float)
    x1 = x+1
    first_arr = (numpy.round(subfield_min + (subfield_max+1-subfield_min)*x/n_process)).astype(int)
    last_arr = (numpy.round(subfield_min + (subfield_max+1-subfield_min)*x1/n_process)-1).astype(int)

    # Loop over PBS scripts, writing them to file and launching them.
    script_names = []
    for i in range(n_process):
        first = first_arr[i]
        last = last_arr[i]
        log("PBS script %d: subfields %d to %d."%(i,first,last))

        # Set up PBS script filename, and open the file.
        script_filename = "mpsg3_%03d.sh"%i
        script_names.append(script_filename)
        with open(script_filename, "w") as f:
            # Write assorted header stuff (system-dependent).
            f.write("#!/bin/sh\n#PBS -l nodes=1:ppn=1\n#PBS -q physics\n#PBS -j oe\n")
            f.write("#PBS -l walltime=24:00:00\n\n")
            f.write("cd $PBS_O_WORKDIR\n")
            # Now, loop over subfields, and write a command for each subfield.
            for subfield_index in range(first, last+1):
                command_str = "./simple_shear.py %d %s %s "%(subfield_index, sim_dir, output_dir)
                if verbose_shear is False: command_str += "--quiet "
                if output_prefix is not None: command_str += "--output_prefix %s "%output_prefix
                if output_type is not None: command_str += "--output_type %s "%output_type
                if clobber is False: command_str += "--no_clobber "
                if sn_weight is False: command_str += "--no_weight "
                if calib_factor is not None: command_str += "--calib_factor %f "%calib_factor
                if coadd: command_str += "--coadd "
                if variable_psf_dir!="": command_str += "--variable_psf_dir %s "%variable_psf_dir
                command_str += "\n"
                f.write(command_str)
        # Finally, submit the job through PBS.
        p = subprocess.Popen("qsub %s"%script_filename, shell=True, close_fds=True)

    # Check whether jobs are done - and wait until that is true.
    log("All jobs submitted, now checking whether they are done...")
    check_done("mpsg3", sleep_time=3)

    # Clean up scripts, etc.
    log("All jobs done, now removing temporary scripts.")
    for name in script_names:
        os.remove(name)

    log("Total time: %f s."%(time.time()-t1))


def main(argv):
    usage = "usage: %prog [options] SIM_DIR WORK_DIR N_PROCESS"
    description = """Drive script to estimate shears for all galaxies in a set of subfields, applying all necessary responsivity and calibration factors.  SIM_DIR is the directory containing GREAT3 images and catalogs for the branch of interest. WORK_DIR is the directory where output files should be placed.  It will be created if it does not exist.  N_PROCESS is the maximum number of processes to simultaneously launch."""
    parser = optparse.OptionParser(usage=usage, description=description)
    parser.add_option("--output_prefix", dest="output_prefix", type=str, default=None,
                      help="Prefix for output file")
    parser.add_option("--output_type", dest="output_type", type=str, default=None,
                      help="Type of output catalog: fits or ascii")
    parser.add_option("--no_clobber", action="store_true",
                      help="Do not clobber pre-existing output files")
    parser.add_option("--no_weight", action="store_true",
                      help="Do not apply S/N-dependent weighting; weight all galaxies equally")
    parser.add_option("--calib_factor", dest="calib_factor", type=float, default=None,
                      help="Multiplicative calibration factor to apply to all shears; if not specified, use default in simple_shear.py")
    parser.add_option("--coadd", dest="coadd", action='store_true', default=False,
                      help="Use to indicate that we are processing coadd_multiepoch.py outputs")
    parser.add_option("--variable_psf_dir", dest="variable_psf_dir", type=str, default="",
                      help="Directory in which to find variable PSF models; only use this option for a variable PSF branch!")
    parser.add_option("--quiet", dest="quiet", action='store_true', default=False,
                      help="Don't print progress statements during mass production")
    parser.add_option("--verbose_shear", dest="verbose_shear", action='store_true', default=False,
                      help="Print progress statements during shear estimation for each subfield")
    parser.add_option("--subfield_min", dest="subfield_min", type=int, default=0,
                      help="Minimum subfield to process (default 0)")
    parser.add_option("--subfield_max", dest="subfield_max", type=int, default=199,
                      help="Maximum subfield to process (default 199)")
    parser.add_option("--pbs", action="store_true",
                      help="Launch jobs via PBS rather than directly in the shell")
    opts, args = parser.parse_args()
    try:
        sim_dir, output_dir, n_process = args
    except ValueError:
        parser.error("exactly three positional arguments are required")
    if not os.path.isdir(sim_dir):
        parser.error("input directory %s does not exist or is not a directory" % sim_dir)
    try:
        n_process = int(n_process)
    except TypeError:
        parser.error("n_process argument '%s' is not an integer" % n_process)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    if opts.output_type is not None:
        if opts.output_type not in ("fits", "ascii"):
            parser.error("output_type '%s' must be one of 'fits' or 'ascii'" % opts.output_type)
    if opts.no_clobber:
        clobber = False
    else:
        clobber = True
    if opts.no_weight:
        sn_weight = False
    else:
        sn_weight = True
    global verbose
    if opts.quiet:
        verbose = False
    else:
        verbose = True
    if opts.pbs:
        MassProducePBS(
            sim_dir, output_dir, n_process,
            output_prefix=opts.output_prefix,
            output_type=opts.output_type,
            clobber=clobber,
            sn_weight=sn_weight,
            calib_factor=opts.calib_factor,
            subfield_min=opts.subfield_min,
            subfield_max=opts.subfield_max,
            verbose_shear=opts.verbose_shear,
            coadd=opts.coadd,
            variable_psf_dir=opts.variable_psf_dir
            )
    else:
        MassProduce(
            sim_dir, output_dir, n_process,
            output_prefix=opts.output_prefix,
            output_type=opts.output_type,
            clobber=clobber,
            sn_weight=sn_weight,
            calib_factor=opts.calib_factor,
            subfield_min=opts.subfield_min,
            subfield_max=opts.subfield_max,
            verbose_shear=opts.verbose_shear,
            coadd=opts.coadd,
            variable_psf_dir=opts.variable_psf_dir
            )

if __name__ == "__main__":
    main(sys.argv)
