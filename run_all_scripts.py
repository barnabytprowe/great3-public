# This script is to be used to launch all the star field jobs for variable_psf and full, so we can
# rerun the buggy ones.  Since the scripts have 'noclobber' in them, we can use the old scripts, and
# they will just do the star fields that we have removed / stashed.

import time
import os
import subprocess
import sys
import shutil
import mass_produce_utils
import getpass

sys.path.append('..')
import great3sims

# Define some basic parameters.  This includes some system-dependent things like directories for
# output.
root = '/physics/rmandelb/great3-v11'
n_config_per_branch = 10 # Number of config files to be run per branch.
sleep_time = 30 # seconds between checks for programs to be done
queue_nicely = 14

# Set which branches to run.
experiments = [
    'variable_psf',
    'full',
]
# Put space before ground, because those run faster.
obs_types = [
    'space',
    'ground',
]
shear_types = [
    'constant',
    'variable',
]

branches = [ (experiment, obs_type, shear_type)
             for experiment in experiments
             for obs_type in obs_types
             for shear_type in shear_types]
n_branches = len(branches)
print "Producing images for ",n_branches," branches"

prefix1 = 'g3_step2_'
all_script_names = []
for experiment, obs_type, shear_type in branches:
    e = experiment[0]
    o = obs_type[0]
    s = shear_type[0]

    # Set up the per-branch prefix
    pbs_name = prefix1+e+o+s
    # Now make the script names by looping over the number of configs
    for i_config in range(n_config_per_branch):
        script_name = pbs_name + '_psf_%03d'%i_config + '.sh'
        print script_name
        all_script_names.append(script_name)
n_scripts = len(all_script_names)
print 'Set up names for ',n_scripts,' scripts'

# Launch processes
t1 = time.time()
for script_name in all_script_names:
    command_str = 'qsub '+script_name
    if queue_nicely:
        mass_produce_utils.check_njobs('g3_', sleep_time=sleep_time, n_jobs=queue_nicely)
    p = subprocess.Popen(command_str, shell=True, close_fds=True)
mass_produce_utils.check_done('g3_', sleep_time=sleep_time)
t2 = time.time()
print
print "Time for generation of images = ",t2-t1
print
