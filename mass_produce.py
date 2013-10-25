# This is a script that can be used for mass-production of sims on a cluster.  Currently it is set
# up for the Warp and Coma clusters at CMU, but hopefully should not be too painful to repurpose
# elsewhere.
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
root = '/home/rmandelb.proj/data-shared/great3-v8'
n_config_per_branch = 10 # Number of config files to be run per branch.
subfield_min = 0
subfield_max = 204 # The total number of subfields is split up into n_config_per_branch config files.
gal_dir = '/home/rmandelb.proj/data-shared/great3_fit_data'
ps_dir = '/home/rmandelb/git/great3-private/inputs/ps/tables'
seed = 2718
delta_seed = 1000 # amount to increment seed for each successive branch
sleep_time = 10 # seconds between checks for programs to be done
package_only = False # only do the packaging and nothing else
do_images = True  # Make images or not?  If False with package_only also False, then just skip the
                  # image-making step: remake catalogs but do not delete files / dirs with images,
                  # then package it all up.
preload = False # preloading for real galaxy branches - irrelevant for others
# Do we want to be nice to others in the queue by rationing the number of image generation scripts
# to be submitted simultaneously?  Or just submit them all, perhaps because there are so few that it
# isn't rude to do them all at once?  If `queue_nicely = False`, then dump them all at once.  If
# `queue_nicely` is set to some number, it means that number is the maximum number to have in the
# queue at once.
queue_nicely = 10

# Set which branches to test.  For now we do the control experiment (all four branches), but nothing
# else.
experiments = [
    #'control',
    'real_galaxy',
    #'variable_psf',
    #'multiepoch',
    'full',
]
obs_types = [
    'ground',
    'space',
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

if not package_only:
    if do_images:
        # Clean up from previous runs
        shutil.rmtree(root, ignore_errors=True)

    # First we set up a process for each branch.  We use a different random seed for each, and we do
    # the following steps: metaparameters, catalogs, config.  For config, there is some additional
    # juggling to do for config file names / dirs.
    prefix1 = 'g3_step1_'
    all_config_names = []
    for experiment, obs_type, shear_type in branches:
        e = experiment[0]
        o = obs_type[0]
        s = shear_type[0]

        pbs_name = prefix1+e+o+s
        pbs_file = pbs_name+'.sh'
        python_file = pbs_name+'.py'

        # Write out some scripts.
        mass_produce_utils.pbs_script_python(pbs_file, pbs_name)
        new_config_names, new_psf_config_names, new_star_test_config_names = \
            mass_produce_utils.python_script(python_file, root, subfield_min, subfield_max,
                                             experiment, obs_type, shear_type, gal_dir, ps_dir,
                                             seed, n_config_per_branch, preload, my_step=1)
        for config_name in new_config_names:
            all_config_names.append(config_name)
        for psf_config_name in new_psf_config_names:
            all_config_names.append(psf_config_name)
        for star_test_config_name in new_star_test_config_names:
            all_config_names.append(star_test_config_name)
        seed += delta_seed

    print "Wrote files necessary to carry out metaparameters, catalogs, and config steps"
    t1 = time.time()
    for experiment, obs_type, shear_type in branches:
        e = experiment[0]
        o = obs_type[0]
        s = shear_type[0]

        pbs_name = prefix1+e+o+s
        pbs_file = pbs_name+'.sh'
        command_str = 'qsub '+pbs_file
        p = subprocess.Popen(command_str,shell=True,close_fds=True)
    # The above command just submitted all the files to the queue.  We have to periodically poll the
    # queue to see if they are still running.
    mass_produce_utils.check_done('g3_step1', sleep_time=sleep_time)
    t2 = time.time()
    # Times are approximate since check_done only checks every N seconds for some N
    print
    print "Time for generation of metaparameters, catalogs, and config files = ",t2-t1
    print

    if do_images:
        # Then we split up into even more processes for images.
        t1 = time.time()
        prefix2 = 'g3_step2_'
        for config_name in all_config_names:
            type, _ = os.path.splitext(config_name)
            pbs_file = prefix2 + type+'.sh'
            mass_produce_utils.pbs_script_yaml(pbs_file, config_name, root)
            command_str = 'qsub '+pbs_file
            if queue_nicely:
                mass_produce_utils.check_njobs('g3_', sleep_time=sleep_time, n_jobs=queue_nicely)
            p = subprocess.Popen(command_str, shell=True, close_fds=True)
        mass_produce_utils.check_done('g3_', sleep_time=sleep_time)
        t2 = time.time()
        # Times are approximate since check_done only checks every N seconds for some N
        print
        print "Time for generation of images = ",t2-t1
        print

# Finally, we go back to a process per branch for the final steps: star_params and packages.
t1 = time.time()
prefix3 = 'g3_step3_'
for experiment, obs_type, shear_type in branches:
    e = experiment[0]
    o = obs_type[0]
    s = shear_type[0]

    pbs_name = prefix3+e+o+s
    pbs_file = pbs_name+'.sh'
    python_file = pbs_name+'.py'

    # Write out some scripts.
    mass_produce_utils.pbs_script_python(pbs_file, pbs_name)
    mass_produce_utils.python_script(python_file, root, subfield_min, subfield_max, experiment,
                                     obs_type, shear_type, gal_dir, ps_dir, seed,
                                     n_config_per_branch, preload, my_step=3)
    # And then submit them
    command_str = 'qsub '+pbs_file
    p = subprocess.Popen(command_str,shell=True, close_fds=True)
# The above command just submitted all the files to the queue.  We have to periodically poll the
# queue to see if they are still running.
mass_produce_utils.check_done('g3_step3', sleep_time=sleep_time)
t2 = time.time()
# Times are approximate since check_done only checks every N seconds for some N
print
print 'Time for great3sims.run star_params, packages = ',t2-t1
print
