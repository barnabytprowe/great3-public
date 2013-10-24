# Utilities needed for mass-production of sims.

def pbs_script_python(filename, jobname):
    """A utility to write a PBS script to filename, with the given jobname for python.
    """
    with open(filename, "w") as f:

        f.write("#!/bin/sh\n")
        f.write("#PBS -l nodes=1:ppn=16\n")
        f.write("#PBS -q physics\n")
        f.write("#PBS -j oe\n")
        f.write("#PBS -N "+jobname+"\n")
        f.write("#PBS -l walltime=18:00:00\n")
        f.write("\n")
        f.write("cd $PBS_O_WORKDIR\n")
        f.write("/opt/python27/bin/python2.7 "+jobname+".py\n")

def pbs_script_yaml(filename, configname, rootname):
    """A utility to write a PBS script to filename, to use the given yaml configname.
    """
    with open(filename, "w") as f:
        import os
        jobname, _ = os.path.splitext(configname)
        jobname = 'g3_'+jobname

        f.write("#!/bin/sh\n")
        f.write("#PBS -l nodes=1:ppn=16\n")
        f.write("#PBS -q physics\n")
        f.write("#PBS -j oe\n")
        f.write("#PBS -N "+jobname+"\n")
        f.write("#PBS -l walltime=45:00:00\n")
        f.write("\n")
        f.write("cd "+rootname+"\n")
        f.write("/home/rmandelb/software/bin/galsim "+configname+" -v1\n")

def python_script(filename, root, subfield_min, subfield_max, experiment, obs_type, shear_type,
                  gal_dir, ps_dir, seed, n_config_per_branch, preload, my_step):
    """A utility to write a python script that just imports GREAT3 and does some steps of simulation
    generation.
    """
    import os
    if my_step == 1:
        with open(filename, "w") as f:
            f.write("import sys\n")
            f.write("import shutil\n")
            f.write("sys.path.append('/home/rmandelb/git/great3-private')\n")
            f.write("import great3sims\n")
            command_str = "great3sims.run('" + root + "', subfield_min=" + str(subfield_min) + \
                ", subfield_max=" + str(subfield_max) + ", experiments=['" + experiment + \
                "'], obs_type=['" + obs_type + "'], shear_type=['" + shear_type + \
                "'], gal_dir='" + gal_dir + "', ps_dir='" + ps_dir + "', seed=" + str(seed) + \
                ", public_dir='" + os.path.join(root, 'public') + \
                "', truth_dir='" + os.path.join(root, 'truth') + \
                "', steps = ['metaparameters', 'catalogs'], preload=" + str(preload) + ")\n"
            f.write(command_str)

            e = experiment[0]
            o = obs_type[0]
            s = shear_type[0]
            dir = os.path.join(root, experiment, obs_type, shear_type)
            config_pref = e + o + s
            psf_config_pref = e + o + s + '_psf'
            star_test_config_pref = e + o + s + '_star_test'

            new_config_names = []
            new_psf_config_names = []
            new_star_test_config_names = []
            for i in range(n_config_per_branch):
                first = subfield_min + (subfield_max-subfield_min+1)/n_config_per_branch * i
                last = subfield_min + (subfield_max-subfield_min+1)/n_config_per_branch * (i+1) - 1
                # Could put nproc setting here.  However, for now we just want to use as many 
                # processors as we have on that node.
                nproc = -1
                command_str = "great3sims.run('" + root + "', subfield_min=" + str(first) + \
                    ", subfield_max=" + str(last) + ", experiments=['" + experiment + \
                    "'], obs_type=['" + obs_type + "'], shear_type=['" + shear_type + \
                    "'], gal_dir='" + gal_dir + "', ps_dir='" + ps_dir + "', seed=" + str(seed) + \
                    ", public_dir='" + os.path.join(root, 'public') + \
                    "', truth_dir='" + os.path.join(root, 'truth') + \
                    "', steps = ['config'], nproc=" + str(nproc) + ", preload=" + str(preload) + \
                    ")\n"
                f.write(command_str)
                new_name = '%s_%02d.yaml'%(config_pref,i)
                new_psf_name = '%s_%02d.yaml'%(psf_config_pref,i)
                new_config_names.append(new_name)
                new_psf_config_names.append(new_psf_name)
                f.write("shutil.move('"+os.path.join(root,config_pref+'.yaml')+"', '"+
                        os.path.join(root,new_name)+"')\n")
                f.write("shutil.move('"+os.path.join(root,psf_config_pref+'.yaml')+"', '"+
                        os.path.join(root,new_psf_name)+"')\n")

            # But don't make multiple star test configs, just one.  We have to rerun great3sims.run
            # for this to get the config file for everything.  It will also remake config files for
            # galaxy and starfield images, but since we won't add their names to our list of config
            # files to run, they just get ignored, and it's all good.
            command_str = "great3sims.run('" + root + "', subfield_min=" + str(subfield_min) + \
                ", subfield_max=" + str(subfield_max) + ", experiments=['" + experiment + \
                "'], obs_type=['" + obs_type + "'], shear_type=['" + shear_type + \
                "'], gal_dir='" + gal_dir + "', ps_dir='" + ps_dir + "', seed=" + str(seed) + \
                ", public_dir='" + os.path.join(root, 'public') + \
                "', truth_dir='" + os.path.join(root, 'truth') + \
                "', steps = ['config'], preload=" + str(preload) + ")\n"
            f.write(command_str)
            new_star_test_name = '%s.yaml'%(star_test_config_pref)
            new_star_test_config_names.append(new_star_test_name)

        return new_config_names, new_psf_config_names, new_star_test_config_names

    elif my_step == 3:
        with open(filename, "w") as f:
            f.write("import sys\n")
            f.write("sys.path.append('/home/rmandelb/git/great3-private')\n")
            f.write("import great3sims\n")
            command_str = "great3sims.run('" + root + "', subfield_min=" + str(subfield_min) + \
                ", subfield_max=" + str(subfield_max) + ", experiments=['" + experiment + \
                "'], obs_type=['" + obs_type + "'], shear_type=['" + shear_type + \
                "'], gal_dir='" + gal_dir + "', ps_dir='" + ps_dir + "', seed=" + str(seed) + \
                ", public_dir='" + os.path.join(root, 'public') + \
                "', truth_dir='" + os.path.join(root, 'truth') + \
                "', steps = ['star_params', 'packages'], preload=" + str(preload) + ")\n"
            f.write(command_str)
    else:
        raise NotImplementedError

def check_done(process_str, sleep_time=60):
    """Check the PBS queue to see if there are any processes that contain process_str in the name.
    """
    import time
    import subprocess
    ind_found = 1000
    while ind_found > 0:
        if ind_found < 1000:
            # only sleep before checking if it's not the first time through this while loop
            time.sleep(sleep_time)
        res = subprocess.check_output('qstat')
        ind_found = res.find(process_str)

def check_njobs(process_str, sleep_time=60, n_jobs=10):
    """Check the PBS queue to see how many processes that contain process_str in the name; wait
    until it is <n_jobs before returning.
    """
    import time
    import subprocess
    import re
    n_found = 1000
    while n_found >= n_jobs:
        if n_found < 1000:
            # only sleep before checking if it's not the first time through this while loop
            time.sleep(sleep_time)
        res = re.findall(process_str, subprocess.check_output('qstat'))
        n_found = len(res)
