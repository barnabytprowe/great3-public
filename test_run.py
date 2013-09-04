#!/sw64/bin/python2.7

import time
import os
import subprocess
import sys
import shutil
import numpy

sys.path.append('..')
import great3sims

# Set which branches to test...
experiments = [
    'variable_psf',
    'control',
    'multiepoch',
    #'real_gal',
    #'full',
]
obs_type = [
    'ground',
    'space',
]
shear_type = [
    'constant',
    'variable',
]

root = 'test_run'
n_config = 3        # Build 3 separate config files to be run separately.
subfield_min = 0
subfield_max = 5    # Total of 3x2 sub-fields  (2 per config file)
data_dir = 'great3_fit_data'    # This should be set up as a sim-link to your Dropbox folder.
ps_dir = '../inputs/ps/tables' 
seed = 12345                    # Whatever.  (But not zero.)
do_config = False   # Do config-related steps?  For now, this is False, because they are not set up
                    # properly for variable PSF branches.

# Clean up possible residue from previous runs.
shutil.rmtree(root, ignore_errors=True)

# Build catalogs, etc.
t1 = time.time()
# Reduce number of galaxies so it won't take so long
great3sims.constants.nrows = 10
great3sims.constants.ncols = 10
# Reduce number of stars so it won't take so long
great3sims.constants.min_star_density = 0.01 # per arcmin^2
great3sims.constants.max_star_density = 0.03
great3sims.run(root, subfield_min=subfield_min, subfield_max=subfield_max,
               experiments=experiments, obs_type=obs_type, shear_type=shear_type,
               gal_dir=data_dir, ps_dir=ps_dir,
               seed=seed, steps=['metaparameters', 'catalogs']
)
t2 = time.time()
print
print 'Time for great3sims.run up to catalogs = ',t2-t1
print

# Make list of config file names and directories:
dirs = []
n_epochs = []
config_names = []
psf_config_names = []
for exp in experiments:
    e = exp[0]
    for obs in obs_type:
        o = obs[0]
        for shear in shear_type:
            s = shear[0]
            dirs.append( os.path.join(root,exp,obs,shear) )
            config_names.append(e + o + s + '.yaml')
            psf_config_names.append(e + o + s + '_psf.yaml')
            if exp == "multiepoch" or exp == "full":
                n_epochs.append(great3sims.constants.n_epochs)
            else:
                n_epochs.append(1)

# Build config files
if do_config:
    t1 = time.time()
    new_config_names = []
    new_psf_config_names = []
    for i in range(n_config):
        first = subfield_min + (subfield_max-subfield_min+1)/n_config * i
        last = subfield_min + (subfield_max-subfield_min+1)/n_config * (i+1) - 1
        great3sims.run(root, subfield_min=first, subfield_max=last,
                       experiments=experiments, obs_type=obs_type, shear_type=shear_type,
                       gal_dir=data_dir, ps_dir=ps_dir,
                       seed=seed, steps=['config']
                       )
        for config_name in config_names:
            name, ext = os.path.splitext(config_name)
            new_name = '%s_%02d.yaml'%(name,i)
            new_config_names.append(new_name)
            shutil.copy(os.path.join(root,config_name),os.path.join(root,new_name))
        for config_name in psf_config_names:
            name, ext = os.path.splitext(config_name)
            new_name = '%s_%02d_psf.yaml'%(name,i)
            new_psf_config_names.append(new_name)
            shutil.copy(os.path.join(root,config_name),os.path.join(root,new_name))
    t2 = time.time()
    print
    print 'Time for great3sims.run config = ',t2-t1
    print

    # build images using galsim executable
    t1 = time.time()
    os.chdir(root)
    for config_name in new_config_names:
        t3 = time.time()
        p = subprocess.Popen(['galsim',config_name,'-v1'])
        p.communicate() # wait until done
        t4 = time.time()
        print 'Time for galsim',config_name,'= ',t4-t3
        print
    for config_name in new_psf_config_names:
        t3 = time.time()
        p = subprocess.Popen(['galsim',config_name,'-v1'])
        p.communicate() # wait until done
        t4 = time.time()
        print 'Time for galsim',config_name,'= ',t4-t3
        print
    os.chdir('..')
    t2 = time.time()
    print 'Total time for galsim = ',t2-t1
    print

    # Move these files to a different name, so we can compare to the ones built in the next step.
    for ind in range(len(dirs)):
        dir = dirs[ind]
        n_epoch = n_epochs[ind]
        for i in range(subfield_min, subfield_max+1):
            for j in range(0, n_epoch):
                f1 = os.path.join(dir,'image-%03d-%1d.fits'%(i,j))
                f2 = os.path.join(dir,'yaml_image-%03d-%1d.fits'%(i,j))
                shutil.move(f1,f2)
                f1 = os.path.join(dir,'starfield_image-%03d-%1d.fits'%(i,j))
                f2 = os.path.join(dir,'yaml_starfield_image-%03d-%1d.fits'%(i,j))
                shutil.move(f1,f2)

# Build images using great3sims.run
t1 = time.time()
great3sims.run(root, subfield_min=subfield_min, subfield_max=subfield_max,
               experiments=experiments, obs_type=obs_type, shear_type=shear_type,
               gal_dir=data_dir, ps_dir=ps_dir,
               seed=seed, steps=['gal_images', 'psf_images']
)
t2 = time.time()
print
print 'Time for great3sims.run images = ',t2-t1
print

# Check that images are the same.
if do_config:
    print 'Checking diffs: (No output means success)'
    sys.stdout.flush()
    for ind in range(len(dirs)):
        dir = dirs[ind]
        n_epoch = n_epochs[ind]
        for i in range(subfield_min, subfield_max+1):
            for j in range(0, n_epoch):
                f1 = os.path.join(dir,'image-%03d-%1d.fits'%(i,j))
                f2 = os.path.join(dir,'yaml_image-%03d-%1d.fits'%(i,j))
                p = subprocess.Popen(['diff',f1,f2],stderr=subprocess.STDOUT)
                p.communicate()
                f1 = os.path.join(dir,'starfield_image-%03d-%1d.fits'%(i,j))
                f2 = os.path.join(dir,'yaml_starfield_image-%03d-%1d.fits'%(i,j))
                p = subprocess.Popen(['diff',f1,f2],stderr=subprocess.STDOUT)
                p.communicate()
    print 'End diffs.'
    print

# Measure star parameters required for metric
t1 = time.time()
great3sims.run(root, subfield_min=subfield_min, subfield_max=subfield_max,
               experiments=experiments, obs_type=obs_type, shear_type=shear_type,
               gal_dir=data_dir, ps_dir=ps_dir,
               seed=seed, steps=['star_params']
)
t2 = time.time()
print
print 'Time for great3sims.run star_params = ',t2-t1
print


# Now package up the data that should be public, and truth tables
t1 = time.time()
great3sims.run(root, subfield_min=subfield_min, subfield_max=subfield_max,
               experiments=experiments, obs_type=obs_type, shear_type=shear_type,
               gal_dir=data_dir, seed=seed, steps=['packages']
)
t2 = time.time()
print
print 'Time for great3sims.run packages = ',t2-t1
print
