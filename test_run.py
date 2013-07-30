#!/sw64/bin/python2.7

import time
import os
import subprocess
import sys
import shutil

sys.path.append('..')
import great3

# Set which branches to test...
experiments = [
    'control',
    #'real_gal',
    #'real_psf',
    #'multiepoch',
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
subfield_min = 199    # Start with the last regular field if we want to check out the deep fields.
subfield_max = 201    # Total of 3 sub-fields
data_dir = 'great3_fit_data'    # This should be set up as a sim-link to your Dropbox folder.
ps_dir = '../inputs/ps/tables' 
seed = 12345                    # Whatever.  (But not zero.)

# Clean up possible residue from previous runs.
shutil.rmtree(root, ignore_errors=True)

# Build catalogs, etc.
t1 = time.time()
great3.constants.nrows = 20
great3.constants.ncols = 20
great3.run(root, subfield_min=subfield_min, subfield_max=subfield_max,
           experiments=experiments, obs_type=obs_type, shear_type=shear_type,
           gal_dir=data_dir, ps_dir=ps_dir,
           seed=seed, steps=['metaparameters', 'catalogs', 'config']
)
t2 = time.time()
print
print 'Time for great3.run up to config = ',t2-t1
print

# build images using galsim_yaml
dirs = []
os.chdir(root)
for exp in experiments:
    e = exp[0]
    if e == 'r': e = exp[5]
    for obs in obs_type:
        o = obs[0]
        for shear in shear_type:
            s = shear[0]
            f = e + o + s + '.yaml'
            dirs.append( os.path.join(root,exp,obs,shear) )
            t1 = time.time()
            p = subprocess.Popen(['galsim_yaml',f,'-v1'])
            p.communicate() # wait until done
            t2 = time.time()
            print
            print 'Time for galsim_yaml',f,'= ',t2-t1
            print
os.chdir('..')

# Build images using great3.run
t1 = time.time()
great3.run(root, subfield_min=subfield_min, subfield_max=subfield_max,
           experiments=experiments, obs_type=obs_type, shear_type=shear_type,
           gal_dir=data_dir, ps_dir=ps_dir,
           seed=seed, steps=['images']
)
t2 = time.time()
print
print 'Time for great3.run images = ',t2-t1
print

# Check that images are the same.
print 'Checking diffs: (No output means success)'
sys.stdout.flush()
for dir in dirs:
    for i in range(subfield_min, subfield_max+1):
        f1 = os.path.join(dir,'image-%03d-0.fits'%i)
        f2 = os.path.join(dir,'yaml_image-%03d-0.fits'%i)
        p = subprocess.Popen(['diff',f1,f2],stderr=subprocess.STDOUT)
        p.communicate()

# Now package up the data that should be public
t1 = time.time()
great3.run(root, subfield_min=subfield_min, subfield_max=subfield_max,
           experiments=experiments, obs_type=obs_type, shear_type=shear_type,
           gal_dir=data_dir, seed=seed, steps=['packages']
)
t2 = time.time()
print
print 'Time for great3.run packages = ',t2-t1
print
