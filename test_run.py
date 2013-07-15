#!/sw64/bin/python2.7

import time
import os
import subprocess
import sys
import shutil

sys.path.append('..')
import great3

root = 'test_run'
n_subfields = 2                 # I keep this small so the test doesn't take too long.
experiments = ['control']       # For now..
obs_type = 'space'              # For now..
shear_type = 'constant'         # For now..
data_dir = 'great3_fit_data'    # This should be set up as a sim-link to your Dropbox folder.
seed = 12345                    # Whatever.


# Clean up possible residue from previous runs.
shutil.rmtree(root, ignore_errors=True)

# Build catalogs, etc.
t1 = time.time()
great3.run(root, subfield_max=n_subfields,
           experiments=experiments, obs_type=obs_type, shear_type=shear_type,
           dir=data_dir, seed=seed, steps=['metaparameters', 'catalogs', 'config']
)
t2 = time.time()
print
print 'Time for great3.run up to config = ',t2-t1
print

# Build images using great3.run
t1 = time.time()
great3.run(root, subfield_max=n_subfields,
           experiments=experiments, obs_type=obs_type, shear_type=shear_type,
           dir=data_dir, seed=seed, steps=['images']
)
t2 = time.time()
print
print 'Time for great3.run images = ',t2-t1
print

# Build images using galsim_yaml
os.chdir(root)
t1 = time.time()
p = subprocess.Popen(['galsim_yaml','csc.yaml','-v1'])
p.communicate() #wait until done
t2 = time.time()
print
print 'Time for galsim_yaml csc.yaml = ',t2-t1
print

# Check that images are the same.
print 'Checking diffs: (No output means success)'
for i in range(n_subfields):
    f1 = 'control/space/constant/image-%03d-0.fits'%i
    f2 = 'control/space/constant/yaml_image-%03d-0.fits'%i
    p = subprocess.Popen(['diff',f1,f2])
    p.communicate()

