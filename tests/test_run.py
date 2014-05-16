#!/sw64/bin/python2.7
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
"""Example script that can be used to run small amounts of GREAT3-like simulations (typically with
fewer galaxies to save time).  This is useful for someone who is trying out the simulation code for
the first time, or making small modifications that they want to test out before mass production."""

import time
import os
import subprocess
import sys
import shutil
import numpy

sys.path.append('..')
import great3sims

# Set which branches to test.  Currently, it does all 20.
experiments = [
    'variable_psf',
    'control',
    'multiepoch',
    'real_galaxy',
    'full',
]
obs_type = [
    'ground',
    'space',
]
shear_type = [
    'constant',
    'variable',
]

# Here are a set of variables that users might need to change depending on what they want to do:
# Define a root directory for the simulations to be produced.
root = 'test_run'
# Decide how many config files to run separately.
n_config = 2
# Decide on a minimum and maximum subfield to generate.  These will be split up into `n_config`
# config files.
subfield_min = 0
subfield_max = 1
# `data_dir` should be set up as a sim-link to the folder containing the catalogs, fits, and actual
# images.
data_dir = './great3_data'
# Choose a random seed.  Can be whatever you want, but if you choose 0 then GalSim will choose a
# seed based on the current time, so outputs will not be deterministic.
seed = 12345

# Now make some decisions about what parts of the simulation process to do:
do_catalogs = True    # Remake the catalogs?
do_images = False     # Make the images in great3 builder code?
do_config = True      # Do config-related steps?
do_final = True       # Do final packaging steps?
preload_real = False  # preload images for RealGalaxy branches?  [always False for parametric branches]

# Choose how many processes to launch simultaneously.  -1 means to let GalSim decide for itself
# (based on the parameters of the system).
nproc = -1

# Note: these definitions have to happen up front.  The purpose is to override the settings in
# great3sims.constants so that the images that are produced are smaller (fewer objects), so that our
# tests run quickly.  These definitions affect image generation in addition to catalog generation.
# Reduce the number of galaxies so it won't take so long.
great3sims.constants.nrows = 10
great3sims.constants.ncols = 10
# Likewise, reduce the number of stars so it won't take so long.
great3sims.constants.min_star_density = 0.01 # per arcmin^2
great3sims.constants.max_star_density = 0.03

# All kwargs for the run command that don't need to change are below.
kwargs = {
    'root' : root,
    'experiments' : experiments,
    'obs_type' : obs_type,
    'shear_type' : shear_type,
    'gal_dir' : data_dir,
    'seed' : seed,
    'preload' : preload_real,
    'nproc' : nproc
}

# Build catalogs, etc.
if do_catalogs:
    # Clean up possible residue from previous runs.
    shutil.rmtree(root, ignore_errors=True)

    t1 = time.time()
    great3sims.run(steps=['metaparameters', 'catalogs'], 
                   subfield_min=subfield_min, subfield_max=subfield_max, **kwargs)
    t2 = time.time()
    print
    print 'Time for great3sims.run up to catalogs = ',t2-t1
    print

# Make list of config file names and directories:
dirs = []
n_epochs = []
config_names = []
psf_config_names = []
star_test_config_names = []
for exp in experiments:
    e = exp[0]
    for obs in obs_type:
        o = obs[0]
        for shear in shear_type:
            s = shear[0]
            dirs.append( os.path.join(root,exp,obs,shear) )
            config_names.append(e + o + s + '.yaml')
            psf_config_names.append(e + o + s + '_psf.yaml')
            star_test_config_names.append(e + o + s + '_star_test.yaml')
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
        great3sims.run(steps=['config'], 
                       subfield_min=first, subfield_max=last, **kwargs)
        for (old_names, new_names) in [ (config_names, new_config_names) ,
                                        (psf_config_names, new_psf_config_names) ]:
            for old_name in old_names:
                base, ext = os.path.splitext(old_name)
                new_name = '%s_%02d.yaml'%(base,i)
                print old_name,'->',new_name
                new_names.append(new_name)
                shutil.copy(os.path.join(root,old_name),os.path.join(root,new_name))
    t2 = time.time()
    print
    print 'Time for great3sims.run config = ',t2-t1
    print

    # build images using galsim executable
    t1 = time.time()
    os.chdir(root)
    for new_names in [ new_config_names, new_psf_config_names, star_test_config_names ]:
        for name in new_names:
            t3 = time.time()
            p = subprocess.Popen(['galsim',name,'-v1'],close_fds=True)
            p.communicate() # wait until done
            t4 = time.time()
            print 'Time for galsim',name,'= ',t4-t3
            print
    os.chdir('..')
    t2 = time.time()
    print 'Total time for galsim = ',t2-t1
    print

    if do_images:
        # Move these files to a different name, so we can compare to the ones built in the next 
        # step.
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
if do_images:
    t1 = time.time()
    great3sims.run(steps=['gal_images', 'psf_images'],
                   subfield_min=subfield_min, subfield_max=subfield_max, **kwargs)
    t2 = time.time()
    print
    print 'Time for great3sims.run images = ',t2-t1
    print

# Helper function to show what the difference is between two image files
def showdiff(f1,f2):
    import pyfits
    import numpy
    a1 = pyfits.open(f1)
    a2 = pyfits.open(f2)
    d1 = a1[0].data
    d2 = a2[0].data
    diff = numpy.where(d1 != d2)
    ndiff = len(diff[0])
    if ndiff > 0:
        print '    The images differ at %d pixel%s.'%(ndiff, ('s' if ndiff>1 else ''))
        maxdiff = numpy.max( d1[diff] - d2[diff] )
        print '    The maximum pixel difference is %e.'%maxdiff
    else:
        print '    Weird.  The images have identical pixel values.'
    a1.close()
    a2.close()

# Check that images are the same.
if do_images and do_config:
    print 'Checking diffs: (No output means success)'
    sys.stdout.flush()
    for ind in range(len(dirs)):
        dir = dirs[ind]
        n_epoch = n_epochs[ind]
        for i in range(subfield_min, subfield_max+1):
            for j in range(0, n_epoch):
                f1 = os.path.join(dir,'image-%03d-%1d.fits'%(i,j))
                f2 = os.path.join(dir,'yaml_image-%03d-%1d.fits'%(i,j))
                p = subprocess.Popen(['diff',f1,f2],stderr=subprocess.STDOUT,close_fds=True)
                p.communicate()
                if p.returncode != 0: showdiff(f1,f2)
                f1 = os.path.join(dir,'starfield_image-%03d-%1d.fits'%(i,j))
                f2 = os.path.join(dir,'yaml_starfield_image-%03d-%1d.fits'%(i,j))
                p = subprocess.Popen(['diff',f1,f2],stderr=subprocess.STDOUT,close_fds=True)
                p.communicate()
                if p.returncode != 0: showdiff(f1,f2)
    print 'End diffs.'
    print

if do_final:
    # Measure star parameters required for metric
    t1 = time.time()
    great3sims.run(steps=['star_params'],
                   subfield_min=subfield_min, subfield_max=subfield_max, **kwargs)
    t2 = time.time()
    print
    print 'Time for great3sims.run star_params = ',t2-t1
    print


    # Now package up the data that should be public, and truth tables
    t1 = time.time()
    great3sims.run(steps=['packages'],
                   subfield_min=subfield_min, subfield_max=subfield_max, **kwargs)
    t2 = time.time()
    print
    print 'Time for great3sims.run packages = ',t2-t1
    print
