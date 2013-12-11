#!/usr/bin/env python
import os
import sys
import numpy as np
import pyfits
import constants
import count_objects
sys.path.append("..")
import great3sims.mapper

EXPERIMENT = "variable_psf"


def make_cats(experiment, obs_type, shear_type, catfile_prefix="./counts/sexcat_"):
    """Make a dict of SExtractor catalogs for all galaxy images in all branches of this experiment
    """
    import subprocess
    mapper = great3sims.mapper.Mapper(
        constants.public_dir, experiment, obs_type, shear_type)
    print "Making catalogs for galaxy images in "+mapper.full_dir
    cats = {}
    for subfield_index in range(200):

        imfile = os.path.join(mapper.full_dir, "image-%03d-0.fits" % subfield_index)
        catfile = catfile_prefix+experiment+"-"+obs_type+"-"+shear_type+\
            ("_image-%03d-0" % subfield_index)+".fits"
        if not os.path.isfile(catfile):
            print "Running SExtractor on "+imfile+" to build "+catfile
            if obs_type == "space":
                config_filename = "sex_space.config"
            else:
                config_filename = "sex.config"
            subprocess.check_call(
                ["sex", imfile, "-c", config_filename, "-CATALOG_NAME", catfile], close_fds=True)
        print "Loading data from "+catfile
        cats[imfile] = pyfits.getdata(catfile)
    # Return the cats dict
    return cats

if __name__ == "__main__":

    import cPickle
    # First of all make the counts output directory and assign the filename for catalog output
    if not os.path.isdir("./counts"): os.mkdir("./counts")
    for obs_type in constants.obs_types:

        for shear_type in constants.shear_types:

            cats = make_cats(EXPERIMENT, obs_type, shear_type)
            for imfile, cat in cats.iteritems():

                snr = cat["FLUX_AUTO"] / cat["FLUXERR_AUTO"]
                plt.hist(snr, bins=100, range=(0, 100), histtype="step")
                plt.xlabel("FLUX_AUTO / FLUXERR_AUTO")
                plt.ylabel("Counts")
                plt.title(EXPERIMENT+"/"+obs_type+"/"+shear_type+"/"+os.path.split(imfile)[-1])
                plt.show()
                1/0


