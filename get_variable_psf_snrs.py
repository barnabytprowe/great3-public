#!/usr/bin/env python
import os
import sys
import numpy as np
import constants
import count_objects
sys.path.append("..")
import great3sims.mapper

EXPERIMENT = "variable_psf"


def make_cat_dict(experiment):
    """Make a dict of SExtractor catalogs for all galaxy images in all branches of this experiment
    """
    print "Beginning generation of catalogs to be stored into "+catfile

    # Setup storage dictionary
    cats = {}
    for obs_type in constants.obs_types:

        for shear_type in constants.shear_types:

            mapper = great3sims.mapper.Mapper(
                constants.public_dir, experiment, obs_type, shear_type)
            for subfield_index in range(200):

                imfile = os.path.join(mapper.full_dir, "image-%03d-0.fits" % subfield_index)
                print "Running SExtractor on "+imfile
                if obs_type == "space":
                    config_filename = "sex_space.config"
                else:
                    config_filename = "sex.config"
                cats[experiment+"-"+obs_type+"-"+shear_type+"/"+(os.path.split(imfile)[-1])] = \
                    count_objects.file_count(imfile, config_filename=config_filename)
    # Return the filled dict
    return cats


if __name__ == "__main__":

    # First of all make the counts output directory and assign the filename for catalog output
    if not os.path.isdir("./counts"): os.mkdir("./counts")
    catfile = os.path.join(".counts", "variable_psf_sexcats.p")
    if not os.path.isfile(catfile):
        make_cat_dict(EXPERIMENT)


