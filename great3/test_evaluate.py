import sys
import os
import logging
import numpy as np
import evaluate
path, module = os.path.split(__file__)
sys.path.append(os.path.join(path, "..", "..", "metrics")) # Appends the great3-private/metrics
                                                           # folder to path
import g3metrics
sys.path.append("..", "..")
import great3sims.mapper


def get_variable_gtrue(experiment, obs_type):
    """Get the full catalog of shears and positions for all fields.

    @return x, y, g1, g2
    """
    mapper = great3sims.mapper.Mapper(truth_dir, experiment, obs_type, "variable")
    g1true = np.empty(
        (evaluate.NGALS_PER_SUBFIELD, evaluate.NSUBFIELDS_PER_FIELD, evaluate.NFIELDS))
    g2true = np.empty(
        (evaluate.NGALS_PER_SUBFIELD, evaluate.NSUBFIELDS_PER_FIELD, evaluate.NFIELDS))
    xx = np.empty(
        (evaluate.NGALS_PER_SUBFIELD, evaluate.NSUBFIELDS_PER_FIELD, evaluate.NFIELDS))
    yy = np.empty(
        (evaluate.NGALS_PER_SUBFIELD, evaluate.NSUBFIELDS_PER_FIELD, evaluate.NFIELDS))
    # Load the offsets
    subfield_indices, offset_deg_x, offset_deg_y = get_generate_variable_offsets(
        experiment, obs_type, storage_dir=evaluate.STORAGE_DIR, truth_dir=evaluate.TRUTH_DIR,
        logger=logger)
    # Then loop over the fields and subfields getting the galaxy catalogues
    import pyfits
    for ifield in range(NFIELDS):

        # Read in all the shears in this field and store
        for jsub in range(NSUBFIELDS_PER_FIELD):

            # Build the x,y grid using the subfield offsets
            isubfield_index = jsub + ifield * NSUBFIELDS_PER_FIELD
            xx[:, jsub, ifield] = xgrid_deg + offset_deg_x[isubfield_index]
            yx[:, jsub, ifield] = ygrid_deg + offset_deg_y[isubfield_index]
            galcatfile = os.path.join(
                mapper.full_dir, ("galaxy_catalog-%03d.fits" % isubfield_index))
            truedata = pyfits.getdata(galcatfile)
            if len(truedata) != evaluate.NGALS_PER_SUBFIELD:
                raise ValueError(
                    "Number of records in "+galcatfile+" (="+str(len(truedata))+") is not "+
                    "equal to NGALS_PER_SUBFIELD (="+str(evaluate.NGALS_PER_SUBFIELD)+")")
            g1true[:, jsub, ifield] = truedata["g1"]
            g2true[:, jsub, ifield] = truedata["g2"]
    # Then return
    return xx, yy, g1true, g2true



if __name__ == "__main__":

    experiment = 'contaol'
    obs_type = 'ground'

    logger = logging.getLogger("test")
    logger.setLevel(logging.DEBUG)

    # Just try getting / building the intermediate products for this branch first
    #sind, g1t, g2t = evaluate.get_generate_const_truth(experiment, obs_type, logger=logger)
    #gdict = evaluate.get_generate_const_subfield_dict(experiment, obs_type, logger=logger)
    #grot = evaluate.get_generate_const_rotations(experiment, obs_type, logger=logger)

    # Try a simple submission, no biases, and see what Q I get
    #label = "sub1"
    #g1sub, g2sub = g3metrics.make_submission_const_shear(
    #   0.001, 0., 0., -0.01, g1t, g2t, 1e4, 0.05, label=label, rotate_cs=grot)
    #subfile = "./g3subs/g3_const_shear_sub."+label+".dat"

    #q, c1, m1, c2, m2, sigc1, sigm1, sigc2, sigm2  = evaluate.q_constant(
    #    subfile, 'control', 'ground')
    #print "Q_c = "+str(q)
    #print "c+ = "+str(c1)+" ("+str(sigc1)+")"
    #print "m+ = "+str(m1)+" ("+str(sigm1)+")"
    #print "cx = "+str(c2)+" ("+str(sigc2)+")"
    #print "mx = "+str(m2)+" ("+str(sigm2)+")"

    #subfield_index, offset_deg_x, offset_deg_y = evaluate.get_generate_variable_offsets(
    #    experiment, obs_type, logger=logger)

    field, theta, map_E, map_B, maperr = evaluate.get_generate_variable_truth(
        experiment, obs_type, logger=logger)

    # Try a basically random variable shear submission from Melanie's code!
    q_v = evaluate.q_variable(
        "../../public-scripts/csv_test.asc", experiment, obs_type, logger=None)
    print "Q_v = "+str(q_v)

    # Then try making
