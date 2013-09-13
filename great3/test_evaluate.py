"""@file test_evaluate.py

This file started off as a simple testbed for the functions in the evaluate.py module, but has
morphed into that plus a calculation of the appropriate normalization for variable shear branches
given the new geometry. 
"""

import sys
import os
import logging
import numpy as np
import evaluate
path, module = os.path.split(__file__)
sys.path.append(os.path.join(path, "..", "..", "metrics")) # Appends the great3-private/metrics
                                                           # folder to path
import g3metrics
sys.path.append(("..", ".."))
import great3sims.mapper


def get_variable_gtrue(experiment, obs_type):
    """Get the full catalog of shears and positions for all fields.

    @return x, y, g1, g2
    """
    mapper = great3sims.mapper.Mapper(evaluate.TRUTH_DIR, experiment, obs_type, "variable")
    g1true = np.empty(
        (evaluate.NGALS_PER_SUBFIELD, evaluate.NSUBFIELDS_PER_FIELD, evaluate.NFIELDS))
    g2true = np.empty(
        (evaluate.NGALS_PER_SUBFIELD, evaluate.NSUBFIELDS_PER_FIELD, evaluate.NFIELDS))
    xx = np.empty(
        (evaluate.NGALS_PER_SUBFIELD, evaluate.NSUBFIELDS_PER_FIELD, evaluate.NFIELDS))
    yy = np.empty(
        (evaluate.NGALS_PER_SUBFIELD, evaluate.NSUBFIELDS_PER_FIELD, evaluate.NFIELDS))
    # Load the offsets
    subfield_indices, offset_deg_x, offset_deg_y = evaluate.get_generate_variable_offsets(
        experiment, obs_type, storage_dir=evaluate.STORAGE_DIR, truth_dir=evaluate.TRUTH_DIR,
        logger=logger)
    # Build basic x and y grids to use for coord positions
    xgrid_deg, ygrid_deg = np.meshgrid(
        np.arange(0., evaluate.XMAX_GRID_DEG, evaluate.DX_GRID_DEG),
        np.arange(0., evaluate.XMAX_GRID_DEG, evaluate.DX_GRID_DEG))
    xgrid_deg = xgrid_deg.flatten() # Flatten these - the default C ordering corresponds to the way
    ygrid_deg = ygrid_deg.flatten() # the true shears are ordered too, which is handy
    if len(xgrid_deg) != evaluate.NGALS_PER_SUBFIELD:
        raise ValueError(
            "Dimensions of xgrid_deg and ygrid_deg do not match NGALS_PER_SUBFIELD.  Please check "+
            "the values of XMAX_GRID_DEG and DX_GRID_DEG in evaluate.py.")
    # Then loop over the fields and subfields getting the galaxy catalogues
    import pyfits
    for ifield in range(evaluate.NFIELDS):

        # Read in all the shears in this field and store
        for jsub in range(evaluate.NSUBFIELDS_PER_FIELD):

            # Build the x,y grid using the subfield offsets
            isubfield_index = jsub + ifield * evaluate.NSUBFIELDS_PER_FIELD
            xx[:, jsub, ifield] = xgrid_deg + offset_deg_x[isubfield_index]
            yy[:, jsub, ifield] = ygrid_deg + offset_deg_y[isubfield_index]
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

def make_variable_submission(x, y, g1true, g2true, c1, c2, m1, m2, outfile, noise_sigma=0.05):
    """Make a fake submission based on input x, y, true shears, bias and noise parameters.

    Saves to outfile.
    """
    # Some sanity checks to start
    if x.shape != (
        evaluate.NGALS_PER_SUBFIELD, evaluate.NSUBFIELDS_PER_FIELD, evaluate.NFIELDS):
        raise ValueError("x.shape does not match required dimensions.")
    if y.shape != x.shape: raise ValueError("y.shape does not match x.shape.")
    if g1true.shape != x.shape: raise ValueError("g1true.shape does not match x.shape.")
    if g2true.shape !=  x.shape: raise ValueError("g2true.shape does not match x.shape.")
    # Then apply the chosen biases 
    g1sub = g1true * (1. + m1) + c1 + np.random.randn(*g1true.shape) * noise_sigma 
    g2sub = g2true * (1. + m2) + c2 + np.random.randn(*g2true.shape) * noise_sigma
    # Define the field array, then theta and map arrays in which we'll store the results
    field = np.arange(evaluate.NBINS_THETA * evaluate.NFIELDS) / evaluate.NBINS_THETA
    theta = np.empty(evaluate.NBINS_THETA * evaluate.NFIELDS)
    map_E = np.empty(evaluate.NBINS_THETA * evaluate.NFIELDS)
    map_B = np.empty(evaluate.NBINS_THETA * evaluate.NFIELDS)
    maperr = np.empty(evaluate.NBINS_THETA * evaluate.NFIELDS)
    for ifield in range(evaluate.NFIELDS):

        # Extracting the x, y and g1, g2 for all the subfields in this field, flatten and use
        # to calculate the map_E
        map_results = g3metrics.run_corr2(
            x[:, :, ifield].flatten(), y[:, :, ifield].flatten(), g1sub[:, :, ifield].flatten(),
            g2sub[:, :, ifield].flatten(), min_sep=evaluate.THETA_MIN_DEG,
            max_sep=evaluate.THETA_MAX_DEG, nbins=evaluate.NBINS_THETA,
            params_file="./corr2.params", xy_units="degrees", sep_units="degrees")
        theta[ifield * evaluate.NBINS_THETA: (ifield + 1) * evaluate.NBINS_THETA] = \
            map_results[:, 0]
        map_E[ifield * evaluate.NBINS_THETA: (ifield + 1) * evaluate.NBINS_THETA] = \
            map_results[:, 1]
        map_B[ifield * evaluate.NBINS_THETA: (ifield + 1) * evaluate.NBINS_THETA] = \
            map_results[:, 2]
        maperr[ifield * evaluate.NBINS_THETA: (ifield + 1) * evaluate.NBINS_THETA] = \
            map_results[:, 5]

    # Finally save in ASCII format
    with open(outfile, "wb") as fout:
        fout.write(
            "# Simulated aperture mass statistics for "+experiment+"-"+obs_type+"-variable\n")
        fout.write("# field_index  theta [deg]  map_E  map_B  maperr\n")
        np.savetxt(
            fout, np.array((field, theta, map_E, map_B, maperr)).T,
            fmt=" %2d %.18e %.18e %.18e %.18e")
    # Job done, return
    return


if __name__ == "__main__":

    experiment = 'control'
    obs_type = 'space'

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
    #q_v = evaluate.q_variable(
    #    "../../public-scripts/csv_test.asc", experiment, obs_type, logger=None)
    #print "Q_v (from presubmission) = "+str(q_v)

    # Then try making a fake submission ourselves
    x, y, g1true, g2true = get_variable_gtrue(experiment, obs_type)
    result = make_variable_submission(x, y, g1true, g2true, 5.e-3, 2.e-4, 0.01, -0.01,
        outfile="junk_test.asc")
    q_biased = evaluate.q_variable("junk_test.asc", experiment, obs_type)
    print "Q_v (from own biased submission simulator) = "+str(q_biased)

    result = make_variable_submission(
        x, y, g1true, g2true, evaluate.CFID, evaluate.CFID, evaluate.MFID, evaluate.MFID, outfile="junk_test.asc")
    q_v2 = evaluate.q_variable("junk_test.asc", experiment, obs_type)
    print "Q_v (from own fiducial submission simulator) = "+str(q_v2)

    qlist = [q_v2]
    for i in range(99):
    
        result = make_variable_submission(
            x, y, g1true, g2true, evaluate.CFID, evaluate.CFID, evaluate.MFID, evaluate.MFID,
            outfile="junk_test.asc")

        q_v2 = evaluate.q_variable("junk_test.asc", experiment, obs_type)
        print "Q_v (from own fiducial submission simulator: "+str(i + 2)+"/100) = "+str(q_v2)
        qlist.append(q_v2)

    qarr = np.asarray(qlist)
    print "Mean of Q_v values = "+str(np.mean(qarr))+"+/-"+str(np.std(qarr) / np.sqrt(len(qarr)))
