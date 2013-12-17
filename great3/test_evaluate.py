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


def get_variable_gtrue(experiment, obs_type, logger=None):
    """Get the full catalog of shears and positions for all fields.

    @return id, x, y, g1, g2
    """
    mapper = great3sims.mapper.Mapper(evaluate.TRUTH_DIR, experiment, obs_type, "variable")
    identifier = np.empty(
        (evaluate.NGALS_PER_SUBFIELD, evaluate.NSUBFIELDS_PER_FIELD, evaluate.NFIELDS), dtype=int)
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
            identifier[:, jsub, ifield] = truedata["ID"]

    # Then return
    return identifier, xx, yy, g1true, g2true

def get_variable_gsuffix(experiment, obs_type, suffix="_intrinsic", file_prefix="galaxy_catalog",
                         logger=None):
    """Get the full catalog of intrinsic "shears" and positions for all fields.

    Gets "g1"+suffix and "g2"+suffix from the subfield_catalog files.

    @return id, x, y, g1, g2
    """
    mapper = great3sims.mapper.Mapper(evaluate.TRUTH_DIR, experiment, obs_type, "variable")
    identifier = np.empty(
        (evaluate.NGALS_PER_SUBFIELD, evaluate.NSUBFIELDS_PER_FIELD, evaluate.NFIELDS), dtype=int)
    g1int = np.empty(
        (evaluate.NGALS_PER_SUBFIELD, evaluate.NSUBFIELDS_PER_FIELD, evaluate.NFIELDS))
    g2int = np.empty(
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
    ygrid_deg = ygrid_deg.flatten() # the shears are ordered too, which is handy
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
                mapper.full_dir, (file_prefix+"-%03d.fits" % isubfield_index))
            truedata = pyfits.getdata(galcatfile)
            if len(truedata) != evaluate.NGALS_PER_SUBFIELD:
                raise ValueError(
                    "Number of records in "+galcatfile+" (="+str(len(truedata))+") is not "+
                    "equal to NGALS_PER_SUBFIELD (="+str(evaluate.NGALS_PER_SUBFIELD)+")")
            g1int[:, jsub, ifield] = truedata["g1"+suffix]
            g2int[:, jsub, ifield] = truedata["g2"+suffix]
            identifier[:, jsub, ifield] = truedata["ID"]

    # Then return
    return identifier, xx, yy, g1int, g2int

def make_variable_submission(x, y, g1true, g2true, g1int, g2int, c1, c2, m1, m2, outfile,
                             noise_sigma=0.05):
    """Make a fake submission based on input x, y, true shears, bias and noise parameters.

    Saves to outfile in the format of Melanie's presubmission.py output.
    """
    # Some sanity checks to start
    if x.shape != (
        evaluate.NGALS_PER_SUBFIELD, evaluate.NSUBFIELDS_PER_FIELD, evaluate.NFIELDS):
        raise ValueError("x.shape does not match required dimensions.")
    if y.shape != x.shape: raise ValueError("y.shape does not match x.shape.")
    if g1true.shape != x.shape: raise ValueError("g1true.shape does not match x.shape.")
    if g2true.shape !=  x.shape: raise ValueError("g2true.shape does not match x.shape.")
    # First apply the chosen biases and some noise, use complex shears for easy addition 
    gbiasc = g1true * (1. + m1) + c1 + np.random.randn(*g1true.shape) * noise_sigma + \
        (g2true * (1. + m2) + c2 + np.random.randn(*g2true.shape) * noise_sigma) * 1j
    gintc = g1int * (1. + m1) + g2int * (1. + m2) * 1j  # m biases also affect intrinsic part
    gsubc = (gintc + gbiasc) / (1. + gbiasc.conj() * gintc)
    g1sub = gsubc.real
    g2sub = gsubc.imag
    # Define the field array, then theta and map arrays in which we'll store the results
    field = np.arange(evaluate.NBINS_THETA * evaluate.NFIELDS) / evaluate.NBINS_THETA
    theta = np.empty(evaluate.NBINS_THETA * evaluate.NFIELDS)
    map_E = np.empty(evaluate.NBINS_THETA * evaluate.NFIELDS)
    map_B = np.empty(evaluate.NBINS_THETA * evaluate.NFIELDS)
    maperr = np.empty(evaluate.NBINS_THETA * evaluate.NFIELDS)
    for ifield in range(evaluate.NFIELDS):

        # Extracting the x, y and g1, g2 for all the subfields in this field, flatten and use
        # to calculate the map_E
        map_results = evaluate.run_corr2(
            x[:, :, ifield].flatten(),
            y[:, :, ifield].flatten(),
            g1sub[:, :, ifield].flatten(),
            g2sub[:, :, ifield].flatten(),
            np.ones_like(x[:, :, ifield]).flatten(),
            min_sep=evaluate.THETA_MIN_DEG,
            max_sep=evaluate.THETA_MAX_DEG,
            nbins=evaluate.NBINS_THETA,
            xy_units="degrees",
            sep_units="degrees")
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

    import os
    import tempfile

    truth_dir = "/Users/browe/great3/beta/truth" # Modify to wherever truth is unpacked
    # Temporarily write this into evaluate.TRUTH_DIR - HACKY!
    evaluate.TRUTH_DIR = truth_dir

    # Set the experiment and observation type to test (both shear_types will be explored)
    experiment = 'control'
    obs_type = 'space'

    # Setup the logger
    logging.basicConfig(stream=sys.stderr)
    logger = logging.getLogger("test")
    logger.setLevel(logging.DEBUG)

    usebins = (evaluate.USEBINS, "subbins") #(evaluate.USEBINS, "subbins")
    poisson = (False, "noweight") 
    fractional = (False, "absdiffs")

    NTEST = 300
    NOISE_SIGMA = 0.10
    cvals = (evaluate.CFID,)# 10. * evaluate.CFID, 100. * evaluate.CFID) 
    mvals = (evaluate.MFID,)# 10. * evaluate.MFID, 100. * evaluate.MFID) 
    qarr = np.empty((NTEST, len(cvals), len(mvals)))

    print usebins[1]
    print poisson[1]
    print fractional[1]

    # Get the x,y, true intrinsic ellips and shears for making fake submissions
    _, x, y, g1true, g2true = get_variable_gtrue(experiment, obs_type)
    _, _, _, g1int, g2int = get_variable_gsuffix(experiment, obs_type)

    for ic, cval in enumerate(cvals):

        for jm, mval in enumerate(mvals):

            print
            print "Running metric simulations for c = %.4f, m = %.4f" % (cval, mval)
            # Perform a simualtion, and do up to NTEST trials, so as to help find an updated
            # normalization factor (do one outside loop so I can check logger action)
            #subfile = tempfile.mktemp(suffix=".dat")
            #result = make_variable_submission(
            #    x, y, g1true, g2true, g1int, g2int, cval, cval, mval, mval,
            #    outfile=subfile, noise_sigma=NOISE_SIGMA)
            #q = evaluate.q_variable(
            #    subfile, experiment, obs_type, logger=logger, usebins=usebins[0],
            #    poisson_weight=poisson[0], fractional_diff=fractional[0])
            #os.remove(subfile)
            #print "%3d/%3d: Q_v (c = %.4f, m = %.4f) = %.5e" % (1, NTEST, cval, mval, q)

            qlist = []
            for i in range(NTEST):
    
                subfile = tempfile.mktemp(suffix=".dat")
                result = make_variable_submission(
                    x, y, g1true, g2true, g1int, g2int, cval, cval, mval, mval,
                    outfile=subfile, noise_sigma=NOISE_SIGMA)
                q = evaluate.q_variable(
                    subfile, experiment, obs_type, logger=None, usebins=usebins[0],
                    poisson_weight=poisson[0], fractional_diff=fractional[0], truth_dir=truth_dir,
                    sigma_min=1.e-6)
                os.remove(subfile)
                print "%3d/%3d: Q_v (c = %.4f, m = %.4f) = %.5e" % (i + 1, NTEST, cval, mval, q)
                qlist.append(q)

            # Collate and print results
            qarr[:, ic, jm] = np.asarray(qlist)
            print "Mean of Q_v values = %.5e +/- %.5e" % (
                np.mean(qarr[:, ic, jm]), np.std(qarr[:, ic, jm]) / np.sqrt(len(qarr[:, ic, jm])))
            print "Std of Q_v values = %.5e +/- %.5e " % (
                np.std(qarr[:, ic, jm]),
                np.std(qarr[:, ic, jm]) / np.sqrt(2 * (len(qarr[:, ic, jm]) - 1)))
            print "Fractional uncertainty on Q_v values = %.5e +/- %.5e " % (
                np.std(qarr[:, ic, jm]) / np.mean(qarr[:, ic, jm]),
                np.std(qarr[:, ic, jm]) / np.sqrt((len(qarr[:, ic, jm]) - 1))
                / np.mean(qarr[:, ic, jm]))
            print "Best estimate of normalization factor = "+str(1000. / np.mean(qarr[:, ic, jm]))

    # Save the arrays
    filename = os.path.join(
        evaluate.STORAGE_DIR,
        "newmetric_sigma_min1.e-6_NOISE_SIGMA"+("%.2f" % NOISE_SIGMA)+"_"+usebins[1]+"_"+poisson[1]+"_"+
        fractional[1]+"_mc_N"+str(NTEST)+".npy")
    print "Saving to "+filename
    np.save(filename, qarr)
