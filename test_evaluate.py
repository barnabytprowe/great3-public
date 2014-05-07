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
"""@file test_evaluate.py

This file started off as a simple testbed for the functions in the evaluate.py module, but has
morphed into that plus a calculation of the appropriate normalization for variable shear branches
given the new geometry. 
"""

import sys
import os
import logging
import numpy as np
path, module = os.path.split(__file__)
sys.path.append(os.path.join(path, "..", "server", "great3")) # Appends the folder
                                                              # great3-private/server/great3 to
                                                              # sys.path

import evaluate
import g3metrics
sys.path.append(("..")) # Add the parent directory to path so as to load up the great3sims.Mapper
import great3sims.mapper


def get_variable_gtrue(experiment, obs_type, logger=None, truth_dir=evaluate.TRUTH_DIR):
    """Get the full catalog of shears and positions for all fields.

    @return id, x, y, g1, g2
    """
    mapper = great3sims.mapper.Mapper(truth_dir, experiment, obs_type, "variable")
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
        experiment, obs_type, storage_dir=evaluate.STORAGE_DIR, truth_dir=truth_dir,
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
                         logger=None, truth_dir=evaluate.TRUTH_DIR):
    """Get the full catalog of intrinsic "shears" and positions for all fields.

    Gets "g1"+suffix and "g2"+suffix from the subfield_catalog files.

    @return id, x, y, g1, g2
    """
    mapper = great3sims.mapper.Mapper(truth_dir, experiment, obs_type, "variable")
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
        experiment, obs_type, storage_dir=evaluate.STORAGE_DIR, truth_dir=truth_dir,
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

        try:
            hstring = \
                "# Simulated aperture mass statistics for "+experiment+"-"+obs_type+"-variable\n"
            fout.write(hstring)
        except NameError: # If called externally then experiment and obs_type aren't defined.
                          # Rather than recode all former uses of this func, simply handle the error
                          # silently and move along...
            pass
        fout.write("# field_index  theta [deg]  map_E  map_B  maperr\n")
        np.savetxt(
            fout, np.array((field, theta, map_E, map_B, maperr)).T,
            fmt=" %2d %.18e %.18e %.18e %.18e")
    # Job done, return
    return

def make_unitc(experiment, obs_type, truth_dir=evaluate.TRUTH_DIR):
    """Make a variable submission to learn what a pure c1 or c2 looks like in map^2
    """
    import tempfile
    # Get the x,y, true intrinsic ellips and shears for making fake submissions
    _, x, y, g1true, g2true = get_variable_gtrue(experiment, obs_type, truth_dir=truth_dir)
    _, _, _, g1int, g2int = get_variable_gsuffix(experiment, obs_type, truth_dir=truth_dir)
    fdtmp, tmpfile = tempfile.mkstemp(suffix=".dat")
    result = make_variable_submission(
        x, y, np.zeros_like(g1true), np.zeros_like(g2true), np.zeros_like(g1true),
        np.zeros_like(g2true), 1., 0., 0., 0., outfile=tmpfile, noise_sigma=0.)
    os.close(fdtmp)
    map_E_c1 = np.loadtxt(tmpfile)[:, 2]
    os.remove(tmpfile)
    fdtmp, tmpfile = tempfile.mkstemp(suffix=".dat")
    result = make_variable_submission(
        x, y, np.zeros_like(g1true), np.zeros_like(g2true), np.zeros_like(g1true),
        np.zeros_like(g2true), 0., 1., 0., 0., outfile=tmpfile, noise_sigma=0.)
    os.close(fdtmp)
    map_E_c2 = np.loadtxt(tmpfile)[:, 2]
    os.remove(tmpfile)
    # Then average these (I checked they are very similar as you would expect) to get our "unit c"
    # term for the modelling.
    map_E_unitc = .5 * (map_E_c1 + map_E_c2)
    return map_E_unitc


if __name__ == "__main__":

    import os
    import tempfile

    truth_dir = "/Users/browe/great3/beta/truth" # Modify to wherever truth is unpacked
    # Temporarily write this into truth_dir - HACKY!
    truth_dir = truth_dir

    # Set the experiment and observation type to test (both shear_types will be explored)
    experiment = 'control'
    obs_type = 'ground'

    # Setup the logger
    logging.basicConfig(stream=sys.stderr)
    logger = logging.getLogger("test")
    logger.setLevel(logging.DEBUG)

    usebins = (evaluate.USEBINS, "subbins") #(evaluate.USEBINS, "subbins")
    poisson = (False, "noweight") 

    NTEST = 300
    NOISE_SIGMA = 0.15
    cvals = (evaluate.CFID, 10. * evaluate.CFID, 100. * evaluate.CFID) 
    mvals = (evaluate.MFID, 10. * evaluate.MFID, 100. * evaluate.MFID)

    # Then define arrays for storing the absolute difference results, the fractional difference
    # results, and the by-m,c results
    qabslarr = np.empty((NTEST, len(cvals), len(mvals)))
    qfracarr = np.empty((NTEST, len(cvals), len(mvals)))
    qsqrdarr = np.empty((NTEST, len(cvals), len(mvals)))
    carr = np.empty((NTEST, len(cvals), len(mvals)))
    marr = np.empty((NTEST, len(cvals), len(mvals)))
    cerrarr = np.empty((NTEST, len(cvals), len(mvals)))
    merrarr = np.empty((NTEST, len(cvals), len(mvals)))
    covcmarr = np.empty((NTEST, len(cvals), len(mvals)))
    qbymcarr = np.empty((NTEST, len(cvals), len(mvals)))

    # Set the sigma2_min for the qabsl and qfrac results:
    sigma2_min = {
        "ground": evaluate.SIGMA2_MIN_VARIABLE_GROUND,
        "space": evaluate.SIGMA2_MIN_VARIABLE_SPACE}[obs_type]
    # And the for the bymc metric, using same values as constant
    sigma2_min_bymc = {
        "ground": evaluate.SIGMA2_MIN_CONSTANT_GROUND,
        "space": evaluate.SIGMA2_MIN_CONSTANT_SPACE}[obs_type]

    print usebins[1]
    print poisson[1]
    print obs_type
    print NOISE_SIGMA

    # Get the x,y, true intrinsic ellips and shears for making fake submissions
    _, x, y, g1true, g2true = get_variable_gtrue(experiment, obs_type, truth_dir=truth_dir)
    _, _, _, g1int, g2int = get_variable_gsuffix(experiment, obs_type, truth_dir=truth_dir)

    # Then make a variable submission to learn what pure c1 and c2 look like in map^2
    map_E_unitc = make_unitc(experiment, obs_type, truth_dir=truth_dir)

    # Then start main loop
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
            for k in range(NTEST):
    
                fdsub, subfile = tempfile.mkstemp(suffix=".dat")
                result = make_variable_submission(
                    x, y, g1true, g2true, g1int, g2int, cval, cval, mval, mval,
                    outfile=subfile, noise_sigma=NOISE_SIGMA)
                os.close(fdsub)
                # First estimate the absolute metric
                qabsl = evaluate.q_variable(
                    subfile, experiment, obs_type, logger=None, usebins=usebins[0],
                    poisson_weight=poisson[0], fractional_diff=False, truth_dir=truth_dir,
                    sigma2_min=sigma2_min, normalization=evaluate.NORMALIZATION_VARIABLE_SPACE)
                print "%3d/%3d: Q_v_absl (c = %.4f, m = %.4f) = %.5e" % (
                    k + 1, NTEST, cval, mval, qabsl)
                # Then the fractional
                qfrac = evaluate.q_variable(
                    subfile, experiment, obs_type, logger=None, usebins=usebins[0],
                    poisson_weight=poisson[0], fractional_diff=True, truth_dir=truth_dir,
                    sigma2_min=sigma2_min_bymc / 10., normalization=300. * 0.6492)
                print "%3d/%3d: Q_v_frac (c = %.4f, m = %.4f) = %.5e" % (
                    k + 1, NTEST, cval, mval, qfrac)
                # Then the sqrd
                qsqrd = evaluate.q_variable(
                    subfile, experiment, obs_type, logger=None, usebins=usebins[0],
                    poisson_weight=poisson[0], fractional_diff=False, squared_diff=True, 
                    truth_dir=truth_dir, sigma2_min=sigma2_min * 3.e-5,
                    normalization=evaluate.NORMALIZATION_VARIABLE_SPACE * 4.687e-6)
                print "%3d/%3d: Q_v_sqrd (c = %.4f, m = %.4f) = %.5e" % (
                    k + 1, NTEST, cval, mval, qsqrd)
                qabslarr[k, ic, jm] = qabsl
                qfracarr[k, ic, jm] = qfrac
                qsqrdarr[k, ic, jm] = qsqrd
                # Then the m,c-derived metric
                qbymc, cest, mest, cerr, merr, covcm = evaluate.q_variable_by_mc(
                    subfile, experiment, obs_type, map_E_unitc, logger=None, truth_dir=truth_dir,
                    cfid=evaluate.CFID, mfid=2.e-2, pretty_print=True, sigma2_min=sigma2_min_bymc,
                    normalization=0.8912)
                print (
                    "%3d/%3d: Q_v_bymc (Estimated c = %.4f +/- %.4f, m = %.4f +/- %.4f) = %.4f" % (
                        k + 1, NTEST, cest, cerr, mest, merr, qbymc))
                carr[k, ic, jm] = cest
                marr[k, ic, jm] = mest
                cerrarr[k, ic, jm] = cerr
                merrarr[k, ic, jm] = merr
                covcmarr[k, ic, jm] = covcm
                qbymcarr[k, ic, jm] = qbymc
                print

            # Collate and print results
            print "Mean of Q_v_absl values = %.5e +/- %.5e" % (
                np.mean(qabslarr[:, ic, jm]),
                np.std(qabslarr[:, ic, jm]) / np.sqrt(len(qabslarr[:, ic, jm])))
            print "Mean of Q_v_frac values = %.5e +/- %.5e" % (
                np.mean(qfracarr[:, ic, jm]),
                np.std(qfracarr[:, ic, jm]) / np.sqrt(len(qfracarr[:, ic, jm])))
            print "Mean of Q_v_sqrd values = %.5e +/- %.5e" % (
                np.mean(qsqrdarr[:, ic, jm]),
                np.std(qsqrdarr[:, ic, jm]) / np.sqrt(len(qsqrdarr[:, ic, jm])))
            print "Mean of Q_v_bymc values = %.5e +/- %.5e" % (
                np.mean(qbymcarr[:, ic, jm]),
                np.std(qbymcarr[:, ic, jm]) / np.sqrt(len(qbymcarr[:, ic, jm])))

            print "Std of Q_v_absl values = %.5e +/- %.5e " % (
                np.std(qabslarr[:, ic, jm]),
                np.std(qabslarr[:, ic, jm]) / np.sqrt(2 * (len(qabslarr[:, ic, jm]) - 1)))
            print "Std of Q_v_frac values = %.5e +/- %.5e " % (
                np.std(qfracarr[:, ic, jm]),
                np.std(qfracarr[:, ic, jm]) / np.sqrt(2 * (len(qfracarr[:, ic, jm]) - 1)))
            print "Std of Q_v_sqrd values = %.5e +/- %.5e " % (
                np.std(qsqrdarr[:, ic, jm]),
                np.std(qsqrdarr[:, ic, jm]) / np.sqrt(2 * (len(qsqrdarr[:, ic, jm]) - 1)))
            print "Std of Q_v_bymc values = %.5e +/- %.5e " % (
                np.std(qbymcarr[:, ic, jm]),
                np.std(qbymcarr[:, ic, jm]) / np.sqrt(2 * (len(qbymcarr[:, ic, jm]) - 1)))

            print "Fractional uncertainty on Q_v_absl values = %.5e +/- %.5e " % (
                np.std(qabslarr[:, ic, jm]) / np.mean(qabslarr[:, ic, jm]),
                np.std(qabslarr[:, ic, jm]) / np.sqrt((len(qabslarr[:, ic, jm]) - 1))
                / np.mean(qabslarr[:, ic, jm]))
            print "Fractional uncertainty on Q_v_frac values = %.5e +/- %.5e " % (
                np.std(qfracarr[:, ic, jm]) / np.mean(qfracarr[:, ic, jm]),
                np.std(qfracarr[:, ic, jm]) / np.sqrt((len(qfracarr[:, ic, jm]) - 1))
                / np.mean(qfracarr[:, ic, jm]))
            print "Fractional uncertainty on Q_v_sqrd values = %.5e +/- %.5e " % (
                np.std(qsqrdarr[:, ic, jm]) / np.mean(qsqrdarr[:, ic, jm]),
                np.std(qsqrdarr[:, ic, jm]) / np.sqrt((len(qsqrdarr[:, ic, jm]) - 1))
                / np.mean(qsqrdarr[:, ic, jm]))
            print "Fractional uncertainty on Q_v_bymc values = %.5e +/- %.5e " % (
                np.std(qbymcarr[:, ic, jm]) / np.mean(qbymcarr[:, ic, jm]),
                np.std(qbymcarr[:, ic, jm]) / np.sqrt((len(qbymcarr[:, ic, jm]) - 1))
                / np.mean(qbymcarr[:, ic, jm]))

            print "Best estimate of Q_v_absl normalization factor = "+str(
                1000. / np.mean(qabslarr[:, ic, jm]))
            print "Best estimate of Q_v_frac normalization factor = "+str(
                1000. / np.mean(qfracarr[:, ic, jm]))
            print "Best estimate of Q_v_sqrd normalization factor = "+str(
                1000. / np.mean(qsqrdarr[:, ic, jm]))
            print "Best estimate of Q_v_bymc normalization factor = "+str(
                1000. / np.mean(qbymcarr[:, ic, jm]))

    # Save the arrays
    filename = os.path.join(
        evaluate.STORAGE_DIR, "testing_qabsl_"+obs_type+"_NOISE_SIGMA"+("%.2f" % NOISE_SIGMA)+"_"+
        usebins[1]+"_"+poisson[1]+"_mc_N"+str(NTEST)+".npy")
    print "Saving to "+filename
    np.save(filename, qabslarr)
    filename = os.path.join(
        evaluate.STORAGE_DIR, "testing_qfrac_"+obs_type+"_NOISE_SIGMA"+("%.2f" % NOISE_SIGMA)+"_"+
        usebins[1]+"_"+poisson[1]+"_mc_N"+str(NTEST)+".npy")
    print "Saving to "+filename
    np.save(filename, qfracarr)
    filename = os.path.join(
        evaluate.STORAGE_DIR, "testing_qsqrd_"+obs_type+"_NOISE_SIGMA"+("%.2f" % NOISE_SIGMA)+"_"+
        usebins[1]+"_"+poisson[1]+"_mc_N"+str(NTEST)+".npy")
    print "Saving to "+filename
    np.save(filename, qsqrdarr)
    filename = os.path.join(
        evaluate.STORAGE_DIR, "testing_qbymc_"+obs_type+"_NOISE_SIGMA"+("%.2f" % NOISE_SIGMA)+"_"+
        usebins[1]+"_"+poisson[1]+"_mc_N"+str(NTEST)+".npy")
    print "Saving to "+filename
    np.save(filename, qbymcarr)
    filename = os.path.join(
        evaluate.STORAGE_DIR, "testing_mcresults_"+obs_type+"_NOISE_SIGMA"+("%.2f" % NOISE_SIGMA)+
        "_"+usebins[1]+"_"+poisson[1]+"_mc_N"+str(NTEST)+".pkl")
    mcresults = {
        "c": carr, "m": marr, "cerr": cerrarr, "merr": merrarr, "covcm": covcmarr, "qbymc": qbymc}
    print "Saving to "+filename
    import cPickle
    with open(filename, "wb") as fout: cPickle.dump(mcresults, fout)
