#!/usr/bin/env python
#
# Copyright 2013 Rachel Mandelbaum & the GREAT3 Team:
# https://github.com/barnabytprowe/great3-public
#
# This file is example code released as part of the GREAT3 lensing competition.
#
# The GREAT3 example code is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This example code is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#

"""
Example module for a simple shear estimator in GREAT3

This Python module provides a simple shear estimation code adapted to work with the GREAT3 data
format.  It is intended to provide a "first-order" shear estimation module for GREAT3 participants
who want to figure out a way to adapt their shear estimation code to work on the GREAT3 simulations.
It is not expected to be competitive with the current state-of-the-art in shear estimation, but it
should not be orders of magnitude worse, either.


DEPENDENCIES

Python 2.6.x-2.7.x
NumPy 1.6+
GalSim 0.5+
PyFITS 1.0+


ALGORITHM / USAGE

The approach in this script is to use a simple shear estimation routine available within GalSim to
estimate a Bernstein & Jarvis (2002)-style *distortion* for each of the 10^4 galaxies in a single
GREAT3 image, then use the ensemble of distortions to determine a "shear responsivity" factor which
can be used to convert from distortion to the mean *shear*.  The results are then stored in an
output catalog in a format required by the GREAT3 presubmission script.

Note that the applied shear responsivity and calibration are only really correct for the ensemble,
and are not entirely meaningful for individual objects.  Likewise, the per-object distortion is
physically meaningful, but the per-object shear is not - it's only the average shear over the
ensemble that is meaningful.

The method of PSF correction that is used by default is the re-Gaussianization method proposed by
Hirata & Seljak (2004).  However, other options are available (see below for more information).

The default shear calibration factor of 0.98 is meant to be applicable for the re-Gaussianization
method and is based on previously published work.  A brief explanation of the source of this factor
is that it comes from two sources:

(1) In Mandelbaum et al. (2012, MNRAS, 420, 1518), section 9.2, a calibration factor of ~1.03 was
found to be necessary due to intrinsic limitations of re-Gaussianization for realistic galaxies and
due to noise bias, in the absence of selection bias.  This was for a COSMOS sample that is
intrinsically brighter than the one used here, for simulated SDSS-like conditions, so the
extrapolation to GREAT3 is somewhat questionable.

(2) In Reyes et al. (2012, MNRAS, 425, 2610), section 4.3.2, it was shown that if we blindly accept
the shape measurement errors coming out of the re-Gaussianization routine without applying some
correction for their known underestimation, then for a galaxy sample that is weighted towards S/N<30
(as in GREAT3), the RMS ellipticity of the sample will be significantly overestimated.  In the
example case given there, the overestimate was 0.41 rather than 0.36.  This would result in the
shear responsivity being underestimated and therefore the shear being overestimated by 5%, giving a
shear calibration factor of 0.95.  When combined with the factor from (1), the net shear calibration
factor is 0.98.


Executing the script looks something like this:

./simple_shear.py --output_prefix my_prefix 27 /sims/control/ground/constant ./output

This will process subfield 27 for the branch located in '/sims/control/ground/constant', placing all
outputs in the './output' directory.  Use of the output_prefix option means the output file will be
named 'my_prefix-027.fits'; the default prefix is 'output_catalog'.  For more help with the
command-line interface, run the following command:

./simple_shear.py --help

This script has been run on all subfields in the 8 branches in the control and real_galaxy
experiments.  On a reasonably up-to-date machine with a single processor, the script takes
approximately 30 s (160 s) to process a single subfield of simulated ground (space) data.  These
correspond to 3 ms (16 ms) per galaxy measured.  Typical success rates for the shear estimation
process are 99.7%, with failures being silently set to large shear values that indicate failure to
presubmission.py.


If the above command is run for all 200 subfields (e.g., using mass_produce_shear.py, a driver
script that calls simple_shear.py once per subfield), then the results could be passed to the GREAT3
presubmission script using

./presubmission.py -o cgc.out --use-fits -b control-ground-constant ./output/my_prefix*fits
-w weight

where the output would be sent to a file called cgc.out; --use-fits tells the presubmission script
that we have made catalogs in FITS format; the -b command is necessary to indicate the branch name
(since our catalog names do not include it); we give it a list of catalogs in
./output/my_prefix*fits; and the -w command is used to tell the script that we will use non-uniform
weighting indicated by the 'weight' field in the FITS catalogs.


KNOWN ISSUES / LIMITATIONS

- The code is set up to work most simply for control and real_galaxy branches.

- For multiepoch branches, the most straightforward usage is to make a set of coadds and process
  them.  The outputs of the coaddition script ('coadd_multiepoch.py') distributed in this repository
  can be processed if they are put in the same directory as the catalogs, and this script is run
  with the --coadd flag.  The --coadd flag results in this script looking for galaxy and starfield
  images with the 'coadd_' prefix that is used by coadd_multiepoch.py, assuming the starfield just
  has a single image of the coadded PSF rather than 9 (like the single epoch starfields).

- For variable PSF branches, this code can be run using the grids of PSF models
  (psf_models-###.fits) output by the publicly distributed psf_models.py script (in this directory).
  To use simple_shear.py in this mode, use the --variable_psf_dir flag to indicate the directory
  containing the psf_models-###.fits files.

- This script has not been set up to easily run on the 'full' branch, for which both coaddition and
  estimation of the variable PSF would be necessary.

- The code assumes that the galsim.hsm routine is being used to estimate a per-object distortion,
  which is true for the (default) re-Gaussianization method and the linear method.  This means that
  under 'default_shear_kwds' it should be possible to set `shear_est="Linear"` without any other
  modification of the code, though the calibration factor that is applied by default will not be
  correct.  However, the KSB method that is part of the galsim.hsm module returns per-object
  shears, so modification of this code (to remove responsivity calculation and application) would be
  necessary in order to get sensible results when running galsim.hsm in KSB mode.

- This code uses a relatively naive RMS ellipticity calculation (after subtraction of a simple
  estimate of measurement noise that is known to be an underestimate).  Furthermore, this scheme
  would ideally use a larger sample than appears in a single subfield; the use of just a single
  subfield introduces some noise into the responsivity and therefore the per-subfield shear.  More
  sophisticated schemes for getting the sample responsivity (e.g., use of equation 5-33 in Bernstein
  & Jarvis, 2002) and for estimating the uncertainty in per-object shear estimates would do a better
  job at overall shear estimation, but are beyond the scope of this module that is intended largely
  for demonstration purposes.  Note that the default calibration factor is in part intended to
  account for biases in the estimated responsivity due to the incorrect RMS ellipticity calculation.

"""

import sys
import time
import os
import optparse
import numpy

try:
    import astropy.io.fits as pyfits
except:
    import pyfits

import galsim

# This dict is used to provide keyword arguments to galsim.hsm.EstimateShear when estimating the
# per-object distortion.  We use the strict=False keyword argument so that the code will not grind
# to a halt if it hits an object for which the shear estimation code raises an exception.  Other
# keyword arguments can be used to select a different shear estimation routine than the default (see
# doctring for more details).
default_shear_kwds = {"strict": False}

verbose = False

def log(msg):
    if verbose:
        print msg

def readData(sim_dir, subfield, coadd, variable_psf_dir=""):
    """Subroutine to read in galaxy and star field images and the galaxy catalog for the specified
    simulation.

    Arguments:

      sim_dir ----------- simulation directory, containing the GREAT3 images and catalogs for a
                          single branch

      subfield ---------- subfield to be processed, as an integer

      coadd ------------- Are we processing the outputs of the coaddition script,
                          coadd_multiepoch.py? If so, set this to True.

      variable_psf_dir--- Directory in which to find PSF model outputs from psf_models.py, if this
                          is a variable PSF branch.  Default value of "" indicates that this is not
                          a variable_psf branch.

    """
    # Construct filename for galaxy image, and read in file.
    if not coadd:
        infile = os.path.join(sim_dir, 'image-%03d-0.fits'%subfield)
    else:
        infile = os.path.join(sim_dir, 'coadd_image-%03d.fits'%subfield)
    try:
        gal_im = galsim.fits.read(infile)
    except:
        raise RuntimeError("Could not read in file %s."%infile)
    log("... Read galaxy image from file "+infile)

    # Construct filename for starfield image, and read in file.  There are three options: a standard
    # control/real_galaxy star field; a variable_psf grid of PSF models, one per galaxy; and a coadd
    # starfield.
    if not coadd and variable_psf_dir=="":
        # This is a standard star field.
        infile = os.path.join(sim_dir, 'starfield_image-%03d-0.fits'%subfield)
    elif variable_psf_dir != "":
        # This is a grid of PSF models for a variable PSF branch.
        infile = os.path.join(variable_psf_dir, 'psf_models-%03d.fits'%subfield)
    else:
        # This is a coadd.
        infile = os.path.join(sim_dir, 'coadd_starfield_image-%03d.fits'%subfield)
    try:
        starfield_im = galsim.fits.read(infile)
    except:
        raise RuntimeError("Could not read in file %s."%infile)
    log("... Read starfield image from file "+infile)

    # Construct filename for galaxy catalog, and read in file.
    infile = os.path.join(sim_dir, 'galaxy_catalog-%03d.fits'%subfield)
    try:
        gal_catalog = pyfits.getdata(infile)
    except:
        raise RuntimeError("Could not read in file %s."%infile)
    log("... Read galaxy catalog from file "+infile)

    return gal_im, starfield_im, gal_catalog

def extractPSF(starfield_im):
    """Subroutine to extract a single PSF image from a 3x3 starfield.

    This routine assumes we are in one of the constant PSF branches, such that the starfield is
    simply a 3x3 grid of stars, for which we wish to extract the lower-left corner (containing a
    star that is perfectly centered in the postage stamp).

    Arguments:

      starfield_im ----- The full starfield image from which we wish to extract the PSF for this
                         galaxy image.

    """
    # Figure out the size of the images based on knowing that the starfield image is a 3x3 grid of
    # stars.
    shape = starfield_im.array.shape
    if shape[0] != shape[1]:
        raise RuntimeError("This starfield image is not square!")
    if shape[0] % 3 != 0:
        raise RuntimeError("Starfield image size is not a multiple of 3!")
    ps_size = shape[0] / 3

    # Cut out the lower-left postage stamp.
    psf_im = starfield_im[galsim.BoundsI(1,ps_size,1,ps_size)]
    return psf_im

def estimateVariance(gal_im):
    """Subroutine to do a fairly naive estimation of the sky variance, using edge pixels.

    This routine uses the fact that the sky variance is the same across the entire large image, so
    we can use a set of pixels at the edge (with minimal contamination from galaxy fluxes) in order
    to estimate the variance.

    Arguments:
      gal_im -------- The full galaxy image to use for estimating the sky variance.

    """
    # Figure out the size of the images
    shape = gal_im.array.shape
    if shape[0] != shape[1]:
        raise RuntimeError("This image is not square!")

    # Choose the 8 rows/columns and make a single array with sky pixels to use
    sky_vec = numpy.concatenate(
        (gal_im.array[0:2, 2:shape[0]-2].flatten(),
         gal_im.array[shape[0]-2:shape[0], 2:shape[0]-2].flatten(),
         gal_im.array[2:shape[0]-2, 0:2].flatten(),
         gal_im.array[2:shape[0]-2, shape[0]-2:shape[0]].flatten()
         ))

    # Estimate and return the variance
    return numpy.var(sky_vec)

def getPS(record, gal_im, ps_size, starfield_image=None):
    """Routine to pull out a galaxy postage stamp based on a record from a catalog.

    Arguments:

      record ---------- A single record from a GREAT3 galaxy catalog for the chosen subfield.

      gal_im ---------- The full galaxy image for that subfield.

      ps_size --------- Total (linear) size of a single postage stamp.

      starfield_image - Grid of PSF models for the galaxies, if this is a variable PSF branch.  In
                        that case, the routine returns not just the galaxy postage stamp, but also
                        the PSF postage stamp.
    """
    # Figure out the galaxy postage stamp bounds based on the record information.
    radius = ps_size / 2
    bounds = galsim.BoundsI(int(numpy.ceil(record['x']) - radius) + 1,
                            int(numpy.ceil(record['x']) + radius),
                            int(numpy.ceil(record['y']) - radius) + 1,
                            int(numpy.ceil(record['y']) + radius))

    # Pull out and return the postage stamp.
    subimage = gal_im[bounds]

    if starfield_image is None:
        return subimage
    else:
        # Then work on the PSF image, if this is a variable PSF branch.
        # This is a general way to find the dimensions of a PSF postage stamp, just in case users
        # use a non-default option for psf_models.py.
        psf_dim = (starfield_image.xmax + 1 - starfield_image.xmin) * ps_size / \
            (gal_im.xmax + 1 - gal_im.xmin)
        # Now find the proper indices
        x_index = (record['x'] + 1 - radius) / ps_size
        y_index = (record['y'] + 1 - radius) / ps_size
        # Take the subimage.
        bounds = galsim.BoundsI(int(x_index*psf_dim), int((x_index+1)*psf_dim)-1,
                                int(y_index*psf_dim), int((y_index+1)*psf_dim)-1)
        psf_subimage = starfield_image[bounds]
        return subimage, psf_subimage

def checkFailures(shear_results):
    """Routine to check for shape measurement failures which should be flagged as such.

    Arguments:

      shear_results --- a list of galsim.hsm.ShapeData objects that comes from measuring shapes of
                        all galaxies in the image.

    """
    n_gal = len(shear_results)

    # Define output structure: a boolean numpy array
    use_shape = numpy.ones(n_gal).astype(bool)

    # Compare resolution factor with galsim.hsm.HSMParams.failed_moments or look for other obvious
    # oddities in shears, quoted errors, etc.
    hsmparams = galsim.hsm.HSMParams()
    for index in range(n_gal):
        test_e = numpy.sqrt(
            shear_results[index].corrected_e1**2 + shear_results[index].corrected_e2**2
            )
        if shear_results[index].resolution_factor == hsmparams.failed_moments or \
                shear_results[index].corrected_shape_err < 0 or test_e > 4. or \
                shear_results[index].corrected_shape_err > 0.5:
            use_shape[index] = False
    return use_shape


def getRMSEllip(shear_results=None, use_shape=None, weight=None, e1=None, e2=None, sigma_e=None):
    """Routine for estimating the RMS ellipticity from a catalog.

    The routine works in one of two ways:

    In the first approach, if we have not done any shape calculations before, we can pass in an
    array of galsim.hsm.ShapeData structures for every object in the image (`shear_results`), a
    boolean NumPy array indicating which ones are useful, and (optionally) a set of weights.  The
    routine will then return the RMS ellipticity (optionally weighted) as well as arrays with e1,
    e2, sigma_e for each object used.

    In the second approach, if we've already gotten e1, e2, sigma_e before, and we just want to
    recalculate RMS ellipticity for some reason (e.g., with weights), then we can pass in e1, e2,
    sigma_e along with an array of weights of the same length, and calculate the weighted RMS
    ellipticity, which is the sole quantity that will be returned.
    """
    mode = None
    if e1 is None and e2 is None and sigma_e is None:
        if shear_results is None or use_shape is None:
            raise RuntimeError("Missing information: need ShapeData and usage flags!")
        else:
            mode = 1

    if shear_results is None and use_shape is None:
        if e1 is None or e2 is None or sigma_e is None:
            raise RuntimeError("Missing information: need e1, e2, sigma_e!")
        else:
            mode = 2

    if mode is None:
        raise RuntimeError("getRMSEllip called without an obvious working mode!")

    if mode == 1:
        n_gal = len(shear_results)
        e1 = []
        e2 = []
        sigma_e = []
        for index in range(n_gal):
            if use_shape[index]:
                e1.append(shear_results[index].corrected_e1)
                e2.append(shear_results[index].corrected_e2)
                # There's a factor of ~2 here, because we want sigma_e and it returns sigma_g.
                # Actually, the exact factor is 2/(1+2*sigma_g^2) according to Bernstein & Jarvis
                # (2002), but the sigma_g values are usually <~0.15 so the denominator is quite
                # close to 1 and therefore ignored for the sake of simplicity.
                sigma_e.append(2.*shear_results[index].corrected_shape_err)
        e1 = numpy.array(e1)
        e2 = numpy.array(e2)
        sigma_e = numpy.array(sigma_e)

    n_shape = len(e1)
    if weight is None:
        weight = numpy.ones(n_shape)
    e_rms_per_component = numpy.sqrt(
        numpy.sum(0.5*(e1**2 + e2**2 - 2*sigma_e**2)*weight) / numpy.sum(weight)
        )
    if mode == 1:
        return e_rms_per_component, e1, e2, sigma_e
    else:
        return e_rms_per_component

def writeOutput(output_dir, output_prefix, subfield, output_type, catalog, clobber = True,
                comment_pref = '#'):
    """Routine for writing outputs to file in some specified format, either fits or ascii.

    Arguments:

      output_dir ------ Output directory.

      output_prefix --- Prefix for output files.

      subfield -------- Subfield that was processed.

      output_type ----- Type of file to output, either FITS or ascii.

      catalog --------- Catalog to write to file.

      clobber --------- Overwrite catalog if it already exists?  Default: true.

      comment_pref ---- String to use to denote comments in ascii outputs.  Default: '#'
    """
    if output_type == "fits":
        output_file = '%s-%03d.fits'%(output_prefix, subfield)
    else:
        output_file = '%s-%03d.dat'%(output_prefix, subfield)
    output_file = os.path.join(output_dir, output_file)
    if not clobber and os.path.exists(output_file):
        raise RuntimeError("Error: file %s already exists!" % output_file)

    if output_type == "fits":
        pyfits.writeto(output_file, catalog, clobber = clobber)
        log("Wrote output catalog to file %s."%output_file)
    else:
        import tempfile
        import shutil
        # First print lines with column names.  This goes into a separate file for now.
        f1, tmp_1 = tempfile.mkstemp(dir=output_dir)
        with open(tmp_1, 'w') as f:
            # Note: this is a stupid way to get it to list the names, but I can't figure out a
            # simpler way.  For example, just `catalog.names` doesn't work.
            name_list = pyfits.new_table(catalog).data.names
            formats = pyfits.new_table(catalog).data.columns.formats
            for name in name_list:
                f.write(comment_pref + " %s\n"%name)

        # Then save catalog itself.
        f2, tmp_2 = tempfile.mkstemp(dir=output_dir)
        # TODO: get format string from column description, so that we can print IDs as ints instead
        # of floats.  This seems to be hard, but I'm not sure why.
        numpy.savetxt(tmp_2, catalog)

        # Finally, concatenate, and remove tmp files
        with open(output_file, 'wb') as destination:
            shutil.copyfileobj(open(tmp_1, 'rb'), destination)
            shutil.copyfileobj(open(tmp_2, 'rb'), destination)
        # Need to close the tempfiles opened by mkstemp.  cf:
        # http://stackoverflow.com/questions/9944135/how-do-i-close-the-files-from-tempfile-mkstemp
        log("Wrote output catalog to file "+output_file)
        os.close(f1)
        os.close(f2)
        os.remove(tmp_1)
        os.remove(tmp_2)

def EstimateAllShears(subfield, sim_dir, output_dir, output_prefix="output_catalog", output_type="fits",
                      clobber=True, sn_weight=True, calib_factor=0.98, coadd=False,
                      variable_psf_dir=""):
    """Main driver for all routines in this file, and the implementation of most of the command-line
    interface.

    This routine has three distinct steps:

      1) Read required inputs from file.

      2) Loop over objects in the catalog and estimate a shear for each.

      3) Collect outputs into the proper format, and write to disk.

    Required arguments:

      subfield ----- subfield to be processed, as an integer

      sim_dir ------ simulation directory, containing the GREAT3 images and catalogs for a single
                     branch

      output_dir --- directory in which outputs should be placed

    Optional arguments:

      output_prefix --- Prefix for output catalog; the subfield (as a 3-digit number) will be
                        appended

      output_type ----- Type for output catalogs: fits (default) or ascii.

      clobber --------- Overwrite pre-existing output files?  Default: true.

      sn_weight ------- Apply S/N-dependent weighting to each object?  Default: true.  If false, all
                        objects get an equal weight.

      calib_factor ---- Multiplicative calibration factor to apply to all shears.  Default value of
                        0.98 is based on estimates from previously published work; see full
                        docstring in simple_shear.py for details.

      coadd ----------- Are we processing the outputs of the coaddition script, coadd_multiepoch.py?
                        If so, set this to True.  Default: false.

      variable_psf_dir- Directory in which to find PSF model outputs from psf_models.py, if this is
                        a variable PSF branch.  Default value of "" indicates that this is not a
                        variable_psf branch.

    """
    if coadd is True and variable_psf_dir!="":
        raise NotImplementedError("Script is not set up to process full experiment.")

    t1 = time.time()
    log("Preparing to read inputs: galaxy and starfield images, and galaxy catalog.")
    # First read inputs.
    gal_im, starfield_im, gal_catalog = readData(sim_dir, subfield, coadd=coadd,
                                                 variable_psf_dir=variable_psf_dir)
    guess_sig = 3.0 # the default guess for PSF sigma
    if variable_psf_dir=="":
        # Since this is a constant PSF branch, make a PSF image from the lower-left (centered) star
        # in the starfield.  We don't have to do this if we already have the output of the
        # coaddition script, which is a single star.
        if coadd:
            psf_im = starfield_im
        else:
            psf_im = extractPSF(starfield_im)
        # Very rarely, the EstimateShear routine can fail for space PSFs that are very aberrated.  A
        # simple fix for this is to change the initial guess for PSF size.  However, that leads to
        # slower convergence for the 99.9% of subfields that don't have this problem, so we don't
        # want to do it all the time.  Instead, we check once for failure of adaptive moments for
        # the PSF, and if it fails, then we adopt a smaller initial guess.
        try:
            galsim.hsm.FindAdaptiveMom(psf_im)
        except:
            guess_sig = 2.0
            try:
                galsim.hsm.FindAdaptiveMom(psf_im, guess_sig=guess_sig)
            except:
                raise RuntimeError("Cannot find a value of PSF size that works for this PSF!")

    # Get number of entries
    n_gal = len(gal_catalog)
    log("Preparing to process %d galaxies"%n_gal)

    # We need to give the shear estimation routines a rough idea of the sky variance, so that it can
    # estimate uncertainty in per-object shears.  For this purpose, we compute the variance from a
    # set of pixels around the outer edge of the image.
    sky_var = estimateVariance(gal_im)
    log("Estimated sky variance: %f"%sky_var)

    shear_results = []
    ps_size = gal_catalog['x'][1] - gal_catalog['x'][0] # size of a postage stamp defined by
                                                        # centroid differences
    # Loop over objects in catalog.
    log("Beginning loop over galaxies: get postage stamp and estimate distortion...")
    t2 = time.time()
    for record in gal_catalog:
        # Pull out the postage stamp corresponding to that record.  If this is a variable PSF
        # branch, then we need to also get the PSF postage stamp.
        if variable_psf_dir=="":
            gal_ps = getPS(record, gal_im, ps_size)
        else:
            starfield_im.setOrigin(0,0)
            gal_ps, psf_im = getPS(record, gal_im, ps_size, starfield_image=starfield_im)
        # Estimate the shear, requiring a silent failure if something goes wrong.  (That flag is in
        # `default_shear_kwds`
        shear_results.append(
            galsim.hsm.EstimateShear(gal_ps, psf_im, sky_var=sky_var,
                                     guess_sig_PSF = guess_sig,  **default_shear_kwds)
            )
    dt = time.time() - t2
    log("...Time per object=%f s, total time for loop=%f s"%(dt/n_gal,dt))

    # First figure out which objects have failures.
    use_shape = checkFailures(shear_results)
    n_success = numpy.round(use_shape.sum())
    log("Number with successful measurements: %d, or %f percent"%
        (n_success,100*float(n_success)/n_gal))

    # Estimate RMS ellipticity.
    e_rms_per_component, e1, e2, sigma_e = getRMSEllip(shear_results=shear_results,
                                                       use_shape=use_shape, weight=None)
    log("Estimated RMS ellipticity: %f"%e_rms_per_component)

    # Use that to define per-object weights.
    if sn_weight:
        use_weight = 1./(e_rms_per_component**2 + sigma_e**2)
    
        # Re-estimate weighted RMS ellipticity and therefore shear responsivity.
        e_rms_per_component = getRMSEllip(e1=e1, e2=e2, sigma_e=sigma_e, weight=use_weight)
        log("Estimated RMS ellipticity with weighting: %f"%e_rms_per_component)

    # The equation below is the simplest possible estimator of the shear responsivity, which we use
    # for the benefit of simplicity.  A better formula is in Bernstein & Jarvis (2002), equation
    # (5-33), and that is the formula that has been used in science papers that use the
    # re-Gaussianization method.
    responsivity = 2.*(1.-e_rms_per_component**2)
    if responsivity > 2. or responsivity < 1.:
        raise RuntimeError("Error: responsivity is %f"%responsivity)
    log("Shear responsivity (2 for round objects, 1 for perfectly flat): %f"%responsivity)

    # Now put outputs into the format that we want.  For default values of g1 and g2, we put 100,
    # to indicate failure.  Then we will change them to some more sensible value if a shape was
    # measurable for that object.
    schema = [("id", int), ("g1", float), ("g2", float), ("weight", float)]
    catalog = numpy.zeros(n_gal, dtype=numpy.dtype(schema))
    catalog["id"] = gal_catalog['ID']
    g1 = numpy.zeros(n_gal).astype(numpy.float64) + 100.
    g2 = numpy.zeros(n_gal).astype(numpy.float64) + 100.
    weight = numpy.zeros(n_gal).astype(numpy.float64) + 1.
    index = 0
    use_index = 0
    for index in range(n_gal):
        if use_shape[index]:
            g1[index] = calib_factor*shear_results[index].corrected_e1 / responsivity
            g2[index] = calib_factor*shear_results[index].corrected_e2 / responsivity
            if sn_weight:
                weight[index] = use_weight[use_index]
            use_index += 1
    catalog["g1"] = g1
    catalog["g2"] = g2
    catalog["weight"] = weight

    # Write to file.
    writeOutput(output_dir, output_prefix, subfield, output_type, catalog, clobber = clobber)
    log("Total time processing this subfield: %f s."%(time.time()-t1))

def main(argv):
    usage = "usage: %prog [options] SUBFIELD SIM_DIR WORK_DIR"
    description = """Estimate shears for all galaxies in the given subfield, applying all necessary responsivity and calibration factors.  SUBFIELD is the number of the subfield to be processed. SIM_DIR is the directory containing GREAT3 images and catalogs for the branch of interest. WORK_DIR is the directory where output files should be placed.  It will be created if it does not exist."""
    parser = optparse.OptionParser(usage=usage, description=description)
    parser.add_option("--output_prefix", dest="output_prefix", type=str, default="output_catalog",
                      help="Prefix for output file")
    parser.add_option("--output_type", dest="output_type", type=str, default="fits",
                      help="Type of output catalog: fits or ascii")
    parser.add_option("--no_clobber", action="store_true",
                      help="Do not clobber pre-existing output files")
    parser.add_option("--no_weight", action="store_true",
                      help="Do not apply S/N-dependent weighting; weight all galaxies equally")
    parser.add_option("--calib_factor", dest="calib_factor", type=float, default=0.98,
                      help="Multiplicative calibration factor to apply to all shears, default 0.98")
    parser.add_option("--coadd", dest="coadd", action='store_true', default=False,
                      help="Use to indicate that we are processing coadd_multiepoch.py outputs")
    parser.add_option("--variable_psf_dir", dest="variable_psf_dir", type=str, default="",
                      help="Directory in which to find variable PSF models; only use this option for a variable PSF branch!")
    parser.add_option("--quiet", dest="quiet", action='store_true', default=False,
                      help="Don't print progress statements")
    opts, args = parser.parse_args()
    try:
        subfield, sim_dir, output_dir = args
    except ValueError:
        parser.error("exactly three positional arguments are required")
    if not os.path.isdir(sim_dir):
        parser.error("input directory %s does not exist or is not a directory" % sim_dir)
    try:
        subfield = int(subfield)
    except TypeError:
        parser.error("subfield argument '%s' is not an integer" % subfield)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    if opts.output_type is not None:
        if opts.output_type not in ("fits", "ascii"):
            parser.error("output_type '%s' must be one of 'fits' or 'ascii'" % opts.output_type)
    if opts.no_clobber:
        clobber = False
    else:
        clobber = True
    if opts.no_weight:
        sn_weight = False
    else:
        sn_weight = True
    global verbose
    if opts.quiet:
        verbose = False
    else:
        verbose = True
    EstimateAllShears(
        subfield, sim_dir, output_dir,
        output_prefix=opts.output_prefix,
        output_type=opts.output_type,
        clobber=clobber,
        sn_weight=sn_weight,
        calib_factor=opts.calib_factor,
        coadd=opts.coadd,
        variable_psf_dir=opts.variable_psf_dir
        )

if __name__ == "__main__":
    main(sys.argv)
