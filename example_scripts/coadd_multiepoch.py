#!/usr/bin/env python
#
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
"""
Example module for combining multiepoch image data in GREAT3.

This Python module provides a simple multiepoch imaging data coaddition code adapted to work with
the GREAT3 data format.  It is intended to provide a "first-order" imaging combination module to
allow GREAT3 participants to compete in branches with multiepoch data without having to write their
own equivalent code, and to provide a pedagogical example for those who might wish to create their 
own image combination module.  It is not expected to be competitive with the current
state-of-the-art in astrophysical image combination, but it should not be orders of magnitude worse,
either.


DEPENDENCIES

Python 2.6.x-2.7.x
NumPy 1.6+
GalSim 1.0
PyFITS 1.0+ / Astropy v0.1+


ALGORITHM / USAGE

Images are combined by interpolative coaddition, resulting in an output "coadd" image of higher
signal-to-noise ratio than the input images taken individually.  For space images this script will
also use interpolation to double the resolution of the output coadd relative to the input pixel
scale.

The coaddition is performed using the GalSim interpolation routines by wrapping the input images
in galsim.InterpolatedImage instances.  (This is in fact done for a set of smaller subimages that
tile the input image region, to reduce memory usage.)  The output is then coadded to an output image
using the galsim.GSObject `.draw()` method having set the  `add_to_image=True` keyword option,
repeating for each input image successively.  For more information about these classes and methods
see the GalSim documentation.

The interpolation used is currently set as the 'galsim.Cubic(tol=1.e-4)' interpolant by default,
corresponding to bicubic interpolation between points on the XY grid.  If you wish to use a
different interpolant, possibly with a larger support, consider increasing the --border keyword to 
increase the overlap between internal subimages.  See the --help message for more details.

Executing the script looks something like this:

./coadd_multiepoch.py --space 42 /sims/multiepoch/space/constant ./msc_coadds

This will process subfield 42 for the GREAT3 branch data located in
/sims/multiepoch/space/constant, generating coadd data from the GREAT3 galaxy simulation images

/sims/multiepoch/space/constant/image-042-0.fits
/sims/multiepoch/space/constant/image-042-1.fits
/sims/multiepoch/space/constant/image-042-2.fits
/sims/multiepoch/space/constant/image-042-3.fits
/sims/multiepoch/space/constant/image-042-4.fits
/sims/multiepoch/space/constant/image-042-5.fits

and the starfield images

/sims/multiepoch/space/constant/starfield_image-042-0.fits
/sims/multiepoch/space/constant/starfield_image-042-1.fits
/sims/multiepoch/space/constant/starfield_image-042-2.fits
/sims/multiepoch/space/constant/starfield_image-042-3.fits
/sims/multiepoch/space/constant/starfield_image-042-4.fits
/sims/multiepoch/space/constant/starfield_image-042-5.fits

These files must be present in the specified directory, along with the corresponding
`epoch_dither-*.txt` files.

For more help with the command-line interface, you can just run:

./coadd_multiepoch.py --help

Two FITS format coadd output files will be placed in the specified directory, `./msc_coadds` in the
example above:

./msc_coadds/coadd_image-042.fits
./msc_coadds/coadd_starfield_image-042.fits

In the starfield image, only the lower-left PSF image (following standard FITS conventions) is
coadded and output.  The other, offset PSF images do have the same offset in each epoch, so they
could be coadded too, but it is not strictly necessary and is left as an exercise for the interested
reader!

If wishing to coadd space branch images, don't forget to use the --space option to double the output
coadd image resolution.  N.B. The output will still suffer from aliasing, see below.

Memory usage should not exceed ~1.3GB for the `multiepoch` branch data, and will be largest when
making coadds from space images.  For the `full` dataset memory usage may be higher for the larger
starfield images: this can be addressed if necessary by experimenting with the `n_split_default`
variable, assigned at module scope.


KNOWN ISSUES / LIMITATIONS

 - Undersampled input images, such as the multiepoch/space/* data, are aliased and will contain
   distortions when interpolated.  As the coaddition here is interpolative, aliasing may be a
   significant factor for the space branches.  Solving this problem to generate oversampled images
   free from aliasing is a non-trivial problem, which we don't attempt to solve in this
   freely-provided example script. 

 - Only the lower-left PSF in coadd output generated from the multiepoch branch starfield images
   is output (see above).

 - Noise in coadded images is necessarily correlated on small spatial scales.  GalSim provides tools
   to estimate the level of noise correlations via the 2D discrete correlation function, and to
   whiten this noise (by adding more noise with suitably chosen anticorrelations) if you wish to
   test the impact of this effect.  See the galsim.CorrelatedNoise class, part of the 
   galsim/correlatednoise.py module in the GalSim repository, for more information.


NOTES

As is usual with NumPy, 2-d arrays representing images are ordered [y, x] (i.e. y=rows, x=cols).

Debug/status printing is controlled by the module-scope "verbose" variable.  Set this to True to
get status reports in long-running routines (also controllable via the command-line interface).
"""

import os
import optparse
import numpy
try: # If astropy.io.fits is available use that
    import astropy.io.fits as pyfits
except:
    import pyfits
import galsim

verbose = False

def log(msg):
    if verbose:
        print msg

n_epochs = 6 # Number of epochs for each subfield when multiepoch or full

# This dict is used to provide the default keyword arguments to the galsim.InterpolatedImage
# instances used for coaddition.  Note that setting calculate_stepk/calculate_maxk would be a Very
# Bad Idea, since a) some of the assumptions that go into these (iterative) routines are broken by
# large images filled with objects; and b) they are unneccesary since all the interpolation used
# to make these coadds happens in real space
default_interp_kwds = {"calculate_stepk": False, "calculate_maxk": False}

# n_split_default - the input images are split into n_split x n_split subimages when coadding, to
# reduce RAM usage as there seems to be a non-trivial overhead in setting up the
# galsim.InterpolatedImage
# N.B. that n_split should divide into the GREAT3 image array dimensions exactly, and the code tests
# for this
n_split_default = 4


def mergeKeywordArgs(defaults, kwds):
    if kwds is None:
        return kwds
    result = defaults.copy()
    result.update(**kwds)
    return result

def main(field, sim_dir, work_dir, upsampling=1, stars=False, n_split=n_split_default,
         border=4, **interp_kwds):
    """Main driver for all routines in this file, and the implementation of most of
    the command-line interface.

    Arguments:

      field -------- subfield in the field to be processed, as an integer

      sim_dir ------ simulation directory, containing the GREAT3 images and catalogs for a single
                     branch

      work_dir ----- work directory; will contain the output coadd image on return

      upsampling --- resolution factor; output resolution is 1. / float(upsampling).  Currently
                     upsampling must be supplied as an integer.  Note that if the upsampling changes
                     then the overall normalization of the output coadd image currently changes too;
                     if you want to change this behaviour use the `normalization` keyword in
                     `interp_kwds`

      stars -------- make coadds for `starfield_image-[subfield_index]-[epoch_index].fits` files if
                     set `True`; otherwise make coadds for the galaxy image
                     `image-[subfield_index]-[epoch_index].fits` files

      n_split ------ n_split x n_split = number of sub "quadrants" to split the input image into
                     when coadding; useful to reduce memory when handling the large input arrays

      border ------- size of the border to use for overlap at the internal edges of the quadrants
                     defined above, which allows us to avoid errors at these quadrant edges due to
                     interpolant edge effects (e.g. galsim.InterpolatedImage extrapolates to zero
                     outside input image bounds).  Note that the size of border needed depends on
                     the x_interpolant being specified (or left as default) in the **interp_kwds
                     sent to galsim.InterpolatedImage

      interp_kwds -- keyword arguments passed unmodified to galsim.InterpolatedImage when shifting
                     the images
 
    """
    # First load up the epoch_index=0 file to get dimensions for coadd image
    if stars:
        inprefix = os.path.join(sim_dir, "starfield_image-")
        outprefix = os.path.join(work_dir, "coadd_starfield_image-")
    else: 
        inprefix = os.path.join(sim_dir, "image-")
        outprefix = os.path.join(work_dir, "coadd_image-")
    interp_kwds = mergeKeywordArgs(default_interp_kwds, interp_kwds)
    log("Reading image from "+inprefix+"%03d-%1d.fits" % (field, 0))
    # Load full image array into an ImageF instance
    image = galsim.ImageViewF(
        pyfits.getdata(inprefix+"%03d-%1d.fits" % (field, 0)).astype(numpy.float32), scale=1.)
    # Test that n_split is appropriately chosen
    if (image.array.shape[0] % n_split != 0) or (image.array.shape[1] % n_split != 0):
        raise ValueError(
            "Module scope variable n_split (="+str(n_split)+") does not divide into the "+
            "dimensions of the image "+inprefix+("03%d-%1d.fits" % (field, 0))+" exactly")
    # Setup the output coadd image with the right dimensions taken from this image, with scale
    # of 1 / upsampling
    coadd = galsim.ImageF(
        image.array.shape[1] * upsampling, image.array.shape[0] * upsampling, scale=1. / upsampling)
    # Loop over the epochs 
    for epoch in range(n_epochs):

        if image is None: # If image not already loaded (i.e. in first step above), load it
            log("Reading image from "+inprefix+"%03d-%1d.fits" % (field, epoch))
            image = galsim.ImageViewF(
                pyfits.getdata(inprefix+"%03d-%1d.fits" % (field, epoch)).astype(numpy.float32),
                scale=1.)
        # Then navigate around the n_split x n_split quadrants, coadding the relevant image parts
        for i in range(n_split):

            # Get the edges of the subimage to cutout / coadd to, and an offset to correct for
            # overlaps of size given by border... Note this, and the stuff for j below, took a bit
            # of testing to get exactly right!
            imin = 1 + i * image.array.shape[0] / n_split
            imax = (1 + i) * image.array.shape[0] / n_split
            imin_coadd = 1 + i * coadd.array.shape[0] / n_split
            imax_coadd = (1 + i) * coadd.array.shape[0] / n_split 
            ioffset_coadd = 0 # Offset to apply to correct for quadrant overlap
            if n_split != 1:
                if i == 0:
                    imax += border
                elif i > 0 and i < n_split - 1:
                    imin -= border
                    imax += border
                    ioffset_coadd = border * upsampling
                elif i == n_split - 1:
                    imin -= border
                    ioffset_coadd = border * upsampling
                else:
                    raise RuntimeError("Quadrant index i out of range")
            # Now to the inner loop...
            for j in range(n_split):

                # Get the bounds of the submiage quadrants for the input image and coadd
                # Note am writing this out explicitly in j as for i so that it's clear
                jmin = 1 + j * image.array.shape[1] / n_split
                jmax = (1 + j) * image.array.shape[1] / n_split
                jmin_coadd = 1 + j * coadd.array.shape[1] / n_split
                jmax_coadd = (1 + j) * coadd.array.shape[1] / n_split
                joffset_coadd = 0 # Offset to apply to correct for quadrant overlap
                if n_split != 1:
                    if j == 0:
                        jmax += border
                    elif j > 0 and j < n_split - 1:
                        jmin -= border
                        jmax += border
                        joffset_coadd = border * upsampling
                    elif j == n_split - 1:
                        jmin -= border
                        joffset_coadd = border * upsampling
                    else:
                        raise RuntimeError("Quadrant index j out of range")
                # Use i,j min/max to get the bounds of the bordered submiage we will extract 
                bounds = galsim.BoundsI(jmin, jmax, imin, imax)
                # Initialize this submimage into an interpim
                interpim = galsim.InterpolatedImage(image[bounds], **interp_kwds)
                # Get the dithers (do not use if stars, since these are differently dithered in
                # of the images not in the lower left of the FITS image)
                if not stars:
                    epoch_dither_x, epoch_dither_y = numpy.loadtxt(
                        os.path.join(sim_dir, "epoch_dither-%03d-%1d.txt" % (field, epoch)),
                        unpack=True)
                    # Apply this dither as a shift - NOTE SIGN HERE DUE TO DEFINITION OF DITHERS!
                    interpim.applyShift(-epoch_dither_x, -epoch_dither_y)
                # Draw this to the same tmp coadd of the right size to cover all of the input image
                # (note the factor of upsampling in this image dimensions)
                coadd_tmp = interpim.draw(
                    galsim.ImageF((jmax - jmin + 1) * upsampling, (imax - imin + 1) * upsampling),
                    dx=coadd.scale, use_true_center=True, add_to_image=True)
                # Set up the bounds in the final coadd image into which to write a subimage of
                # coadd_tmp...
                bounds_coadd = galsim.BoundsI(jmin_coadd, jmax_coadd, imin_coadd, imax_coadd)
                # Assign the correct pixels in coadd_tmp to the right region in coadd
                coadd[bounds_coadd] += coadd_tmp[
                    galsim.BoundsI(
                        1 + joffset_coadd, 1 + jmax_coadd - jmin_coadd + joffset_coadd,
                        1 + ioffset_coadd, 1 + imax_coadd - imin_coadd + ioffset_coadd)] 

        # Having done these n_split x n_split coadds, set the image to None for it to be read in
        # from the next epoch
        image = None

    # If the input image was a starfield, only output the lower left of the 3 x 3 array
    # (note it would have been quicker not to coadd the whole image first, but this is only a small
    # cost relative to galaxy field coadds for these much smaller starfield images...)
    # TODO: This assumption only holds for multiepoch, needs changing for full!
    if stars:
        coadd = coadd[galsim.BoundsI(1, coadd.array.shape[1] / 3, 1, coadd.array.shape[0] / 3)]

    # Normalize the coadd correctly so that flux is conserved
    coadd /= float(n_epochs)
    # Write the output
    log("Writing coadd output to "+outprefix+("%03d.fits" % field))
    coadd.setOrigin(0,0)
    coadd.write(outprefix+("%03d.fits" % field), clobber=True)
    return coadd

 
if __name__ == "__main__":
    usage = "usage: ./coadd_multiepoch.py [options] FIELD SIM_DIR WORK_DIR"
    description = """Coadd multiple epoch images using a specified galsim.Interpolant. 

FIELD is the number of the subfield to process
(starting at zero as in the filenames).  SIM_DIR is the directory
containing GREAT3 images and catalogs for the branch of interest. 
WORK_DIR is the directory where output files should be placed: it
will be created if it does not exist.
"""
    parser = optparse.OptionParser(usage=usage, description=description)
    parser.add_option(
        "--interpolant", dest="interpolant", type=str, default="galsim.Cubic(tol=1.e-4)",
        metavar="INTERPOLANT", help="String representation of the galsim.Interpolant that will "+
        "be used to perform the coaddition.  In Python, eval(INTERPOLANT) should yield a valid "+
        "galsim.Interpolant instance [default = 'galsim.Cubic(tol=1.e-4)']")
    parser.add_option(
        "--border", dest="border", type=int, metavar="BORDER", default=3,
        help="Size of overlap border in input pixels used to avoid galsim.InterpolatedImage edge "+
        "effects when splitting each input image into smaller sub-images.  The "+
        "galsim.InterpolatedImage class extrapolates to zero outside input image bounds, so not "+
        "including this border leads to visible artifacts.  Note that the size of the border "+
        "required depends on the support of the INTERPOLANT being specified [default = 3]")
    parser.add_option(
        "--space", dest="space", action="store_true", default=False,
        help="Images being input are 'space'-type.  If set, will increase "+
        "resolution in output coadds by factor of 2 along each dimension")
    parser.add_option(
        "--quiet", dest="quiet", action='store_true', default=False,
        help="Don't print progress statements")
    opts, args = parser.parse_args()
    try:
        field, sim_dir, work_dir = args
    except ValueError:
        parser.error("Exactly three positional arguments are required")
    if not os.path.isdir(sim_dir):
        parser.error("Input directory %s does not exist or is not a directory" % sim_dir)
    try:
        field = int(field)
    except TypeError:
        parser.error("Field argument '%s' is not an integer" % field)
    if not os.path.isdir(work_dir):
        os.makedirs(work_dir)
    if opts.space:
        upsampling = 2
    else:
        upsampling = 1
    if opts.quiet:
        verbose = False
    else:
        verbose = True
    interp_kwds = {}
    log("Using interpolant: "+opts.interpolant)
    interp_kwds['x_interpolant'] = eval(opts.interpolant) 
    main(field, sim_dir, work_dir, upsampling=upsampling, border=opts.border, **interp_kwds)
    main(
        field, sim_dir, work_dir, upsampling=upsampling, border=opts.border, stars=True,
        **interp_kwds)
