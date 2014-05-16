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
Example module for PSF modeling in GREAT3

This Python module provides a simple PSF modeling code adapted to work with the GREAT3 data
format.  It is intended to provide a "first-order" PSF module to allow GREAT3 participants to
compete in branches with variable PSFs without writing their own PSF estimation module, and to
provide an example for those who do wish to create their own PSF estimation module.  It is not
expected to be competitive with the current state-of-the-art in PSF estimation, but it should not
be orders of magnitude worse, either.


DEPENDENCIES

Python 2.6.x-2.7.x
NumPy 1.6+
GalSim 1.0
PyFITS 1.0+
Matplotlib 0.10+ (optional; used for diagnostic plots)


ALGORITHM / USAGE

The algorithmic approach is a simple PCA of star images, interpolated using Chebyshev polynomials;
this is similar to the PSF modeling approach used in the SDSS Photo reduction pipeline and in
Jee & Tyson (2011).  It should not be confused with the spatial PCA of Jarvis and Jain (2004).
Overall, the approach is:

 1) Measure the sub-pixel centroids and fluxes of stars.  This is done using the HSM adaptive
    moments code in GalSim, which is overkill (as it measures shapes too, which we don't need)
    but effective nonetheless.  At the same time, we reorganize the from subfield-based files
    to tile-based files.  This produces 'star_meas' files.  See measureStars().

 2) Resample and rescale the images of input stars such that they have a common center and unit
    flux, using the measurements from the previous step.  In this step, we loop over each
    tile-based star_meas file and then over each subfield-based starfield_image file, extracting
    only the stars associated with the current tile.  This results in more I/O but lower memory
    consumption.  The outputs of this are the `star_data` files, which are FITS binary tables
    that contain the centered, normalized postage stamps as a 2-d array field.  See
    buildDataTables().

 3) Compute an orthonormal basis for the star images using Principle Components Analysis
    (see computePCA).  The PCA is computed after subtracting the mean star image, and we
    keep only the first few basis images, which represent most of the variance in the sample.
    A separate basis is computed for each tile.

 4) Project each star onto the PCA basis, and fit a 2-d Chebyshev polynomial to the resulting
    coefficients as a function of position on the tile (LinearChebyshevModel.fit).  Note that
    there is a different Chebyshev polynomial for each basis image, so the model image at a
    point is computed by evaluating the Chebyshev polynomials at that point, and using the
    those values as the weights in a weighted sum of the basis images (which also includes the
    mean star image, unweighted).  A separate model is computed for each tile, but these can be
    be gathered into a FieldModelSuite object to represent the entire field.

 5) (optional) Evaluate the PSF model at the position of each galaxy.  This can be done en-mass
    using FieldModelSuite.makeImageGrid, which creates a grid of images for a subfield galaxy
    catalog.  It can also be done one galaxy at a time, using FieldModelSuite.evaluate().

Steps 3 and 4 can be done for a full field by calling FieldModelSuite.fit, and all 5 steps can
be performed by calling the main() function or executing this file as a script.

The LinearChebyshevModel class also has a few inspect* methods that can be used to visualize
the goodness of fit of the model given the per-tile star_data catalog it was constructed from.

Executing the script looks something like this:

./psf_models.py --model-file psf_models.p 40 /sims/variable_psf/ground/constant output

This will process subfields 40-59 (all of the subfields in a single field) for the branch located in
/sims/variable_psf/ground/constant, placing all outputs in the 'output' directory.  The
FieldModelSuite object will be saved to output/psf_models.p for future use (in addition to creating
PSF image grids that correspond to the galaxy image grids).  For more help with the command-line
interface, you can just run:

./psf_models.py --help


KNOWN ISSUES / LIMITATIONS

 - The code has only been tested on 'variable_psf' branches, and may require some modification to
   work on 'full' branches.  Some of these may be modestly backwards-incompatible.

 - Constructing the PSF models for an entire field while attempting to use all the stars can consume
   a significant amount of memory (nearly 16 GB) and take ~3 hours on a single CPU.  Using a smaller
   image dimension for the PSF model images or using the max_stars parameter can be used to adress
   this; just setting max_stars=1000 (the default) reduces the memory footprint to <2GB and the
   single-CPU runtime to ~22 minutes.  Once constructed, a PSF model suite is generally quite
   lightweight unless a very large number of basis functions and/or high Chebyshev degree is used.

 - In the GREAT3 multi-epoch space branches the individual images are not Nyquist sampled, which
   invalidates the simple centroiding and interpolation algorithms use to construct the centered,
   normalized star data tables here.  For correct handling of this branch, these should be
   replaced (after which code should be valid: the PCA and Chebyshev spatial fit make no
   assumptions about the sampling, and indeed should be valid even if the star images are
   upsampled relative to the original data).

 - The PSFs in GREAT3 all have more than 4 degrees of freedom, which is the number of PCA
   basis images we use by default here.  At least some of these additional degrees of freedom
   are likely important, but it is not clear whether using PCA for dimensionality reduction
   is a valid way to obtain these additional basis functions, as the basis images typically get
   much noisier after the first four.

 - While polynomial spatial interpolation is a common choice today, the true spatial variation
   of realistic PSFs (including those simulated in GREAT3) cannot be fit exactly using low-order
   polynomials.

The classes in this file are designed to be modular; you should be able to replace the spatial
fitting component while retaining the shifting/normalization and PCA basis code, for instance,
or use basis images derived some other way with the Chebyshev spatial model.


NOTES

As is usual with NumPy, 2-d arrays representing images are ordered [y,x] (i.e. y=rows, x=cols).

NumPy conventions for 2-d polynomials (which we also adopt) reverse this, however, so an array of
2-d polynomial coefficients will be ordered [x,y].

Debug/status printing is controlled by the module-scope "verbose" variable.  Set this to True to
get status reports in long-running routines (also controllable via the command-line interface).
"""

import os
import re
import optparse
import sys
import numpy

try:
    import astropy.io.fits as pyfits
except:
    import pyfits

import galsim

try:
    import cPickle as pickle
except ImportError:
    import pickle

try:
    from numpy.polynomial.chebyshev import chebval2d
except ImportError:
    # chebval2d only available in 1.7+; we implement
    # a more limited replacement (scalar x,y only)
    # if that's not available
    import numpy.polynomial.chebyshev
    def chebval2d(x, y, c):
        vy = numpy.zeros(c.shape[0], dtype=float)
        for j in range(c.shape[0]):
            vy[j] = numpy.polynomial.chebyshev.chebval(y, c[j,:])
        return numpy.polynomial.chebyshev.chebval(x, vy)

nsubfields = 20  # There are 20 subfields, but this lets you easily change how many
                 # get processed by this program.

verbose = False

def log(msg):
    if verbose:
        print msg

obs_type_bounds = {"space": (0.0, 0.5, 0.0, 0.5),
                   "ground": (0.0, 2.0, 0.0, 2.0)}

obs_type_basis_size = {"space": 8, "ground": 4}

# This dict is used to provide the default keyword arguments to galsim.InterpolatedImage
# when resampling stars to center them.  It can be modified here to avoid having to
# set keywords repeatedly when calling buildDataTables() or main().
default_interp_kwds = {"calculate_stepk": False, "calculate_maxk": False}

def mergeKeywordArgs(defaults, kwds):
    if kwds is None:
        return defaults.copy()
    result = defaults.copy()
    result.update(**kwds)
    return result

def makeDataSchema(nx, ny):
    """A numpy.dtype object to be used for the structured array "data" argument to
    LinearChebyshevModel.fit() and computePCA().

    Arguments:

      nx ---------- size of star images (and PSF model images) in the x direction

      ny ---------- size of star images (and PSF model images) in the y direction

    Both nx and ny must be odd (and are usually the same).

    The fields of the returned dtype object are:

      image ------- float64, shape (ny, nx): centered star image, normalized to unit flux

      weight ------ float64, weight for the star relative to others in the data table

      x,y --------- float64, position of the star, in the coordinates the spatial interpolation
                    will be done; usually tile_[x,y]_pos_deg
      g1,g2 ------- float64, observed_shape.g1 and observed_shape.g2 adaptive moments ellipticity
                    estimates for the star as determined by HSM's FindAdaptiveMom(), taken from the
                    star_meas catalog
      sigma ------- float64, the moments_sigma determined by HSM's FindAdaptiveMom(), taken from the
                    star_meas catalog, in pixels
    """
    if nx % 2 == 0 or ny % 2 == 0:
        raise ValueError("Image dimensions (%d, %d) must be odd" % (nx, ny))
    return numpy.dtype([("image", float, (ny, nx)),
                        ("weight", float), ("x", float), ("y", float),
                        ("g1", float), ("g2", float), ("sigma", float), ])

class LinearChebyshevModel(object):
    """A spatial interpolation class for PSF models: it assumes the model is a linear combination
    of basis images at each point, with the coefficients of the basis images varying as Chebyshev
    polynomials of the first kind.
    """

    @staticmethod
    def _get_transform_parameters(lower, upper):
        offset = (lower + upper) / (lower - upper)
        scale = 2.0 / (upper - lower)
        return offset, scale

    def __init__(self, basis, coeff, bounds, image0=None, normalize=True, clip=True):
        """Construct the spatial model.

        Usually this constructor will not be called directly, as the fit() static method will
        generally be used instead to compute the Chebyshev coefficients from a table of star data.

        Arguments:

          basis --------- a 3-d array containing the basis images, with shape
                          (basis_size, ny, nx), where basis_size is the number of basis images,
                          and (ny, nx) is the shape of each image.

          coeff --------- a 3-d array of Chebyshev polynomial coefficients, with shape
                          (basis_size, cx, cy), where basis_size is the number of basis images,
                          and (cx, cy) is a matrix of Chebyshev coefficients.

          bounds -------- tuple of (xmin, xmax, ymin, ymax) indicating the region the spatial
                          model is valid.

          image0 -------- A static image to add to the model (useful for including the mean image
                          when basis functions represent deviations from the mean).

          normalize ----- Sets the default for the normalize argument to __call__: whether to
                          rescale returned images so that they sum to exactly 1.

          clip ---------- Sets the default for the clip argument to __call__: whether to reset
                          negative pixels to zero before returning model images (clipping is done
                          before normalization, if applicable).
        """
        self.basis = basis
        self.coeff = coeff
        if self.basis.shape[0] != self.coeff.shape[0]:
            raise ValueError(("First dimension of basis (%d) does not match first dimension of "
                              "coeff (%d)") % (self.basis.shape[0], self.coeff.shape[0]))
        self.xmin = float(bounds[0])
        self.xmax = float(bounds[1])
        self.ymin = float(bounds[2])
        self.ymax = float(bounds[3])
        self.x_offset, self.x_scale = self._get_transform_parameters(self.xmin, self.xmax)
        self.y_offset, self.y_scale = self._get_transform_parameters(self.ymin, self.ymax)
        self.image0 = image0
        if self.image0 is not None and self.image0.shape != self.basis.shape[1:]:
            raise ValueError("Shape of image0 %s does not match shape of basis images %s"
                             % (self.image0.shape, self.basis.shape[1:]))
        self.clip = bool(clip)
        self.normalize = bool(normalize)

    def __getinitargs__(self):
        return (self.basis, self.coeff, (self.xmin, self.xmax, self.ymin, self.ymax),
                self.image0, self.clip, self.normalize)

    def evaluate(self, x, y, normalize=None, clip=None):
        """Evaluate the spatial model at the given point, returning a spatially-weighted linear
        combination of the basis images.

        Arguments:

          x ------------- X position at which to evaluate the model

          x ------------- Y position at which to evaluate the model

          normalize ----- Whether to rescale returned images so that they sum to exactly 1.
                          Passing None defers to the default set in __init__.

          clip ---------- Whether to reset negative pixels to zero before returning model images
                          (clipping is done before normalization, if applicable).  None defers to
                          the default set in __init__.

        """
        if x < self.xmin or x > self.xmax:
            raise ValueError("x coordinate %s out of range: should be between %s and %s"
                             % (x, self.xmin, self.xmax))
        if y < self.ymin or y > self.ymax:
            raise ValueError("y coordinate %s out of range: should be between %s and %s"
                             % (y, self.ymin, self.ymax))
        xt = self.x_scale * x + self.x_offset
        yt = self.y_scale * y + self.y_offset
        result = numpy.zeros(self.basis.shape[1:], dtype=float)
        for i in range(self.basis.shape[0]):
            result += self.basis[i] * chebval2d(xt, yt, self.coeff[i])
        if self.image0 is not None:
            result += self.image0
        if clip or clip is None and self.clip:
            result[result < 0.0] = 0.0
        if normalize or normalize is None and self.normalize:
            result /= result.sum()
        return result

    def inspectStar(self, record):
        """Display model and residual images for a single star.

        Given a record from a data table (i.e. the one passed to the fit method), that 'data'
        image will be displayed, along with model and residual pairs for each of:

          image0 ----- Constant component of the model; usually the mean obtained before computing
                       a mean-subtracted PCA basis.

          proj ------- Direct projection of the basis onto the star image, with no concern for
                       the spatial interpolation (also includes the image0 term).  This will
                       generally be the best fit to the data, as it does not include the spatial
                       smoothing from the Chebyshev polynomials.

          interp ----- The full PSF model, obtained by evaluating the basis coefficients from the
                       Chebyshev spatial polynomial at the position of this star (also includes
                       the image0 term).  This is what is returned by the evaluate() method.

        Note that all images will be centered and normalized to unit flux, because
        the data is constructed after shifting and rescaling the stars.
        """
        try:
            import matplotlib
        except:
            raise ImportError('Could not find matplotlib, which is needed to display PSF '
                              'model diagnostics!')
        interp_image = self.evaluate(record['x'], record['y'])
        interp_res = record['image'] - interp_image
        if self.image0 is None:
            n_rows = 3
            zero_mean_data = record['image']
            worst_res = interp_res
        if self.image0 is not None:
            n_rows = 4
            image0_res = record['image'] - self.image0
            zero_mean_data = image0_res
            worst_res = image0_res
        proj_coeff = self.basis.reshape(self.basis.shape[0],-1).dot(zero_mean_data.flatten())
        proj_image = numpy.tensordot(proj_coeff, self.basis, axes=1)
        proj_res = record['image'] - proj_image
        if self.image0 is not None:
            proj_image += self.image0
            proj_res -= self.image0
        kwds1 = dict(interpolation='nearest', origin='lower', cmap=matplotlib.cm.Blues,
                     vmin=record['image'].min()-1E-8, vmax=record['image'].max()+1E-8)
        kwds2 = dict(interpolation='nearest', origin='lower', cmap=matplotlib.cm.coolwarm_r,
                     vmin=worst_res.min(), vmax=worst_res.max())
        def addPlot(image, title, idx, **kw):
            ax = fig.add_subplot(2, n_rows, idx)
            mappable = ax.imshow(image, **kw)
            ax.set_title(title)
            ax.axis("off")
            return ax, mappable
        fig = matplotlib.pyplot.figure(figsize=(10,5))
        addPlot(record['image'], 'data', 1, **kwds1)
        index = 2
        if self.image0 is not None:
            addPlot(self.image0, 'image0', 2, **kwds1)
            addPlot(image0_res, 'image0 res', 2 + n_rows, **kwds2)
            index = 3
        addPlot(proj_image, 'proj', index, **kwds1)
        addPlot(proj_res, 'proj res', index + n_rows, **kwds2)
        index += 1
        ax1, mappable1 = addPlot(interp_image, 'interp', index, **kwds1)
        ax2, mappable2 = addPlot(interp_res, 'interp res', index + n_rows, **kwds2)
        fig.subplots_adjust(top=0.925, bottom=0.025, left=0.025, right=0.85,
                            wspace=0.025, hspace=0.15)
        box1 = ax1.get_position()
        box2 = ax2.get_position()
        cax1 = fig.add_axes([box1.x1+0.025, box1.y0, 0.025, box1.height])
        cax2 = fig.add_axes([box2.x1+0.025, box2.y0, 0.025, box2.height])
        fig.colorbar(mappable1, cax=cax1)
        fig.colorbar(mappable2, cax=cax2)
        fig.canvas.draw()
        return fig

    def inspectSpatial(self, data, func=lambda r: (r**2).mean(), **kwds):
        """Create a colormapped scatter plot showing some function of the residuals between the
        star data and the interpolated PSF model at the same points.

        Arguments:

          data ------- a structure NumPy array data table with dtype as provided by
                       makeDataSchema().  Usually the same one passed to the fit() method when
                       constructing the model.

          func ------- a function that takes a single residual image (NumPy array) and returns a
                       single value.  The default computes the RMS.

          **kwds ----- passed unmodified to matplotlib's 'scatter' routine.

        """
        try:
            from matplotlib import pyplot
        except:
            raise ImportError('Could not find matplotlib, which is needed to display PSF '
                              'model diagnostics!')
        v = numpy.zeros(data.size, dtype=float)
        for i, record in enumerate(data):
            r = self.evaluate(record['x'], record['y'])
            r -= record['image']
            r *= record['weight']
            v[i] = func(r)
        fig = pyplot.figure()
        ax = fig.add_subplot(1,1,1)
        ax.scatter(data['x'], data['y'], c=v, **kwds)
        fig.canvas.draw()
        return fig

    def inspectSigmaSpatial(self, data, **kwds):
        """Create three colormapped scatter plot showing HSM-measured moments_sigma values for data,
        model and residuals.

        These are evaluated at the location of each star in the input data table.

        Arguments:

          data ------- a structure NumPy array data table with dtype as provided by
                       makeDataSchema().  Usually the same one passed to the fit() method when
                       constructing the model.

          **kwds ----- passed unmodified to matplotlib's 'scatter' routine.

        """
        try:
            from matplotlib import pyplot
        except:
            raise ImportError('Could not find matplotlib, which is needed to display PSF '
                              'model diagnostics!')
        v = numpy.zeros(data.size, dtype=float)
        fig = pyplot.figure()
        ax = fig.add_subplot(1,1,1)
        sc = ax.scatter(data['x'], data['y'], c=data['sigma'], **kwds)
        pyplot.colorbar(sc)
        fig.canvas.draw()
        return fig

    def inspectWhisker(self, data, label="", **kwds):
        """Create a whisker plot showing the HSM-measured observed_shape.g1 and observed_shape.g2
        at the location of each star.

        Arguments:

          data ------- a structure NumPy array data table with dtype as provided by
                       makeDataSchema().  Usually the same one passed to the fit() method when
                       constructing the model.

          label ------ a string with which to label the plots, which will be given the title
                       label+" Data", label+" Model", and label+" (Data - Model) Residuals"

          **kwds ----- passed unmodified to matplotlib's 'quiver' routine.

        """
        try:
            from matplotlib import pyplot
        except:
            raise ImportError('Could not find matplotlib, which is needed to display PSF '
                              'model diagnostics!')
        # Get the data shapes into polars
        gmag = numpy.sqrt(data['g1']**2 + data['g2']**2)
        gphi = .5 * numpy.arctan2(data['g2'], data['g1'])
        fig = pyplot.figure(figsize=(10, 11))
        ax = fig.add_subplot(2,2,1)
        # Set the necessary default kwargs for making a nice whisker plot
        kwds["headwidth"] = 0.
        kwds["headlength"] = 0.
        kwds["headaxislength"] = 0.
        kwds["pivot"] = "mid"
        # First do the Data plot
        # Call the quiver routine, and plot a key
        qq = pyplot.quiver(
            data['x'], data['y'], gmag * numpy.cos(gphi), gmag * numpy.sin(gphi), **kwds)
        xran = max(data['x']) - min(data['x']) # Get the x and y ranges
        yran = max(data['y']) - min(data['y']) #
        qqkey = pyplot.quiverkey( # Put this centrally above the whisker plot, a bit above whiskers
            qq, numpy.mean(data['x']), max(data['y']) + 0.10 * yran, 0.05, "|g| = 0.05",
            coordinates="data")
        pyplot.xlabel('x [deg]')
        pyplot.ylabel('y [deg]')
        pyplot.title(label+' Data')
        # Sort out the limits so that we don't have too-early quiver clipping by borders
        pyplot.xlim(min(data['x']) - 0.05 * xran, max(data['x']) + 0.05 * xran)
        pyplot.ylim(min(data['y']) - 0.05 * yran, max(data['y']) + 0.20 * yran)
        # Then do the Model plot, requires calculation
        model_shapes = numpy.zeros( # Set up storage table
            len(data['x']), dtype=numpy.dtype([("g1", float), ("g2", float)]))
        for model_record, data_record in zip(model_shapes, data):

            hsm_results = (
                galsim.ImageViewD(self.evaluate(data_record['x'], data_record['y']))
                ).FindAdaptiveMom()
            model_record['g1'] = hsm_results.observed_shape.g1
            model_record['g2'] = hsm_results.observed_shape.g2

        gmag = numpy.sqrt(model_shapes['g1']**2 + model_shapes['g2']**2)
        gphi = .5 * numpy.arctan2(model_shapes['g2'], model_shapes['g1'])
        ax = fig.add_subplot(2,2,2)
        # Call the quiver routine, and plot a key
        qq = pyplot.quiver(
            data['x'], data['y'], gmag * numpy.cos(gphi), gmag * numpy.sin(gphi),
            **kwds)
        qqkey = pyplot.quiverkey( # Put this centrally above the whisker plot, a bit above whiskers
            qq, numpy.mean(data['x']), max(data['y']) + 0.10 * yran, 0.05, "|g| = 0.05",
            coordinates="data")
        pyplot.xlabel('x [deg]')
        pyplot.ylabel('y [deg]')
        pyplot.title(label+' Model')
        # Sort out the limits so that we don't have too-early quiver clipping by borders
        pyplot.xlim(min(data['x']) - 0.05 * xran, max(data['x']) + 0.05 * xran)
        pyplot.ylim(min(data['y']) - 0.05 * yran, max(data['y']) + 0.20 * yran)
        # Then do the residuals plot
        gmag = numpy.sqrt(
            (data['g1'] - model_shapes['g1'])**2 + (data['g2'] - model_shapes['g2'])**2)
        gphi = .5 * numpy.arctan2(data['g2'] - model_shapes['g2'], data['g1'] - model_shapes['g1'])
        ax = fig.add_subplot(2,2,3)
        qq = pyplot.quiver(
            data['x'], data['y'], gmag * numpy.cos(gphi), gmag * numpy.sin(gphi),
            **kwds)
        qqkey = pyplot.quiverkey( # Put this centrally above the whisker plot, a bit above whiskers
            qq, numpy.mean(data['x']), max(data['y']) + 0.10 * yran, 0.01, "|g| = 0.01",
            coordinates="data")
        pyplot.xlabel('x [deg]')
        pyplot.ylabel('y [deg]')
        pyplot.title(label+' (Data - Model) Residuals')
        # Sort out the limits so that we don't have too-early quiver clipping by borders
        pyplot.xlim(min(data['x']) - 0.05 * xran, max(data['x']) + 0.05 * xran)
        pyplot.ylim(min(data['y']) - 0.05 * yran, max(data['y']) + 0.20 * yran)
        fig.canvas.draw()
        return fig

    @classmethod
    def fit(cls, basis, data, dim, degree, bounds, **kwds):
        """Fit a orthonormal image basis with Chebyshev spatial variation from a star data table.

        Arguments:

          basis --------- a 3-d array containing the basis images, with shape (basis_size, ny, nx),
                          where basis_size is the number of basis images, and (ny, nx) is the shape
                          of each image.

          data ---------- a structured NumPy array with dtype of the type returned by
                          makeDataSchema: fields include 'image' (2-d float), 'weight', 'x', 'y'.

          dim ----------- size of the PSF model images (on a side).  Must be odd.

          degree -------- maximum combined degree of the 2-d Chebyshev polynomial
                          (degree_x + degree_y <= degree)

          bounds -------- tuple of (xmin, xmax, ymin, ymax) indicating the region the spatial model
                          is valid.

          **kwds -------- passed unmodified to __init__
        """
        if basis.shape[1:] != (dim, dim):
            raise ValueError("Basis shape %s does not match provided dimension %s"
                             % (basis.shape[1:], dim))

        # We construct the model object we'll return first, and fill in its coefficients later.
        coeff = numpy.zeros((basis.shape[0], degree+1, degree+1), dtype=float)
        model = cls(basis, coeff, bounds, **kwds)

        # Start by projecting data images onto the basis (which we assume to be orthonormal).
        # At the same time we'll map the data x,y coordinates to the (-1,1) window appropriate
        # for Chebyshevs
        proj_dtype = numpy.dtype([("coeff", float, basis.shape[0]), ("weight", float),
                                  ("xt", float), ("yt", float)])
        proj = numpy.zeros(data.size, dtype=proj_dtype)
        proj['xt'] = model.x_scale * data['x'] + model.x_offset
        proj['yt'] = model.y_scale * data['y'] + model.y_offset
        proj['weight'] = data['weight']
        proj_matrix = basis.reshape(basis.shape[0], -1)
        for n in range(data.size):
            # The image shape should already be (dim, dim), but older versions of pyfits
            # can screw it up on the round trip through writing and reading. 
            # So reshape it just to be sure.
            zero_mean_image = data[n]['image'].reshape((dim,dim)).copy()
            if model.image0 is not None:
                zero_mean_image -= model.image0
            proj[n]['coeff'][:] = proj_matrix.dot(zero_mean_image.flatten())

        # Now we fit the projected basis coefficients with the spatial Chebyshev.  We only
        # include terms in the lower triange (i.e. degree_x+degree_y <= degree), and pack
        # x and y into the second dimension of the cheby_m matrix.  Note that we can reuse
        # the same matrix for all basis functions, as these are completely independent
        # (thanks to the orthonormality of the basis).
        cheby_m = numpy.zeros((data.size, (degree+1)*(degree+2)/2), dtype=float)
        cheby_x = numpy.polynomial.chebyshev.chebvander(proj['xt'], degree)
        cheby_y = numpy.polynomial.chebyshev.chebvander(proj['yt'], degree)
        im = 0
        for ix in range(1+degree):
            for iy in range(0, 1+degree-ix):
                cheby_m[:,im] = cheby_x[:,ix] * cheby_y[:,iy]
                im += 1

        # Apply the weights to the matrix and the data vectors
        cheby_m *= proj['weight'][:,numpy.newaxis]
        proj['coeff'] *= proj['weight'][:,numpy.newaxis]

        # Linear least squares (with a 2-d matrix on the rhs, columns are the basis elements)
        fit_coeff = numpy.zeros((cheby_m.shape[1], basis.shape[0]), dtype=float)
        for i in range(basis.shape[0]):
            fit_coeff[:,i], chisq, _, _ = numpy.linalg.lstsq(cheby_m, proj['coeff'][:,i])

        # Now we set the lower triangle of model.coeff attribute by unpacking the best-fit coeffs
        im = 0
        for ix in range(1+degree):
            for iy in range(0, 1+degree-ix):
                model.coeff[:,ix,iy] = fit_coeff[im,:]
                im += 1

        return model

def computePCA(data, dim, basis_size=4, weighted_mean=True, weighted_pca=True):
    """Create an image basis appropriate for use with LinearChebyshevModel using
    Principle Components Analysis.

    Arguments:

      data ------------ a structured NumPy array with dtype of the type returned by
                        makeDataSchema: fields include 'image' (2-d float), 'weight', 'x', 'y'.

      dim ------------- size of the PSF model images (on a side).  Must be odd.

      basis_size ------ number of basis functions to keep

      weighted_mean --- use the weights in the data array when computing the mean image

      weighted_pca ---- use the weights in the data array when computing the PCA

    Returns a tuple of (basis, image0), which can be used directly as inputs to
    LinearChebyshevModel's constructor or fit method.
    """
    data_flat = data['image'].reshape(data.size, -1).copy()
    weight_sum = data['weight'].sum()
    if weighted_mean:
        mean_flat = (data_flat * data['weight'][:,numpy.newaxis]).sum(axis=0) / weight_sum
    else:
        mean_flat = data_flat.mean(axis=0)
    data_flat -= mean_flat
    if weighted_pca:
        data_flat *= (data['weight'] / weight_sum)[:,numpy.newaxis]
    u, s, vt = numpy.linalg.svd(data_flat, full_matrices=False)
    basis_flat = vt[:basis_size,:]
    #image_shape = data['image'].shape[1:]   # This doesn't always work!
    image_shape = (dim, dim)
    basis = basis_flat.reshape(basis_size, *image_shape)
    image0 = mean_flat.reshape(*image_shape)
    return basis, image0

def measureStars(field, sim_dir, work_dir):
    """Measure the sub-pixel centroids and fluxes of star images, and reorganize
    the outputs by tile instead of subfield.

    We use galsim.hsm to do the measurement (this does more than we need, as it measures
    shapes as well as centroids and fluxes, but it's convenient).

    Writes "star_meas-[x_tile_index]-[y_tile_index].fits" FITS binary tables to the work directory.

    Arguments:

      field -------- first subfield in the field to be processed, as an integer (should be a
                     multiple of 20)

      sim_dir ------ simulation directory, containing the GREAT3 images and catalogs for a single
                     branch

      work_dir ----- work directory; will contain the output star_meas catalogs on return

    """
    log("Measuring star images and reorganizing into tiles")
    schema = numpy.dtype([("subfield", int), ("x", float), ("y", float),
                          ("tile_x", float), ("tile_y", float),
                          ("dx", float), ("dy", float), ("flux", float),
                          ("g1", float), ("g2", float), ("sigma", float), ])
    tiles = {}
    for subfield in range(field, field + nsubfields):
        log("  Measuring subfield %d" % subfield)
        sim_cat = pyfits.getdata(os.path.join(sim_dir, "star_catalog-%03d.fits" % subfield))
        image = galsim.fits.read(os.path.join(sim_dir, "starfield_image-%03d-0.fits" % subfield))
        meas_cat = numpy.zeros(sim_cat.size, dtype=schema)
        # it's a little unfortunate it's this difficult to compute the bounding box of a star
        # postage stamp given the record; apparently we don't record the size of the boxes
        # directly anywhere in the challenge data itself
        dx0 = 0 - int(sim_cat['x'].min())  # add this to x to get min x coordinate
        dy0 = 0 - int(sim_cat['y'].min())  # add this to y to get min y coordinate
        dx1 = 1 - dx0   # add this to x to get max x coordinate
        dy1 = 1 - dy0   # add this to y to get max y coordinate
        skipped = {}
        for sim_record, meas_record in zip(sim_cat, meas_cat):
            bounds = galsim.BoundsI(int(sim_record['x'] + dx0), int(sim_record['x'] + dx1),
                                    int(sim_record['y'] + dy0), int(sim_record['y'] + dy1))
            subimage = image[bounds]
            meas_record['subfield'] = subfield
            meas_record['x'] = sim_record['x']
            meas_record['y'] = sim_record['y']
            meas_record['tile_x'] = sim_record['tile_x_pos_deg']
            meas_record['tile_y'] = sim_record['tile_y_pos_deg']
            try:
                # Use HSM algorithms in GalSim to compute centroid and flux
                meas = subimage.FindAdaptiveMom(guess_x_centroid=float(sim_record['x']),
                                                guess_y_centroid=float(sim_record['y']))
                meas_record['dx'] = meas.moments_centroid.x - sim_record['x']
                meas_record['dy'] = meas.moments_centroid.y - sim_record['y']
                meas_record['flux'] = meas.moments_amp
                meas_record['g1'] = meas.observed_shape.g1
                meas_record['g2'] = meas.observed_shape.g2
                meas_record['sigma'] = meas.moments_sigma
                # Add meas_record to a tile-indexed dict if successful
                tile_index = (int(sim_record['x_tile_index']), int(sim_record['y_tile_index']))
                tiles.setdefault(tile_index, []).append(meas_record)
            except RuntimeError as err:
                # sometimes HSM fails; maybe other things too?
                msg = str(err)
                skipped[msg] = skipped.get(msg, 0) + 1
        if skipped:
            log("    skipped %d of %d records due to the following errors:" %
                (sum(skipped.itervalues()), sim_cat.size))
            for msg, count in skipped.iteritems():
                log("      %d records: %s" % (count, msg.strip()))
    log("  Reorganizing into tiles")
    if not os.path.isdir(work_dir):
        os.makedirs(work_dir)
    for tile_index, records in sorted(tiles.items()):
        log("  Sorting and writing catalog for tile %d,%d" % tile_index)
        records.sort(key=lambda r: r['flux'], reverse=True)
        data_array = numpy.array(records, dtype=schema)
        out_file = os.path.join(work_dir, "star_meas-%03d-%02d-%02d.fits" % ((field,) + tile_index))
        pyfits.writeto(out_file, data_array, clobber=True)


def buildDataTables(field, sim_dir, work_dir, dim, max_stars=1000, interp_kwds=None):
    """Loop over all subfields in a field, building a tile-indexed dictionary of star data tables
    appropriate for use with computePCA and LinearChebyshevModel.fit.  The main work done by this
    routine is the shifting and rescaling necessary to center the star images and normalize their
    fluxes to unity.

    Writes "star_data-[x_tile_index]-[y_tile_index].fits" FITS binary tables to the work directory.

    Requires measureStars() to have been run on the same directories previously.

    Arguments:

      field -------- first subfield in the field to be processed, as an integer (should be a
                     multiple of 20)

      sim_dir ------ simulation directory, containing the GREAT3 images and catalogs for a single
                     branch

      work_dir ----- work directory; contains the star_meas catalogs used as inputs and will
                     contain the output star_data catalogs on return

      dim ---------- size of the PSF model images (on a side).  Must be odd.

      max_stars ---- maximum number of stars to use per tile (None = no limit); the highest SNR
                     stars will be kept without regard for their spatial distribution.

      interp_kwds -- keyword arguments passed unmodified to galsim.InterpolatedImage when shifting
                     the images.

    """
    log("Shifting star images and building data tables")
    interp_kwds = mergeKeywordArgs(default_interp_kwds, interp_kwds)
    radius = dim // 2
    assert radius * 2 + 1 == dim
    schema = makeDataSchema(dim, dim)
    regex = re.compile("star_meas-%03d-(\d\d)-(\d\d).fits" % field)

    tile_index_list = []
    meas_cat_list = []
    data_cat_list = []

    for meas_file in sorted(os.listdir(work_dir)):
        m = regex.match(meas_file)
        if not m: continue
        tile_index = (int(m.group(1)), int(m.group(2)))

        meas_cat = pyfits.getdata(os.path.join(work_dir, meas_file))

        # truncate the catalog if desired
        if max_stars is not None:
            stop = min(max_stars, meas_cat.size)
            meas_cat = meas_cat[:stop].copy()

        # Allocate the output array, and create an iterator to it
        # so we can iterate over it while processing multiple subfields.
        # This means meas_cat[i] won't correspond to sim_cat[i].
        data_cat = numpy.zeros(meas_cat.size, dtype=schema)

        tile_index_list.append(tile_index)
        meas_cat_list.append(meas_cat)
        data_cat_list.append(data_cat)
    n_tiles = len(tile_index_list)

    # The subfield images are pretty huge (especially for space) so we want to make sure
    # we only need to read the once each.  So make the subfield loop the outer loop.
    for subfield in range(field, field + nsubfields):
        log("  Shifting images for subfield %d" % subfield)

        image = galsim.fits.read(
            os.path.join(sim_dir, "starfield_image-%03d-0.fits" % subfield))
        pix_scale_deg = (image.scale * galsim.arcsec) / galsim.degrees

        for i_tile in range(n_tiles):

            tile_index = tile_index_list[i_tile]
            meas_cat = meas_cat_list[i_tile]
            data_cat = data_cat_list[i_tile]

            data_cat_iter = iter(data_cat)
            sub_meas_cat = meas_cat[meas_cat['subfield'] == subfield]
            if sub_meas_cat.size == 0:
                log("    skipping tile %d,%d; no records" % tile_index)
                continue

            for meas_record, data_record in zip(sub_meas_cat, data_cat_iter):
                bounds = galsim.BoundsI(int(numpy.floor(meas_record['x'])) - radius,
                                        int(numpy.ceil(meas_record['x'])) + radius,
                                        int(numpy.floor(meas_record['y'])) - radius,
                                        int(numpy.ceil(meas_record['y'])) + radius)
                subimage = image[bounds]

                # Use InterpolatedImage to shift so we're centered on the center pixel
                # and have unit flux.
                interp = galsim.InterpolatedImage(subimage, flux=meas_record['flux'], **interp_kwds)
                interp.applyShift(dx=-meas_record['dx']*image.scale,
                                  dy=-meas_record['dy']*image.scale)
                interp.setFlux(1.0)
                out_image = galsim.ImageView[numpy.float64](data_record['image'],
                                                            xmin=-radius, ymin=-radius,
                                                            scale=image.scale)
                interp.draw(image=out_image)

                # Compute the true position of the star in degrees from the origin of the tile.
                # The measured sub-pixel offset almost certainly doesn't matter, but we'll include
                # it for completeness.
                data_record['x'] = meas_record['tile_x'] + meas_record['dx'] * pix_scale_deg
                data_record['y'] = meas_record['tile_y'] + meas_record['dy'] * pix_scale_deg
                # Copy across the HSM moments estimates for the star stored in meas_data too
                data_record['g1'] = meas_record['g1']
                data_record['g2'] = meas_record['g2']
                data_record['sigma'] = meas_record['sigma']
                # within a field, we have the same noise level, so SNR ~ flux
                data_record['weight'] = meas_record['flux']

    for i_tile in range(n_tiles):
        tile_index = tile_index_list[i_tile]
        data_cat = data_cat_list[i_tile]

        data_cat['weight'] /= data_cat['weight'].sum()
        filename = os.path.join(work_dir, "star_data-%03d-%02d-%02d.fits" % ((field,) + tile_index))
        pyfits.writeto(filename, data_cat, clobber=True)


class FieldModelSuite(object):
    """Represents a suite of PSF models (e.g. LinearChebyshevModel instances) that represent
    the PSF of a GREAT3 field, split up into tiles.

    It provides dict-like access to the models themselves.  For instance:

      s = FieldModelSuite(...)
      s[1,2].inspectStar(...)    # inspect a star's fit in tile (1,2)

    It also provides an evaluate() method that looks up the appropriate model and calls its
    evaluate() method:

      s.evaluate(x_tile_index=1, y_tile_index=2, tile_x_pos_deg=0.321, tile_y_pos_=0.456)

    Most often, FieldModelSuite objects will be constructed via the fit() static method
    (or the main() function, which delegates to that), and then either makeImages() is
    called or the suite object is saved for later use.
    """

    def __init__(self, models, dim):
        self._models = models
        self.dim = dim

    def __getitem__(self, k):
        x_tile, y_tile = k
        return self._models[int(x_tile), int(y_tile)]

    def __setitem__(self, k, v):
        x_tile, y_tile = k
        self._models[int(x_tile), int(y_tile)] = v

    def save(self, filename):
        """Save the model suite to disk

        This is a workaround for the fact that you can't pickle an object from a script run via
        __main__, as __main__ doesn't know what module it's in.  Instead we pickle the __init__
        arguments of all of the constituent models, and call them directly when loading.
        """
        m = {}
        for k, v in self._models.iteritems():
            m[k] = v.__getinitargs__()
        with open(filename, 'w') as f:
            pickle.dump((m, self.dim), f, protocol=2)

    @staticmethod
    def load(filename):
        """Load a model suite from disk

        See save() for an explanation of why we can't just rely on pickle
        """
        with open(filename, 'r') as f:
            m, dim = pickle.load(f)
        models = {}
        for k, v in m.iteritems():
            models[k] = LinearChebyshevModel(*v)
        return FieldModelSuite(models, dim)

    def evaluate(self, x_tile_index, y_tile_index, tile_x_pos_deg, tile_y_pos_deg, **kwds):
        """Return the PSF model at the given tile index and position

        Arguments match the GREAT3 galaxy catalog columns with the same names; **kwds are
        passed unmodified to the PSF model class' evaluate method (e.g.
        LinearChebyshevModel.evaluate).
        """
        return self[x_tile_index, y_tile_index].evaluate(tile_x_pos_deg, tile_y_pos_deg, **kwds)

    def makeImageGrid(self, galaxy_cat):
        """Create a grid of PSF images that corresponds to the PSF model evaluated at the
        positions of galaxies.

        The 100x100 grid of PSF images will directly correspond to the grid of galaxy images
        the given galaxy catalog refers to, but the size of the full PSF image will be different,
        as GREAT3 galaxy images are typically 48x48 (in most branches) but PSF models must have
        odd dimensions.

        Arguments:

          galaxy_cat ---- GREAT3 galaxy catalog to read positions from.  Must be ordered
                          by row, then column (as GREAT3 catalogs are, by default).

        """
        full_image = numpy.zeros((100*self.dim, 100*self.dim), dtype=numpy.float32)
        n = 0
        for i in range(100):

            for j in range(100):

                record = galaxy_cat[n]
                image = self.evaluate(record['x_tile_index'], record['y_tile_index'],
                                      record['tile_x_pos_deg'], record['tile_y_pos_deg'])
                full_image[i*self.dim:(i+1)*self.dim,j*self.dim:(j+1)*self.dim] = image
                n += 1

        return full_image

    @staticmethod
    def fit(field, work_dir, dim, bounds, degree=4, basis_size=4, pca_kwds=None, spatial_kwds=None):
        """
        Create a new FieldModelSuite using data tables created by buildDataTables.

        Arguments:

          field -------- first subfield in the field to be processed, as an integer (should be a
                         multiple of 20)

          work_dir ----- work directory, containing the star_data catalogs used as inputs

          dim ---------- size of the PSF model images (on a side).  Must be odd.

          degree ------- degree of the Chebyshev polynomial used for interpolation between stars

          basis_size --- number of basis functions to keep

          bounds ------- tuple of (xmin, xmax, ymin, ymax) that sets the bounding box for each
                         Model object.

          pca_kwds ----- keyword arguments passed unmodified to computePCA

          spatial_kwds - keyword arguments passed unmodified to LinearChebyshevModel.fit

        """
        log("Fitting models")

        if pca_kwds is None: pca_kwds = {}
        if spatial_kwds is None: spatial_kwds = {}

        models = {}
        regex = re.compile("star_data-%03d-(\d\d)-(\d\d).fits" % field)
        for data_file in sorted(os.listdir(work_dir)):
            m = regex.match(data_file)
            if not m: continue
            tile_index = (int(m.group(1)), int(m.group(2)))
            data = pyfits.getdata(os.path.join(work_dir, data_file))
            # The shape doesn't always survive the roundtrip through pyfits, so reshape it
            # and use that to check against the provided dimension.  (If the read in shape
            # is incommensurate with dim,dim, then the reshape will fail.)
            image = data['image'].reshape(data['image'].shape[0], dim, dim)
            assert dim == image.shape[1]
            assert dim == image.shape[2]
            log("  Fitting model for tile %d,%d" % tile_index)
            basis, image0 = computePCA(data, dim, basis_size, **pca_kwds)
            models[tile_index] = LinearChebyshevModel.fit(
                basis, data, dim, degree, bounds, image0=image0, **spatial_kwds)
        return FieldModelSuite(models, dim)

def main(field, sim_dir, work_dir, obs_type, dim=47, max_stars=1000, degree=4, basis_size=None,
         bounds=None, use_old_meas=False, use_old_data=False, model_file=None,
         use_saved_models=False, make_psf_images=True, interp_kwds=None, pca_kwds=None,
         spatial_kwds=None):
    """Main driver for all routines in this file, and the implementation of most of
    the command-line interface.

    This routine has four distinct steps, each with a (possibly optional) intermediate output:

      1) Measure the star images and reorganize into tiles, using measureStars(), writing star_meas
         catalogs to the work directory.  Skipped if any of the use_* arguments are True.

      2) Resample and rescale the star images using buildDataTables(), writing star_data catalogs to
         the work directory.  Skipped if use_old_data or use_saved_models is True.

      3) Compute a PCA basis and fit Chebyshev spatial polynomials for each tile, building a
         FieldModelSuite object using FieldModelSuite.fit().  Saves this object if (and only if)
         model_file is set to the filename to write to.  Skipped (and the FieldModelSuite loaded
         from model_file) if use_saved_models is True.

      4) Create PSF model images that correspond to the positions of all galaxies in the field,
         using FieldModelSuite.makeImageGrid(), writing psf_models FITS images to the work
         directory.  Skipped if make_psf_images is False.

    Arguments:

      field -------- first subfield in the field to be processed, as an integer (should be a
                     multiple of 20)

      sim_dir ------ simulation directory, containing the GREAT3 images and catalogs for a single
                     branch

      work_dir ----- work directory; contains the star_meas catalogs used as inputs and will
                     contain the output star_data catalogs on return

      obs_type ----- one of 'ground' or 'space'; used to set the bounding boxes for the models
                     (space branches have 0.5 degree tiles, ground has 2.0 degree tiles) and
                     the default value for basis_size.

      dim ---------- size of the PSF model images (on a side).  Must be odd.

      max_stars ---- maximum number of stars to use per tile (None = no limit); the highest SNR
                     stars will be kept without regard for their spatial distribution.

      degree ------- degree of the Chebyshev polynomial used for interpolation between stars

      basis_size --- number of basis functions to keep

      use_old_meas - whether to skip step 1 above, and use existing star_meas catalogs.

      use_old_data - whether to skip steps 1 and 2 above, and use existing star_data catalogs.

      model_file --- name of the file to save the FieldModelSuite to, or (if use_saved_models)
                     to load it from

      use_saved_models -- whether skip steps 1-3 above, and instead load a FieldModelSuite by
                          unpickling model_file.

      make_psf_images --- whether to perform step (4) above

      interp_kwds -- keyword arguments passed unmodified to galsim.InterpolatedImage when shifting
                     the star images

      pca_kwds ----- keyword arguments passed unmodified to computePCA

      spatial_kwds - keyword arguments passed unmodified to LinearChebyshevModel.fit

    """
    bounds = obs_type_bounds[obs_type]
    if basis_size is None:
        basis_size = obs_type_basis_size[obs_type]
    if use_saved_models:
        suite = FieldModelSuite.load(os.path.join(work_dir, model_file))
    else:
        if not use_old_data:
            if not use_old_meas:
                measureStars(field, sim_dir, work_dir)
            buildDataTables(field, sim_dir, work_dir, dim=dim, max_stars=max_stars,
                            interp_kwds=interp_kwds)
        suite = FieldModelSuite.fit(
            field, work_dir, dim, 
            bounds=bounds, degree=degree, basis_size=basis_size, pca_kwds=pca_kwds,
            spatial_kwds=spatial_kwds)
        if model_file is not None:
            suite.save(os.path.join(work_dir, model_file))
    if make_psf_images:
        log("Creating PSF model images")
        for subfield in range(field, field + nsubfields):
            log(" creating model images for subfield %d" % subfield)
            galaxy_cat_file = os.path.join(sim_dir, "galaxy_catalog-%03d.fits" % subfield)
            galaxy_cat = pyfits.getdata(galaxy_cat_file)
            out_file = os.path.join(work_dir, "psf_models-%03d.fits" % subfield)
            image_grid = suite.makeImageGrid(galaxy_cat)
            pyfits.writeto(out_file, image_grid, clobber=True)

if __name__ == "__main__":
    usage = "usage: %prog [options] FIELD SIM_DIR WORK_DIR OBS_TYPE"
    description = """Create PSF models for the given field, and use
them to create images containing grids of PSF images for the galaxies
in that field.  FIELD is the number of the first subfield in the field
to process (a multiple of 20); all subfields in that field will be
processed together.  SIM_DIR is the directory containing GREAT3 images
and catalogs for the branch of interest.  WORK_DIR is the directory
where output files should be placed.  It will be created if it does
not exist.  There will be one file for each subfield, with PSF images
on a grid that corresponds to the same grid for galaxies (though the
grid spacing will be different, as defined by the --dim argument).
OBS_TYPE is one of "ground" or "space", and is used to set the
expected tile size and the default value for --basis-size.

By default, all steps of reorganizing the data and building the PSF
models will be carried out, with intermediate outputs written to disk
as well as the final PSF model images that correspond to the galaxy
images.  The --use* options can be used to restart the job from these
intermediate outputs instead of having to repeat them.
"""
    parser = optparse.OptionParser(usage=usage, description=description)
    parser.add_option("--dim", dest="dim", type=int, default=41, metavar="N",
                      help="Width and height of PSF model images")
    parser.add_option("--basis-size", dest="basis_size", type=int, default=None, metavar="N",
                      help="Number of PCA basis images to use (default is 4 for ground and 8 "+
                      "for space)")
    parser.add_option("--degree", dest="degree", type=int, default=4, metavar="N",
                      help="Maximum order of the 2-d Chebyshev polynomial spatial functions")
    parser.add_option("--max-stars", dest="max_stars", type=int, default=1000, metavar="N",
                      help="Maximum number of stars to include per-tile (clip by SNR)")
    parser.add_option("--no-max-stars", dest="max_stars", action='store_const', const=None,
                      help="Make the number of stars to use per tile unlimited")
    parser.add_option("--model-file", dest="model_file", type=str, default="psf_model.p", 
                      metavar="FILE", help="Save (or load) the model suite object with this "+
                      "filename in the WORK_DIR directory [default='psf_model.p']")
    parser.add_option("--no-psf-images", dest="make_psf_images", action="store_false", default=True,
                      help="Don't create PSF model images (should be used with --model-file)")
    parser.add_option("--use-old-meas", dest="use_old_meas", action="store_true", default=False,
                      help="Reuse star_meas files from a previous run")
    parser.add_option("--use-old-data", dest="use_old_data", action="store_true", default=False,
                      help="Reuse star_data files from a previous run")
    parser.add_option("--use-saved-models", dest="use_saved_models", action="store_true",
                      default=False, metavar="FILE",
                      help="Load models from disk; --model-file must be set")
    parser.add_option("--quiet", dest="quiet", action='store_true', default=False,
                      help="Don't print progress statements")
    opts, args = parser.parse_args()
    try:
        field, sim_dir, work_dir, obs_type = args
    except ValueError:
        parser.error("exactly four positional arguments are required")
    if not os.path.isdir(sim_dir):
        parser.error("input directory %s does not exist or is not a directory" % sim_dir)
    try:
        field = int(field)
    except TypeError:
        parser.error("field argument '%s' is not an integer" % field)
    obs_type = obs_type.strip().lower()
    if obs_type not in ("space", "ground"):
        parser.error("obs_type '%s' must be one of 'space' or 'ground'" % obs_type)
    if opts.use_saved_models and not opts.model_file:
        parser.error("--use-pickled-models requires --model-file")
    if not opts.make_psf_images and not opts.model_file:
        sys.stderr.write("WARNING: not making PSF images or saving the PSF model\n")
    if not os.path.isdir(work_dir):
        os.makedirs(work_dir)
    if opts.quiet:
        verbose = False
    else:
        verbose = True
    main(
        field, sim_dir, work_dir, obs_type,
        dim=opts.dim, basis_size=opts.basis_size,
        degree=opts.degree, use_old_meas=opts.use_old_meas, use_old_data=opts.use_old_data,
        use_saved_models=opts.use_saved_models, max_stars=opts.max_stars,
        model_file=opts.model_file, make_psf_images=opts.make_psf_images)
