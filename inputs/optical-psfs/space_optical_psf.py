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
"""This file contains the classes needed to represent space-based optical PSFs for an early design
for WFIRST-AFTA, namely the OpticalPSFModel class, which reads in tabulated information about an
optical PSF model.

This class is imported and used by the GREAT3 simulation software."""

import pyfits
import galsim
import numpy as np
import itertools
import exceptions

class OpticalPSFModel:
    """This class is used for obtaining a space-based optical PSF at an arbitrary position (x,y).

    One should prepare a fits file that contains field positions in second HDU and Zernike
    coefficients in third HDU. For details look at an example file
    "afta_wfirst_example_psf_exaggerated.fields_and_coefs.fits"
    This fits file comes from an original fits file we got by a courtesy of
    Alden Jurling and Dave Content (NASA), which contains WFIRST optics information.
    Since the original files is about 250 MB, I cut out field and Zernike-coefficient information.

    This class interpolates based the Zernike coefficients across the FOV.
    One can add field-by-field rms error of Zernike coefficients by specifying total rms in meters.
    Also one can add additional errors for each Zernike coefficients by numpy array when generating
    PSF by get_psf(). One should specify the error in units of meters.
    One can specify wavelength in nm, diameter in m, obscuration in ratio between mirror diameter
    and obscuration diameter, number of support struts, and thickness of struts as a fraction of
    pupil diameter.
    One can get minimum/maximum x/y of the data used for the interpolation. They are stored as
    variables xmin, xmax, ymin, and ymax.

    A recommended setup based on the WFIRST optics design with a realistic value of the total rms
    is adopted as a default of this class. Since this is based on an actual telescope, the results
    could get unrealistic if the diameter is changed by much. Since the optical PSFs given here is
    "ideal" case and does not contain an error budget they allow, Alden and Dave recommend us to 
    add the total rms corresponds to diffraction limit which they define as
    rms = lambda/13.3
    .
    """
    def __init__(self, filename = "afta_wfirst_example_psf_exaggerated.fields_and_coefs.fits",
                 lam = 1000., diameter = 2.4, obscuration = 0.28,
                 nstruts = 6, strut_thick = 0.01, strut_angle = 0.*galsim.degrees,
                 pad_factor = None,
                 rms = 0.075, interpolant2d = None, seed = None):
        """
        Inputs
        - filename: filename of fits file with information of optics.
        - lam: wavelength [nm]
        - diameter: diameter of telescope [m]
        - obscuration: central obscuration [ratio between mirror diameter and obscuration diameter]
        - nstruts: number of radial support struts to add to the central obscuration
        - strut_thick: thickness of support struts as a fraction of pupil diameter
        - strut_angle: angle made between the vertical and the strut starting closest to it,
                       defined to be positive in the counter-clockwise direction; must be a
                       galsim.Angle instance
        - pad_factor: optional padding specification if 1.5 is not good enough
        - rms: total rms of the random Zernike coefficients [wavelength]
        - interpolant2d: galsim._galsim.InterpolantXY
                         If None, galsim.InterpolantXY(galsim.Quintic())
        - seed: random seed to use for numpy routines that make random additional aberrations (if
                None, then let numpy seed routines based on time)
        """

        self.lam = lam*1e-9 # meters
        self.lam_over_diam = self.lam/diameter*206265 # arcsec
        self.obscuration = obscuration
        self.nstruts = nstruts
        self.strut_thick = strut_thick
        self.strut_angle = strut_angle
        self.pad_factor = pad_factor

        # read file
        hdulist = pyfits.open(filename)
        primary_hdu, fields, coefs = hdulist

        # note that fields.data['x'] is actually y-axis and fields.data['y'] is x-axis
        fields_1d = np.array(zip(fields.data['x'],fields.data['y']))

        n = int(np.sqrt(fields_1d.shape[0]))
        self.dy = sorted(fields.data['x'])[n] - sorted(fields.data['x'])[n-1]
        self.dx = sorted(fields.data['y'])[n] - sorted(fields.data['y'])[n-1]

        self.xmin = fields.data['y'].min()
        self.xmax = fields.data['y'].max()
        self.ymin = fields.data['x'].min()
        self.ymax = fields.data['x'].max()

        sorted_coordinates_raster = np.array(
            [row[2] for row in sorted([(r[0], r[1], i) for i, r in enumerate(fields_1d)])])
        field_mapping_index = sorted_coordinates_raster.reshape(n, n)
        # Zernike coefficients in the input file is 10 times larger than actual ones
        mapped_coefs = coefs.data[field_mapping_index]/10.

        # interpolate coefficients
        self.n_coefs = 8 # defocus, a1, a2, c1, c2, t1, t2, spher
        self.interpolated_coefficients = list()
        for i_coefs in range(self.n_coefs):
            im_coef = galsim.ImageViewD(np.ascontiguousarray(mapped_coefs[:, :, i_coefs+3]))
            im_coef.setScale(1.)
            if interpolant2d == None:
                interpolant2d = galsim.InterpolantXY(galsim.Quintic())
            self.interpolated_coefficients.append(galsim.InterpolatedImage(im_coef,
                                                                   x_interpolant = interpolant2d,
                                                                   normalization = "sb",
                                                                   calculate_stepk = False,
                                                                   calculate_maxk = False,
                                                                   ))
        # generate aberration errors
        if rms != 0.:
            if seed is not None:
                np.random.seed(seed)
            self.aberration_errors = np.random.normal(0., rms/np.sqrt(self.n_coefs), self.n_coefs)
        else:
            self.aberration_errors = np.zeros(self.n_coefs)

        hdulist.close()

    def get_zernike_coefficients(self, x, y):
        """
        Inputs
        - x: x position on FOV [deg]
        - y: y position on FOV [deg]

        Outputs
        - delta_coefs: ndarray, Zernike coefficients at (x, y) including effects by misalignments
                       [wavelength].
        """
        if x < self.xmin or x > self.xmax or y < self.ymin or y > self.ymax:
            import warnings
            warnings.warn(
                    "Warning: position (%f,%f) not within the bounds "%(x,y) +
                    "of the gridded values.  Min, max x: (%f,%f) "%(self.xmin,self.xmax) +
                    "and min, max y: (%f, %f) "%(self.ymin,self.ymax))

        coefs = list()
        for interpolated_coefficient, aberration_error in zip(self.interpolated_coefficients,
                                                      self.aberration_errors):
            coefs.append(interpolated_coefficient.xValue(galsim.PositionD(x/self.dx, y/self.dy))
                         /self.lam + aberration_error)
        return np.array(coefs)

    def get_psf(self, x, y, additional_coefs = None):
        """
        - x: x position on FOV [deg]
        - y: y position on FOV [deg]
        - additional_coefs: ndarray (Noll ordering: defocus, astig1, astig2, coma1, coma2,
                            trefoil1, trefoil2, spher), additional Zernike coefficients.
                            If None, no additional errors are added [wavelength]

        Outputs
        - optics: galsim.optics.OpticalPSF
        """
        if x < self.xmin or x > self.xmax or y < self.ymin or y > self.ymax:
            import warnings
            warnings.warn(
                    "Warning: position (%f,%f) not within the bounds "%(x,y) +
                    "of the gridded values.  Min, max x: (%f,%f) "%(self.xmin,self.xmax) +
                    "and min, max y: (%f, %f) "%(self.ymin,self.ymax))

        coefs = self.get_zernike_coefficients(x, y)
        if additional_coefs is not None:
            coefs += additional_coefs
        optics = galsim.OpticalPSF(lam_over_diam = self.lam_over_diam, 
                                   defocus = coefs[0],
                                   astig1 = coefs[1],
                                   astig2 = coefs[2],
                                   coma1 = coefs[3],
                                   coma2 = coefs[4],
                                   trefoil1 = coefs[5],
                                   trefoil2 = coefs[6],
                                   spher = coefs[7],
                                   obscuration = self.obscuration,
                                   nstruts = self.nstruts,
                                   strut_thick = self.strut_thick,
                                   strut_angle = self.strut_angle,
                                   pad_factor = self.pad_factor,
                                   suppress_warning = True)

        return optics
