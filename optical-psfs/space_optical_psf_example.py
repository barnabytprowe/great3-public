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
"""This file contains a very simple example of usage for the classes in space_optical_psf.py"""
import space_optical_psf
import galsim
import numpy as np

rms = 1./13. # wavelength
tel_diam = 2.4 # meters
f_number = 8.1
dx = 1.e-6/tel_diam/f_number*206265 # arcsec

# read wfirst fits file
optical_psf_model = space_optical_psf.OpticalPSFModel(rms = rms)
xgrid = [optical_psf_model.xmin, 0.5*optical_psf_model.xmin,
         0., 0.5*optical_psf_model.xmax, optical_psf_model.xmax]
ygrid = [optical_psf_model.ymin, 0.5*optical_psf_model.ymin,
         0., 0.5*optical_psf_model.ymax, optical_psf_model.ymax]


# make image
stamp_size = 128
n_grid = len(xgrid)
image = galsim.ImageF(stamp_size * n_grid, stamp_size * n_grid)
additional_errors = np.array([1e-7, 0., 0., 0., 0., 0., 0., 0.])
for iy, y in enumerate(ygrid):
    for ix, x in enumerate(xgrid):
        print "generating PSF at (%f, %f)" % (x, y)
        optics = optical_psf_model.get_psf(x, y)

        b = galsim.BoundsI(ix*stamp_size+1 , (ix+1)*stamp_size-1, 
                           iy*stamp_size+1 , (iy+1)*stamp_size-1)
        sub_image = image[b]
        optics.draw(image = sub_image, dx = dx)
filename = "space_optical_psf_example.fits"
image.write(filename)
print "output is written in:", filename
