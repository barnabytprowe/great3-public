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
"""File containing the classes that generate parameters for pixel noise in the images."""
import galsim.noise
import numpy

from . import constants

def makeBuilder(obs_type, multiepoch, variable_psf):
    """Return a NoiseBuilder appropriate for the given options.

    @param[in] obs_type      Observation type: either "ground" or "space"
    @param[in] multiepoch    If True, this is a multiepoch simulation, and this is just the noise of
                             a single exposure.  If False, the noise should be that of a coadd with
                             constants.n_epochs epochs
    @param[in] variable_psf  If True, this is a variable PSF branch, which affects the seeing within
                             a given image.
    """
    return PlaceholderNoiseBuilder(
        obs_type=obs_type, multiepoch=multiepoch, variable_psf=variable_psf)

class NoiseBuilder(object):
    """A NoiseBuilder is a class that can carry out the steps necessary to define the pixel noise
    model for GREAT3."""

    def generateEpochParameters(self, rng, subfield_index, epoch_index, seeing, noise_mult):
        """Return a dict of metaparameters for the given subfield and epoch.  These will be passed
        to generateCatalog() when it is called.

        @param[in] rng             A galsim.UniformDeviate to be used for any random numbers.
        @param[in] subfield_index  Index of the simulated patch of sky we are "observing" with this
                                   PSF.
        @param[in] epoch_index     Index of this "epoch" among those of the same subfield.
        @param[in] seeing          PSF FWHM in arcsec for ground-based data.  This quantity is
                                   ignored for space-based data. In poor seeing, less noise is
                                   necessary to achieve a given galaxy S/N.
        @param[in] noise_mult      A factor by which to multiply the noise variance when actually
                                   generating the noise.  The other quantities stored internally in
                                   the noise builder are not modified.  Do not use this to change
                                   the noise level in exposures for multiepoch sims; noise builder
                                   will take care of that internally based on the value of
                                   multiepoch.
        """
        raise NotImplementedError("NoiseBuilder is abstract.")

    def addNoise(self, rng, parameters, image, current_var):
        """Add noise to a postage stamp image.

        @param[in] rng         A galsim.UniformDeviate to be used for any random numbers.
        @param[in] parameters  A dict of metaparameters, as returned by the
                               generateEpochParameters() method.
        @param[in,out] image   Postage stamp image to add noise to.
        @param[in] current_var The current variance of noise already applied to the image, if any.
        """
        raise NotImplementedError("NoiseBuilder is abstract.")

    def addStarImageNoise(self, rng, parameters, snr, image):
        """Add noise to a star postage stamp in a variable PSF branch.
        """
        raise NotImplementedError("NoiseBuilder is abstract.")

class PlaceholderNoiseBuilder(NoiseBuilder):
    """
    A PlaceholderNoiseBuilder includes a simple noise model: just Gaussian noise with the same
    variance throughout the entire image.  While this is overly simplistic, it is what was used for
    GREAT3, so we never made a more complex NoiseBuilder.
    """
    def __init__(self, obs_type, multiepoch, variable_psf):
        self.obs_type = obs_type
        self.multiepoch = multiepoch
        self.variable_psf = variable_psf
        if obs_type == "ground":
            self.fwhm_vals = [0.5, 0.65, 0.8, 0.95]
            max_var_vals = [0.0088, 0.0060, 0.0046, 0.0040]
            self.max_var_tab = galsim.LookupTable(self.fwhm_vals, max_var_vals, f_log=True)

    def generateEpochParameters(self, rng, subfield_index, epoch_index, seeing, noise_mult):
        """Generate the noise parameters for this epoch, given its seeing."""
        if self.multiepoch:
            n_epochs = constants.n_epochs
            # The next line accounts for the fact that we've defined values of noise variance for
            # single-epoch data, so if we change the pixel scale in the multiepoch data (as we do
            # for space) then the noise variance in the multiepoch sims should be made 4x as large
            # (corresponding to collection of 4x as many sky photons in a given pixel).
            noise_mult *= \
                (constants.pixel_scale[self.obs_type][True]/constants.pixel_scale[self.obs_type][False])**2
        else:
            n_epochs = 1
        # Variance gets decreased by n_epochs, because we're assuming it's like taking the same
        # total exposure time but splitting it up.  If the exposure time is 1/n_epochs shorter, then
        # the sky level is lower by that factor, and so is the effective noise variance.
        self.noise_mult = noise_mult / n_epochs
        if self.obs_type == "space":
            self.min_variance = 1.35e-3
            self.max_variance = 1.40e-3
        else:
            # Interpolate between tabulated values.  Note that changing the value of noise variance
            # based on the per-epoch seeing is not what we want; we want to preserve the fact that
            # S/N differs for good- and bad-seeing images.  So for single epoch, we use the seeing
            # in that epoch.  For multi-epoch or variable PSF (where there are multiple values per
            # tile), we use some effective seeing, the harmonic mean.
            if not self.multiepoch and not self.variable_psf:
                effective_seeing = seeing
            else:
                effective_seeing = 1. / numpy.mean(1./seeing)
            if effective_seeing < min(self.fwhm_vals):
                effective_seeing = min(self.fwhm_vals)
            if effective_seeing > max(self.fwhm_vals):
                effective_seeing = max(self.fwhm_vals)
            max_var = self.max_var_tab(effective_seeing)
            self.min_variance = 0.95*max_var
            self.max_variance = 1.05*max_var
        # Note that typical variance was determined for a regular (not deep) field, assuming a
        # single combined image at the single epoch pixel scale.  Hence it does not include (a)
        # factors of 1/n_epochs, (b) factors related to the deep field variance reduction, or (c)
        # differences in variance to account for different pixel scales in single vs. multi-epoch
        # sims.  The only considerations that determine it are what noise variance
        # we want for a given obs_type (ground/space) and seeing (if ground-based imaging).
        self.typical_variance = 0.5*(self.min_variance + self.max_variance)

        variance = rng() * (self.max_variance - self.min_variance) + self.min_variance
        # However, the dict return must include all factors that multiply the nominal noise
        # variance, to account for whether this is a regular or deep field, whether the imaging is
        # single vs. multi-epoch.  This number will be used in the image generation process.
        return dict(variance=variance*self.noise_mult)

    def addNoise(self, rng, parameters, image, current_var):
        """Actually add noise according to our noise model and the given variance to a postage stamp
        image.  This includes whitening the noise if it's a real galaxy image."""
        # Depending on the path we've taken here, this could be an array of length 1, or a scalar.
        # Make sure it's just a scalar (not a length-1 NumPy array) so that we don't run into issues
        # later on.
        variance = parameters['variance']
        if isinstance(variance, numpy.ndarray):
            variance = variance[0]
        variance -= current_var
        if variance < 0.0:
            raise ValueError("After whitening, desired noise level cannot be achieved")
        new_noise = galsim.GaussianNoise(rng, sigma=variance**0.5)
        new_noise.applyTo(image)

    def addStarImageNoise(self, rng, parameters, snr, image):
        """Add noise to a star postage stamp image in a variable PSF branch.
        """
        if self.variable_psf:
            # Add noise with the same variance as in the galaxy image for this epoch, and adjust
            # star fluxes to achieve the desired SNR.
            noise = galsim.GaussianNoise(rng,
                                         sigma = numpy.sqrt(parameters['variance']))
            image.addNoiseSNR(noise, snr, preserve_flux=False)
        else:
            pass
