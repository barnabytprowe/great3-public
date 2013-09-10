import galsim.noise
import numpy

from . import constants

def makeBuilder(obs_type, multiepoch, variable_psf):
    """Return a NoiseBuilder appropriate for the given options.

    @param[in] obs_type      Observation type: either "ground" or "space"
    @param[in] multiepoch    If True, this is a multiepoch simulation, and this
                             is just the noise of a single exposure.  If False,
                             the noise should be that of a coadd with
                             constants.n_epochs epochs
    @param[in] variable_psf  If True, this is a variable PSF branch, which affects the seeing
                             within a given image.
    """
    return PlaceholderNoiseBuilder(
        obs_type=obs_type, multiepoch=multiepoch, variable_psf=variable_psf)

class NoiseBuilder(object):

    def generateEpochParameters(self, rng, subfield_index, epoch_index, seeing, noise_mult):
        """Return a dict of metaparameters for the given subfield and
        epoch.  These will be passed to generateCatalog when it
        is called.

        @param[in] rng             A galsim.UniformDeviate to be used for any random numbers.
        @param[in] subfield_index  Index of the simulated patch of sky we are "observing" with this
                                   PSF.
        @param[in] epoch_index     Index of this "epoch" among those of the same subfield.
        @param[in] seeing          PSF FWHM in arcsec for ground-based data.  This quantity is
                                   ignored for space-based data.
        @param[in] noise_mult      A factor by which to multiply the noise variance when actually
                                   generating the noise.  The other quantities stored internally in
                                   the noise builder are not modified.  Do not use this to change
                                   the noise level in exposures for multiepoch sims; noise builder
                                   will take care of that internally based on the value of
                                   multiepoch.
        """
        raise NotImplementedError("NoiseBuilder is abstract.")

    def addNoise(self, rng, parameters, image, noise):
        """Add noise to a postage stamp image.

        @param[in] rng         A galsim.UniformDeviate to be used for any random numbers.
        @param[in] parameters  A dict of metaparameters, as returned by the
                               generateEpochParameters() method.
        @param[in,out] image   Postage stamp image to add noise to.
        @param[in] noise       A galsim.CorrelatedNoise instance that specifies the noise already in
                               the image, and can be used to whiten it or otherwise modify it before
                               adding the desired noise. May be None to indicate there is no noise
                               already in the image.
        """
        raise NotImplementedError("NoiseBuilder is abstract.")

    def addStarImageNoise(self, rng, parameters, snr, image):
        """Add noise to a real-PSF star postage stamp.
        """
        raise NotImplementedError("NoiseBuilder is abstract.")

class PlaceholderNoiseBuilder(NoiseBuilder):
    """
    Some notes from MJ about what the real thing (non-placeholder) ought to look like:

    I think it should just have the parameters for CCDNoise: sky_level, gain, read_noise.

    For ground, the sky_level should definitely vary. We could keep the other two constant probably.

    For space, I guess sky_level is either 0 or very low. And probably everything stays constant.

    Also, we discussed having the galaxies have a minimum S/N, so we need to somehow couple the
    parameters here to the range of flux values we choose for the galaxies. We can do it by hand if
    we have to, but it might be nice to have the "maximum variance" be accessible to the galaxy
    builder. So perhaps a getMaxVariance() function here. And then we could pass the noise_builder
    to galaxy_builder.generateSubfieldParameters.
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
        if self.multiepoch:
            n_epochs = constants.n_epochs
        else:
            n_epochs = 1
        # Variance gets decreased by n_epochs, because we're assuming it's like taking the same
        # total exposure time but splitting it up.  If the exposure time is 1/n_epochs shorter, then
        # the sky level is lower by that factor, and so is the effective noise variance (in the
        # limit that sky noise dominates, which is what we're assuming).
        self.noise_mult = noise_mult / n_epochs
        if self.obs_type == "space":
            self.min_variance = 1.35e-3
            self.max_variance = 1.40e-3
        else:
            # Interpolate between tabulated values.  Note that (I believe) changing the value of
            # noise variance based on the per-epoch seeing is not what we want; we want to preserve
            # the fact that S/N differs for good- and bad-seeing images.  So for single epoch, we
            # use the seeing in that epoch.  For multi-epoch or variable PSF (where there are
            # multiple values per tile), we use some effective seeing, the harmonic mean.
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
        self.typical_variance = 0.5*(self.min_variance + self.max_variance)

        variance = rng() * (self.max_variance - self.min_variance) + self.min_variance
        return dict(variance=variance*self.noise_mult)

    def addNoise(self, rng, parameters, image, noise):
        # Depending on the path we've taken here, this could be an array of length 1, or a scalar.
        # Make sure it's just a scalar so galsim doesn't barf later on.
        variance = parameters['variance']
        if isinstance(variance, numpy.ndarray):
            variance = variance[0]
        if noise is not None:
            variance -= noise.applyWhiteningTo(image)
        if variance < 0.0:
            raise ValueError("After whitening, desired noise level cannot be achieved")
        new_noise = galsim.GaussianNoise(rng, sigma=variance**0.5)
        new_noise.applyTo(image)

    def addStarImageNoise(self, rng, parameters, snr, image):
        """Add noise to a star postage stamp image.
        """
        if self.variable_psf:
            # Add noise with the same variance as in the galaxy image for this epoch, and adjust
            # star fluxes to achieve the desired SNR.
            noise = galsim.GaussianNoise(rng,
                                         sigma = numpy.sqrt(parameters['variance']))
            image.addNoiseSNR(noise, snr, preserve_flux=False)
        else:
            pass
