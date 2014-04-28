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
Global constants, some of which are dictionaries that depend on the experiment, observation, or
shear type.
"""
import numpy

image_size_deg = 10. # linear extent of one of our images.

nrows = 100   # Number of rows, columns of galaxies per subfield
ncols = 100

# Define the size of each galaxy postage stamp in pixels.  We want to enable different values
# depending on the observation type ("space" vs. "ground") and/or on single- vs. multi-epoch, since
# we've found that for the space-based data, the tiny pixels used for single-epoch sims mean that
# some galaxies are rather extended compared to the postage stamp size.  Make this a two-level dict,
# indexed first by observation type and then by single vs. multi-epoch.
xsize = {
    "space": { # index second level of dictionary based on bool multiepoch
        True: 48,
        False : 96,
        },
    # For ground, we use the same pixel scale for both, but it's neater to
    # maintain the same 2-level structure and it leaves open the
    # possibility of having different pixel scales in the two cases if
    # we decide it is necessary later on.
    "ground": {
        True: 48,
        False: 48,
        },
    }
# Let's just require the PS to always be square.
ysize = xsize

xcells = 15   # When placing stars for variable_psf images, we split the image into cells and assign
ycells = 15   # one star per cell.  This avoids blending and rare pathological clustering.

# When making star fields for constant PSF images, we make a grid of PSFs.  One of them is perfectly
# centered, the others have random offsets that are unknown to the participants, in the range
# +/-centroid_shift_max.
nx_constpsf = 3
ny_constpsf = 3

n_subfields = 220   # Number of galaxy/shear catalogs per experiment.  This must be an integer
# multiple of n_subfields_per_field, and should include the number of deep subfields (for which more
# parameters are defined below).

centroid_shift_max = 1.0  # Shift galaxy centroids by uniform random in [-, +] in PIXELS.

epoch_shift_max = 1.0    # Shift epochs relative to the subfield by uniform random in [-, +]
                         # in PIXELS.

n_epochs = 6  # Number of epochs for each subfield when multiepoch

pixel_scale = {
    "space": { # index second level of dictionary based on bool multiepoch
        True: 0.10, # not Nyquist sampled when there are multiple epochs
        False : 0.05, # Nyquist sampled if just a single epoch
        },
    # For ground, we use the same pixel scale for both, but it's neater to
    # maintain the same 2-level structure and it leaves open the
    # possibility of having different pixel scales in the two cases if
    # we decide it is necessary later on.
    "ground": {
        True: 0.2,
        False: 0.2,
        },
    }

# Number of subfields per field: this is a two-level dict that is indexed first based on the shear
# type, and then based on the boolean variable_psf, which is True for the variable_psf and full experiments
# and False for the others.  The idea is that branches with variable *anything* should have 20
# subfields per field, and branches with constant shear AND PSF should have 1 subfield per field.

n_subfields_per_field = {
    # first layer of indexing based on shear type
    "variable": { # second layer of indexing based on bool variable_psf
        True: 20,
        False: 20,
        },
    "constant": {
        True: 20,
        False: 1,
        },
}

# Subsampling parameter for subfields within a field.  i.e., we allow this many offset positions
# (per dimension) for a subfield with respect to the overall field position.
subfield_grid_subsampling = 7

# Parameters related to the "deep" subfields, which are always generated at the end of the branch,
# after all the regular subfields.
n_deep_subfields = 20 # Number of deep subfields to generate per branch (must be multiple of
                      # n_subfields_per_field).
deep_frac = 0.025     # Relative size of the deep dataset compared to the regular one, when it comes
                      # to what is actually distributed with the challenge.  It must be the case
                      # that deep_frac <= n_deep_subfields / (n_subfields - n_deep_subfields).
deep_variance_mult = 0.0625 # Factor by which the noise variance should be reduced compared to that
                            # in the regular fields when creating the deep fields.  This comes from
                            # the idea that the depth is supposed to differ by ~1 mag, which would
                            # mean increasing fluxes by 4 or keeping fluxes fixed while lowering the
                            # noise by 1/4 (or noise variance by 1/4^2).

# For variable PSF case, set a stellar density and number of stars per image
min_star_density = 1. # per arcmin^2
max_star_density = 3.
