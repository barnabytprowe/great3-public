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

#truth_dir = "/Users/browe/great3/truth-alpha-release-2"
truth_dir = "/great3/beta/truth"
#public_dir = "/Users/browe/great3/great3-private/tests/test_run"
#public_dir = "/Users/browe/great3/public"
public_dir = "/great3/tmp/unpacked/beta/public"
#public_dir = /Volumes/Data/GREAT3/beta/public

# Experiments, obs_types and shear_types to validate
experiments = [
    "control",
    "real_galaxy",
    "variable_psf",
    "multiepoch",
    "full"
]

obs_types = [
    "ground",
    "space"
]

shear_types = [
    "constant",
    "variable"
]

# Information about pipeline location, and tags and column layouts
pipeline_dir = "/Users/browe/great3/validation"  # Directory in which the pipeline output is stored
pipeline_columns = { 
    "imcat-1": { # First gen. alpha data
        "id": 0, "g1": 3, "g2": 4},
    "hscpipe_shape_hsm_regauss-1": { # First gen. alpha data
        "id": 0, "g1": 3, "g2": 4},
    "im3shape-1": { # First gen. alpha data
        "id": 0, "g1": 3, "g2": 4},
    "im3shape-great3-beta": { # First gen. beta data
        "id": 0, "g1": 3, "g2": 4}, 
    "hscpipe_shape_hsm_regauss-2": { # Second gen. alpha data with weight
        "id": 0, "g1": 3, "g2": 4, "w": 15}, 
}

# Directory into which to put submissions
submission_dir = "./submissions"
