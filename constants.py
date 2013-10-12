truth_dir = "/Users/browe/great3/truth-alpha-release-2"
#public_dir = "/Users/browe/great3/great3-private/tests/test_run"
public_dir = "/Users/browe/great3/public"

# Experiments, obs_types and shear_types to validate
experiments = [
    "control",
    "real_gal",
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
    "hscpipe_shape_hsm_regauss-2": { # Second gen. alpha data with weight
        "id": 0, "g1": 3, "g2": 4, "w": 5}, 
}

# Directory into which to put submissions
submission_dir = "./submissions"
