GREAT3 simulation scripts
=========================

The scripts in this directory were used to generate the simulations for the
GREAT3 challenge.

In particular, the following files form the core of the 'great3sims' package:

1. `__init__.py` includes the top-level driver.

2. `builder.py` contains the SimBuilder class and the routines that do the heavy
   lifting: figuring out how to generate parameters, catalogs, config files,
   images, and packages for each branch.

3. `constants.py` contains definitions for a number of constants related to the
   simulations that are produced, such as the number of galaxies per field and
   the number of fields per branch.

4. `galaxies.py`, `shear.py`, `psf.py`, and `noise.py` contain the detailed
   implementation of the galaxy populations, shear fields, PSFs, and noise
   models, respectively.

5. `mapper.py` contains functionality related to i/o.

The scripts `mass_produce.py` and utilities in `mass_produce_utils.py` can be
used to generate large sets of simulations; they were used to drive the
production of the GREAT3 simulations on a large cluster with many cores and lots
of storage space.

## How do these scripts work?

The scripts are heavily commented with detailed docstrings.  For example, you
can import the `great3sims` package and do `help(great3sims.run)` to get a
top-level docstring.  We recommend going through these directly; for those who
wish to reproduce simulations similar to those in GREAT3, you may simply want to
use `mass_produce.py` or `../tests/test_run.py` as an example of how to carry
out all steps of the simulation generation process.  For those who wish to try
something very simple, the docstrong for the `great3sims.run()` command gives an
example of basic usage.

## What software is required?

The `great3sims` package uses
[GalSim](https://github.com/GalSim-developers/GalSim) to generate the
simulations.  In particular, it is set up to work with GalSim v1.0.1, which has
a few necessary bug fixes beyond v1.0.0.  To get GalSim v1.0.1, either download
and install the code package from the GalSim [tags
page](https://github.com/GalSim-developers/GalSim/releases), or if you have
cloned the GalSim repository, you can do `git checkout v1.0.1` to get and
install that tagged version.  Note that due to a number of API changes, the
`great3sims` package is not compatible with the latest version of GalSim on the
master branch.

Other requirements are NumPy, PyFITS (both of which are GalSim dependencies in
any case), and PyYAML.

## What data is required?

Some of the GREAT3 branches require additional data in order to generate the
simulations.  You need to obtain that data separately and pass the path to the
data to the `great3sims.run()` command via keyword arguments in order to
generate simulated data for those branches:

- All GREAT3 branches require catalogs of COSMOS galaxies, various selection
  criteria, and parametric fits to the galaxies in order to select which
  galaxies to use.  The 'real_galaxy' and 'full' experiments also require the
  galaxy images from HST.  A tarball with all required files is available from
  the [GREAT3 download page](http://great3.projects.phys.ucl.ac.uk/leaderboard/data).
- The variable shear branches require tabulated shear power spectra as
  inputs.  The ones used for GREAT3 are in `../inputs/shear-ps/tables/`, so if
  you have cloned or downloaded the contents of the great3-public repository,
  you should already have these.  You will just have to give the path to that
  directory to the `ps_dir` keyword argument when you call `great3sims.run()`.
- The variable PSF, ground-based simulations require atmospheric PSF anisotropy
  power spectra, which are also tabulated in this repository, in
  `../inputs/atmospsf/pk_math`.  You should give the path to that directory to
  the `atmos_ps_dir` keyword argument when you call `great3sims.run()`.
- The variable PSF simulations (ground and space) require optical PSF models.
  The relevant files are in this repositry, in ../inputs/optical-psfs.  You
  should give the path to that directory to the `optical_psf_dir` keyword
  argument when you call `great3sims.run()`.

## Getting help

Questions can be posted to the [GREAT3 issue
page](https://github.com/barnabytprowe/great3-public/issues?state=open) or
e-mailed to the [GREAT3 leaders](http://www.great3challenge.info/?q=contacts).

## Future work

Anyone who wishes to use these scripts in a significant way may wish to consider
the following updates.  (If you do them and wish to share them with others,
please send a pull request to the great3-public repository so we can check them
out!)

- Compatibility with the latest version of GalSim.
- More convenient way to package the data than giant tarballs.
- Make it possible to run the 'catalogs' step of processing on its own.
- Centralize the calculation of field positions, which currently lives in
  `galaxies.py`, `shear.py`, and `psf.py`.
- Split up some of the large functions.
- Make it possible to generate star test images in python `builder.py` rather
  than just via config.
- Propagate info from precomputation script into `galaxies.py` (e.g., pixel
  scales and FWHM array) rather than hard-coding everything.
