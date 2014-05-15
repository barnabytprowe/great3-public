# GREAT3 test scripts
=====================

The scripts in this directory were used to test out the GREAT3 simulation
scripts, and simulation design, in several ways that users may wish to try.

The contents are as follows:

1. `test_run.py` is a script that can be used to go through all the steps of the
GREAT3 simulation process, but for a limited number of subfields and with images
that have fewer galaxies than a real GREAT3 simulation.  Because of the limited
amount of data, this can be run quickly at the command line without the need to
use a large cluster or any multiprocessing.  Users will need to make sure they
have all the dependencies for the great3sims package (see
../great3sims/README.md for details) and all the necessary data described there.

2. `flux_test.py` is a script that was used to flag galaxies that cannot be used
for simulation purposes, because of issues like blends or masking trouble.

3. `galaxy_props.py` is a script that can be used to check the effect on the
galaxy population of the cuts imposed on galaxies in GREAT3.  It allows for
successive imposition of cuts, and plots histograms of properties like
ellipticity, size, etc. after each cut is imposed.

4. `test_bmode_noise.py` is a script that investigated the use of B-mode only
galaxy intrinsic shape noise in the GREAT3 simulations, designed to improve the
sensitivity of metrics based on the E-mode aperture mass dispersion signal.

All of these are documented, and users should be able to run them if they have
the galaxy data used for GREAT3 simulations.
