# GREAT3 galaxy sample
======================

This directory contains scripts that were used to produce the galaxy catalogs
and selection criteria for GREAT3.  This process can be roughly divided into two
steps:

1. Making a COSMOS galaxy sample to I<23.5 with associated files for parametric
files.  The outputs of this process were packaged and made available on the
GREAT3 download page, http://great3.projects.phys.ucl.ac.uk/leaderboard/data

2. Calculating additional quantities that were needed to cut the sample in
various ways, for example, to ensure the galaxy images were within our target
S/N range, were sufficiently resolved in the ground-based simulated images, were
not too contaminated by blends or image defects, etc.

Below we describe the steps in and files associated with both of these
processes.


## Files related to the basic COSMOS I<23.5 sample
--------------------------------------------------

The first step in the processing is based on legacy code from the SHERA software
package.  It consists of five IDL scripts, the inputs for which are not
available here, so they are primarily made available for users to see what was
done.

1. `make_shera_cat.pro` was used to put together a catalog listing all galaxies
in the F814W<23.5 sample, in a format that can be used by SHERA.

2. `run_many.pro` was used to drive the three scripts below, in sequence, to
generate the postage stamps and a GalSim-style catalog representing the galaxy
sample.  It was used with the following calling sequence in IDL:

    run_many,'_23.5',0,56062,1000,10

This calling sequence tells the script about the size of the galaxy sample, how
many postage stamps to include in a given image file (1000), and a random seed
to use when randomly reordering the sample (10).

3. `makecatalog_many.pro` is the first script driven by `run_many.pro`.  It
takes the SHERA catalog, and puts it into the format needed by GalSim.

4. `makestamps_many.pro` is the second script driven by `run_many.pro`.  It
creates postage stamps for all the galaxies in the proper format, and calculates
their noise properties (variance).

5. `makecatalog_many_var.pro` is the third script driven by `run_many.pro`.  It
adds the noise variances from step (4) to the catalogs from step (3).



## Files related to the additional selection criteria
-----------------------------------------------------
