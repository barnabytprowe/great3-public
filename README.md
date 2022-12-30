great3-public
=============

This is the public repository for the Third Gravitational Lensing Accuracy Testing
Challenge (GREAT3), a community challenge by the weak lensing community to test methods of inferring weak lensing shear from images.

## Key references ##

* GREAT3 challenge [handbook](https://ui.adsabs.harvard.edu/abs/2014ApJS..212....5M/abstract) describing the purpose of the challenge, the challenge dataset, the overall challenge process, and more

* GREAT3 challenge [results paper](https://ui.adsabs.harvard.edu/abs/2015MNRAS.450.2963M/abstract) describing the main results

* [GalSim](https://github.com/GalSim-developers/GalSim) community image simulation software developed initially for the challenge, and now widely used for other purpose

## How to get the data ##

The data may be downloaded from the Pittsburgh Supercomputing Center via Globus or direct download from a web server:

* [Globus link](https://app.globus.org/file-manager?origin_id=ab25ca84-b8ee-11eb-9d92-5f1f6f07872f&origin_path=%2Fhildafs%2Fdatasets%2Fpublic%2FGreat3%2F)

* [Web server](https://cosmo.vera.psc.edu/Great3/public/)

## Public code and other information in this repository ##

For a full description of the contents of great3-public,
please go into each sub-directory and take a look at the `README.md` file.

great3-public includes:

* The code used to generate the GREAT3 challenge data (in `great3sims/`)

* The scripts used to generate and explore the inputs to the GREAT3 simulation
  code (in `inputs/`)

* The metric evaluation code, and scripts used to simulate the performace of a
  variety of metrics (in `metrics/`)

* The example data processing scripts that were public throughout the challenge
  (in `example_scripts/`)

* The code provided to assist in generating and checking submissions (in
  `presubmission_script`)

Please also see the [wiki](https://github.com/barnabytprowe/great3-public/wiki)
for more information, an FAQ, and links to downloads for supporting data.
