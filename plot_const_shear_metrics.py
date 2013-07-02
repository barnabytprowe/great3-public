import os
import cPickle
import numpy as np
import g3metrics

NIMS = 200               # Number of images per set, always 200 for G10
NGALS_PER_IM = 10000     # In GREAT08/GREAT10 there were 10000 galaxies per image
TRUE_MIN = 0.01          # Range of true input |shear| for random selection from an annulus around
TRUE_MAX = 0.05          # the origin
NFIELDS = 10             # Don't necessarily need to have NIMS input shears, but rather NFIELDS
                         # where NFIELDS * NSUBFIELDS = NIMS

CFID = 2.e-4 # Fiducial, "target" m and c values
MFID = 2.e-3 #

# Plotting ranges of interest
CMIN = CFID
CMAX = 2.e-2
MMIN = MFID
MMAX = 2.e-1

NBINS = 5     # Number of bins to plot in the ranges above
NMONTE = 100  # Number of montecarlo samples
NOISE_SIGMA = 0.05  # Noise due to pixel shot noist on a shear estimate, per galaxy

# Generate arrays of values
cvals = CMIN * (CMAX / CMIN)**(np.arange(NBINS) / float(NBINS - 1.)) # geometric series
mvals = MMIN * (MMAX / MMIN)**(np.arange(NBINS) / float(NBINS - 1.))
cgrid, mgrid = np.meshgrid(cvals, mvals) # 2D arrays covering full space if needed

# Create empty storage arrays
Q08 = np.empty((NBINS, NBINS, NMONTE))
QZ1 = np.empty((NBINS, NBINS, NMONTE))
QZ2 = np.empty((NBINS, NBINS, NMONTE))

# File for storing pickled output
OUTFILE = os.path.join('results', 'const_shear_metrics.pkl')

if not os.path.isfile(OUTFILE):
    for krepeat in range(NMONTE):

        # Generate the truth tables for this realization
        g1true, g2true = g3metrics.make_const_truth_uniform_annular(
            NFIELDS, NIMS, range_min=TRUE_MIN, range_max=TRUE_MAX)
    
        # Loop over each c, m combination
        for i in range(NBINS):
        
            for j in range(NBINS):

                # Make the submissions
                g1sub, g2sub = g3metrics.make_submission_const_shear(
                    cvals[i], cvals[i], mvals[j], mvals[j], g1true, g2true,
                    ngals_per_im=NGALS_PER_IM, noise_sigma=NOISE_SIGMA)
                # Calculate the metrics and store
                Q08[i, j, krepeat] = g3metrics.metricQ08_const_shear(
                    g1sub, g2sub, g1true, g2true, nfields=NFIELDS)
                QZ1[i, j, krepeat] = g3metrics.metricQZ1_const_shear(
                    g1sub, g2sub, g1true, g2true, cfid=CFID, mfid=MFID)[0]
                QZ2[i, j, krepeat] = g3metrics.metricQZ2_const_shear(
                    g1sub, g2sub, g1true, g2true, cfid=CFID, mfid=MFID)[0]

        print "Calculated const shear metrics for "+str(krepeat + 1)+"/"+str(NMONTE)+" realizations"

    # Save the results as a tuple of NumPy arrays
    cPickle.dump((Q08, QZ1, QZ2), open(OUTFILE, 'wb'))

else:
    Q08, QZ1, QZ2 = cPickle.load(open(OUTFILE, 'rb'))

