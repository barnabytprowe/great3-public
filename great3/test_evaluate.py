import sys
import os
import logging
import numpy as np
import evaluate
path, module = os.path.split(__file__)
sys.path.append(os.path.join(path, "..", "..", "metrics")) # Appends the great3-private/metrics
                                                           # folder to path
import g3metrics

experiment = 'control'
obs_type = 'ground'

logger = logging.getLogger("test")
logger.setLevel(logging.DEBUG)

# Just try getting / building the intermediate products for this branch first
sind, g1t, g2t = evaluate.get_generate_const_truth(experiment, obs_type, logger=logger)
gdict = evaluate.get_generate_const_subfield_dict(experiment, obs_type, logger=logger)
grot = evaluate.get_generate_const_rotations(experiment, obs_type, logger=logger)

# Try a simple submission, no biases, and see what Q I get
label = "sub1"
g1sub, g2sub = g3metrics.make_submission_const_shear(
   0.001, 0., 0., -0.01, g1t, g2t, 1e4, 0.05, label=label, rotate_cs=grot)
subfile = "./g3subs/g3_const_shear_sub."+label+".dat"

q, c1, m1, c2, m2, sigc1, sigm1, sigc2, sigm2  = evaluate.Q_const(subfile, 'control', 'ground')

print "Q_c = "+str(q)
print "c+ = "+str(c1)+" ("+str(sigc1)+")"
print "m+ = "+str(m1)+" ("+str(sigm1)+")"
print "cx = "+str(c2)+" ("+str(sigc2)+")"
print "mx = "+str(m2)+" ("+str(sigm2)+")"

subfield_index, offset_deg_x, offset_deg_y = evaluate.get_generate_variable_offsets(
    experiment, obs_type, logger=logger)

field, theta, map_E, map_B, maperr = evaluate.get_generate_variable_truth(
    experiment, obs_type, logger=logger)



