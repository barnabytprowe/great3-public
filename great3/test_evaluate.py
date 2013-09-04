import sys
import os
import logging
import evaluate
path, module = os.path.split(__file__)
sys.path.append(os.path.join(path, "..", "..", "metrics")) # Appends the great3-private/metrics
                                                           # folder to path
import g3metrics


logger = logging.getLogger("test")
logger.setLevel(logging.DEBUG)

sind, g1t, g2t = evaluate.get_generate_const_truth('control', 'ground', logger=logger)
gdict = evaluate.get_generate_const_subfield_dict('control', 'ground', logger=logger)
grot = evaluate.get_generate_const_rotations('control', 'ground', logger=logger)

label = "sub1"
g1sub, g2sub = g3metrics.make_submission_const_shear(
    0., 0., 0., 0., g1t, g2t, 1e4, 0.05, label=label)
subfile = "./g3subs/g3_const_shear_sub."+label+".dat"

qc = evaluate.Q_const(subfile, 'control', 'ground')

