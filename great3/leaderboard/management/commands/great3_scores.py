from django.core.management.base import BaseCommand, CommandError
from leaderboard.models import Entry, PLACEHOLDER_SCORE, recompute_scoring
from great3.settings import installation_base
import logging
import numpy as np
import os
import datetime
import pytz
import codecs
import sys
sys.path.append("/home/jaz/great3-private")
import evaluate

EXPERIMENTS_DICT = { # A dictionary for translating between Joe's experiment names and Rachel's
    "Control": "control",
    "Realistic Galaxy": "real_galaxy",
    "Realistic PSF": "variable_psf",
    "Multi-epoch": "multiepoch",
    "Everything": "full"}

SCORE_LOG_FILE = os.path.join(installation_base, "results", "results.log")

# At least for now for debugging purposes, would be good to have a logfile into which the
# logger supplied to the evaluate module can store a record of the pitfalls it falls into
EVALUATE_LOG_FILE = os.path.join(installation_base, "results", "evaluation.log")
LOG_LEVEL = logging.INFO # Logging level - currently info, could set to WARN once testing is done

# Directory in which the truth catalogues are stored (and assumed unpacked!)
TRUTH_DIR = os.path.join("/great3", "truth") 
# Directory into which to store intermediate products (truth catalogues, etc.) as a hardcopy cache
# that speeds up the metric calculation
STORAGE_DIR = os.path.join(installation_base, "results", "intermediates") 

# Some things required by corr2:
CORR2_EXEC = "/usr/local/bin/corr2" # Path to executable for MJ's corr2 correlation function code
                                    # on the server

# Note the Q scores for both the constant and variable shear branches are affected by
# empricially-determined normalization factor, and (in the case of the constant branches) fiducial,
# target values of m & c biases.
#
# These numbers are stored in the great3-private/server/great3/evaluate.py module, and should be
# changed or updated there if necessary.


class Command(BaseCommand):
    args = ''
    help = 'Checks for entries that do not have a score yet and processes them'

    def handle(self, *args, **options):
        outfile = codecs.open(SCORE_LOG_FILE, mode='a', encoding='utf-8')
        entries = Entry.objects.filter(score=PLACEHOLDER_SCORE)
        # Define a logging.Logger
        logging.basicConfig(filename=EVALUATE_LOG_FILE, level=LOG_LEVEL)
        logger = logging.getLogger(__name__)
        for entry in entries:

            print "Processing entry %s" % entry.name.encode('utf-8')
            filename = entry.get_filename()
            print "Computing score, loading entry filename "+str(filename)
            if entry.board.experiment not in EXPERIMENTS_DICT:
                raise ValueError("Experiment not recognised.")
            experiment = EXPERIMENTS_DICT[entry.board.experiment]
            if entry.board.space:
                obs_type = "space"
            else:
                obs_type = "ground"
            if entry.board.varying:
                shear_type = "variable"
                entry.score = evaluate.q_variable(
                    filename, experiment, obs_type, 
                    normalization=evaluate.NORMALIZATION_VARIABLE, truth_dir=TRUTH_DIR,
                    storage_dir=STORAGE_DIR, logger=logger, corr2_exec=CORR2_EXEC)
            else:
                shear_type = "constant"
                entry.score = evaluate.q_constant(
                    filename, experiment, obs_type,  cfid=evaluate.CFID, mfid=evaluate.MFID,
                    normalization=evaluate.NORMALIZATION_CONSTANT, just_q=True, truth_dir=TRUTH_DIR,
                    storage_dir=STORAGE_DIR, logger=logger)
            datestamp = datetime.datetime.utcnow().replace(tzinfo=pytz.utc).isoformat()
            outfile.write('%s\t%s\t%r\t%s\n' % (entry.name,filename,entry.score,datestamp))
            entry.save()
            outfile.close()
            recompute_scoring()
