from django.core.management.base import BaseCommand, CommandError
from leaderboard.models import Entry, PLACEHOLDER_SCORE, recompute_scoring
from great3.settings import installation_base
import evaluate
import numpy as np
import os
import datetime
import pytz
import codecs

EXPERIMENTS = (
    "control",
    "real_galaxy",
    "variable_psf",
    "multiepoch",
    "full",
    )


SCORE_LOG_FILE = os.path.join(installation_base, "results", "results.log")

CORR2_EXEC = "/usr/local/bin/corr2"
CORR2_PARAMS = os.path.join(installation_base, "server", "corr2.params")


class Command(BaseCommand):
	args = ''
	help = 'Checks for entries that do not have a score yet and processes them'

	def handle(self, *args, **options):
		outfile = codecs.open(SCORE_LOG_FILE,mode='a',encoding='utf-8')
		entries = Entry.objects.filter(score=PLACEHOLDER_SCORE)
		for entry in entries:
			print "Processing entry %s" % entry.name.encode('utf-8')
			filename = entry.get_filename()
                        if entry.board.experiment is not in experiments:
                            raise ValueError("Experiment not recognised.")
                        if entry.board.space:
                            obs_type = "space"
                        else:
                            obs_type = "ground"
                        if entry.board.varying:
                            shear_type = "variable"
                        else:
                            shear_type = "constant"
			print "should compute score here, from loading entry filename ", filename
			print "Instead using random number from 1-1000"
			entry.score = evaluate.getq(
                            filename, experiment, obs_type, shear_type, truth_dir="/great3/truth",
                            storage_dir=STORAGE_DIR, logger=logger)
			datestamp = datetime.datetime.utcnow().replace(tzinfo=pytz.utc).isoformat()
			outfile.write('%s\t%s\t%r\t%s\n' % (entry.name,filename,entry.score,datestamp))
			entry.save()
		outfile.close()
		recompute_scoring()
