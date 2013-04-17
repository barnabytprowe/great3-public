from django.core.management.base import BaseCommand, CommandError
from leaderboard.models import Entry, PLACEHOLDER_SCORE, recompute_scoring
from great3.settings import installation_base
import numpy as np
import os
import datetime
import pytz

SCORE_LOG_FILE = os.path.join(installation_base, "results", "results.log")

class Command(BaseCommand):
	args = ''
	help = 'Checks for entries that do not have a score yet and processes them'

	def handle(self, *args, **options):
		outfile = open(SCORE_LOG_FILE,'a')
		entries = Entry.objects.filter(score=PLACEHOLDER_SCORE)
		for entry in entries:
			print "Processing entry %s" % entry.name
			filename = entry.get_filename()
			print "should compute score here, from loading entry filename ", filename
			print "Instead using random number from 1-1000"
			entry.score = np.random.uniform(1, 1000)
			datestamp = datetime.datetime.utcnow().replace(tzinfo=pytz.utc).isoformat()
			outfile.write('%s\t%s\t%r\t%s\n' % (entry.name,filename,entry.score,datestamp))
			entry.save()
		outfile.close()
		recompute_scoring()
