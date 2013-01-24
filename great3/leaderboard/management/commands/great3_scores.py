from django.core.management.base import BaseCommand, CommandError
from leaderboard.models import Entry, PLACEHOLDER_SCORE, recompute_scoring
import numpy as np

class Command(BaseCommand):
	args = ''
	help = 'Checks for entries that do not have a score yet and processes them'

	def handle(self, *args, **options):
		entries = Entry.objects.filter(score=PLACEHOLDER_SCORE)
		for entry in entries:
			print "Processing entry %s" % entry.name
			filename = entry.get_filename()
			print "should compute score here, from loading entry filename ", filename
			print "Instead using random number from 1-1000"
			entry.score = np.random.uniform(1, 1000)
			entry.save()
		recompute_scoring()
