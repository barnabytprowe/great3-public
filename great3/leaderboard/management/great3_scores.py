from django.core.management.base import BaseCommand, CommandError
from leaderboard.models import Entry, PLACEHOLDER_SCORE
import numpy as np

class Command(BaseCommand):
	args = ''
	help = 'Checks for entries that do not have a score yet and processes them'

	def handle(self, *args, **options):
		entries = Entry.objects.filter(score==PLACEHOLDER_SCORE)
		for entry in entries:
			filename = entry.get_filename()
			print "should compute score here, from loading entry filename"
			entry.score = np.random.randn()
			entry.save()
