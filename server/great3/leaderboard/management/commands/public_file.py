from django.core.management.base import BaseCommand, CommandError
from leaderboard.models import PublicDataFile
import os
import sys

class Command(BaseCommand):
	args = 'filename1 [filename2 ...]'
	help = 'Create miscellaneous PublicDataFiles for the given files'

	def handle(self, *args, **options):
		for filename in args:
			abspath=os.path.abspath(filename)
			filename=os.path.basename(abspath)
			dataFile=PublicDataFile(abspath=abspath, filename=filename, miscellaneous=True)
			dataFile.save()
			sys.stdout.write("Saved %s  %s\n"%(dataFile.filename, dataFile.abspath))
