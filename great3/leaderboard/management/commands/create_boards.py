from django.core.management.base import BaseCommand, CommandError
import django.db.utils
from leaderboard.models import Board, EXPERIMENT_CHOICES
from great3.settings import installation_base
import numpy as np
import os
import datetime
import pytz
import codecs
import sys

SCORE_LOG_FILE = os.path.join(installation_base, "results", "results.log")
texts = {
    'Control': 'The baseline branch with parametric galaxy and PSF models, with a single image per galaxy. ',
    'Realistic Galaxy': 'Galaxies models are drawn from real deep data. ',
    'Realistic PSF': 'PSF models are simulations of realistic experiments. ',
    'Multi-epoch': 'Each galaxy is observed multiple times with different PSFs and noise realizations. ',
    'Everything': 'All the effects from the different experiments are present: galaxies and PSFs are realistic, and there are multiple images of each galaxy. ',
}



class Command(BaseCommand):    
    args = ''
    help = 'Checks for entries that do not have a score yet and processes them'

    def handle(self, *args, **options):
        for experiment in EXPERIMENT_CHOICES:
            for variable in [True,False]:
                for space in [True,False]:
                    name = experiment[0]
                    notes = texts[experiment[0]]
                    if variable:
                        name += '-variable'
                        notes += 'Shear and magnification vary across the image according to a power spectrum. '
                    else:
                        name += '-constant'
                        notes += 'Shear and magnification are the same for every galaxy. '
                    if space:
                        name += '-space'
                        notes += ' Observing conditions simulate those of a next-generation space telescope.'
                    else:
                        name += '-ground'
                        notes += ' Observing conditions simulate those of a ground-based survey.'
                    try:
                        board = Board(experiment=experiment[0], varying=variable, space=space, notes=notes, name=name)
                        board.save()
                    except django.db.utils.IntegrityError:
                        print 'You cannot run this command a second time without deleting the database'
                        sys.exit(1)
