from django.core.management.base import BaseCommand, CommandError
from leaderboard.models import Entry, PLACEHOLDER_SCORE
from optparse import make_option
import sys

class Command(BaseCommand):
    args = ''
    help = 'Checks for entries that do not have a score yet and processes them'
    option_list = BaseCommand.option_list + (
        make_option('--i-am-sure',
            action='store_true',
            dest='sure',
            default=False,
            help='Confirm that you are really sure that you want to do this'),
        )

    def handle(self, *args, **options):
        sure = options['sure']
        if not sure:
            self.stdout.write('This command resets all great3 scores to their placeholder value.\n')
            self.stdout.write('Run this command with the --i-am-sure flag to confirm you really want to.\n')
            return
        entries = Entry.objects.all()
        for entry in entries:
            sys.stdout.write('Resetting%s\n'%entry.name)
            entry.score = PLACEHOLDER_SCORE
            entry.save()
        sys.stdout.write('\nAll scores reset\n\n')

