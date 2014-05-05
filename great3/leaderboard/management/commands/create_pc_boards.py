from django.core.management.base import BaseCommand, CommandError
import django.db.utils
from leaderboard.models import Board, Entry, PLACEHOLDER_SCORE
from great3.settings import installation_base
import numpy as np
import os
import datetime
import pytz
import codecs
import sys



class Command(BaseCommand):    
    args = ''
    help = 'Clone the challenge leaderboards into the post-challenge ones'

    def handle(self, *args, **options):

        #First clone the boards
        boards = Board.objects.all()
        for b in boards:
            # The new board should not be blind or frozen
            b.frozen = False
            b.blind = False

            #Give it a new name
            b.name += '-post'

            #Unset the primary key to save to a new object
            b.pk = None
            b.save()

        #Now the entries on those boards
        entries = Entry.objects.all()
        for e in entries:
            #This entry was originally blind, so reflect that.
            #Actually this is probably true already
            e.blind = True

            #Add something to the name
            e.name += '_pc'

            #Tell the system that this was imported
            e.imported = True

            # Set the new leaderboard to which it will apply
            e.board = Board.objects.get(name=e.board.name+'-post')

            #Reset the score so it will be re-analyzed.
            #This will also set m and c
            e.score = PLACEHOLDER_SCORE

            #Unset the primary key to save to a new object
            e.pk = None
            e.save()
