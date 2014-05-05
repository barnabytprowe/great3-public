from django.db import models, transaction
from django.contrib.auth.models import User
from django.db.models.signals import post_save
from collections import defaultdict
import datetime
import os
import hashlib
import random
import pytz
import unicodedata
import re
import itertools

#Setting this to true will suspend all submissions!
SUBMISSION_SUSPENDED = False


# This function taken from django and slightly modified
def slugify(value):
    """
    Normalizes string, converts to lowercase, removes non-alpha characters,
    and converts spaces to hyphens.
    """
    value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore')
    value = unicode(re.sub('[^\w\s-]', '', value).strip().lower())
    return re.sub('[-\s]+', '-', value)


# Create your models here.
SUBMISSION_SAVE_PATH=os.path.join(os.path.split(__file__)[0], '..','..','results')
PLACEHOLDER_SCORE = -1.0
PLACEHOLDER_ERROR = 1.0
PLACEHOLDER_RANK = 1000
MAXIMUM_ENTRIES_PER_DAY = 4
MAX_BOARDS_FOR_SCORING = 5
NO_TIEBREAK = 0
MIN_SCORE_FOR_POINTS = 1.0
TIEBREAK_ALL_SCORES = 1
TIEBREAK_TIMESTAMP = 2
ENTRIES_WITHOUT_SUBMISSION_RATE_RESTRICTION = 5

EXPERIMENT_CHOICES = [
	('Control','Control'),
	('Realistic Galaxy','Realistic Galaxy'),
	('Realistic PSF', 'Realistic PSF'),
	('Multi-epoch', 'Multi-epoch'),
	('Everything', 'Everything')
	]

class UserProfile(models.Model):
    user = models.OneToOneField(User)
    teams = models.ManyToManyField('Team', related_name='users', null=True, blank=True)
    def __unicode__(self):
    	return self.user.username


def create_user_profile(sender, instance, created, **kwargs):
    if created:
        UserProfile.objects.create(user=instance)

post_save.connect(create_user_profile, sender=User)


class Method(models.Model):
	name = models.CharField(max_length=128, unique=True, error_messages={"unique":"Another team has already used this method name"})
	team = models.ForeignKey('Team')
	def __unicode__(self):
		return self.name

class Team(models.Model):
	name = models.CharField(max_length=128, unique=True)
	notes = models.CharField(max_length=512)
	score = models.IntegerField(default=0)
	rank = models.IntegerField(default=PLACEHOLDER_RANK)
	tainted = models.BooleanField(default=False)

	def __unicode__(self):
		return self.name

	def score_text(self):
		if self.tainted:
			return '*'
		else:
			return str(self.score)

	def earliest_ranked_entry_time(self):
		entries = [(entry.get_points(), entry) for entry in self.entry_set.all()]
		entries =sorted(entries)[::-1] #sorted by score (first tuple element), ascending
		entries = [e[1] for e in entries if e[0]>0] #all the entries with any points
		timestamps = sorted([entry.date for entry in entries])
		if timestamps:
			return timestamps[0]
		else:
			return None

	def calculate_tiebreak_score(self):
		scores = [entry.get_points() for entry in self.entry_set.all()]
		total_score = sum(scores)
		return total_score

	def calculate_score(self):
		scores = [entry.get_points() for entry in self.entry_set.all()]
		scores = sorted(scores)[::-1]
		new_score = sum(scores[:MAX_BOARDS_FOR_SCORING])
		self.score = new_score
		self.save()
		return new_score

	@classmethod
	@transaction.commit_manually()
	def break_ties(cls):
		# go through the scores and break any ties by adding
		# single points to the winner of them.
		teams = cls.objects.order_by('-score')
		scores = [team.score for team in teams]
		#now we need to iterate through groups
		#of teams with the same score and break any ties
		first_element = lambda x: x[0]
		pairs = zip(scores, teams)
		#Okay, this is a bit complicated.
		#First go through all the teams, grouped by their score
		for (score, teams_with_score) in itertools.groupby(pairs, first_element):
			teams_with_score = list(teams_with_score)
			# If there is only one team with that score then nothing to worry about;
			# we just carry on - no change to the score.  This is the usual case
			if len(teams_with_score)==1:
				continue
			#We also ignore teams with a score of zero.
			#one the grounds that it is not important what happens there
			if score==0:
				continue

			#If there is more than one team with that score we will need
			#the tiebreak scores, in sorted order
			tiebreak_scores = [(team.calculate_tiebreak_score(), team) for (_,team) in teams_with_score]
			tiebreak_scores.sort()

			# But now there can be a double tiebreak!
			# So we need to group up these scores again.
			# I'm sure this could be done much more elegantly.
			# We will get an ordering for all the teams with the same base
			#(non-tiebroken) score
			tiebroken_order = []
			# and now the grouping again:
			for tiebreak_score, doubly_tied_teams in itertools.groupby(tiebreak_scores, first_element):
				doubly_tied_teams = list(doubly_tied_teams)
				#Once again, it is most likely that there is only one
				#team with each tiebreak score. In that case the ordering
				#is clear - just the tiebreak score itself.
				#So the element in tiebroken_score is the first and only team.
				if len(doubly_tied_teams)==1:
					tiebroken_order.append(doubly_tied_teams[0][1])
					continue
				#On the other hand it is possible there is more than one
				#team with the same tiebreak score.  In that case we will need
				#an ordering within them.  We use the earliest ranked entry time.
				timestamps = [(team.earliest_ranked_entry_time(), team) for (_,team) in doubly_tied_teams]
				timestamps.sort()
				#Now we have a sorted list of the timestamps.
				#So the ordering in the tiebroken_order is the same as in this list.
				#the "_" is the timestamp itself that we used to sort
				for _, team in timestamps:
					tiebroken_order.append(team)
			#Now we have a tiebroken order.  We give the teams extra
			#individual points to break the tie.  Then we are done.
			for i,team in enumerate(tiebroken_order):
				team.score += i
				team.save()
		#Save!  We do this manually because there are a bunch of commits here
		#all to be done at once.
		transaction.commit()





	def top_entries_by_rank(self, n=MAX_BOARDS_FOR_SCORING):
		return self.entry_set.order_by('rank','date')[:n]

	@classmethod
	@transaction.commit_manually()
	def update_ranks(cls):
		teams = cls.objects.order_by('-score')
		for rank,team in enumerate(teams):
			team.rank = rank+1
			team.save()
		transaction.commit()

	def number_entries(self):
		return len(self.entry_set.all())

	def rank_text(self):
		if self.tainted:
			return "*"
		if self.rank==PLACEHOLDER_RANK:
			return "-"
		else:
			return str(self.rank)

	@classmethod
	def update_scores_and_ranks(cls):
		for team in cls.objects.all():
			team.calculate_score()
		cls.break_ties()
		cls.update_ranks()

	@classmethod
	def winning_teams(cls):
		teams = cls.objects.order_by('-score')
		if len(teams)==0:
			return [], PLACEHOLDER_SCORE, NO_TIEBREAK
		top_team = teams[0]
		winners = [top_team]
		best_score = top_team.score
		for team in teams[1:]:
			if team.score==best_score:
				winners.append(team)
			else:
				break
		tiebreak=0
		if len(winners)>1:
			scores = [(winner.calculate_tiebreak_score(), winner) for winner in winners]
			scores.sort()  #python trick.  sorts by the first element of the tuple
			best_tiebreak_score = scores[-1][0] #last score should be highest
			winners = [w[1] for w in scores if w[0]==best_tiebreak_score]
			tiebreak=1
		# Second tie-break!  Time-stamp on the earliest entry
		if len(winners)>1:
			timestamps = [(winner.earliest_ranked_entry_time(), winner) for winner in winners]
			timestamps.sort()
			winners = [timestamps[0][1]]
			tiebreak=2

		return winners, best_score, tiebreak


def score_for_rank(rank):
	""" The score that a team gets if their top-ranked
		entry into a board is at the given rank.
	"""
	if   rank==1: return 16000
	elif rank==2: return 8000
	elif rank==3: return 4000
	elif rank==4: return 2000
	elif rank==5: return 1000
	else: return 0



class Board(models.Model):
	name = models.CharField(max_length=128, unique=True)
	experiment = models.CharField(max_length=20, choices=EXPERIMENT_CHOICES)
	notes = models.CharField(max_length=512)
	space = models.BooleanField()
	varying = models.BooleanField()
	datafile = models.ForeignKey('PublicDataFile', null=True, blank=True)
	enabled = models.BooleanField()
	frozen = models.BooleanField()
	blind =  models.BooleanField()

	@transaction.commit_manually()
	def assign_ranks(self):
		entries = self.entry_set.order_by('-score', 'date')
		ranked_teams = []
		rank=1
		for entry in entries:
			if entry.team in ranked_teams or entry.team.tainted:
				entry.rank=PLACEHOLDER_RANK
			else:
				ranked_teams.append(entry.team)
				entry.rank = rank
				rank += 1
			entry.save()
		transaction.commit()

	@classmethod
	def assign_all_ranks(cls):
		for board in cls.objects.filter(enabled=True):
			board.assign_ranks()

	def number_entries(self):
		return len(self.entry_set.all())

	def __unicode__(self):
		return self.name

	def get_entry_at_rank(self, rank):
		try:
			return self.entry_set.filter(rank=rank).get()
		except Entry.DoesNotExist:
			return None


	def winner(self):
		return self.get_entry_at_rank(1)


class Entry(models.Model):
	team = models.ForeignKey('Team')
	name = models.CharField(max_length=128, unique=True)
	method = models.ForeignKey('Method')
	notes = models.CharField(max_length=512)
	user = models.ForeignKey(User)
	board = models.ForeignKey('Board')
	score = models.FloatField(default=PLACEHOLDER_SCORE)
	date = models.DateTimeField(auto_now_add=True)
	rank = models.IntegerField(default=PLACEHOLDER_RANK)
	points_text = models.CharField(max_length=24)
	m = models.FloatField(default=PLACEHOLDER_ERROR)
	c = models.FloatField(default=PLACEHOLDER_ERROR)
	delta_m = models.FloatField(default=PLACEHOLDER_ERROR)
	delta_c = models.FloatField(default=PLACEHOLDER_ERROR)
	blind = models.BooleanField(default=False)
	
	def __unicode__(self):
		return self.name

	def rank_text(self):
		if self.team.tainted:
			return '*'
		r = self.rank
		if r==PLACEHOLDER_RANK:
			return ""
		return str(r)

	def score_text(self):
		if self.score == PLACEHOLDER_SCORE:
			return "...Calculating..."
		else:
			score = "%.1f" % self.score
			if self.team.tainted:
				score += ' *'
			return score

	def get_points(self):
		if self.score<MIN_SCORE_FOR_POINTS:
			return 0
		return score_for_rank(self.rank)


	def compute_points_text(self):
		p = self.get_points()
		if p==0:
			self.points_text = ""
			self.save()
			return
		top_entries = self.team.top_entries_by_rank()
		in_top = self in top_entries
		if in_top:
			self.points_text = str(p)
		else:
			self.points_text = str(p) + ' [*]'
		self.save()


	def get_filename(self):
		team_name = slugify(self.team.name)
		method_name = slugify(self.method.name)
		entry_name = slugify(self.name)
		date_text = self.date.isoformat()
		filename = '%s-%s-%s-%s.g3' % (team_name, method_name, entry_name, date_text)
		return os.path.join(SUBMISSION_SAVE_PATH, filename)

	@classmethod
	def assign_all_points_texts(cls):
		for entry in cls.objects.all():
			entry.compute_points_text()


def recompute_scoring(*boards):
	if not boards:
		boards = Board.objects.filter(enabled=True)
	for board in boards:
		board.assign_ranks()
	Team.update_scores_and_ranks()
	Entry.assign_all_points_texts()

#Some initial teams and people
# Boards

def create_data():
	# JAZ I used this once to generate some initial data and then 
	# dumped it to YAML, from which it can be regenerated.
	# Preserved just in case.
	vanilla = Board(name="Vanilla", space=False, varying=False)
	space = Board(name="Space", space=True, varying=False)
	des = Board(name="DES", space=False, varying=True)
	euclid = Board(name="Euclid", space=True, varying=True)
	wfirst = Board(name="WFIRST", space=True, varying=False)

	vanilla.save()
	space.save()
	des.save()
	euclid.save()
	wfirst.save()

	lensfit = Team(name="Lensfit", members="Lance et al")
	ucl = Team(name="UCL-Manchester-Alliance", members="Sarah et al")
	libertarians = Team(name="RonPaul4Eva", members="Mike et al")
	nasa = Team(name="NASA", members="Jason et al")

	lensfit.save()
	ucl.save()
	libertarians.save()
	nasa.save()

	barney = User(username="barney", email="barney@example.com")
	rachel = User(username="rachel", email="rachel@example.com")
	mike = User(username="mike", email="mike@example.com")
	jaz = User.objects.get(username="jaz")

	barney.set_password('barney')
	rachel.set_password('rachel')
	mike.set_password('mike')

	barney.save()
	rachel.save()
	mike.save()
	# jaz.save()

	def add_teams(user, *teams):
		profile = user.get_profile()
		for team in teams:
			profile.teams.add(team)
		profile.save()
		user.save()

	add_teams(jaz, ucl)
	add_teams(barney, ucl, lensfit)
	add_teams(rachel, nasa)
	add_teams(mike, libertarians, nasa)

	Entry(name='nbc1', team=ucl, user=barney, board=vanilla, score=14.5).save()
	Entry(name='nbc2', team=ucl, user=jaz, board=vanilla, score=20.4).save()
	Entry(name='nbc3', team=ucl, user=barney, board=vanilla, score=884.4).save()

	Entry(name='shapelets1', team=libertarians, user=mike, board=vanilla, score=510.2).save()
	Entry(name='shapelets2', team=libertarians, user=mike, board=vanilla, score=456.4).save()

	Entry(name='deimos1', team=nasa, user=mike, board=space, score=101.2).save()
	Entry(name='deimos2', team=nasa, user=mike, board=space, score=141.4).save()

	recompute_scoring()


def save_submission_file(submission, name, notes, method, user, team, board):
	print "Sanity check the file size here"
	entry = Entry(team=team, name=name, notes=notes, user=user, method=method, board=board)
	entry.save()
	try:
		with open(entry.get_filename(), 'wb+') as destination:
			for chunk in submission.chunks():
				destination.write(chunk)
	except Exception as E:
		print "Could not save: %r!" % E
		entry.delete()
		return False
	return True




class MembershipRequest(models.Model):
	user = models.ForeignKey(User)
	team = models.ForeignKey(Team)
	token = models.CharField(max_length=40, unique=True)

	def generate_token(self):
		salt = hashlib.sha1(str(random.random())).hexdigest()[:5]
		username = unicode(self.user.username).encode('utf-8')
		teamname = unicode(self.team)
		self.token = hashlib.sha1(salt+username+teamname).hexdigest()
		return self.token



def user_is_member_of_team(user, team):
	teams = user.get_profile().teams.all()
	return team in teams


def too_many_entries_in_last_day(team, board):
	recent_entries = Entry.objects.filter(team=team, board=board).order_by('-date')
	if len(recent_entries)<ENTRIES_WITHOUT_SUBMISSION_RATE_RESTRICTION:
		return False
	try:
		test_entry = recent_entries[MAXIMUM_ENTRIES_PER_DAY-1]
	except IndexError:
		return False
	one_day_ago = datetime.datetime.utcnow().replace(tzinfo=pytz.utc) - datetime.timedelta(days=1.0)
	#I know this looks odd - it does not mean "more than one day ago."
	return test_entry.date > one_day_ago

class PublicDataFile(models.Model):
	filename = models.CharField(max_length=128)
	abspath = models.CharField(max_length=512)
	info = models.CharField(max_length=512)
	#miscellaneous files are files that are not
	#one of the data set tarballs, they are extra stuff
	miscellaneous = models.BooleanField(default=True)
	mirror1 = models.URLField(blank=True, null=True)
	mirror2 = models.URLField(blank=True, null=True)  # Just in case
	mirror3 = models.URLField(blank=True, null=True)
	def __unicode__(self):
		return self.filename

class AdminDataFile(models.Model):
	filename = models.CharField(max_length=128)
	abspath = models.CharField(max_length=512)
	def __unicode__(self):
		return self.filename
