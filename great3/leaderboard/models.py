from django.db import models
from django.contrib.auth.models import User
from django.db.models.signals import post_save
import os

# Create your models here.
SUBMISSION_SAVE_PATH="/Users/jaz/src/great3-server/results/"


class UserProfile(models.Model):
    user = models.OneToOneField(User)
    teams = models.ManyToManyField('Team', related_name='users')


def create_user_profile(sender, instance, created, **kwargs):
    if created:
        UserProfile.objects.create(user=instance)

post_save.connect(create_user_profile, sender=User)


class Team(models.Model):
	name = models.CharField(max_length=128)
	members = models.CharField(max_length=512)
	# users = models.ManyToManyField(User, related_name='teams')

	def __unicode__(self):
		return self.name


class Board(models.Model):
	name = models.CharField(max_length=128)
	space = models.BooleanField()
	varying = models.BooleanField()

	def __unicode__(self):
		return self.name

	def winner(self):
		if self.entry_set.all():
			return self.entry_set.order_by('score')[0]
		else:
			return " "

PLACEHOLDER_SCORE = -1.0

class Entry(models.Model):
	team = models.ForeignKey('Team')
	name = models.CharField(max_length=128)
	user = models.ForeignKey(User)
	board = models.ForeignKey('Board')
	score = models.FloatField(default=PLACEHOLDER_SCORE)
	date = models.DateTimeField(auto_now_add=True)

	def __unicode__(self):
		return self.name

	def get_filename(self):
		return os.path.join(SUBMISSION_SAVE_PATH, str(self.id))


#Some initial teams and people
# Boards

def create_data():
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
	Entry(name='nbc1', team=ucl, user=jaz, board=vanilla, score=20.4).save()
	Entry(name='nbc1', team=ucl, user=barney, board=vanilla, score=884.4).save()

	Entry(name='shapelets1', team=libertarians, user=mike, board=vanilla, score=510.2).save()
	Entry(name='shapelets2', team=libertarians, user=mike, board=vanilla, score=456.4).save()

	Entry(name='deimos1', team=nasa, user=mike, board=space, score=101.2).save()
	Entry(name='deimos2', team=nasa, user=mike, board=space, score=141.4).save()



# class Entry(models.Model):
# 	team = models.ForeignKey('Team')
# 	user = models.ForeignKey('user')
# 	leaderboard = models.ForeignKey('Board')
# 	score = models.FloatField()
# 	date = models.DateTimeField(auto_now_add=True)

	



