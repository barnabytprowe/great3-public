from django.core.management.base import BaseCommand, CommandError
from leaderboard.models import User, create_user_profile


class Command(BaseCommand):
	args = ''
	help = 'Rebuild user profiles after importing users from backup'

	def handle(self, *args, **options):
		users=User.objects.all()
		for user in users:
			create_user_profile(user, user, True)
