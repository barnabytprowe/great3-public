from leaderboard.models import Team
from django.http import HttpResponse
from django.shortcuts import render
from django.contrib.auth.decorators import login_required


def index(request):
	teams = Team.objects.all().order_by('name')
	data = dict(teams=teams)
	return render(request, 'leaderboard/team_list.html', data)


def detail(request, team_id):
	team = Team.objects.get(id=team_id)
	data = dict(team=team)
	return render(request, 'leaderboard/team_detail.html', data)

@login_required
def setup(request):
	profile = request.user.get_profile()
	teams = profile.teams.all()
	data=dict(teams=teams)
	return render(request, 'leaderboard/team_setup.html', data)