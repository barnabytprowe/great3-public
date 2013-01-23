from leaderboard.models import Team
from django.http import HttpResponse
from django.shortcuts import render


def index(request):
	teams = Team.objects.all().order_by('name')
	data = dict(team=team)
	return render(request, 'leaderboard/team_list.html', data)


def detail(request, team_id):
	team = Team.objects.get(id=team_id)
	data = dict(team=team)
	return render(request, 'leaderboard/team_detail.html', data)

