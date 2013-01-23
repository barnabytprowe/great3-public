from django.shortcuts import render
from django.http import HttpResponse
from models import Team, Board, create_data


def front_page(request):
	# create_data()
	# print "CREATED DATA"
	winners, score = Team.winning_teams()
	data = dict(winners=winners, score=score)
	if request.user.is_authenticated():
		return render(request, 'leaderboard/front_page_logged_in.html', data)
	else:
		return render(request, 'leaderboard/front_page.html', data)


