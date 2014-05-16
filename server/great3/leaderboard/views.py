from django.shortcuts import render
from django.http import HttpResponse
import models


def front_page(request):
	# print "CREATING DATA ... "
	# models.create_data()
	# print "CREATED DATA"
	winners, score, tiebreak = models.Team.winning_teams()
	if score%100==0:
		tiebreak_text=""
	else:
		tiebreak_text=" (on a tiebreak)"
	data = dict(winners=winners, score=score, tiebreak=tiebreak_text)
	if request.user.is_authenticated():
		return render(request, 'leaderboard/front_page_logged_in.html', data)
	else:
		return render(request, 'leaderboard/front_page.html', data)


