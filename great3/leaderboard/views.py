from django.shortcuts import render
from django.http import HttpResponse
import models


def front_page(request):
	# print "CREATING DATA ... "
	# models.create_data()
	# print "CREATED DATA"
	winners, score, tiebreak = models.Team.winning_teams()
	print tiebreak
	if tiebreak==models.NO_TIEBREAK:
		tiebreak_text=""
	elif tiebreak==models.TIEBREAK_ALL_SCORES:
		tiebreak_text=" (on a tiebreak by all entry scores) "
	elif tiebreak==models.TIEBREAK_TIMESTAMP:
		tiebreak_text=" (on a double-tiebreak by earliest ranked entry)"
	else:
		tiebreak_text="" #this should not happen
	data = dict(winners=winners, score=score, tiebreak=tiebreak_text)
	if request.user.is_authenticated():
		return render(request, 'leaderboard/front_page_logged_in.html', data)
	else:
		return render(request, 'leaderboard/front_page.html', data)


