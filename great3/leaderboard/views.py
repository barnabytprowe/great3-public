from django.shortcuts import render
from django.http import HttpResponse
from django.contrib.auth.decorators import login_required
from models import Team, Board, create_data


def front_page(request):
	# create_data()
	# print "CREATED DATA"
	if request.user.is_authenticated():
		return render(request, 'leaderboard/front_page_logged_in.html')
	else:
		return render(request, 'leaderboard/front_page.html')


@login_required
def submit(request):
	try:
		user_teams = request.user.get_profile().teams.all()
	except Team.DoesNotExist:
		return render(request, 'leaderboard/noteam.html')
	if len(user_teams)==1:
		request.session['active_team'] = user_teams[0]
	if 'active_team' in request.session:
		return render(request, 'leaderboard/submit.html', {'active_team':request.session['active_team']})
	else:
		return render(request,'leaderboard/active_team.html', )
