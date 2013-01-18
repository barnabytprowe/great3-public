from leaderboard.models import Team
from django.http import HttpResponse
from django.views.generic import DetailView, ListView

def team_view_index(request):
	teams = Team.objects.all().order_by('name')
	output = ""
	for team in teams:
		users = ", ".join(user.user.username for user in team.users.all())
		output += '<BR> %s  [%s]' % (team.name,users)
	return HttpResponse(output)


class TeamListView(ListView):
	model=Team


class TeamDetailView(DetailView):
	model=Team
