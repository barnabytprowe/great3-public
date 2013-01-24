from leaderboard.models import Board
from django.shortcuts import render

def index(request):
	boards = Board.objects.all()
	data = dict(boards=boards)
	return render(request, "leaderboard/data_list.html", data)
