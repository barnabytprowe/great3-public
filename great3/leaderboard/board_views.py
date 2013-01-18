from leaderboard.models import Board, Entry
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import render
from django import forms




class SubmissionForm(forms.Form):
	title = forms.CharField(max_length=128)
	file_upload = forms.FileField()
	def __init__(self, *args, **kwargs):
		teams = kwargs.pop("teams",None)
		super(SubmissionForm, self).__init__(*args, **kwargs)
		if teams is not None:
			self.fields['team'] = forms.ModelChoiceField(teams, empty_label="--Select Team--")


def index(request):
	""" List of all leaderboards """
	boards = Board.objects.all().order_by('name')
	data = dict(boards=boards)
	return render(request,'leaderboard/board_list.html',data)


def detail(request, board_id):
	""" Detail view of a single leaderboard """
	board = Board.objects.get(id=board_id)
	entries = board.entry_set.all()
	data = dict(board=board, entries=entries)
	return render(request,'leaderboard/board_detail.html',data)


def save_submission(submission, name, user, team, board):
	print "Sanity check the file size here"
	entry = Entry(team=team, name=name, user=user, board=board)
	entry.save()
	try:
		with open(entry.get_filename(), 'wb+') as destination:
			for chunk in submission.chunks():
				destination.write(chunk)
	except IOError:
		entry.delete()
		return False
	return True



def submit(request, board_id):
	teams = request.user.get_profile().teams.all()
	board = Board.objects.get(id=board_id)

	if request.method == 'POST':
		if len(teams)==1:
			team = teams[0]
			form = SubmissionForm(request.POST, request.FILES)
		else:
			team=None
			form = SubmissionForm(request.POST, request.FILES, teams=teams)
		if form.is_valid():
			if team is None: team = form.cleaned_data['team']
			save_submission(request.FILES['file_upload'], form.cleaned_data['title'], request.user, team, board)
			return HttpResponseRedirect('/leaderboard/board/submitted/')
		else: 
			print "Errors ", form.errors
			data= dict(form=form, board=board)
			return render(request, 'leaderboard/submit.html', data)
	else:
		if len(teams)==1:
			form = SubmissionForm()
		else:
			form = SubmissionForm(teams=teams)
		data = dict(form=form, board=board)
		return render(request, 'leaderboard/submit.html', data)

def submitted(request):
	return render(request, 'leaderboard/submitted.html')


def test(request):
	""" Make a new submission to a leaderboard """
	return HttpResponse('<BR>Submit page</BR>')