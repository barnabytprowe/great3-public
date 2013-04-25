from leaderboard.models import Board, Entry, save_submission_file
from django.http import HttpResponse, HttpResponseRedirect, Http404
from django.shortcuts import render
from django import forms
from django.contrib.auth.decorators import login_required


def check_unique_submission_name(name):
	if Entry.objects.filter(name=name).exists():
		raise forms.ValidationError("This name has already been used for a submission.")
		
class SubmissionForm(forms.Form):
	title = forms.CharField(max_length=128, validators=[check_unique_submission_name])
	file_upload = forms.FileField()
	notes = forms.CharField(max_length=512, required=False, widget = forms.Textarea)
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
	try:
		board = Board.objects.get(id=board_id)
	except Board.DoesNotExist:
		raise Http404
	entries = board.entry_set.order_by('-score')
	data = dict(board=board, entries=entries)
	return render(request,'leaderboard/board_detail.html',data)



@login_required
def submit(request, board_id):
	teams = request.user.get_profile().teams.all()
	try:
		board = Board.objects.get(id=board_id)
	except Board.DoesNotExist:
		raise Http404

	if request.method == 'POST':
		if len(teams)==1:
			team = teams[0]
			form = SubmissionForm(request.POST, request.FILES)
		else:
			team=None
			form = SubmissionForm(request.POST, request.FILES, teams=teams)

		if form.is_valid():
			if team is None: team = form.cleaned_data['team']
			ok = save_submission_file(request.FILES['file_upload'], form.cleaned_data['title'], form.cleaned_data['notes'], request.user, team, board)
			if ok:
				return HttpResponseRedirect('/leaderboard/board/%s/submitted/'%board_id)
			else:
				raise ValueError("Submission Failed")
		else: 
			data= dict(form=form, board=board, teams=teams)
			return render(request, 'leaderboard/submit.html', data)
	else:
		if len(teams)==1:
			form = SubmissionForm()
		else:
			form = SubmissionForm(teams=teams)
		data = dict(form=form, board=board, teams=teams)
		return render(request, 'leaderboard/submit.html', data)

def submitted(request, board_id):
	data = dict(board_id=board_id)
	return render(request, 'leaderboard/submitted.html', data)

