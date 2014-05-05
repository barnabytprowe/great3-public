from leaderboard.models import Board, Entry, Team, save_submission_file, too_many_entries_in_last_day, MAXIMUM_ENTRIES_PER_DAY, Method, SUBMISSION_SUSPENDED
from django.http import HttpResponse, HttpResponseRedirect, Http404
from django.shortcuts import render
from django import forms
from django.contrib.auth.decorators import login_required

def check_unique_submission_name(name):
	if Entry.objects.filter(name=name).exists():
		raise forms.ValidationError("This name has already been used for a submission.")

def check_non_empty_name(name):
	if not name.strip():
		raise forms.ValidationError("Please fill in this field")


		
class SubmissionForm(forms.Form):
	title = forms.CharField(max_length=128, validators=[check_unique_submission_name, check_non_empty_name], required=True)
	file_upload = forms.FileField()
	notes = forms.CharField(max_length=512, required=False, widget = forms.Textarea)
	class Media:
		js = ('js/entry_form.js',)

	def __init__(self, *args, **kwargs):
		teams = kwargs.pop("teams",None)
		if teams is None:
			self.only_valid_team = kwargs.pop("team")
		methods = kwargs.pop("methods", None)
		super(SubmissionForm, self).__init__(*args, **kwargs)
		if teams is not None:
			teams = [('None','- Select team -')] + [(team.id,team.name) for team in teams]
			self.fields['team'] = forms.ChoiceField(choices=teams)
		if methods is None:
			self.fields['method'] = forms.CharField(max_length=128,validators=[check_non_empty_name], required=True)
		else:
			method_names = [(method,method) for method in methods]
			method_names.append(('None','New method...'))
			self.fields['method'] = forms.ChoiceField(choices=method_names)
			self.fields['new_method_name'] = forms.CharField(max_length=128,required=False)

	def clean(self):
		if 'team' in self.fields:
			team_id = self.cleaned_data.get('team')
			if team_id == 'None':
				raise forms.ValidationError("Please select a team.")
			try:
				team = Team.objects.get(id=team_id)
			except Team.DoesNotExist:
				raise forms.ValidationError("Please select a valid team.")
			self.cleaned_data['team'] = team
		else:
			team = self.only_valid_team

		method_name=self.cleaned_data.get('method')
		if method_name=='None':
			method_name=self.cleaned_data.get('new_method_name').strip()
		if method_name=='':
			raise forms.ValidationError("Please select a method name from the drop-down list or enter a new name.")
		try:
			method=Method.objects.get(name=method_name)
			if method.team != team:
				raise forms.ValidationError("A different team has already used this method name")
		except Method.DoesNotExist:
			method=Method(name=method_name, team=team)
			method.save()
		self.cleaned_data['method']=method

		return self.cleaned_data

def index(request):
	""" List of all leaderboards """
	boards = Board.objects.filter(enabled=True, blind=True).order_by('name')
	data = dict(boards=boards)
	return render(request,'leaderboard/board_list.html',data)


def post_challenge_index(request):
	""" List of all leaderboards """
	boards = Board.objects.filter(enabled=True, blind=False).order_by('name')
	data = dict(boards=boards)
	return render(request,'leaderboard/post_challenge_board_list.html',data)



def detail(request, board_id):
	""" Detail view of a single leaderboard """
	try:
		board = Board.objects.get(id=board_id)
	except Board.DoesNotExist:
		raise Http404
	if not board.enabled:
		raise Http404
	entries = board.entry_set.order_by('-score')
	data = dict(board=board, entries=entries)
	return render(request,'leaderboard/board_detail.html',data)

def methods_for_teams(*teams):
	methods = []
	for team in teams:
		methods += sorted([method.name for method in team.method_set.all()])
	return methods


@login_required
def submit(request, board_id):
	if SUBMISSION_SUSPENDED:
		return render(request, 'leaderboard/contest_finished.html')
	try:
		board = Board.objects.get(id=board_id)
	except Board.DoesNotExist:
		raise Http404
	if not board.enabled:
		raise Http404
	if board.frozen:
		return render(request, 'leaderboard/contest_finished.html')
	all_teams = request.user.get_profile().teams.all()

	if not all_teams:
		data = dict(user=request.user, limit=MAXIMUM_ENTRIES_PER_DAY)
		return render(request, 'leaderboard/noteam.html', data)

	valid_teams = [team for team in all_teams if not too_many_entries_in_last_day(team, board)]
	excluded_teams = [team for team in all_teams if team not in valid_teams]
	if not valid_teams:
		data=dict(excluded_teams=excluded_teams, entry_limit=MAXIMUM_ENTRIES_PER_DAY)
		return render(request, 'leaderboard/toosoon.html',data)

	if request.method == 'POST':
		if len(valid_teams)==1:
			team = valid_teams[0]
			methods = methods_for_teams(team)
			form = SubmissionForm(request.POST, request.FILES, team=team, methods=methods)
		else:
			team=None
			methods = methods_for_teams(*valid_teams)
			form = SubmissionForm(request.POST, request.FILES, teams=valid_teams, methods=methods)

		if form.is_valid():
			if team is None: team = form.cleaned_data['team']
			#check file size here, before upload!
			ok = save_submission_file(request.FILES['file_upload'], 
				form.cleaned_data['title'], 
				form.cleaned_data['notes'], 
				form.cleaned_data['method'], 
				request.user, team, board)
			if ok:
				return HttpResponseRedirect('/leaderboard/board/%s/submitted/'%board_id)
			else:
				raise ValueError("Submission Failed")
		else: 
			data= dict(form=form, board=board, teams=valid_teams, excluded_teams=excluded_teams, entry_limit=MAXIMUM_ENTRIES_PER_DAY)
			return render(request, 'leaderboard/submit.html', data)
	else:
		methods=methods_for_teams(*valid_teams)
		if len(valid_teams)==1:
			form = SubmissionForm(methods=methods, team=valid_teams[0])
		else:
			form = SubmissionForm(teams=valid_teams, methods=methods)
		data = dict(form=form, board=board, teams=valid_teams, excluded_teams=excluded_teams, entry_limit=MAXIMUM_ENTRIES_PER_DAY)
		return render(request, 'leaderboard/submit.html', data)

def submitted(request, board_id):
	data = dict(board_id=board_id)
	return render(request, 'leaderboard/submitted.html', data)

