from leaderboard.models import Board, Entry
from django.http import HttpResponse, HttpResponseRedirect, Http404
from django.shortcuts import render
from django import forms


def detail(request, entry_id):
	try:
		entry = Entry.objects.get(id=entry_id)
	except Entry.DoesNotExist:
		raise Http404
	data=dict(entry=entry)
	return render(request, "leaderboard/entry_detail.html", data)
