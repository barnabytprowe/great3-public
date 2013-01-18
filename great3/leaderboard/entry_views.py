from leaderboard.models import Board, Entry
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import render
from django import forms


def detail(request, entry_id):
	entry = Entry.objects.get(id=entry_id)
	data=dict(entry=entry)
	return render(request, "leaderboard/entry_detail.html", data)
