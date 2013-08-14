from leaderboard.models import Board, AdminDataFile, PublicDataFile
from django.shortcuts import render
from django.contrib.auth.decorators import login_required
import django.http


def index(request):
	boards = Board.objects.all()
	misc_public_files = PublicDataFile.objects.filter(miscellaneous=True)
	data = dict(boards=boards, public_files=misc_public_files)
	if request.user and request.user.is_staff:
		admin_files = AdminDataFile.objects.all()
		data['admin_files'] = admin_files
		return render(request, "leaderboard/data_list_admin.html", data)
	return render(request, "leaderboard/data_list.html", data)



@login_required
def admin_file(request, filename):
	if not request.user.is_staff:
		raise django.http.Http404		
	dataFile = AdminDataFile.objects.get(filename=filename)
	abs_filename = dataFile.abspath
	response = django.http.HttpResponse()
	del response['content-type']
	response['X-Sendfile'] = abs_filename
	return response


def public_file(request, filename):
	dataFile = PublicDataFile.objects.get(filename=filename)
	abs_filename = dataFile.abspath
	response = django.http.HttpResponse()
	del response['content-type']
	response['X-Sendfile'] = abs_filename
	return response

