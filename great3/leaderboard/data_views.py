from leaderboard.models import Board, AdminDataFile
from django.shortcuts import render
from django.contrib.auth.decorators import login_required

def index(request):
	boards = Board.objects.all()
	data = dict(boards=boards)
	if request.user and request.user.is_staff:
		admin_files = AdminDataFile.objects.all()
		data['admin_files'] = admin_files
		return render(request, "leaderboard/data_list_admin.html", data)
	return render(request, "leaderboard/data_list.html", data)



@login_required
def admin_file(request, index):
	if not request.user.is_staff:
		raise Http404		
	abs_filename = AdminDataFile.objects.get(id=index).abspath
	response = django.http.HttpResponse()
	del response['content-type']
	response['X-Sendfile'] = abs_filename
	return response

