from django.conf.urls import patterns, include, url
from team_views import TeamListView, TeamDetailView, team_view_index
from django.contrib import databrowse
import models
databrowse.site.register(models.Team, models.Board)



urlpatterns = patterns('leaderboard.board_views',
    url(r'^board/(\d+)/?$', 'detail'),
    url(r'^board/?$', 'index'),
	url(r'^board/(\d+)/submit/?$', 'submit'),
	url(r'^board/submit/?$', 'test'),
	url(r'^board/submitted/?$', 'submitted'),


)

urlpatterns += patterns('leaderboard.team_views',
    url(r'^team$', TeamListView.as_view()),
    url(r'^team/?$', team_view_index),
    url(r'^browse/(.*)', databrowse.site.root),

)

urlpatterns += patterns('leaderboard.views',
	url(r'^submit/?$', 'submit'),
	url(r'^$', 'front_page'),
	)