from django.conf.urls import patterns, include, url
from team_views import TeamListView, TeamDetailView
import models



urlpatterns = patterns('leaderboard.board_views',
    url(r'^board/(\d+)/?$', 'detail'),
    url(r'^board/?$', 'index'),
	url(r'^board/(\d+)/submitted/?$', 'submitted'),
	url(r'^board/(\d+)/submit/?$', 'submit'),

)

urlpatterns += patterns('leaderboard.team_views',
    url(r'^team$', TeamListView.as_view()),
    url(r'^team/?$', 'team_view_index'),
)

urlpatterns += patterns('leaderboard.entry_views',
    url(r'^entry/(\d+)/?$', 'detail'),

)
urlpatterns += patterns('leaderboard.views',
	url(r'^$', 'front_page'),
	)