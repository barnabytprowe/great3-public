from django.conf.urls import patterns, include, url
import models



urlpatterns = patterns('leaderboard.board_views',
    url(r'^board/(\d+)/?$', 'detail'),
    url(r'^board/?$', 'index'),
	url(r'^board/(\d+)/submitted/?$', 'submitted'),
	url(r'^board/(\d+)/submit/?$', 'submit'),

)

urlpatterns += patterns('leaderboard.team_views',
    url(r'^team/(\d+)/?$', 'detail'),
    url(r'^team/?$', 'index'),
    url(r'^team/setup?$', 'setup'),
)

urlpatterns += patterns('leaderboard.entry_views',
    url(r'^entry/(\d+)/?$', 'detail'),

)
urlpatterns += patterns('leaderboard.views',
	url(r'^$', 'front_page'),
	)


urlpatterns += patterns('leaderboard.data_views',
	url(r'^data/?$', 'index'),
	)
