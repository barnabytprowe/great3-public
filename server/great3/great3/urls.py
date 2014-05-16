from django.conf.urls import patterns, include, url
from  django.contrib.auth.views import login
from django.views.generic.simple import redirect_to
import registration

# This enables the admin site to find
# all the models it controls
from django.contrib import admin
admin.autodiscover()

# These are the main URLs provided by the system.
# Most of them delegate to the leaderboard/urls.py
# patterns
urlpatterns = patterns('',
    url(r'^admin/', include(admin.site.urls)),
    url(r'^/?$', redirect_to, {'url': 'leaderboard/'}),    
    url(r'^leaderboard/', include('leaderboard.urls')),
    # url(r'^login$','django.contrib.auth.views.login'),
    # url(r'^logout$','django.contrib.auth.views.logout', {'next_page': '/'}),
	url(r'^accounts/', include('registration.backends.default.urls')),
)

