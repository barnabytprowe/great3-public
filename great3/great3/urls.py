from django.conf.urls import patterns, include, url
from  django.contrib.auth.views import login

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',
    # Examples:
    # url(r'^$', 'great3.views.home', name='home'),
    # url(r'^great3/', include('great3.foo.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    url(r'^admin/', include(admin.site.urls)),
    url(r'^leaderboard/', include('leaderboard.urls')),
    url(r'^login$','django.contrib.auth.views.login'),
)

