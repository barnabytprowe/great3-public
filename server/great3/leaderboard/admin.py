from leaderboard.models import Team, Board, Entry, UserProfile, MembershipRequest, AdminDataFile, PublicDataFile
from django.contrib import admin
from django.contrib.auth.admin import UserAdmin
from django.contrib.auth.models import User

admin.site.register(Team)
admin.site.register(Board)
admin.site.register(Entry)
admin.site.register(MembershipRequest)

admin.site.register(AdminDataFile)
admin.site.register(PublicDataFile)


class UserProfileInline(admin.StackedInline):
    model = UserProfile
    can_delete = False
    verbose_name_plural = 'profile'

# Define a new User admin
class UserAdmin(UserAdmin):
    inlines = (UserProfileInline, )

# Re-register UserAdmin
admin.site.unregister(User)
admin.site.register(User, UserAdmin)
