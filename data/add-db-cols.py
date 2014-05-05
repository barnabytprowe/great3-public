import sqlite3

con = sqlite3.connect('great3.sqlite3')
c = con.cursor()
c.execute("alter table leaderboard_board add column 'frozen' 'bool' NOT NULL default 1")
c.execute("alter table leaderboard_board add column 'blind' 'bool' NOT NULL  default 1")

c.execute("alter table leaderboard_entry add column 'm' 'float' NOT NULL default 1.0")
c.execute("alter table leaderboard_entry add column 'c' 'float' NOT NULL default 1.0")
c.execute("alter table leaderboard_entry add column 'delta_m' 'float' NOT NULL default 1.0")
c.execute("alter table leaderboard_entry add column 'delta_c' 'float' NOT NULL default 1.0")
c.execute("alter table leaderboard_entry add column 'blind' 'bool' NOT NULL default 1")
