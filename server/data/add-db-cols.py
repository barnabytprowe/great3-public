import sqlite3

con = sqlite3.connect('great3.sqlite3')
c = con.cursor()
c.execute("alter table leaderboard_board add column 'frozen' 'bool' NOT NULL default 1")
c.execute("alter table leaderboard_board add column 'blind' 'bool' NOT NULL  default 1")

c.execute("alter table leaderboard_entry add column 'm_plus' 'float' NOT NULL default 1.0")
c.execute("alter table leaderboard_entry add column 'c_plus' 'float' NOT NULL default 1.0")
c.execute("alter table leaderboard_entry add column 'm_cross' 'float' NOT NULL default 1.0")
c.execute("alter table leaderboard_entry add column 'c_cross' 'float' NOT NULL default 1.0")
c.execute("alter table leaderboard_entry add column 'delta_m_plus' 'float' NOT NULL default 1.0")
c.execute("alter table leaderboard_entry add column 'delta_m_cross' 'float' NOT NULL default 1.0")
c.execute("alter table leaderboard_entry add column 'delta_c_plus' 'float' NOT NULL default 1.0")
c.execute("alter table leaderboard_entry add column 'delta_c_cross' 'float' NOT NULL default 1.0")
c.execute("alter table leaderboard_entry add column 'blind' 'bool' NOT NULL default 0")
c.execute("alter table leaderboard_entry add column 'imported' 'bool' NOT NULL default 0")
