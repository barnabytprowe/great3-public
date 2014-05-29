import matplotlib
matplotlib.use("Agg")
import astropy.table
import iso8601
import time
import pylab
import numpy as np
import collections
from matplotlib.patches import Ellipse
import os


def parse_filename_date(filename):
	date_text = filename[-35:-3]
	d = iso8601.parse_date(date_text)
	return d

def timestamp(filename):
	d = parse_filename_date(filename)
	return time.mktime(d.timetuple())


class TeamColorsIndices(object):
	_colors = ['red', 'blue', 'ForestGreen', 'orange', 'black', 'yellow', 'palegreen', 
	'PaleVioletRed','cyan', 'grey','maroon', 'lavender','lightblue', 
	'magenta', 'olive', 'orangered']
	def __init__(self):
		self.teams = {}
		self.index = 0
	def __call__(self, team):
		if team in self.teams:
			return self.teams[team]
		else:
			c = self._colors[self.index%len(self._colors)], self.index
			self.index+=1
			self.teams[team] = c
			return c
	def name(self, index):
		for t, x in self.teams.items():
			if x[1]==index:
				return t
		return ""


def make_video(board):
	try:
		os.mkdir("plots1/%s"%board)
	except:
		pass
	table = astropy.table.Table.read("./export.fits")
	teams = TeamColorsIndices()
	rows_by_day = collections.defaultdict(list)
	for i,row in enumerate(table):
		d = parse_filename_date(row['filename']).date()
		if i==0:
			start_day = d
		day_count = (d-start_day).days
		rows_by_day[day_count].append(row)
	best_by_team = {}
	team_colors = TeamColorsIndices()
	count = 0
	for d in xrange(day_count+1):
		for row in rows_by_day[d]:
			if row['board']!=board+'-post': continue
			team = row['team']
			score = row['score']
			if team not in best_by_team or score>best_by_team[team]['score']:
				best_by_team[team] = row
				count += 1
				#print team, float(best_by_team[team]['score'])
		pylab.figure()
		ax = pylab.gca()
		for team,row in best_by_team.items():
			color,_ = team_colors(team)
			mp = row['m_plus']
			mm = row['m_cross']
			dmp = row['delta_m_plus']
			dmm = row['delta_m_cross']

			ell = Ellipse((mp,mm), width=dmp, height=dmm)
			ell.set_color(color)
			ell.set_alpha(0.2)
			ax.add_artist(ell)
		ax.set_xlim(-8e-2, 8e-2)
		ax.set_ylim(-8e-2, 8e-2)
		pylab.savefig("plots1/%s/%.5d.png"%(board,d))
		pylab.close()
	print 
	print "COUNT = ", count
	print
	cmd = "ffmpeg -pattern_type glob -i 'plots1/%s/*.png' plots1/%s.mp4" % (board, board)
	print cmd
	os.system(cmd)


boards = ['%s-%s-%s' % (a,b,c) 
	for a in ['control', 'multiepoch', 'variable_psf', 'real_galaxy']
	for b in ['ground', 'space']
	for c in ['constant', 'variable']
]


for board in boards:
	try:
		os.mkdir("plots1")
	except:
		pass
	if not board.endswith("variable"):
		print board
		make_video(board)
