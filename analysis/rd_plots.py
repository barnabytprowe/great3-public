'''Code for creating 'Football Manager'-style radar charts for comparing
performance of methods across different branches.

Ian Harrison
ian.harrison-2@manchester.ac.uk
@itrharrison
'''
import numpy as np
import pyfits as pf
from radar_factory import *

# get data, column names
hdulist = pf.open('./export.fits')
cols = hdulist[1].columns
tbdata = hdulist[1].data

# split the 'board' field into useable parts
branches_ns = tbdata.field('board')
branches = branches_ns.split('-')
branches, locations, shears, challenges = [list(block) for
                                           block in zip(*branches)]
branches = np.asarray(branches)
locations = np.asarray(locations)
shears = np.asarray(shears)
challenges = np.asarray(challenges)

teams = tbdata.field('team')
scores = np.asarray(tbdata.field('score'), dtype=float)

# find out all of the unique entries
team_list = np.unique(teams)
shear_list = np.unique(shears)
branch_list = np.unique(branches)
loc_list = np.unique(locations)
shear_list = np.unique(shears)
chall_list = np.unique(challenges)

# find how many branches each team entered
team_n_entries = np.asarray(np.zeros_like(team_list), dtype=int)
for i, team in enumerate(team_list):
  team_n_entries[i] = len(np.unique(branches_ns[(teams==team)]))

# for plotting nicety, only plot teams which entered n > n_limit branches
n_limit = 19
n_plots = sum(team_n_entries > n_limit)

plt.close('all')
for shear in shear_list:
  N = len(branch_list)
  theta = radar_factory(N, frame='polygon')
  spoke_labels = branch_list
  colors = ['r', 'g']
  fig = plt.figure()
  fig.suptitle(shear+' shear')
  
  data = np.ones_like(branch_list, dtype=float)
  i_pl = 1
  for i_t, team in enumerate(team_list):
    if team_n_entries[i_t] > n_limit:
      for i_c, loc in enumerate(loc_list):
        for i_b, branch in enumerate(branch_list):
          cut = (teams==team)*(locations==loc)*(shears==shear)*(branches==branch)
          if (sum(cut)==0):
            data[i_b] = 0.e0
          else:
            # team's greatest score for this branch
            data[i_b] = scores[cut].max()
        
        ax = fig.add_subplot(n_plots,1,i_pl, projection='radar')
        ax.plot(theta, data, color=colors[i_c])
        ax.fill(theta, data, facecolor=colors[i_c], alpha=0.25)
        ax.set_title(team, position=(0.5, 1.1),
                     horizontalalignment='center',
                     verticalalignment='center')
        ax.set_rmax(scores.max())
        ax.set_varlabels(spoke_labels)
        if i_pl==1:
          plt.legend(('space', 'ground'), loc='upper right')
      i_pl += 1

plt.show()