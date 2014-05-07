# Copyright (c) 2014, the GREAT3 executive committee (http://www.great3challenge.info/?q=contacts)
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted
# provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions
# and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of
# conditions and the following disclaimer in the documentation and/or other materials provided with
# the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to
# endorse or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""@file plot_const_shear_metrics_2D.py

Plot some constant shear results as nice 2D surfaces.
"""

import os
import numpy as np
import g3metrics

# Currently setup for doing NOISE FREE simulations!

NIMS = 200               # Number of images per set, always 200 for G10
NGALS_PER_IM = 1     # In GREAT08/GREAT10 there were 10000 galaxies, but here there's no noise
TRUE_SIGMA = 0.04        # Standard deviation of true input shears for normal distribution
TRUE_RANGE = 0.08        # Range of true input shears for a uniform distribution
NFIELDS = 10             # Don't necessarily need to have NIMS input shears. But easiest if
                         # NFIELDS is an integral fraction of NIMS..

CFID = 2.e-4 # Fiducial, "target" m and c values
MFID = 2.e-3 #

# Plotting ranges of interest
CMIN = CFID 
CMAX = 1.e-1
MMIN = MFID
MMAX = 1.e0

NBINS = 7 # Number of bins to plot in the ranges above
NMONTE = 30  # Number of montecarlo samples
NOISE_SIGMA = 0.  # Noise due to pixel shot noist on a shear estimate, per galaxy: noise free now!

# Generate arrays of values
cvals = CMIN * (CMAX / CMIN)**(np.arange(NBINS) / float(NBINS - 1.)) # geometric series
mvals = MMIN * (MMAX / MMIN)**(np.arange(NBINS) / float(NBINS - 1.))

cgrid, mgrid = np.meshgrid(cvals, mvals) # 2D arrays covering full space

# Generate the truth tables
g1true, g2true = g3metrics.make_const_truth_normal_dist(NFIELDS, NIMS, true_sigma=TRUE_SIGMA)

# Create empty storage arrays in which to put
QZ1_mcboth = np.empty((NBINS, NBINS))
QZ2_mcboth = np.empty((NBINS, NBINS))

# Loop over mvalues making independent submissions at each c, m combination
for i in range(NBINS):

    for j in range(NBINS):
        
        g1sub, g2sub = g3metrics.make_submission_const_shear(cgrid[i, j], cgrid[i, j],
                                                           mgrid[i, j], mgrid[i, j],
                                                           g1true, g2true,
                                                           ngals_per_im=NGALS_PER_IM,
                                                           noise_sigma=NOISE_SIGMA)
        QZ1_mcboth[i, j] = g3metrics.metricQZ1_const_shear(g1sub, g2sub, g1true, g2true,
                                                         cfid=CFID, mfid=MFID)[0]
        QZ2_mcboth[i, j] = g3metrics.metricQZ2_const_shear(g1sub, g2sub, g1true, g2true,
                                                         cfid=CFID, mfid=MFID)[0]

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt

if not os.path.isdir('./plots'):
    os.mkdir('./plots')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(np.log10(mgrid), np.log10(cgrid), np.log10(QZ1_mcboth), rstride=1, cstride=1,
                alpha=0.3, color='r')
ax.view_init(25, 60)
plt.xlabel(r'log$_{10}$(m = m1 = m2)')
plt.ylabel(r'log$_{10}$(c = c1 = c2)')
plt.title('QZ1 metric')
ax.set_zlabel(r'log$_{10}$(Q)')
ax.set_zlim(0, 3)
outfile = './plots/QZ1_2D_mcboth.png'
plt.savefig(outfile)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(np.log10(mgrid), np.log10(cgrid), np.log10(QZ2_mcboth), rstride=1, cstride=1,
                alpha=0.3, color='g')
ax.view_init(25, 60)
plt.xlabel(r'log$_{10}$(m = m$_1$ = m$_2$)')
plt.ylabel(r'log$_{10}$(c = c$_1$ = c$_2$)')
plt.title('QZ2 metric')
ax.set_zlabel(r'log$_{10}$(Q)')
ax.set_zlim(0, 3)
outfile = './plots/QZ2_2D_mcboth.png'
plt.savefig(outfile)


# Create empty storage arrays in which to put
QZ1_m1c1 = np.empty((NBINS, NBINS))
QZ2_m1c1 = np.empty((NBINS, NBINS))

# Loop over mvalues making independent submissions at each c, m combination
for i in range(NBINS):

    for j in range(NBINS):
        
        g1sub, g2sub = g3metrics.make_submission_const_shear(cgrid[i, j], CFID,
                                                             mgrid[i, j], MFID,
                                                             g1true, g2true,
                                                             ngals_per_im=NGALS_PER_IM,
                                                             noise_sigma=NOISE_SIGMA)
        QZ1_m1c1[i, j] = g3metrics.metricQZ1_const_shear(g1sub, g2sub, g1true, g2true,
                                                         cfid=CFID, mfid=MFID)[0]
        QZ2_m1c1[i, j] = g3metrics.metricQZ2_const_shear(g1sub, g2sub, g1true, g2true,
                                                         cfid=CFID, mfid=MFID)[0]

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(np.log10(mgrid), np.log10(cgrid), np.log10(QZ1_m1c1), rstride=1, cstride=1,
                alpha=0.3, color='r')
ax.view_init(25, 60)
plt.xlabel(r'log$_{10}$(m$_1$)')
plt.ylabel(r'log$_{10}$(c$_1$)')
plt.title('QZ1 metric [m$_2$ = '+str(MFID)+', c$_2$ = '+str(CFID)+']')
ax.set_zlabel(r'log$_{10}$(Q)')
ax.set_zlim(0, 3)
outfile = './plots/QZ1_2D_m1c1.png'
plt.savefig(outfile)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(np.log10(mgrid), np.log10(cgrid), np.log10(QZ2_m1c1), rstride=1, cstride=1,
                alpha=0.3, color='g')
ax.view_init(25, 60)
plt.xlabel(r'log$_{10}$(m$_1$)')
plt.ylabel(r'log$_{10}$(c$_1$)')
plt.title('QZ2 metric [m$_2$ = '+str(MFID)+', c$_2$ = '+str(CFID)+']')
ax.set_zlabel(r'log$_{10}$(Q)')
ax.set_zlim(0, 3)
outfile = './plots/QZ2_2D_m1c1.png'
plt.savefig(outfile)


m1grid, m2grid = np.meshgrid(mvals, mvals) # 2D arrays covering full space for m1 m2 case

# Create empty storage arrays in which to put
QZ1_m1m2 = np.empty((NBINS, NBINS))
QZ2_m1m2 = np.empty((NBINS, NBINS))

# Loop over mvalues making independent submissions at each c, m combination
for i in range(NBINS):

    for j in range(NBINS):
        
        g1sub, g2sub = g3metrics.make_submission_const_shear(CFID, CFID,
                                                             m1grid[i, j], m2grid[i, j],
                                                             g1true, g2true,
                                                             ngals_per_im=NGALS_PER_IM,
                                                             noise_sigma=NOISE_SIGMA)
        QZ1_m1m2[i, j] = g3metrics.metricQZ1_const_shear(g1sub, g2sub, g1true, g2true,
                                                         cfid=CFID, mfid=MFID)[0]
        QZ2_m1m2[i, j] = g3metrics.metricQZ2_const_shear(g1sub, g2sub, g1true, g2true,
                                                         cfid=CFID, mfid=MFID)[0]

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(np.log10(m1grid), np.log10(m2grid), np.log10(QZ1_m1m2), rstride=1, cstride=1,
                alpha=0.3, color='r')
ax.view_init(25, 60)
plt.xlabel(r'log$_{10}$(m$_1$)')
plt.ylabel(r'log$_{10}$(m$_2$)')
plt.title('QZ1 metric [c$_1$ = c$_2$ = '+str(CFID)+']')
ax.set_zlabel(r'log$_{10}$(Q)')
ax.set_zlim(0, 3)
outfile = './plots/QZ1_2D_m1m2.png'
plt.savefig(outfile)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(np.log10(m1grid), np.log10(m2grid), np.log10(QZ2_m1m2), rstride=1, cstride=1,
                alpha=0.3, color='g')
ax.view_init(25, 60)
plt.xlabel(r'log$_{10}$(m$_1$)')
plt.ylabel(r'log$_{10}$(m$_2$)')
plt.title('QZ2 metric [c$_1$ = c$_2$ = '+str(CFID)+']')
ax.set_zlabel(r'log$_{10}$(Q)')
ax.set_zlim(0, 3)
outfile = './plots/QZ2_2D_m1m2.png'
plt.savefig(outfile)
