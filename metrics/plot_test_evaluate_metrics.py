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
"""@file plot_test_evaluate_metrics.py

Plotting routines for outputs of the test_evaluate.py script (made by modifying module scope
constants appropriately in that module).
"""
import os
import sys
import cPickle
import numpy as np
import matplotlib.pyplot as plt
path, module = os.path.split(__file__)
sys.path.append(os.path.join(path, "..", "server", "great3")) # Appends the folder
                                                              # great3-private/server/great3 to
                                                              # sys.path
import evaluate
import g3metrics

NTEST = 300 # Number of realizations


if __name__ == "__main__":
    # Setup the storage dicts
    qabsl = {}
    qfrac = {}    
    qsqrd = {}
    qbymc = {}
    mcresults = {}
    plt.figure(figsize=(8, 8))
    # Loop over ground, space
    for obs_type in ("space", "ground"):

        noise_sigma = {"ground": 0.15, "space": 0.10}[obs_type]
        cvals = (evaluate.CFID, 10. * evaluate.CFID, 100. * evaluate.CFID) 
        mvals = (evaluate.MFID, 10. * evaluate.MFID, 100. * evaluate.MFID)

        qabslfile = os.path.join(
            evaluate.STORAGE_DIR,
            "testing_qabsl_"+obs_type+"_NOISE_SIGMA"+("%.2f" % noise_sigma)+"_subbins_"+
            "noweight_mc_N"+str(NTEST)+".npy")
        qfracfile = os.path.join(
            evaluate.STORAGE_DIR,
            "testing_qfrac_"+obs_type+"_NOISE_SIGMA"+("%.2f" % noise_sigma)+"_subbins_"+
            "noweight_mc_N"+str(NTEST)+".npy")
        qsqrdfile = os.path.join(
            evaluate.STORAGE_DIR,
            "testing_qsqrd_"+obs_type+"_NOISE_SIGMA"+("%.2f" % noise_sigma)+"_subbins_"+
            "noweight_mc_N"+str(NTEST)+".npy")
        qbymcfile = os.path.join(
            evaluate.STORAGE_DIR,
            "testing_qbymc_"+obs_type+"_NOISE_SIGMA"+("%.2f" % noise_sigma)+"_subbins_"+
            "noweight_mc_N"+str(NTEST)+".npy")
        mcfile = os.path.join(
            evaluate.STORAGE_DIR,
            "testing_mcresults_"+obs_type+"_NOISE_SIGMA"+("%.2f" % noise_sigma)+"_subbins_"+
            "noweight_mc_N"+str(NTEST)+".pkl")
        # Load these data
        print "Loading "+qabslfile
        qabsl[obs_type] = np.load(qabslfile)
        print "Loading "+qfracfile
        qfrac[obs_type] = np.load(qfracfile)
        print "Loading "+qabslfile
        qsqrd[obs_type] = np.load(qsqrdfile)
        print "Loading "+qbymcfile
        qbymc[obs_type] = np.load(qbymcfile)
        print "Loading "+mcfile
        with open(mcfile, "rb") as fin: mcresults[obs_type] = cPickle.load(fin)
        # Normalize data which needs it to Q=1000 for unbiased case in *space* (leave tricky qbymc
        # to later)
        qfrac[obs_type] *= 1000. / qfrac["space"][:, 0, 0].mean()
        qsqrd[obs_type] *= 1000. / qsqrd["space"][:, 0, 0].mean()
        # QABSL results first...
        # Let's do the by c then m
        plt.clf()
        plt.subplot(2,1,1)
        for j, mval in enumerate(mvals):

            plt.errorbar(
                np.asarray(cvals) * (1.01)**j,
                qabsl[obs_type][:, :, j].mean(axis=0),
                yerr=qabsl[obs_type][:, :, j].std(axis=0),
                label=r'm$_i$ = %3.e' % mval)

        plt.xscale('log')
        plt.xlim(1.e-1, 1.e-5)
        plt.ylim(0., 1400.)
        plt.xlabel(r'Input c$_i$ = c$_1$ = c$_2$')
        plt.ylabel(r'Q$_V$ (abs diffs)')
        plt.title(
            r'Test results for Q$_V$ with absolute differences ('+obs_type+', NTEST='+
            str(NTEST)+')')
        plt.legend()
        plt.subplot(2,1,2)
        for i, cval in enumerate(cvals):

            plt.errorbar(
                np.asarray(mvals) * (1.01)**j,
                qabsl[obs_type][:, i, :].mean(axis=0),
                yerr=qabsl[obs_type][:, i, :].std(axis=0),
                label=r'c$_i$ = %3.e' % cval)

        plt.xscale('log')
        plt.xlim(1.e0, 1.e-4)
        plt.ylim(0., 1400.)
        plt.xlabel(r'Input m$_i$ = m$_1$ = m$_2$')
        plt.ylabel(r'Q$_V$ (abs diffs)')
        plt.legend()
        qabslplotfile = os.path.join("plots", os.path.split(qabslfile)[-1].rsplit('.npy')[0]+'.png')
        print "Saving plot to "+qabslplotfile
        plt.savefig(qabslplotfile)

        # QFRAC results next...
        # Let's do the by c then m
        plt.clf()
        plt.subplot(2,1,1)
        for j, mval in enumerate(mvals):

            plt.errorbar(
                np.asarray(cvals) * (1.01)**j,
                qfrac[obs_type][:, :, j].mean(axis=0),
                yerr=qfrac[obs_type][:, :, j].std(axis=0),
                label=r'm$_i$ = %3.e' % mval)

        plt.xscale('log')
        plt.xlim(1.e-1, 1.e-5)
        plt.ylim(0., 1400.)
        plt.xlabel(r'Input c$_i$ = c$_1$ = c$_2$')
        plt.ylabel(r'Q$_V$ (abs frac diffs)')
        plt.title(
            r'Test results for Q$_V$ with absolute fractional differences ('+obs_type+', NTEST='+
            str(NTEST)+')')
        plt.legend()
        plt.subplot(2,1,2)
        for i, cval in enumerate(cvals):

            plt.errorbar(
                np.asarray(mvals) * (1.01)**j,
                qfrac[obs_type][:, i, :].mean(axis=0),
                yerr=qfrac[obs_type][:, i, :].std(axis=0),
                label=r'c$_i$ = %3.e' % cval)

        plt.xscale('log')
        plt.xlim(1.e0, 1.e-4)
        plt.ylim(0., 1400.)
        plt.xlabel(r'Input m$_i$ = m$_1$ = m$_2$')
        plt.ylabel(r'Q$_V$ (abs frac diffs)')
        plt.legend()
        qfracplotfile = os.path.join("plots", os.path.split(qfracfile)[-1].rsplit('.npy')[0]+'.png')
        print "Saving plot to "+qfracplotfile
        plt.savefig(qfracplotfile)

        # QSQRD results next...
        # Let's do the by c then m
        plt.clf()
        plt.subplot(2,1,1)
        for j, mval in enumerate(mvals):

            plt.errorbar(
                np.asarray(cvals) * (1.01)**j,
                qsqrd[obs_type][:, :, j].mean(axis=0),
                yerr=qsqrd[obs_type][:, :, j].std(axis=0),
                label=r'm$_i$ = %3.e' % mval)

        plt.xscale('log')
        plt.xlim(1.e-1, 1.e-5)
        plt.ylim(0., 1400.)
        plt.xlabel(r'Input c$_i$ = c$_1$ = c$_2$')
        plt.ylabel(r'Q$_V$ (squared diffs)')
        plt.title(
            r'Test results for Q$_V$ with squared differences ('+obs_type+', NTEST='+
            str(NTEST)+')')
        plt.legend()
        plt.subplot(2,1,2)
        for i, cval in enumerate(cvals):

            plt.errorbar(
                np.asarray(mvals) * (1.01)**j,
                qsqrd[obs_type][:, i, :].mean(axis=0),
                yerr=qsqrd[obs_type][:, i, :].std(axis=0),
                label=r'c$_i$ = %3.e' % cval)

        plt.xscale('log')
        plt.xlim(1.e0, 1.e-4)
        plt.ylim(0., 1400.)
        plt.xlabel(r'Input m$_i$ = m$_1$ = m$_2$')
        plt.ylabel(r'Q$_V$ (squared diffs)')
        plt.legend()
        qsqrdplotfile = os.path.join("plots", os.path.split(qsqrdfile)[-1].rsplit('.npy')[0]+'.png')
        print "Saving plot to "+qsqrdplotfile
        plt.savefig(qsqrdplotfile)

        # Then do the mc histograms
        plt.clf()
        for j, mval in enumerate(mvals):

            plt.subplot(3, 1, j + 1)
            if j==0: plt.title("Error on best-fitting c parameters for "+obs_type+" data")
            for i, cval in enumerate(cvals):

                plt.hist(
                    mcresults[obs_type]["cerr"][:, i, j], bins=np.logspace(-5, 2, 140),
                    label="c=%.1e, m=%.1e" % (cval, mvals[j]))

            plt.xscale('log')
            if j == 2: plt.xlabel(r"$\sigma_c$", fontsize="x-large")
            plt.ylabel("Counts")
            plt.legend()

        cerrhistfile = os.path.join("plots", "testing_cerr_hists_"+obs_type+".png")
        print "Saving plot to "+cerrhistfile 
        plt.savefig(cerrhistfile)

        plt.clf()
        for j, mval in enumerate(mvals):

            plt.subplot(3, 1, j + 1)
            if j==0: plt.title("Error on best-fitting m parameters for "+obs_type+" data")
            for i, cval in enumerate(cvals):

                plt.hist(
                    mcresults[obs_type]["merr"][:, i, j], bins=np.logspace(-3, 3, 120),
                    label="c=%.1e, m=%.1e" % (cval, mvals[j]))

            plt.xscale('log')
            if j == 2: plt.xlabel(r"$\sigma_m$", fontsize="x-large")
            plt.ylabel("Counts")
            plt.legend()

        merrhistfile = os.path.join("plots", "testing_merr_hists_"+obs_type+".png")
        print "Saving plot to "+merrhistfile 
        plt.savefig(merrhistfile)

        plt.clf()
        for j, mval in enumerate(mvals):

            plt.subplot(3, 1, j + 1)
            if j==0: plt.title("Best-fitting c parameters for "+obs_type+" data")
            for i, cval in enumerate(cvals):

                plt.hist(
                    mcresults[obs_type]["c"][:, i, j], bins=370, range=(-0.002, 0.058),
                    label="c=%.1e, m=%.1e" % (cval, mvals[j]))
                plt.axvline(cval, ls='--', color=["b", "g", "r"][i])

            plt.xlim(-0.002, 0.06)
            if j == 2: plt.xlabel("Best-fitting c (dashed lines show input)")
            plt.ylabel("Counts")
            plt.legend()

        chistfile = os.path.join("plots", "testing_c_hists_"+obs_type+".png")
        print "Saving plot to "+chistfile 
        plt.savefig(chistfile)

        plt.clf()
        for j, mval in enumerate(mvals):

            plt.subplot(3, 1, j + 1)
            if j==0: plt.title("Best-fitting m parameters for "+obs_type+" data")
            for i, cval in enumerate(cvals):

                plt.hist(
                    mcresults[obs_type]["m"][:, i, j], bins=101, range=(-1, 1),
                    label="c=%.1e, m=%.1e" % (cval, mval))

            plt.xlim(-1, 1)
            plt.axvline(mval, ls='--', color='k', lw=2)
            if j == 2: plt.xlabel("Best-fitting m (black dashed line shows input)")
            plt.ylabel("Counts")
            plt.legend(loc=2)

        mhistfile = os.path.join("plots", "testing_m_hists_"+obs_type+".png")
        print "Saving plot to "+mhistfile 
        plt.savefig(mhistfile)


        
