#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
from matplotlib import pyplot as plt
from matplotlib import cm

import pickle
import netCDF4 as nc

import climtools_lib as ctl
import climdiags as cd

from matplotlib.colors import LogNorm
from datetime import datetime

from scipy import stats
import pandas as pd

plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 28
plt.rcParams['axes.labelsize'] = 22

#############################################################################
cart_out_orig = '/home/fabiano/Research/lavori/CMIP6/Check_ece_vs_mpi/FINAL/'
ctl.mkdir(cart_out_orig)

cart_in = '/data-hobbes/fabiano/WR_CMIP6/'
file_hist_tot = cart_in + 'out_NEW_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refCLUS_dtr_light.p'
#
# fil_ece_hist = cart_in + 'out_eceens_hist_NDJFM_{}_4clus_4pcs_1964-2014_refCLUS_dtr.p'
# fil_ece = cart_in + 'out_eceens_ssp585_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr_reb.p'
# fil_ece_r4 = cart_in + 'out_eceens_ssp585_r4_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr_reb.p'
# fil_mpi = cart_in + 'mpiens_ssp585/out_mpiens_ssp585_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr_reb.p'

fil_hist = cart_in + '{}ens_hist_rbtot/out_{}ens_hist_rbtot_NDJFM_{}_4clus_4pcs_1964-2014_refCLUS_dtr_reb.p'
fil_ssp = cart_in + '{}ens_ssp585_rbtot/out_{}ens_ssp585_rbtot_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr_reb.p'

tips = 'ece mpi uk'.split()
tips = ['mpi', 'uk']

ssp = 'ssp585'
numclus = 4
reg_names_area = dict()
reg_names_area['EAT'] = ['NAO+', 'SBL', 'AR', 'NAO-']
reg_names_area['PNA'] = ['PT', 'PNA+', 'PNA-', 'BR']

clatlo = dict()
clatlo['EAT'] = (70., -20.)
clatlo['PNA'] = (70., -120.)

#allssps = 'ssp119 ssp126 ssp245 ssp370 ssp585'.split()
#allssps = 'ssp126 ssp245 ssp370 ssp585'.split()
allssps = ['ssp585']

#oknam = ['EC-Earth3_r1i1p1f1', 'EC-Earth3_r4i1p1f1', 'MPI-ESM1-2-LR_r1i1p1f1']

ttests = dict()
for area in ['EAT', 'PNA']:
    cart_out = cart_out_orig + '{}_NDJFM/'.format(area)
    ctl.mkdir(cart_out)

    results_hist, results_ref = ctl.load_wrtool(file_hist_tot.format(area))
    histme = np.mean([results_hist[ke]['freq_clus'] for ke in results_hist.keys()], axis = 0)
    del results_hist

    rescoso = dict()
    for tip in tips:
        rescoso[(tip, 'hist')], _ = ctl.load_wrtool(fil_hist.format(tip, tip, area))
        rescoso[(tip, 'ssp585')], _ = ctl.load_wrtool(fil_ssp.format(tip, tip, area))

    # appiccico hist e ssp
    #for cosone, histmem in zip([ece_ssp, ece_ssp_r4, mpi_ssp], ['EC-Earth3_r1i1p1f1', 'EC-Earth3_r4i1p1f1', 'MPI-ESM1-2-LR_r1i1p1f1']):
    resdict = dict()
    for tip in tips:
        cosone = dict()
        print(rescoso[(tip, 'ssp585')].keys())
        for mem in rescoso[(tip, 'ssp585')].keys():
            cosone[mem] = dict()
            labs = np.concatenate([rescoso[(tip, 'hist')][mem]['labels'], rescoso[(tip, 'ssp585')][mem]['labels']])
            dats = np.concatenate([rescoso[(tip, 'hist')][mem]['dates'], rescoso[(tip, 'ssp585')][mem]['dates']])
            pcs = np.concatenate([rescoso[(tip, 'hist')][mem]['pcs'], rescoso[(tip, 'ssp585')][mem]['pcs']], axis = 0)

            cosone[mem]['labels'] = labs
            cosone[mem]['dates'] = dats
            cosone[mem]['pcs'] = pcs

        resdict[tip] = cosone

    alltips = tuple(resdict.keys())
    colorz = ctl.color_set(len(alltips))

    runfreq = dict()
    seasfreq = dict()

    figs = []
    for tip, col in zip(alltips, colorz):
        okmods = resdict[tip].keys()
        fig = plt.figure(figsize = (16,12))
        for reg in range(4):
            ax = fig.add_subplot(2, 2, reg+1)
            cosi = []
            for mem in okmods:
                seasfr, yr = ctl.calc_seasonal_clus_freq(resdict[tip][mem]['labels'], resdict[tip][mem]['dates'], numclus)
                seasfreq[(tip, mem, reg)] = seasfr[reg, :]
                seas1 = np.array(ctl.running_mean(seasfr[reg, :], 20))
                runfreq[(tip, mem, reg)] = seas1
                ax.plot(yr, seas1, label = tip, color = col)
            ax.set_title(reg_names_area[area][reg])
            #ax.legend()
        fig.suptitle(tip)
        figs.append(fig)

    ctl.plot_pdfpages(cart_out + 'check_single_members_{}.pdf'.format(area), figs)

    fig = plt.figure(figsize = (16,12))
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg+1)
        for tip, col in zip(alltips, colorz):
            okmods = resdict[tip].keys()
            cosi = [runfreq[(tip, mem, reg)] for mem in okmods]
            coso = np.mean(cosi, axis = 0)
            runfreq[(tip, reg)] = coso
            coserr = np.std(cosi, axis = 0)/np.sqrt(len(okmods)-1)
            ax.fill_between(yr, coso-coserr, coso+coserr, color = col, alpha = 0.3)
            ax.plot(yr, coso, label = tip, color = col)

        ax.set_title(reg_names_area[area][reg])
        #ax.legend()

    ctl.custom_legend(fig, colorz, ['MPI-ESM1-2-LR', 'UKESM1-0-LL'], ncol = 2, add_space_below = 0.12)
    fig.savefig(cart_out + 'check_ece_vs_mpi_{}.pdf'.format(area))

    trend_ssp = dict()

    for tip in alltips:
        okmods = resdict[tip].keys()
        for reg in range(4):
            for mem in okmods:
                seasfr = seasfreq[(tip, mem, reg)]
                seas10 = np.array(ctl.running_mean(seasfr, 10))

                okcose = ~(np.isnan(seas10)) & (yr >= 2015)
                m, c, err_m, err_c = ctl.linear_regre_witherr(np.array(yr[okcose]), np.array(seas10[okcose]))
                trend_ssp[(mem, tip, reg)] = m
                trend_ssp[(mem, tip, reg, 'err')] = err_m

            trend_ssp[(tip, reg)] = np.array([trend_ssp[(mem, tip, reg)] for mem in okmods])
            trend_ssp[(tip, reg, 'err')] = np.array([trend_ssp[(mem, tip, reg, 'err')] for mem in okmods])

    freqs = dict()
    for tip in alltips:
        okmods = resdict[tip].keys()
        for mem in okmods:
            dat1 = pd.Timestamp('09-01-1964').to_pydatetime()
            dat2 = pd.Timestamp('04-01-2014').to_pydatetime()
            labs, dats = ctl.sel_time_range(resdict[tip][mem]['labels'], resdict[tip][mem]['dates'], (dat1, dat2))
            freqs[('hist', mem, tip)] = ctl.calc_clus_freq(labs, numclus)

            dat1 = pd.Timestamp('09-01-2050').to_pydatetime()
            dat2 = pd.Timestamp('04-01-2100').to_pydatetime()
            labs, dats = ctl.sel_time_range(resdict[tip][mem]['labels'], resdict[tip][mem]['dates'], (dat1, dat2))

            # restim, _, _ = ctl.calc_regime_residtimes(labs, dats)
            freqs[(ssp, mem, tip)] = ctl.calc_clus_freq(labs, numclus)

        freqs[('hist', tip)] = np.stack([freqs[('hist', mem, tip)] for mem in okmods])
        freqs[(ssp, tip)] = np.stack([freqs[(ssp, mem, tip)] for mem in okmods])

    pickle.dump([seasfreq, runfreq, freqs, trend_ssp], open(cart_out + 'check_refit_{}.p'.format(area), 'wb'))

    reg_names = reg_names_area[area]

    figall = plt.figure(figsize = (16,12))
    axes = []
    for reg in range(4):
        ax = figall.add_subplot(2, 2, reg + 1)
        axes.append(ax)

        allpercs = dict()
        for nu in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(nu)] = [np.percentile(freqs[(ssp, tip)][:, reg], nu) - histme[reg] for tip in alltips]

        ctl.boxplot_on_ax(ax, allpercs, alltips, colorz, plot_mean = False, plot_ensmeans = False, plot_minmax = False)
        ax.axhline(0, color = 'gray', linewidth = 0.5)
        ax.set_xticks([])
        ax.set_title(reg_names[reg])

        if reg == 0 or reg == 2: ax.set_ylabel('Regime frequency anomaly')

    ctl.adjust_ax_scale(axes)

    #ctl.custom_legend(figall, colorz, alltips, ncol = 2)
    ctl.custom_legend(figall, colorz, ['MPI-ESM1-2-LR', 'UKESM1-0-LL'], ncol = 2, add_space_below = 0.1)
    figall.savefig(cart_out + 'check_ece_vs_mpi_WRfreq_{}.pdf'.format(area))


    figall = plt.figure(figsize = (16,12))
    axes = []
    for reg in range(4):
        ax = figall.add_subplot(2, 2, reg + 1)

        axes.append(ax)
        histmean = dict()
        for tip in alltips:
            histmean[tip] = np.mean(freqs[('hist', tip)][:, reg])

        allpercs = dict()
        for nu in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(nu)] = [np.percentile(freqs[(ssp, tip)][:, reg]-freqs[('hist', tip)][:, reg], nu) for tip in alltips]

        ctl.boxplot_on_ax(ax, allpercs, alltips, colorz, plot_mean = False, plot_ensmeans = False, plot_minmax = False)
        ax.axhline(0, color = 'gray', linewidth = 0.5)
        ax.set_xticks([])
        ax.set_title(reg_names[reg])

        if reg == 0 or reg == 2: ax.set_ylabel('Regime frequency anomaly')

    ctl.adjust_ax_scale(axes)

    ctl.custom_legend(figall, colorz, alltips, ncol = 2)
    figall.savefig(cart_out + 'check_ece_vs_mpi_WRfreqREL_{}.pdf'.format(area))


    figall = plt.figure(figsize = (16,12))
    axes = []
    for reg in range(4):
        ax = figall.add_subplot(2, 2, reg + 1)

        allpercs = dict()
        for nu in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(nu)] = [np.percentile(trend_ssp[(tip, reg)], nu) for tip in alltips]

        ctl.boxplot_on_ax(ax, allpercs, alltips, colorz, plot_mean = False, plot_ensmeans = False, plot_minmax = False)

        ax.axhline(0, color = 'gray', linewidth = 0.5)
        ax.set_xticks([])
        if reg == 0 or reg == 2: ax.set_ylabel('Trend in regime frequency (1/yr)')
        axes.append(ax)

    ctl.adjust_ax_scale(axes)

    ctl.custom_legend(figall, colorz, alltips, ncol = 3)
    figall.savefig(cart_out + 'check_ece_vs_mpi_trends_{}.pdf'.format(area))


    figall = plt.figure(figsize = (16,12))
    axes = []
    for reg in range(4):
        ax = figall.add_subplot(2, 2, reg + 1)
        axes.append(ax)

        histmean = dict()
        for tip in alltips:
            histmean[tip] = np.mean(freqs[('hist', tip)][:, reg])

        data = [freqs[(ssp, tip)][:, reg]-histmean[tip] for tip in alltips]

        parts = ax.violinplot(data, positions = None, widths=0.4, showmeans=True, showextrema=True, showmedians=True, bw_method=0.5)#, quantiles=[0.25, 0.75]
        for pc, col in zip(parts['bodies'], colorz):
            pc.set_facecolor(col)
            pc.set_edgecolor(col)
            pc.set_alpha(0.5)

        ax.axhline(0, color = 'gray', linewidth = 0.5)
        ax.set_xticks([])
        ax.set_title(reg_names[reg])

        if reg == 0 or reg == 2: ax.set_ylabel('Regime frequency anomaly')

    ctl.adjust_ax_scale(axes)

    ctl.custom_legend(figall, colorz, alltips, ncol = 2)
    figall.savefig(cart_out + 'check_ece_vs_mpi_WRfreq_violin_{}.pdf'.format(area))
