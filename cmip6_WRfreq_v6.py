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
cart_cmip5 = '/home/fabiano/Research/lavori/CMIP6/Results_cmip5/{}_NDJFM/'

yr10 = 10 # length of running mean
cart_in = '/data-hobbes/fabiano/WR_CMIP6/'

tip = 'r1_rebase'
if tip == 'r1_rebase':
    cart_out_orig = '/home/fabiano/Research/lavori/CMIP6/Results_v6_rebase/'
    cart_v5 = '/home/fabiano/Research/lavori/CMIP6/Results_v5_rebase/{}_NDJFM/'
    file_hist = cart_in + 'out_NEW_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refCLUS_dtr_light.p'
elif tip == 'ensrebase':
    cart_out_orig = '/home/fabiano/Research/lavori/CMIP6/Results_v6_ensrebase/'
    cart_v5 = '/home/fabiano/Research/lavori/CMIP6/Results_v5_ensrebase/{}_NDJFM/'
    file_hist = cart_in + 'out_cmip6_ensrebase_hist_NDJFM_{}_4clus_4pcs_1964-2014_refCLUS_dtr_reb.p'

ctl.mkdir(cart_out_orig)

numclus = 4
reg_names_area = dict()
reg_names_area['EAT'] = ['NAO+', 'SBL', 'AR', 'NAO-']
reg_names_area['PNA'] = ['PT', 'PNA+', 'PNA-', 'BR']

clatlo = dict()
clatlo['EAT'] = (70., -20.)
clatlo['PNA'] = (70., -120.)

#allssps = 'ssp119 ssp126 ssp245 ssp370 ssp585'.split()
allssps = 'ssp126 ssp245 ssp370 ssp585'.split()

ttests = dict()
cosette = dict()

area = 'EAT'
for area in ['EAT', 'PNA']:
    cart_out = cart_out_orig + '{}_NDJFM/'.format(area)
    ctl.mkdir(cart_out)
    reg_names = reg_names_area[area]

    results_hist, results_ref = ctl.load_wrtool(file_hist.format(area))

    trend_ssp, residtime_ssp = pickle.load(open(cart_v5.format(area) + 'trends_wrfreq_e_restime_{}.p'.format(area), 'rb'))
    seasfreq, runfreq = pickle.load(open(cart_v5.format(area) + 'seasfreqs_{}_v4.p'.format(area), 'rb'))
    freqs, residtimes, patterns, num_event = pickle.load(open(cart_v5.format(area) + 'allresults_dicts_{}_v3.p'.format(area), 'rb'))

    freqs_cmip5, trend_ssp_cmip5, residtimes_cmip5, num_event_cmip5 = pickle.load(open(cart_cmip5.format(area) + 'freqs_cmip5_{}.p'.format(area), 'rb'))
    seasfreq_cmip5, runfreq_cmip5 = pickle.load(open(cart_cmip5.format(area) + 'seasfreqs_cmip5_{}.p'.format(area), 'rb'))
    freqs.update(freqs_cmip5)
    trend_ssp.update(trend_ssp_cmip5)
    residtimes.update(residtimes_cmip5)
    seasfreq.update(seasfreq_cmip5)
    runfreq.update(runfreq_cmip5)
    num_event.update(num_event_cmip5)

    allssps = ['ssp126', 'ssp245', 'ssp370', 'ssp585', 'rcp85_cmip5']
    allsims = ['hist', 'ssp126', 'ssp245', 'ssp370', 'ssp585', 'rcp85_cmip5']
    allsims_wcmip5 = ['hist', 'ssp126', 'ssp245', 'ssp370', 'ssp585', 'hist_cmip5', 'rcp85_cmip5']
    allsimcol = ['hist', 'ssp126', 'ssp245', 'hist_cmip5', 'ssp370', 'ssp585', 'rcp85_cmip5']

    colorzi = ctl.color_set(7)
    allsimcol = ['ssp370', 'ssp245', 'ssp126', 'bau', 'rcp85_cmip5', 'ssp585', 'bau2']
    coldic = dict(zip(allsimcol, colorzi))
    coldic['hist'] = 'darkslategray'
    coldic['hist_cmip5'] = 'saddlebrown'
    colsim = [coldic[ssp] for ssp in allsims]
    colssp = [coldic[ssp] for ssp in allssps]
    colsim_wcmip5 = [coldic[ssp] for ssp in allsims_wcmip5]

    yr = np.arange(1965, 2100)

    fig = plt.figure(figsize = (16,12))
    axes = []
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg+1)
        for ssp in allssps:
            col = coldic[ssp]
            coso = runfreq[(ssp, 'run20', reg)]-np.mean(freqs[('hist', 'all', 'tot50')][:, reg])
            coserr = runfreq[(ssp, 'run20err', reg)]
            ax.fill_between(yr, coso-coserr, coso+coserr, color = col, alpha = 0.15)
            ax.plot(yr, coso, label = ssp, color = col, linewidth = 2)
        ax.set_title(reg_names_area[area][reg])
        ax.axvline(2015, color = 'lightslategray', linewidth = 0.2, linestyle = '--')
        ax.axhline(0., color = 'lightslategray', linewidth = 0.2)
        axes.append(ax)

    ctl.adjust_ax_scale(axes)
    ctl.custom_legend(fig, colssp, allssps, ncol = 3)
    fig.savefig(cart_out + 'allssps_freq20_{}_anom_wcmip5.pdf'.format(area))

    fig = plt.figure(figsize = (27,6))
    axes = []
    for reg in range(4):
        ax = fig.add_subplot(1, 4, reg+1)
        for ssp in allssps:
            col = coldic[ssp]
            coso = runfreq[(ssp, 'run20', reg)]-np.mean(freqs[('hist', 'all', 'tot50')][:, reg])
            coserr = runfreq[(ssp, 'run20err', reg)]
            ax.fill_between(yr, coso-coserr, coso+coserr, color = col, alpha = 0.15)
            ax.plot(yr, coso, label = ssp, color = col, linewidth = 2)
        ax.set_title(reg_names_area[area][reg])
        ax.axvline(2015, color = 'lightslategray', linewidth = 0.2, linestyle = '--')
        ax.axhline(0., color = 'lightslategray', linewidth = 0.2)
        axes.append(ax)
        if reg == 0:
            ax.set_ylabel('Frequency anomaly')

    areana = area
    if area == 'PNA': areana = 'PAC'
    axes[0].text(0.05, 0.5, areana, horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 35)

    ctl.adjust_ax_scale(axes)
    ctl.custom_legend(fig, colssp, allssps, ncol = 3, add_space_below = 0.15)
    fig.savefig(cart_out + 'allssps_freq20_{}_anom_wcmip5_line.pdf'.format(area))

    figall = plt.figure(figsize = (28,12))
    axes = []
    cos = 'tot50'
    for reg in range(4):
        ax = figall.add_subplot(2, 4, reg + 1)
        axes.append(ax)

        allpercs = dict()
        for nu in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(nu)] = [np.percentile(freqs[(ssp, 'all', cos)][:, reg], nu) for ssp in allsims]

        ax.axhline(allpercs['p50'][0], color = 'gray', linewidth = 0.5)

        positions = list(np.arange(len(allsims)-1)*0.7)
        positions.append(positions[-1]+0.3+0.7)
        ctl.boxplot_on_ax(ax, allpercs, allsims, colsim, plot_mean = False, plot_ensmeans = False, plot_minmax = False, positions = positions)
        ax.set_xticks([])
        ax.set_title(reg_names[reg])
        if reg == 0: ax.set_ylabel('Regime frequency')

        ax.scatter(0, results_ref['freq_clus'][reg], color = 'black', marker = '*', s = 5)
        ax.axvline(np.mean([positions[-1], positions[-2]]), color = 'lightslategray', linewidth = 0.2, linestyle = '--')
        xli = ax.get_xlim()

    #ctl.adjust_ax_scale(axes)
    axes = []
    # qui ci metto il trend
    for reg in range(4):
        ax = figall.add_subplot(2, 4, 4 + reg + 1)
        allpercs = dict()
        for nu in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(nu)] = [np.percentile(trend_ssp[(ssp, 'all', 'trend', 'freq10', reg)], nu) for ssp in allssps]

        ctl.boxplot_on_ax(ax, allpercs, allssps, colssp, plot_mean = False, plot_ensmeans = False, plot_minmax = False, positions = positions[1:])

        ax.axhline(0, color = 'gray', linewidth = 0.5)
        ax.set_xticks([])
        if reg == 0: ax.set_ylabel('Trend in regime frequency (1/yr)')
        axes.append(ax)
        ax.set_xlim(xli)
        ax.axvline(np.mean([positions[-1], positions[-2]]), color = 'lightslategray', linewidth = 0.2, linestyle = '--')

    ctl.adjust_ax_scale(axes)

    ctl.custom_legend(figall, colsim, allsims, ncol = 3)
    figall.savefig(cart_out + 'WRfreq_allssp_{}_8box_wtrend_wcmip5.pdf'.format(area, cos))

    print('T-TEST for {}'.format(area))
    for reg in range(4):
        print('REGIME',reg)
        for ssp in allssps:
            if ssp == 'rcp85_cmip5':
                a = freqs[('hist_cmip5', 'all', cos)][:, reg]
            else:
                a = freqs[('hist', 'all', cos)][:, reg]
            print(ssp)
            b = freqs[(ssp, 'all', cos)][:, reg]
            ttests[('freq', area, reg, ssp)] = stats.ttest_ind(a, b, equal_var = False)
            print(ttests[('freq', area, reg, ssp)])

    figall = plt.figure(figsize = (16,12))
    axes = []
    cos = 'tot50'
    for reg in range(4):
        ax = figall.add_subplot(2, 2, reg + 1)
        axes.append(ax)

        histmean = np.mean(freqs[('hist', 'all', 'tot50')][:, reg])
        allpercs = dict()
        for nu in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(nu)] = [np.percentile(freqs[(ssp, 'all', cos)][:, reg], nu) - histmean for ssp in allsims]

        ax.axhline(allpercs['p50'][0], color = 'gray', linewidth = 0.5)

        positions = list(np.arange(len(allsims)-1)*0.7)
        positions.append(positions[-1]+0.3+0.7)
        ctl.boxplot_on_ax(ax, allpercs, allsims, colsim, plot_mean = False, plot_ensmeans = False, plot_minmax = False, positions = positions)
        ax.set_xticks([])
        ax.set_title(reg_names[reg])
        #ax.text(1.0, 1.0, na, horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 20)
        if reg == 0 or reg == 2: ax.set_ylabel('Regime frequency anomaly')

        #ax.scatter(0, results_ref['freq_clus'][reg], color = 'black', marker = '*', s = 5)
        for pos, ssp in zip(positions[1:], allsims[1:]):
            if ttests[('freq', area, reg, ssp)].pvalue < 0.05:
                ax.scatter(pos, -10, color = 'black', marker = '*', s = 30)

        ax.axvline(np.mean([positions[-1], positions[-2]]), color = 'lightslategray', linewidth = 0.2, linestyle = '--')
        xli = ax.get_xlim()

    ctl.adjust_ax_scale(axes)

    ctl.custom_legend(figall, colsim, allsims, ncol = 3)
    figall.savefig(cart_out + 'WRfreq_{}_{}_FINAL.pdf'.format(area, cos))

    figall = plt.figure(figsize = (16,12))
    axes = []
    cos = 'tot50'
    for reg in range(4):
        ax = figall.add_subplot(2, 2, reg + 1)
        axes.append(ax)

        histmean = np.mean(freqs[('hist', 'all', 'tot50')][:, reg])
        ax.axhline(0, color = 'gray', linewidth = 0.5)

        positions = list(np.arange(len(allsims)-1)*0.7)
        positions.append(positions[-1]+0.3+0.7)

        data = [freqs[(ssp, 'all', cos)][:, reg] - histmean for ssp in allsims]

        parts = ax.violinplot(data, positions = positions, widths=0.4, showmeans=False, showextrema=True, showmedians=True)#, quantiles=[0.25, 0.75]
        for pc, col in zip(parts['bodies'], colsim):
            pc.set_facecolor(col)
            pc.set_edgecolor(col)
            pc.set_alpha(0.5)

        ax.set_xticks([])
        ax.set_title(reg_names[reg])
        #ax.text(1.0, 1.0, na, horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 20)
        if reg == 0 or reg == 2: ax.set_ylabel('Regime frequency anomaly')

        #ax.scatter(0, results_ref['freq_clus'][reg], color = 'black', marker = '*', s = 5)
        for pos, ssp in zip(positions[1:], allsims[1:]):
            if ttests[('freq', area, reg, ssp)].pvalue < 0.05:
                ax.scatter(pos, -10, color = 'black', marker = '*', s = 30)

        ax.axvline(np.mean([positions[-1], positions[-2]]), color = 'lightslategray', linewidth = 0.2, linestyle = '--')
        xli = ax.get_xlim()

    ctl.adjust_ax_scale(axes)

    ctl.custom_legend(figall, colsim, allsims, ncol = 3)
    figall.savefig(cart_out + 'WRfreq_{}_{}_violin.pdf'.format(area, cos))


    figall = plt.figure(figsize = (16,12))
    axes = []
    cos = 'tot50'
    for reg in range(4):
        ax = figall.add_subplot(2, 2, reg + 1)
        axes.append(ax)

        histmean = np.mean(freqs[('hist', 'all', 'tot50')][:, reg])
        allpercs = dict()
        for nu in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(nu)] = [np.percentile(freqs[(ssp, 'all', cos)][:, reg], nu) - histmean for ssp in allsims_wcmip5]

        ax.axhline(allpercs['p50'][0], color = 'gray', linewidth = 0.5)

        positions = list(np.arange(len(allsims)-1)*0.7)
        positions.append(positions[-1]+0.3+0.7)
        positions.append(positions[-1]+0.7)
        ctl.boxplot_on_ax(ax, allpercs, allsims_wcmip5, colsim_wcmip5, plot_mean = False, plot_ensmeans = False, plot_minmax = False, positions = positions)
        ax.set_xticks([])
        ax.set_title(reg_names[reg])
        #ax.text(1.0, 1.0, na, horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 20)
        if reg == 0 or reg == 2: ax.set_ylabel('Regime frequency anomaly')
        okpos = [pos for pos, ssp in zip(positions, allsims_wcmip5) if 'hist' not in ssp]
        okssps = [ssp for ssp in allsims_wcmip5 if 'hist' not in ssp]
        for pos, ssp in zip(okpos, okssps):
            if ttests[('freq', area, reg, ssp)].pvalue < 0.05:
                ax.scatter(pos, -10, color = 'black', marker = '*', s = 30)

        ax.axvline(np.mean([positions[-3], positions[-2]]), color = 'lightslategray', linewidth = 0.2, linestyle = '--')
        xli = ax.get_xlim()

    ctl.adjust_ax_scale(axes)

    ctl.custom_legend(figall, colsim_wcmip5, allsims_wcmip5, ncol = 3)
    figall.savefig(cart_out + 'WRfreq_{}_{}_FINAL_whistcmip5.pdf'.format(area, cos))

    positions = list(np.arange(len(allsims)-1)*0.7)
    positions.append(positions[-1]+0.3+0.7)

    figall = plt.figure(figsize = (16,12))
    axes = []
    for reg in range(4):
        ax = figall.add_subplot(2, 2, reg + 1)
        allpercs = dict()
        for nu in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(nu)] = [np.percentile(trend_ssp[(ssp, 'all', 'trend', 'freq10', reg)], nu) for ssp in allssps]

        ctl.boxplot_on_ax(ax, allpercs, allssps, colssp, plot_mean = False, plot_ensmeans = False, plot_minmax = False, positions = positions[1:])

        ax.axhline(0, color = 'gray', linewidth = 0.5)
        ax.set_xticks([])
        if reg == 0 or reg == 2: ax.set_ylabel('Trend in regime frequency (1/yr)')
        axes.append(ax)
        ax.set_xlim(xli)
        ax.axvline(np.mean([positions[-1], positions[-2]]), color = 'lightslategray', linewidth = 0.2, linestyle = '--')

    ctl.adjust_ax_scale(axes)

    ctl.custom_legend(figall, colsim, allsims, ncol = 3)
    figall.savefig(cart_out + 'Trends_{}_FINAL.pdf'.format(area))

    cos = 'mean'
    print('T-TEST for {}'.format(area))
    for reg in range(4):
        print('REGIME',reg)
        for ssp in allssps:
            if ssp == 'rcp85_cmip5':
                a = residtimes[('hist_cmip5', 'all', cos, reg)]
            else:
                a = residtimes[('hist', 'all', cos, reg)]
            print(ssp)
            b = residtimes[(ssp, 'all', cos, reg)]
            ttests[('residtimes', area, reg, ssp)] = stats.ttest_ind(a, b, equal_var = False)
            print(ttests[('residtimes', area, reg, ssp)])

    for reg in range(4):
        print('REGIME',reg)
        for ssp in allssps:
            if ssp == 'rcp85_cmip5':
                a = num_event[('hist_cmip5', 'all', reg)]
            else:
                a = num_event[('hist', 'all', reg)]
            print(ssp)
            b = num_event[(ssp, 'all', reg)]
            ttests[('num_event', area, reg, ssp)] = stats.ttest_ind(a, b, equal_var = False)

            print(ttests[('num_event', area, reg, ssp)])

    fig = plt.figure(figsize = (16,12))
    axes = []
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg + 1)
        axes.append(ax)

        allpercs = dict()

        for nu in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(nu)] = [np.percentile(residtimes[(ssp, 'all', cos, reg)], nu) for ssp in allsims]

        ax.axhline(allpercs['p50'][0], color = 'gray', linewidth = 0.5)
        ctl.boxplot_on_ax(ax, allpercs, allsims, colsim, plot_mean = False, plot_ensmeans = False, plot_minmax = False, positions = positions)
        # ax.axhline(0, color = 'gray', linewidth = 0.5)
        ax.set_xticks([])
        ax.set_title(reg_names[reg])
        ax.axvline(np.mean([positions[-1], positions[-2]]), color = 'lightslategray', linewidth = 0.2, linestyle = '--')
        for pos, ssp in zip(positions[1:], allsims[1:]):
            if ttests[('residtimes', area, reg, ssp)].pvalue < 0.05:
                ax.scatter(pos, 3.5, color = 'black', marker = '*', s = 30)

        if reg == 0 or reg == 2:
            ax.set_ylabel('Av. persistence (days)')

    ctl.adjust_ax_scale(axes)
    ctl.custom_legend(fig, colsim, allsims, ncol = 3)
    fig.savefig(cart_out + 'Restime_allssp_{}_{}_wcmip5.pdf'.format(area, cos))


    fig = plt.figure(figsize = (16,12))
    axes = []
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg + 1)
        axes.append(ax)

        allpercs = dict()

        for nu in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(nu)] = [np.percentile(num_event[(ssp, 'all', reg)], nu) for ssp in allsims]

        ax.axhline(allpercs['p50'][0], color = 'gray', linewidth = 0.5)
        ctl.boxplot_on_ax(ax, allpercs, allsims, colsim, plot_mean = False, plot_ensmeans = False, plot_minmax = False, positions = positions)
        # ax.axhline(0, color = 'gray', linewidth = 0.5)
        ax.set_xticks([])
        ax.set_title(reg_names[reg])
        ax.axvline(np.mean([positions[-1], positions[-2]]), color = 'lightslategray', linewidth = 0.2, linestyle = '--')
        for pos, ssp in zip(positions[1:], allsims[1:]):
            if ttests[('num_event', area, reg, ssp)].pvalue < 0.05:
                ax.scatter(pos, 2.7, color = 'black', marker = '*', s = 30)

        if reg == 0 or reg == 2:
            ax.set_ylabel('Num. events per 100 days')

    ctl.adjust_ax_scale(axes)
    ctl.custom_legend(fig, colsim, allsims, ncol = 3)
    fig.savefig(cart_out + 'Num_event_allssp_{}_{}_wcmip5.pdf'.format(area, cos))


    fig = plt.figure(figsize = (27,6))
    axes = []
    for reg in range(4):
        ax = fig.add_subplot(1, 4, reg + 1)
        axes.append(ax)

        allpercs = dict()

        for nu in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(nu)] = [np.percentile(residtimes[(ssp, 'all', cos, reg)], nu) for ssp in allsims]

        cosette[(area, 'residtimes', reg)] = allpercs

        ax.axhline(allpercs['p50'][0], color = 'gray', linewidth = 0.5)
        ctl.boxplot_on_ax(ax, allpercs, allsims, colsim, plot_mean = False, plot_ensmeans = False, plot_minmax = False, positions = positions)
        # ax.axhline(0, color = 'gray', linewidth = 0.5)
        ax.set_xticks([])
        ax.set_title(reg_names[reg])
        ax.axvline(np.mean([positions[-1], positions[-2]]), color = 'lightslategray', linewidth = 0.2, linestyle = '--')
        if reg == 0 or reg == 2:
            ax.set_ylabel('Num. events per 100 days')
        for pos, ssp in zip(positions[1:], allsims[1:]):
            if ttests[('residtimes', area, reg, ssp)].pvalue < 0.05:
                ax.scatter(pos, 3.5, color = 'black', marker = '*', s = 30)

    ctl.adjust_ax_scale(axes)
    ctl.custom_legend(fig, colsim, allsims, ncol = 3)
    fig.savefig(cart_out + 'Restime_allssp_{}_{}_wcmip5_line.pdf'.format(area, cos))


    fig = plt.figure(figsize = (27,6))
    axes = []
    for reg in range(4):
        ax = fig.add_subplot(1, 4, reg + 1)
        axes.append(ax)

        allpercs = dict()

        for nu in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(nu)] = [np.percentile(num_event[(ssp, 'all', reg)], nu) for ssp in allsims]

        cosette[(area, 'num_event', reg)] = allpercs

        ax.axhline(allpercs['p50'][0], color = 'gray', linewidth = 0.5)
        ctl.boxplot_on_ax(ax, allpercs, allsims, colsim, plot_mean = False, plot_ensmeans = False, plot_minmax = False, positions = positions)
        # ax.axhline(0, color = 'gray', linewidth = 0.5)
        ax.set_xticks([])
        ax.set_title(reg_names[reg])
        ax.axvline(np.mean([positions[-1], positions[-2]]), color = 'lightslategray', linewidth = 0.2, linestyle = '--')
        if reg == 0 or reg == 2:
            ax.set_ylabel('Av. persistence (days)')
        for pos, ssp in zip(positions[1:], allsims[1:]):
            if ttests[('residtimes', area, reg, ssp)].pvalue < 0.05:
                ax.scatter(pos, 3.5, color = 'black', marker = '*', s = 30)

    ctl.adjust_ax_scale(axes)
    ctl.custom_legend(fig, colsim, allsims, ncol = 3)
    fig.savefig(cart_out + 'Num_event_allssp_{}_{}_wcmip5_line.pdf'.format(area, cos))


for nometti in ['residtimes', 'num_event']:
    fig = plt.figure(figsize = (27,11))
    axes = []
    i = 0
    for area in ['EAT', 'PNA']:
        reg_names = reg_names_area[area]
        for reg in range(4):
            ax = fig.add_subplot(2, 4, 4*i + reg + 1)
            axes.append(ax)

            allpercs = cosette[(area, nometti, reg)]

            ax.axhline(allpercs['p50'][0], color = 'gray', linewidth = 0.5)
            ctl.boxplot_on_ax(ax, allpercs, allsims, colsim, plot_mean = False, plot_ensmeans = False, plot_minmax = False, positions = positions)
            # ax.axhline(0, color = 'gray', linewidth = 0.5)
            ax.set_xticks([])
            ax.set_title(reg_names[reg])
            ax.axvline(np.mean([positions[-1], positions[-2]]), color = 'lightslategray', linewidth = 0.2, linestyle = '--')
            if reg == 0:
                if nometti == 'residtimes':
                    ax.set_ylabel('Av. persistence (days)')
                    poscoso = 3.5
                else:
                    ax.set_ylabel('Num. events per 100 days')
                    poscoso = 2.7

            for pos, ssp in zip(positions[1:], allsims[1:]):
                if ttests[(nometti, area, reg, ssp)].pvalue < 0.05:
                    ax.scatter(pos, poscoso, color = 'black', marker = '*', s = 30)

        i += 1

    ctl.adjust_ax_scale(axes)
    ctl.custom_legend(fig, colsim, allsims, ncol = 3)

    ax.text(0.05, 0.75, 'EAT', horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 35)
    ax.text(0.05, 0.25, 'PAC', horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 35)

    fig.savefig(cart_out_orig + '{}_allssp_{}_wcmip5_FINAL.pdf'.format(nometti, cos))



with open(cart_out + 'ttests.p', 'wb') as fillo:
    pickle.dump(ttests, fillo)
