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

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 18

#############################################################################
cart_v5 = '/home/fabiano/Research/lavori/CMIP6/Results_v5_rebase/{}_NDJFM/'
cart_cmip5 = '/home/fabiano/Research/lavori/CMIP6/Results_cmip5/{}_NDJFM/'

yr10 = 10 # length of running mean
#dtrtyp = 'light'
dtrtyp = 'histrebase'
cart_out_orig = '/home/fabiano/Research/lavori/CMIP6/Results_v6_rebase/'
ctl.mkdir(cart_out_orig)

numclus = 4
reg_names_area = dict()
reg_names_area['EAT'] = ['NAO+', 'SBL', 'AR', 'NAO-']
reg_names_area['PNA'] = ['PT', 'PNA+', 'PNA-', 'AR']

clatlo = dict()
clatlo['EAT'] = (70., -20.)
clatlo['PNA'] = (70., -120.)

#allssps = 'ssp119 ssp126 ssp245 ssp370 ssp585'.split()
allssps = 'ssp126 ssp245 ssp370 ssp585'.split()

area = 'EAT'
for area in ['EAT']:#, 'PNA']:
    cart_out = cart_out_orig + '{}_NDJFM/'.format(area)
    ctl.mkdir(cart_out)

    trend_ssp, residtime_ssp = pickle.load(open(cart_v5.format(area) + 'trends_wrfreq_e_restime_{}.p'.format(area), 'rb'))
    seasfreq, runfreq = pickle.load(open(cart_v5.format(area) + 'seasfreqs_{}_v4.p'.format(area), 'rb'))
    freqs, residtimes, patterns = pickle.load(open(cart_v5.format(area) + 'allresults_dicts_{}_v3.p'.format(area), 'rb'))

    freqs_cmip5, trend_ssp_cmip5, residtimes_cmip5 = pickle.load(open(cart_cmip5.format(area) + 'freqs_cmip5_{}.p'.format(area), 'rb'))
    freqs.update(freqs_cmip5)
    trend_ssp.update(trend_ssp_cmip5)
    residtimes.update(residtimes_cmip5)

    allssps = ['rcp85_cmip5', 'ssp126', 'ssp245', 'ssp370', 'ssp585']
    allsims = ['hist_cmip5', 'rcp85_cmip5', 'hist', 'ssp126', 'ssp245', 'ssp370', 'ssp585']
    coldic = dict(zip(ctl.color_set(7), allsims))
    colsim = ctl.color_set(7)
    colssp = [coldic[ssp] for ssp in allssps]

    yr = np.arange(1965, 2100)

    fig = plt.figure(figsize = (16,12))
    axes = []
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg+1)
        for ssp in allssps:
            col = colsim[ssp]
            coso = runfreq[(ssp, 'run20', reg)]-np.mean(freqs[('hist', 'all', 'tot50')][:, reg])
            coserr = runfreq[(ssp, 'run20err', reg)]
            ax.fill_between(yr, coso-coserr, coso+coserr, color = col, alpha = 0.15)
            ax.plot(yr, coso, label = ssp, color = col, linewidth = 2)
        ax.set_title(reg_names_area[area][reg])
        ax.axvline(2015, color = 'lightslategray', linewidth = 0.2, linestyle = '--')
        ax.axhline(0., color = 'lightslategray', linewidth = 0.2)
        axes.append(ax)

    ctl.adjust_ax_scale(axes)
    ctl.custom_legend(fig, colssp, allssps, ncol = 5)
    fig.savefig(cart_out + 'allssps_freq20_{}_anom_wcmip5.pdf'.format(area))


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

        ctl.boxplot_on_ax(ax, allpercs, allsims, colsim, plot_mean = False, plot_ensmeans = False, plot_minmax = False)
        ax.set_xticks([])
        ax.set_title(reg_names[reg])
        if reg == 0: ax.set_ylabel('Regime frequency')

        ax.scatter(0, results_ref['freq_clus'][reg], color = 'black', marker = '*', s = 5)

    #ctl.adjust_ax_scale(axes)
    axes = []
    # qui ci metto il trend
    for reg in range(4):
        ax = figall.add_subplot(2, 4, 4 + reg + 1)
        allpercs = dict()
        for nu in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(nu)] = [np.percentile(trend_ssp[(ssp, 'all', 'trend', 'freq10', reg)], nu) for ssp in allssps]

        ctl.boxplot_on_ax(ax, allpercs, allssps, colssp, plot_mean = False, plot_ensmeans = False, plot_minmax = False)

        ax.axhline(0, color = 'gray', linewidth = 0.5)
        ax.set_xticks([])
        if reg == 0: ax.set_ylabel('Trend in regime frequency (1/yr)')
        axes.append(ax)

    ctl.adjust_ax_scale(axes)

    ctl.custom_legend(figall, colsim, allsims, ncol = 4)
    figall.savefig(cart_out + 'WRfreq_allssp_{}_8box_wtrend_wcmip5.pdf'.format(area, cos))


    fig = plt.figure(figsize = (16,12))
    axes = []
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg + 1)
        axes.append(ax)

        allpercs = dict()

        for nu in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(nu)] = [np.percentile(residtimes[(ssp, 'all', cos, reg)], nu) for ssp in allsims]

        ax.axhline(allpercs['p50'][0], color = 'gray', linewidth = 0.5)
        ctl.boxplot_on_ax(ax, allpercs, allsims, colsim, plot_mean = False, plot_ensmeans = False, plot_minmax = False)
        # ax.axhline(0, color = 'gray', linewidth = 0.5)
        ax.set_xticks([])
        ax.set_title(reg_names[reg])

    ctl.custom_legend(fig, colsim, allsims, ncol = 4)
    fig.savefig(cart_out + 'Restime_allssp_{}_{}_wcmip5.pdf'.format(area, cos))
