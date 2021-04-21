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

import xarray as xr
from scipy import stats

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 28
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.axisbelow'] = True

#############################################################################

# leggi cose
#cart_in = '/home/fabiano/Research/lavori/WeatherRegimes/tipes_nnetau_proj/'
cart_in = '/home/fedef/Research/lavori/tipes/'


patnames = ['NAO+', 'SBL', 'AR', 'NAO-']

#for tip in ['refCLUS', 'refEOF', 'refCLUS_dtr_reb']:
tip = 'refCLUS_dtr_reb'

keall = ['pi', 'ho03', 'c3r5', 'eta']
colorz = ['steelblue', 'indianred', 'forestgreen', 'orange']
# read output
resu, resu_ref = ctl.load_wrtool(cart_in + 'out_tipes_hosing_DJFM_EAT_4clus_4pcs_allyrs_{}.p'.format(tip))

# calculate frequency in running windows
# 30-yr windows, distribution of seasonal frequency
timeseries = dict()
for ke in resu.keys():
    timeseries[ke] = xr.DataArray(resu[ke]['freq_clus_seasonal'], dims = ('reg', 'time'), coords = {'reg': [0, 1, 2, 3], 'time': resu[ke]['freq_clus_seasonal_years']})

resu2, resu_ref2 = ctl.load_wrtool(cart_in + 'out_tipes_nnetau_proj_NDJFM_EAT_4clus_4pcs_allyrs_{}.p'.format(tip))
timeseries['eta'] = xr.DataArray(resu2['eta']['freq_clus_seasonal'], dims = ('reg', 'time'), coords = {'reg': [0, 1, 2, 3], 'time': resu2['eta']['freq_clus_seasonal_years']})

if 'freq_clus_seasonal' not in resu_ref:
    resu_ref['freq_clus_seasonal'], resu_ref['freq_clus_seasonal_years'] = ctl.calc_seasonal_clus_freq(resu_ref['labels'], resu_ref['dates'], 4)

gigi_obs = xr.DataArray(resu_ref['freq_clus_seasonal'], dims = ('reg', 'time'), coords = {'reg': [0, 1, 2, 3], 'time': resu_ref['freq_clus_seasonal_years']})

fig = plt.figure(figsize = (16,12))
axes = []
for num, patt in enumerate(patnames):
    ax = plt.subplot(2, 2, num+1)

    for ke, col in zip(keall, colorz):
        timeseries[ke].sel(reg = num).plot(ax = ax, color = col, linewidth = 0.2)
        rupi = ctl.running_mean(timeseries[ke].sel(reg = num).values, 10)
        ax.plot(timeseries[ke].time, rupi, color = col, linewidth = 3, label = ke)

    ax.set_title(patt, fontsize = 16)
    axes.append(ax)
    ax.grid()
    if num in [2,3]:
        ax.set_xlabel('year')
    else:
        ax.set_xlabel(None)
    if num == 1:
        ax.legend()
    if num in [0,2]: ax.set_ylabel('Regime frequency')

ctl.adjust_ax_scale(axes)
fig.savefig(cart_in + 'hosing_freq_timeseries_{}.pdf'.format(tip))

fig = plt.figure(figsize = (16,12))
axes = []
for num, patt in enumerate(patnames):
    ax = plt.subplot(2, 2, num+1)

    dist_pi = timeseries['pi'].sel(reg = num).values
    dist_ho = timeseries['ho03'].sel(reg = num, time = slice(1860, 1880)).values
    dist_rec = timeseries['c3r5'].sel(reg = num).values
    dist_eta = timeseries['eta'].sel(reg = num, time = slice(1860, 1880)).values
    dist_obs = gigi_obs.sel(reg = num).values

    dist_cose = [dist_pi, dist_ho, dist_rec, dist_eta]

    ttests01 = []
    ttests05 = []
    for co in dist_cose:
        ttests01.append(stats.ttest_ind(dist_pi, co, equal_var = False).pvalue < 0.01)
        ttests05.append(stats.ttest_ind(dist_pi, co, equal_var = False).pvalue < 0.05)

    allpercs = dict()
    for nu in [10, 25, 50, 75, 90]:
        allpercs['p{}'.format(nu)] = [np.percentile(cose, nu) for cose in dist_cose]
    allpercs['mean'] = [np.mean(cose) for cose in dist_cose]

    obsperc = dict()
    for nu in [10, 25, 50, 75, 90]:
        obsperc['p{}'.format(nu)] = np.percentile(dist_obs, nu)
    obsperc['mean'] = np.mean(dist_obs)

    #nomi = ['pi (1850-1999)', 'ho03 (1850-1879)', 'c3r5 (1900-1999)', 'eta (1850-1879)']

    positions = [0., 0.7, 1.4, 2.1, 3.2]
    #colors = ['steelblue', 'indianred', 'steelblue', 'indianred']
    ctl.boxplot_on_ax(ax, allpercs, keall, colorz, plot_mean = True, plot_ensmeans = False, obsperc = obsperc, obs_color = 'black', obs_name = 'ERA', plot_minmax = False, positions = positions)

    for pos, tte1, tte5 in zip(positions, ttests01, ttests05):
        if tte1:
            ax.scatter(pos, 5, color = 'black', marker = '*', s = 50)
        elif tte5:
            ax.scatter(pos, 5, color = 'black', marker = '*', s = 50, facecolors='none')

    ax.set_title(patt, fontsize = 16)
    ax.set_xticklabels(keall+['ERA'], fontsize = 14)
    axes.append(ax)
    #ax.grid()

    if num in [0,2]: ax.set_ylabel('Seasonal regime frequency')

ctl.adjust_ax_scale(axes)
fig.savefig(cart_in + 'hosing_freq_boxplot_{}.pdf'.format(tip))

# plots. with boxplots of variability.

# end

#### rebase to historical
