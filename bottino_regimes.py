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
#cart_in = '/home/fedef/Research/lavori/tipes/'
cart_in = '/home/fedef/Research/lavori/BOTTINO/midlat/'

allru = ['pi', 'b025', 'b050', 'b100']
allnams = ['piControl', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']

colorz = ['black', 'forestgreen', 'orange', 'violet']

patnames = dict()
patnames['DJFM'] = ['NAO+', 'SBL', 'NAO-', 'AR']
patnames['JJAS'] = ['SBL', 'NAO+', 'NAO-/AL', 'AR']

tip = 'refCLUS_dtr_reb'

# read output
for seas in ['JJAS', 'DJFM']:
    resu, resu_ref = ctl.load_wrtool(cart_in + 'out_bottino_{}_EAT_4clus_4pcs_allyrs_{}.p'.format(seas, tip))

    # calculate frequency in running windows
    # 30-yr windows, distribution of seasonal frequency
    timeseries = dict()
    for ke in resu.keys():
        if 'freq_clus_seasonal' not in resu[ke]:
            resu[ke]['freq_clus_seasonal'], resu[ke]['freq_clus_seasonal_years'] = ctl.calc_seasonal_clus_freq(resu[ke]['labels'], resu[ke]['dates'], 4)
        timeseries[ke] = xr.DataArray(resu[ke]['freq_clus_seasonal'], dims = ('reg', 'time'), coords = {'reg': [0, 1, 2, 3], 'time': resu[ke]['freq_clus_seasonal_years']})

    # resu2, resu_ref2 = ctl.load_wrtool(cart_in + 'out_tipes_nnetau_proj_NDJFM_EAT_4clus_4pcs_allyrs_{}.p'.format(tip))
    # timeseries['eta'] = xr.DataArray(resu2['eta']['freq_clus_seasonal'], dims = ('reg', 'time'), coords = {'reg': [0, 1, 2, 3], 'time': resu2['eta']['freq_clus_seasonal_years']})

    if 'freq_clus_seasonal' not in resu_ref:
        resu_ref['freq_clus_seasonal'], resu_ref['freq_clus_seasonal_years'] = ctl.calc_seasonal_clus_freq(resu_ref['labels'], resu_ref['dates'], 4)

    #gigi_obs = xr.DataArray(resu_ref['freq_clus_seasonal'], dims = ('reg', 'time'), coords = {'reg': [0, 1, 2, 3], 'time': resu_ref['freq_clus_seasonal_years']})

    fig = plt.figure(figsize = (16,12))
    axes = []
    for num, patt in enumerate(patnames[seas]):
        ax = plt.subplot(2, 2, num+1)

        for ke, col in zip(allru, colorz):
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
    fig.savefig(cart_in + 'bottino_freqts_{}_{}.pdf'.format(seas, tip))

    fig = plt.figure(figsize = (16,12))
    axes = []
    for num, patt in enumerate(patnames[seas]):
        ax = plt.subplot(2, 2, num+1)

        dist_cose = []
        for ru in allru:
            dist_cose.append(timeseries[ru].sel(reg = num).values)

        ttests01 = []
        ttests05 = []
        for co in dist_cose:
            ttests01.append(stats.ttest_ind(dist_cose[0], co, equal_var = False).pvalue < 0.01)
            ttests05.append(stats.ttest_ind(dist_cose[0], co, equal_var = False).pvalue < 0.05)

        allpercs = dict()
        for nu in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(nu)] = [np.percentile(cose, nu) for cose in dist_cose]
        allpercs['mean'] = [np.mean(cose) for cose in dist_cose]

        # obsperc = dict()
        # for nu in [10, 25, 50, 75, 90]:
        #     obsperc['p{}'.format(nu)] = np.percentile(dist_obs, nu)
        # obsperc['mean'] = np.mean(dist_obs)

        positions = [0., 0.7, 1.4, 2.1]#, 3.2]
        #ctl.boxplot_on_ax(ax, allpercs, keall, colorz, plot_mean = True, plot_ensmeans = False, obsperc = obsperc, obs_color = 'black', obs_name = 'ERA', plot_minmax = False, positions = positions)
        ctl.boxplot_on_ax(ax, allpercs, allru, colorz, plot_mean = True, plot_ensmeans = False, obsperc = None, obs_color = 'black', obs_name = None, plot_minmax = False, positions = positions, wi = 0.4)

        for pos, tte1, tte5 in zip(positions, ttests01, ttests05):
            if tte1:
                ax.scatter(pos, 0, color = 'black', marker = '*', s = 80)
            elif tte5:
                ax.scatter(pos, 0, color = 'black', marker = '*', s = 80, facecolors='none')

        ax.set_title(patt, fontsize = 16)
        ax.set_xticklabels(allru, fontsize = 14)
        axes.append(ax)
        #ax.grid()

        if num in [0,2]: ax.set_ylabel('Seasonal regime frequency')

    ctl.adjust_ax_scale(axes)
    fig.savefig(cart_in + 'bottino_fbox_{}_{}.pdf'.format(seas, tip))

# plots. with boxplots of variability.

# end

#### rebase to historical
