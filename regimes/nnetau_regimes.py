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

for tip in ['refCLUS', 'refEOF', 'refCLUS_dtr_reb']:
    # read output
    resu, resu_ref = ctl.load_wrtool(cart_in + 'out_tipes_nnetau_proj_NDJFM_EAT_4clus_4pcs_allyrs_{}.p'.format(tip))

    # calculate frequency in running windows
    # 30-yr windows, distribution of seasonal frequency
    years = resu['pi']['freq_clus_seasonal_years']
    seas_freq = resu['pi']['freq_clus_seasonal']
    gigi_pi = xr.DataArray(seas_freq, dims = ('reg', 'time'), coords = {'reg': [0, 1, 2, 3], 'time': years})
    seas_freq = resu['eta']['freq_clus_seasonal']
    gigi_eta = xr.DataArray(seas_freq, dims = ('reg', 'time'), coords = {'reg': [0, 1, 2, 3], 'time': years})

    if 'freq_clus_seasonal' not in resu_ref:
        resu_ref['freq_clus_seasonal'], resu_ref['freq_clus_seasonal_years'] = ctl.calc_seasonal_clus_freq(resu_ref['labels'], resu_ref['dates'], 4)

    gigi_obs = xr.DataArray(resu_ref['freq_clus_seasonal'], dims = ('reg', 'time'), coords = {'reg': [0, 1, 2, 3], 'time': resu_ref['freq_clus_seasonal_years']})

    fig = plt.figure(figsize = (16,12))
    axes = []
    for num, patt in enumerate(patnames):
        ax = plt.subplot(2, 2, num+1)

        gigi_pi.sel(reg = num).plot(ax = ax, color = 'steelblue', linewidth = 0.2)
        gigi_eta.sel(reg = num).plot(ax = ax, color = 'indianred', linewidth = 0.2)

        rupi = ctl.running_mean(gigi_pi.sel(reg = num).values, 10)
        rueta = ctl.running_mean(gigi_eta.sel(reg = num).values, 10)
        ax.plot(years, rupi, color = 'steelblue', linewidth = 3, label = 'pi')
        ax.plot(years, rueta, color = 'indianred', linewidth = 3, label = 'eta')

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
    fig.savefig(cart_in + 'freq_timeseries_{}.pdf'.format(tip))

    fig = plt.figure(figsize = (16,12))
    axes = []
    for num, patt in enumerate(patnames):
        ax = plt.subplot(2, 2, num+1)

        dist_pi = gigi_pi.sel(reg = num).values
        dist_eta = gigi_eta.sel(reg = num).values
        dist_obs = gigi_obs.sel(reg = num).values

        dist_cose = [dist_pi[:50], dist_eta[:50], dist_pi[50:], dist_eta[50:]]

        allpercs = dict()
        for nu in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(nu)] = [np.percentile(cose, nu) for cose in dist_cose]
        allpercs['mean'] = [np.mean(cose) for cose in dist_cose]

        obsperc = dict()
        for nu in [10, 25, 50, 75, 90]:
            obsperc['p{}'.format(nu)] = np.percentile(dist_obs, nu)
        obsperc['mean'] = np.mean(dist_obs)

        nomi = ['pi1', 'eta1', 'pi2', 'eta2']

        positions = [0., 0.7, 1.8, 2.5, 3.6]
        colors = ['steelblue', 'indianred', 'steelblue', 'indianred']
        ctl.boxplot_on_ax(ax, allpercs, nomi, colors, plot_mean = True, plot_ensmeans = False, obsperc = obsperc, obs_color = 'black', obs_name = 'ERA', plot_minmax = False, positions = positions)

        ax.set_title(patt, fontsize = 16)
        ax.set_xticklabels(nomi+['ERA (1964-2014)'], fontsize = 14)
        axes.append(ax)
        #ax.grid()

        if num in [0,2]: ax.set_ylabel('Seasonal regime frequency')

    ctl.adjust_ax_scale(axes)
    fig.savefig(cart_in + 'freq_boxplot_{}.pdf'.format(tip))

# plots. with boxplots of variability.

# end

#### rebase to historical
