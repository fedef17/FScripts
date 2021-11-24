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
#from tunlib import gregplot_on_ax

from matplotlib.colors import LogNorm
from datetime import datetime

from scipy import stats
import xarray as xr
import glob
import xclim

plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
titlefont = 22
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.axisbelow'] = True
plt.rcParams['legend.fontsize'] = 20

#############################################################################

if os.uname()[1] == 'hobbes':
    cart_or = '/home/fabiano/Research/lavori/BOTTINO/'
    cartind_or = '/nas/BOTTINO/indices/{}500_xr/'
elif os.uname()[1] == 'xaru':
    cart_or = '/home/fedef/Research/lavori/BOTTINO/'
    cartind_or = '/home/fedef/Research/lavori/BOTTINO/indices/data/{}500_xr/'

cart_or = cart_or + 'indices/'
ctl.mkdir(cart_or)


allru0 = ['pi', 'b025', 'b050', 'b100']
allnams0 = ['piControl', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']

colors0 = ['black', 'forestgreen', 'orange', 'violet']

allrup1 = ['pi', 'b990', 'b025', 'b050', 'b100']
allnamsp1 = ['piControl', 'stabilization-hist-1990', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']

colorsp1 = ['black', 'steelblue', 'forestgreen', 'orange', 'violet']

colors_vtr = ['black', 'lightgreen', 'forestgreen', 'moccasin', 'orange', 'thistle', 'violet']

####################################################################################################

allinds = dict()

for ind in ['enso', 'amv', 'pdo', 'nam', 'sam']:
    if ind in ['enso', 'amv']:
        varn = 'tos'
    else:
        varn = ind

    if ind == 'enso':
        allru = allrup1
        colors = colorsp1
        allnams = allnamsp1
    else:
        allru = allru0
        colors = colors0
        allnams = allnams0

    cart_out = cart_or + '{}/'.format(ind)
    ctl.mkdir(cart_out)

    cartind = cartind_or.format(ind)

    enso = dict()

    for ru in allru:
        if ru == 'pi':
            enso[ru] = xr.load_dataset(cartind + 'piControl_{}_360day.nc'.format(ind), use_cftime = True)[varn]
        elif ru == 'b990':
            enso[ru] = xr.load_dataset(cartind + 'Nino34_ORCA1.L75-b990.nc', use_cftime = True)[varn]

            base = ctl.lowpass_butter(enso[ru], 30*12)
            enso[ru] = enso[ru] - base

            enso[ru] = enso[ru].assign_coords({'time': xr.cftime_range('1990-01-01', '2208-12-31', freq = 'MS')})

            enso[ru] = enso[ru].groupby('time.month') - enso[ru].groupby('time.month').mean()

        else:
            enso[ru] = xr.load_dataset(cartind + '{}_{}_360day.nc'.format(ru, ind), use_cftime = True)[varn]


    if ind == 'enso':
        obsname = 'HadISST'
        obsnino = ctl.read_from_txt(cartind + 'nino34.long.data')
        obsni = obsnino[1].flatten()[:-3]
        base = ctl.lowpass_butter(obsni, 30*12)
        obsni = obsni-base
        date = xr.cftime_range('1870', '2021-09', freq='MS')

        enso['obs'] = xr.DataArray(data = obsni, dims = ['time'], coords = [date], name = 'tos')
        enso['obs'] = enso['obs'].groupby('time.month') - enso['obs'].groupby('time.month').mean()
    elif ind == 'amv':
        obsname = 'Kaplan SST'
        obsnino = ctl.read_from_txt(cartind + 'amon.us.long.data')
        obsni = obsnino[1].flatten()[:-3]
        date = xr.cftime_range(str(obsnino[0][0]), '2021-09', freq='MS')

        enso['obs'] = xr.DataArray(data = obsni, dims = ['time'], coords = [date], name = 'tos')
    elif ind == 'pdo':
        obsname = 'ersst'
        obsnino = ctl.read_from_txt(cartind + 'ersst.v5.pdo.dat')
        obsni = obsnino[1].flatten()[:-3]
        date = xr.cftime_range(str(obsnino[0][0]), '2021-09', freq='MS')

        enso['obs'] = xr.DataArray(data = obsni, dims = ['time'], coords = [date], name = 'pdo')


    # for ru in allru[1:]:
    #     firstye = enso[ru].time.values[0].year
    #     lastye = enso[ru].time.values[-1].year
    #     enso[ru+'_st'] = enso[ru].sel(time = slice('{}-01-01'.format(lastye-200), '{}-12-30'.format(lastye)))
    #     enso[ru+'_tr'] = enso[ru].sel(time = slice('{}-01-01'.format(firstye), '{}-12-30'.format(firstye+50)))

    if ind in ['enso', 'nam', 'sam']:
        shi = 20
        delta = 50
    else:
        shi = 40
        delta = 100

    if 'obs' in enso:
        allruK = allru + ['obs']
    else:
        allruK = allru

    enso_std50 = dict()
    enso_abs50 = dict()
    enso_yr = dict()
    for ru in allruK:
        if ind not in ['nam', 'sam']:
            # if ru == 'b990':
            #     years = np.reshape(enso[ru].time.values, (-1, 12))[:, 0].astype(int)
            #     piuz = np.reshape(enso[ru]['tos'].values, (-1, 12)).mean(axis = 1)
            #     piuz = xr.DataArray(data = piuz, dims = ['year'], coords = [years])
            # else:
            piuz = enso[ru].groupby('time.year').mean()
        elif ind == 'nam':
            piuz = ctl.seasonal_set(enso[ru], season = 'NDJFM', seasonal_stat = 'mean')
        elif ind == 'sam':
            piuz = ctl.seasonal_set(enso[ru], season = 'MJJAS', seasonal_stat = 'mean')

        # if ind == 'amv' and ru == 'pi':
        #     # remove 200yr oscillation
        #     piuzlow = ctl.lowpass_butter(piuz, 150)
        #     piuz = piuz-piuzlow

        enso_yr[ru] = piuz

        enso_std50[ru] = []
        enso_abs50[ru] = []

        y0 = 0
        y1 = y0+delta
        while y1+shi < len(piuz):
            y0 += shi
            y1 += shi
            #piuz_ch = piuz.isel(year = slice(y0, y1))
            piuz_ch = piuz[y0:y1]
            enso_std50[ru].append(piuz_ch.std())
            enso_abs50[ru].append(piuz_ch.mean())

        if ind == 'enso':
            enso_std50[(ru, 'mon')] = []

            y0 = enso[ru]['time.year'][0].values
            y1 = y0+delta
            while y1+shi < enso[ru]['time.year'][-1].values:
                y0 += shi
                y1 += shi
                piuz_ch = enso[ru].sel(time = slice(str(y0), str(y1)))
                enso_std50[(ru, 'mon')].append(piuz_ch.std())

    allinds[(ind, 'mon')] = enso
    allinds[(ind, 'yr')] = enso_yr
    allinds[(ind, 'std50')] = enso_std50
    allinds[(ind, 'abs50')] = enso_abs50

    ###
    # 1) plot index with red/blue shades (multipanels with single mems)
    fig, axs = plt.subplots(len(allru), 1, figsize = (24,16))

    for ru, ax in zip(allru, axs):
        years = enso_yr[ru].year
        dat = enso_yr[ru].squeeze()

        if ind != 'enso':
            #dat = ctl.lowpass_butter(dat, 10)

            ax.fill_between(years, np.zeros(len(years)), dat, where = dat >= 0., interpolate = True, color = 'indianred', lw = 0.)
            ax.fill_between(years, dat, np.zeros(len(years)), where = dat <= 0., interpolate = True, color = 'steelblue', lw = 0.)
            ax.plot(years, dat, color = 'grey', lw = 0.1)
        else:
            dat = enso[ru].squeeze()

            ax.fill_between(np.arange(len(dat)), np.zeros(len(dat)), dat, where = dat >= 0., interpolate = True, color = 'indianred')
            ax.fill_between(np.arange(len(dat)), dat, np.zeros(len(dat)), where = dat <= 0., interpolate = True, color = 'steelblue')
            ax.plot(np.arange(len(dat)), dat, color = 'grey', lw = 0.1)

        if ru != allru[-1]:
            ax.set_xticks([])
        else:
            ax.set_xticks(np.arange(2100, 2601, 100))
            ax.set_xticklabels(np.arange(0, 501, 100))

    ctl.adjust_ax_scale(axs, sel_axis = 'y')

    fig.savefig(cart_out + 'allru_redblue_{}.pdf'.format(ind))


    # 2) plot variance: a) running mean all-in-one, b) boxplot (with dot for transient) (as function of eq. temperature? no)
    fig, ax = plt.subplots(figsize = (12,8))

    allpercs = dict()
    for nu in [10, 25, 50, 75, 90]:
        allpercs['p{}'.format(nu)] = [np.percentile(enso_abs50[ru], nu) for ru in allru]
    allpercs['mean'] = [np.mean(enso_abs50[ru]) for ru in allru]
    allpercs['min'] = [np.min(enso_abs50[ru]) for ru in allru]
    allpercs['max'] = [np.max(enso_abs50[ru]) for ru in allru]

    nams = ['pi'] + allru[1:]
    edgecol = colors

    positions = 0.7*np.arange(len(allru))
    posticks = 0.7*np.arange(len(allru))

    ax.axhline(0., color = 'grey', lw = 0.1)
    ctl.boxplot_on_ax(ax, allpercs, nams, colors, positions = positions, edge_colors = edgecol, plot_mean = False, plot_minmax = False, plot_ensmeans = False)#, obsperc = obsperc, obs_color = 'black', obs_name = 'pi')

    # for ru, col, pos in zip(allru[1:], colors[1:], positions[1:]):
    #     ax.scatter(pos, enso_abs50[ru][0], color = col, marker = 'D', s = 100)

    # ax.axhline(0, color = 'gray', linewidth = 0.5)
    ax.set_xticks(posticks)
    ax.set_xticklabels(allru)
    #ax.set_title(tit)

    #ctl.custom_legend(fig, colors_vtr, ['pi'] + nams, ncol = 4, add_space_below = 0.1)
    ax.set_ylabel('{} index'.format(ind))

    fig.savefig(cart_out + '{}_boxplot.pdf'.format(ind))



    fig, ax = plt.subplots(figsize = (12,8))

    allpercs = dict()
    for nu in [10, 25, 50, 75, 90]:
        allpercs['p{}'.format(nu)] = [np.percentile(enso_std50[ru], nu) for ru in allru]
    allpercs['mean'] = [np.mean(enso_std50[ru]) for ru in allru]
    allpercs['min'] = [np.min(enso_std50[ru]) for ru in allru]
    allpercs['max'] = [np.max(enso_std50[ru]) for ru in allru]

    nams = ['pi'] + allru[1:]
    edgecol = colors

    positions = 0.7*np.arange(len(allru))
    posticks = 0.7*np.arange(len(allru))

    # if ind == 'enso':
    if 'obs' in enso_std50:
        obsperc = dict()
        for nu in [10, 25, 50, 75, 90]:
            obsperc['p{}'.format(nu)] = np.percentile(enso_std50['obs'], nu)
            obsperc['mean'] = np.mean(enso_std50['obs'])
            obsperc['max'] = np.max(enso_std50['obs'])
            obsperc['min'] = np.min(enso_std50['obs'])

        positions = np.append(positions, positions[-1] + 1.1)
        posticks = positions
        ctl.boxplot_on_ax(ax, allpercs, nams, colors, positions = positions, edge_colors = edgecol, plot_mean = False, plot_minmax = False, plot_ensmeans = False, obsperc = obsperc, obs_color = 'grey', obs_name = obsname)
        # ax.axhline(0, color = 'gray', linewidth = 0.5)
        ax.set_xticks(posticks)
        ax.set_xticklabels(allru + [obsname + '({}-2021)'.format(enso_yr['obs'].year[0])])

        ax.scatter(posticks[-1], enso_std50['obs'][-1], color = 'grey', marker = 'X', s = 100)
        ax.scatter(posticks[-1], enso_std50['obs'][0], color = 'grey', marker = 'D', s = 100)
    else:
        ctl.boxplot_on_ax(ax, allpercs, nams, colors, positions = positions, edge_colors = edgecol, plot_mean = False, plot_minmax = False, plot_ensmeans = False)#, obsperc = obsperc, obs_color = 'black', obs_name = 'pi')
        ax.set_xticks(posticks)
        ax.set_xticklabels(allru)

    for ru, col, pos in zip(allru[1:], colors[1:], positions[1:]):
        ax.scatter(pos, enso_std50[ru][0], color = col, marker = 'D', s = 100)

    # ax.axhline(0, color = 'gray', linewidth = 0.5)
    #ax.set_title(tit)

    #ctl.custom_legend(fig, colors_vtr, ['pi'] + nams, ncol = 4, add_space_below = 0.1)
    ax.set_ylabel('std dev of {} index'.format(ind))

    fig.savefig(cart_out + '{}_std_boxplot.pdf'.format(ind))


    if ind == 'enso':
        fig, ax = plt.subplots(figsize = (12,8))

        allpercs = dict()
        for nu in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(nu)] = [np.percentile(enso_std50[(ru, 'mon')], nu) for ru in allru]
        allpercs['mean'] = [np.mean(enso_std50[(ru, 'mon')]) for ru in allru]
        allpercs['min'] = [np.min(enso_std50[(ru, 'mon')]) for ru in allru]
        allpercs['max'] = [np.max(enso_std50[(ru, 'mon')]) for ru in allru]

        obsperc = dict()
        for nu in [10, 25, 50, 75, 90]:
            obsperc['p{}'.format(nu)] = np.percentile(enso_std50[('obs', 'mon')], nu)
        obsperc['mean'] = np.mean(enso_std50[('obs', 'mon')])
        obsperc['max'] = np.max(enso_std50[('obs', 'mon')])
        obsperc['min'] = np.min(enso_std50[('obs', 'mon')])


        nams = ['pi'] + allru[1:]
        edgecol = colors

        positions = 0.7*np.arange(len(allru))
        positions = np.append(positions, positions[-1] + 1.1)
        posticks = positions

        ctl.boxplot_on_ax(ax, allpercs, nams, colors, positions = positions, edge_colors = edgecol, plot_mean = False, plot_minmax = False, plot_ensmeans = False, obsperc = obsperc, obs_color = 'grey', obs_name = 'HadISST')

        for ru, col, pos in zip(allru[1:], colors[1:], positions[1:]):
            ax.scatter(pos, enso_std50[(ru, 'mon')][0], color = col, marker = 'D', s = 100)

        ax.scatter(posticks[-1], enso_std50[('obs', 'mon')][-1], color = 'grey', marker = 'X', s = 100)
        ax.scatter(posticks[-1], enso_std50[('obs', 'mon')][0], color = 'grey', marker = 'D', s = 100)

        # ax.axhline(0, color = 'gray', linewidth = 0.5)
        ax.set_xticks(posticks)
        ax.set_xticklabels(allru + ['HadISST (1870-2021)'])
        #ax.set_title(tit)

        #ctl.custom_legend(fig, colors_vtr, ['pi'] + nams, ncol = 4, add_space_below = 0.1)
        ax.set_ylabel('std dev of {} index'.format(ind))

        fig.savefig(cart_out + '{}_std_boxplot_monthly.pdf'.format(ind))


    # 3) plot low-res spectra all-in-one
    fig, ax = plt.subplots(figsize = (16,12))
    #allshi = [-0.3, -0.1, 0.1, 0.3]

    ax.set_xscale('log')
    ax.set_yscale('log')

    if 'obs' in enso:
        allruK = allru + ['obs']
        colorsK = colors + ['grey']
    else:
        allruK = allru
        colorsK = colors

    for ru, col in zip(allruK, colorsK):
        # piuz = enso[ru]['tos'].groupby('time.year').mean()
        # data_all = piuz.values.flatten()
        piuz = enso_yr[ru]

        if ind == 'amv' and ru == 'pi':
            print('removing amv pi oscillation')
            # remove 200yr oscillation for spectrum
            piuzlow = ctl.lowpass_butter(piuz, 150)
            piuz = piuz-piuzlow

        data_all = piuz.squeeze()

        if ind == 'enso':
            delta = 60
            shi = 15
            frme = np.arange(1, 10, 0.5)
            frmain = [2,3,4,5,6,8,10,20]
        elif ind in ['nam', 'sam']:
            delta = 100
            shi = 25
            # frbins = [0, 2, 3, 4, 5, 6, 7, 8, 10, 15, 20]
            # frme = [np.mean([fr1,fr2]) for fr1, fr2 in zip(frbins[:-1], frbins[1:])]
            frme = np.arange(1, 10, 0.5)
            frmain = [2,4,6,8,10,20]
        elif ind in ['amv', 'pdo']:
            delta = 200
            shi = 50
            frme = np.arange(0, 50, 2)
            frmain = [3, 5, 10, 20, 30, 50]
            # frbins = [0, 5, 10, 15, 20, 25, 30, 40, 50, 75, 100]
            # frme = [np.mean([fr1,fr2]) for fr1, fr2 in zip(frbins[:-1], frbins[1:])]

        allcose = []

        y0 = 0
        y1 = y0+delta

        data = data_all[y0:y1]
        ps = np.abs(np.fft.rfft(data, norm='forward'))**2
        freqs = np.fft.rfftfreq(data.size, 1)
        invfr = 1/freqs
        allcose.append(ps)

        while y1+shi < len(piuz):
            y0 += shi
            y1 += shi

            data = data_all[y0:y1]
            ps = np.abs(np.fft.rfft(data, norm='forward'))**2

            freqs = np.fft.rfftfreq(data.size, 1)

            invfr = 1/freqs

            # barz = []
            # xba = []
            # #for bi0, bi1 in zip(frbins[:-1], frbins[1:]):
            # frstep = (frme[1]-frme[0])/2
            # for fr in frme:
            #     bi0 = fr - frstep
            #     bi1 = fr + frstep
            #     xba.append('{} - {}'.format(bi0, bi1))
            #     okke = (bi0 <= invfr) & (invfr < bi1)
            #     gig = np.sum(ps[okke])
            #     barz.append(gig)
            # allcose.append(barz)
            allcose.append(ps)

        allcose = np.stack(allcose)
        nboxs = allcose.shape[1]

        spect_mean = np.mean(allcose, axis = 0)
        spect_q1 = np.percentile(allcose, 25, axis = 0)
        spect_q3 = np.percentile(allcose, 75, axis = 0)

        if ind == 'enso':
            zuc = np.sum(spect_mean)
            spect_mean /= zuc
            spect_q1 /= zuc
            spect_q3 /= zuc

        if ind == 'enso':
            ruwi = 3
        else:
            ruwi = 5

        spect_mean = ctl.running_mean(spect_mean, ruwi)
        spect_q1 = ctl.running_mean(spect_q1, ruwi)
        spect_q3 = ctl.running_mean(spect_q3, ruwi)

        #due giri di smoothing
        spect_mean = ctl.running_mean(spect_mean, ruwi)
        spect_q1 = ctl.running_mean(spect_q1, ruwi)
        spect_q3 = ctl.running_mean(spect_q3, ruwi)

        # if ind == 'enso':
        #     thres = 10
        # else:
        #     thres = 50
        # ax.fill_between(frme, spect_q1, spect_q3, color = col, alpha = 0.3)
        # ax.plot(frme, spect_mean, color = col, lw = 2)

        # ax.fill_between(invfr[invfr < thres], spect_q1[invfr < thres], spect_q3[invfr < thres], color = col, alpha = 0.1)
        # ax.plot(invfr[invfr < thres], spect_mean[invfr < thres], color = col, lw = 2)

        #ax.fill_between(freqs, spect_q1, spect_q3, color = col, alpha = 0.1)
        if ru == 'obs':
            ax.plot(invfr, spect_mean, color = col, lw = 2, ls = '--', label = obsname)
        else:
            ax.plot(invfr, spect_mean, color = col, lw = 2, label = ru)

        #if ind != 'enso':

        # allpercs = dict()
        # for nu, realnu in zip([10, 25, 50, 75, 90], [0, 20, 50, 80, 100]):
        #     allpercs['p{}'.format(nu)] = [np.percentile(allcose[:, iup], realnu) for iup in range(nboxs)]
        #
        # positions = np.arange(nboxs) + shi
        #
        # ctl.boxplot_on_ax(ax, allpercs, xba, nboxs*[col], positions = positions, edge_colors = nboxs*[col], plot_mean = False, plot_minmax = False, plot_ensmeans = False, wi = 0.1)

    # if ind == 'enso':
    #     ax.set_ylim(1.e-3, 0.2)
    # elif ind == 'amv':
    #     ax.set_ylim(5.e-6, 3e-3)

    ax.legend()
    ax.set_xticks(frmain)
    ax.set_xticklabels([str(int(fr)) for fr in frmain])
    # for ii in np.arange(nboxs-1) + 0.5:
    #     ax.axvline(ii, color = 'grey', linestyle = ':', linewidth = 0.1)
    # ax.set_xticklabels(xba, rotation = 30)
    ax.set_xlabel('Period (yr)')
    #ax.set_xlabel('Frequency (yr-1)')
    if ind == 'enso':
        ax.set_ylabel(r'Spectral power (normalized)')
    elif ind in ['amv', 'pdo']:
        ax.set_ylabel(r'Spectral power ($K^2$)')
    else:
        ax.set_ylabel(r'Spectral power')


    fig.savefig(cart_out + '{}_spectrum_all.pdf'.format(ind))
