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

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 22
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['axes.axisbelow'] = True

#############################################################################

if os.uname()[1] == 'hobbes':
    cart_or = '/home/fabiano/Research/lavori/BOTTINO/'
    cartind = '/nas/BOTTINO/indices/amv500_xr/'
elif os.uname()[1] == 'xaru':
    cart_or = '/home/fedef/Research/lavori/BOTTINO/'
    cartind = '/home/fedef/Research/lavori/BOTTINO/amv500_xr/'
elif os.uname()[1] == 'tintin':
    cart_or = '/home/fabiano/work/lavori/BOTTINO/'

cart_out = cart_or + 'indices/'
ctl.mkdir(cart_out)

cart_out = cart_out + 'amv/'
ctl.mkdir(cart_out)

allru = ['pi', 'b025', 'b050', 'b100']
allnams = ['piControl', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']

colors = ['black', 'forestgreen', 'orange', 'violet']
colors_vtr = ['black', 'lightgreen', 'forestgreen', 'moccasin', 'orange', 'thistle', 'violet']

####################################################################################################

enso = dict()

for ru in allru:
    if ru == 'pi':
        enso[ru] = xr.load_dataset(cartind + 'piControl_amv_360day.nc', use_cftime = True)
    else:
        enso[ru] = xr.load_dataset(cartind + '{}_amv_360day.nc'.format(ru), use_cftime = True)

        firstye = enso[ru].time.values[0].year
        lastye = enso[ru].time.values[-1].year
        enso[ru+'_st'] = enso[ru].sel(time = slice('{}-01-01'.format(lastye-200), '{}-12-30'.format(lastye)))
        enso[ru+'_tr'] = enso[ru].sel(time = slice('{}-01-01'.format(firstye), '{}-12-30'.format(firstye+50)))


fig, ax = plt.subplots(figsize = (12,8))

allpercs = dict()
for nu in [10, 25, 50, 75, 90]:
    allpercs['p{}'.format(nu)] = [np.percentile(enso['pi']['tos'], nu)] + [np.percentile(enso[ru+'_'+vers]['tos'], nu) for ru in allru[1:] for vers in ['tr', 'st']]
allpercs['mean'] = [np.mean(enso['pi']['tos']).values] + [np.mean(enso[ru+'_'+vers]['tos']).values for ru in allru[1:] for vers in ['tr', 'st']]
allpercs['min'] = [np.min(enso['pi']['tos']).values] + [np.min(enso[ru+'_'+vers]['tos']).values for ru in allru[1:] for vers in ['tr', 'st']]
allpercs['max'] = [np.max(enso['pi']['tos']).values] + [np.max(enso[ru+'_'+vers]['tos']).values for ru in allru[1:] for vers in ['tr', 'st']]

# ru = 'pi'
# obsperc = dict()
# for nu in [10, 25, 50, 75, 90]:
#     obsperc['p{}'.format(nu)] = np.percentile(enso[ru]['tos'], nu)
# obsperc['mean'] = np.mean(enso[ru]['tos']).values
# obsperc['min'] = np.min(enso[ru]['tos']).values
# obsperc['max'] = np.max(enso[ru]['tos']).values

nams = ['pi'] + [ru+'_'+vers for ru in allru[1:] for vers in ['tr', 'st']]
edgecol = np.append(['black'], np.concatenate([(col, col) for col in colors[1:]]))

positions = [0.]
posticks = [0.]
for i in range(len(allru[1:])):
    positions.append(positions[-1]+0.7+0.4)
    positions.append(positions[-1]+0.7)
    posticks.append(np.mean(positions[-2:]))

ctl.boxplot_on_ax(ax, allpercs, nams, colors_vtr, positions = positions, edge_colors = edgecol, plot_mean = False, plot_minmax = True, plot_ensmeans = False)#, obsperc = obsperc, obs_color = 'black', obs_name = 'pi')
# ax.axhline(0, color = 'gray', linewidth = 0.5)
ax.set_xticks(posticks)
ax.set_xticklabels(allru)
#ax.set_title(tit)

#ctl.custom_legend(fig, colors_vtr, ['pi'] + nams, ncol = 4, add_space_below = 0.1)
ax.set_ylabel('AMV index (K)')

fig.savefig(cart_out + 'amv_boxplot.pdf')


############# enso plot with tas/net_toa

glomeans, pimean = pickle.load(open(cart_or + 'seasmean/bottino_glomeans.p', 'rb'))

fig, axs = plt.subplots(1,2, figsize = (16,12))

for ru, col in zip(allru, colors):
    piuz = enso[ru]['tos'].groupby('time.year').mean()
    tass = glomeans[(ru, 'tas')][1]
    toaz = glomeans[(ru, 'net_toa')][1]
    # if ru == 'pi':
    #     tass = tass[:-1]
    #     toaz = toaz[:-1]
    axs[0].scatter(tass, piuz, c = toaz)
    axs[1].scatter(toaz, piuz, c = tass)

axs[0].set_title('amv vs tas')
axs[1].set_title('amv vs net_toa')
fig.savefig(cart_out + 'amv_vs_tasytoa.pdf')


fig, ax = plt.subplots(figsize = (16,12))

nye = 50
nchu = 10

markers = ['x', 'o', 's', 'D']
for ru, col, mark in zip(allru, colors, markers):
    print(ru)
    piuz = enso[ru]['tos'].groupby('time.year').mean()
    tass = glomeans[(ru, 'tas')][1]
    toaz = glomeans[(ru, 'net_toa')][1]
    # if ru == 'pi':
    #     tass = tass[:-1]
    #     toaz = toaz[:-1]
    #     ### PROBLEM WITH TOANET FOR PI!!!!
    #     toaz = np.zeros(len(toaz))

    enso_std = []
    tass50 = []
    toaz50 = []
    for ich in range(nchu):
        #DIVIDERE IN CHUNKS DI 50 anni
        #Calcolare STD di enso e plottare scatter con quello
        piuz_ch = piuz.isel(year = slice(ich*nye, (ich+1)*nye))
        enso_std.append(piuz_ch.std())
        tass50.append(tass[ich*nye:(ich+1)*nye].mean())
        toaz50.append(toaz[ich*nye:(ich+1)*nye].mean())

    ax.plot(tass50, enso_std, color = 'grey', linewidth = 0.3)
    sc = ax.scatter(tass50, enso_std, c = toaz50, vmin = 0., vmax = 1., edgecolor = 'grey', s = 100, marker = mark)

#ax.set_title('enso vs tas')
ax.set_ylabel('Amv std (K)')
ax.set_xlabel('Global tas (K)')

cb = plt.colorbar(sc)
cb.set_label('Global net TOA (W/m2)')

fig.savefig(cart_out + 'amvstd_vs_tasytoa_50.pdf')

### Trying separately with p90 e p10

fig, ax = plt.subplots(figsize = (16,12))

nye = 50
nchu = 10

markers = ['x', 'o', 's', 'D']
for ru, col, mark in zip(allru, colors, markers):
    print(ru)
    piuz = enso[ru]['tos'].groupby('time.year').mean()
    tass = glomeans[(ru, 'tas')][1]
    toaz = glomeans[(ru, 'net_toa')][1]
    # if ru == 'pi':
    #     tass = tass[:-1]
    #     toaz = toaz[:-1]
    #     ### PROBLEM WITH TOANET FOR PI!!!!
    #     toaz = np.zeros(len(toaz))

    enso_p10 = []
    enso_p90 = []
    tass50 = []
    toaz50 = []
    for ich in range(nchu):
        #DIVIDERE IN CHUNKS DI 50 anni
        #Calcolare STD di enso e plottare scatter con quello
        piuz_ch = piuz.isel(year = slice(ich*nye, (ich+1)*nye))
        enso_p90.append(np.percentile(piuz_ch, 90))
        enso_p10.append(np.percentile(piuz_ch, 10))
        tass50.append(tass[ich*nye:(ich+1)*nye].mean())
        toaz50.append(toaz[ich*nye:(ich+1)*nye].mean())

    ax.plot(tass50, enso_p90, color = 'grey', linewidth = 0.3)
    sc = ax.scatter(tass50, enso_p90, c = toaz50, vmin = 0., vmax = 1., edgecolor = 'grey', s = 100, marker = mark)

    ax.plot(tass50, enso_p10, color = 'grey', linewidth = 0.3)
    sc = ax.scatter(tass50, enso_p10, c = toaz50, vmin = 0., vmax = 1., edgecolor = 'grey', s = 100, marker = mark)

#ax.set_title('enso vs tas')
ax.set_ylabel('AMV p10/p90 (K)')
ax.set_xlabel('Global tas (K)')
ax.grid()

cb = plt.colorbar(sc)
cb.set_label('Global net TOA (W/m2)')

fig.savefig(cart_out + 'amvperc_vs_tasytoa_50.pdf')


fig, axs = plt.subplots(2, 2, figsize = (16,12))

for ru, ax in zip(allru, axs.flatten()):
    piuz = enso[ru]['tos'].groupby('time.year').mean()
    data = piuz.values.flatten()
    ps = np.abs(np.fft.fft(data))**2

    freqs = np.fft.fftfreq(data.size, 1)
    idx = freqs > 0
    ax.plot(1/freqs[idx], ps[idx], linewidth = 0.5)
    #ps_low = ctl.butter_filter(ps, 20)
    ps_low = ctl.running_mean(ps, 20)
    ax.plot(1/freqs[idx], ps_low[idx], linewidth = 2)
    #ps_low = ctl.butter_filter(ps, 50)
    ps_low = ctl.running_mean(ps, 50)
    ax.plot(1/freqs[idx], ps_low[idx], linewidth = 2)

    ax.set_title(ru)
    if ax in axs[1, :]:
        ax.set_xlabel('period (yr)')
        #ax.set_ylabel('amplitude')
    ax.set_xlim(2, 50)

#ctl.adjust_ax_scale(axs.flatten())
fig.savefig(cart_out + 'amv_spectra_all500.pdf')


#frbins = [0, 5, 10, 20, 30, 50, 100, 200]
frbins = [1, 5, 20, 50, 100, 200]

fig, axs = plt.subplots(2, 2, figsize = (16,12))

for ru, ax, col in zip(allru, axs.flatten(), colors):
    piuz = enso[ru]['tos'].groupby('time.year').mean()
    data = piuz.values.flatten()
    ps = np.abs(np.fft.rfft(data))**2

    freqs = np.fft.rfftfreq(data.size, 1)
    invfr = 1/freqs

    barz = []
    xba = []
    for bi0, bi1 in zip(frbins[:-1], frbins[1:]):
        xba.append('{} - {}'.format(bi0, bi1))
        okke = (bi0 <= invfr) & (invfr < bi1)
        gig = np.sum(ps[okke])
        barz.append(gig)
    ax.bar(np.arange(len(barz)), barz, color = col)
    ax.set_xticks(np.arange(len(barz)))
    ax.set_xticklabels(xba, rotation = 30)

    ax.grid()
    ax.set_title(ru)
    if ax in axs[1, :]:
        ax.set_xlabel('period (yr)')

ctl.adjust_ax_scale(axs.flatten())
fig.savefig(cart_out + 'amv_spectra_bins.pdf')

allshi = [-0.3, -0.1, 0.1, 0.3]
fig, ax = plt.subplots(figsize = (16,12))

for ru, col, shi in zip(allru, colors, allshi):
    piuz = enso[ru]['tos'].groupby('time.year').mean()
    data = piuz.values.flatten()
    ps = np.abs(np.fft.rfft(data))**2

    freqs = np.fft.rfftfreq(data.size, 1)
    invfr = 1/freqs

    barz = []
    xba = []
    for bi0, bi1 in zip(frbins[:-1], frbins[1:]):
        xba.append('{} - {}'.format(bi0, bi1))
        okke = (bi0 <= invfr) & (invfr < bi1)
        gig = np.sum(ps[okke])
        barz.append(gig)
    ax.bar(np.arange(len(barz))+shi, barz, color = col, width = 0.2)

ax.set_xticks(np.arange(len(barz)))
ax.set_xticklabels(xba, rotation = 30)
ax.set_xlabel('period (yr)')

fig.savefig(cart_out + 'amv_spectra_bins11.pdf')


fig, ax = plt.subplots(figsize = (16,12))

for ru, col, shi in zip(allru, colors, allshi):
    piuz = enso[ru]['tos'].groupby('time.year').mean()
    data = piuz.values.flatten()
    ps = np.abs(np.fft.rfft(data))**2

    freqs = np.fft.rfftfreq(data.size, 1)
    invfr = 1/freqs

    normco = np.sum(ps)

    barz = []
    xba = []
    for bi0, bi1 in zip(frbins[:-1], frbins[1:]):
        xba.append('{} - {}'.format(bi0, bi1))
        okke = (bi0 <= invfr) & (invfr < bi1)
        gig = np.sum(ps[okke])
        barz.append(gig)
    ax.bar(np.arange(len(barz))+shi, barz, color = col, width = 0.2)

ax.set_xticks(np.arange(len(barz)))
ax.set_xticklabels(xba, rotation = 30)
ax.set_xlabel('period (yr)')

fig.savefig(cart_out + 'amv_spectra_bins11_rel.pdf')


frbins = [1, 5, 20, 50, 100]
nchu = 3
nye = 200
yesta = [0, 150, 300]
colzz = ctl.color_set(3, cmap = 'viridis', use_seaborn=False)

ps_max = dict()
ps_ext = dict()

fig, axs = plt.subplots(2, 2, figsize = (16,12))
allshi = [-0.2, 0., 0.2]

for ru, ax in zip(allru, axs.flatten()):
    piuz = enso[ru]['tos'].groupby('time.year').mean()
    data_all = piuz.values.flatten()

    for ich, ye1, col, shi in zip(range(nchu), yesta, colzz, allshi):
        data = data_all[ye1:ye1+nye]
        ps = np.abs(np.fft.rfft(data))**2

        freqs = np.fft.rfftfreq(data.size, 1)

        invfr = 1/freqs

        barz = []
        xba = []
        for bi0, bi1 in zip(frbins[:-1], frbins[1:]):
            xba.append('{} - {}'.format(bi0, bi1))
            okke = (bi0 <= invfr) & (invfr < bi1)
            gig = np.sum(ps[okke])
            barz.append(gig)
        ax.bar(np.arange(len(barz))+shi, barz, color = col, width = 0.2)

    ax.set_xticks(np.arange(len(barz)))
    ax.set_xticklabels(xba, rotation = 30)

    ax.set_title(ru)
    if ax in axs[1, :]:
        ax.set_xticks(np.arange(len(barz)))
        ax.set_xticklabels(xba, rotation = 30)
        ax.set_xlabel('period (yr)')


ctl.adjust_ax_scale(axs.flatten())
fig.savefig(cart_out + 'amv_spectra_varying.pdf')
