#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
from matplotlib import pyplot as plt
from matplotlib import cm

import pickle

import climtools_lib as ctl
import xarray as xr
import glob

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 22
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['axes.axisbelow'] = True

#############################################################################

if os.uname()[1] == 'hobbes':
    cart_out = '/home/fabiano/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'xaru':
    cart_out = '/home/fedef/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'tintin':
    cart_out = '/home/fabiano/work/lavori/BOTTINO/'

cart_out = cart_out + 'amoc/'

allru = ['b990', 'b025', 'b050', 'b065', 'b080', 'b100']
allsyear = [1990, 2025, 2050, 2065, 2080, 2100]
allnams = ['stabilization-hist-1990', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2065', 'stabilization-ssp585-2080', 'stabilization-ssp585-2100']

#colors = ['teal', 'forestgreen', 'orange', 'chocolate', 'maroon', 'violet']
colors = ['lightslategray', 'forestgreen', 'orange', 'chocolate', 'maroon', 'violet']


allru_wI = ['b065', 'b65I', 'b080', 'b80I', 'b100', 'b00I']
colors_wI = ['chocolate', 'peru', 'maroon', 'firebrick', 'violet', 'plum']

amoc_all = pickle.load(open('/home/fabiano/Research/lavori/BOTTINO/amoc/amoc_all_1000.p', 'rb'))

amoc_pi = dict()

pino = xr.load_dataset('/nas/TIPES/Bottino/piControl/msftyz_Omon_EC-Earth3_piControl_r1i1p1f1_gn_1850-2350.nc', use_cftime = True, decode_times=False)['msftyz']
amoc_max = pino.sel(basin = 1, rlat = slice(30, 50), lev = slice(500., 2000.)).max(['rlat', 'lev']).values
amoc_pi[('pi', 'amoc_max')] = amoc_max
gino = pino.sel(basin = 1, rlat = slice(30, 50))
zino = np.all(gino.values < amoc_max[:, np.newaxis, np.newaxis]/2., axis = 2)
amoc_wid = []
for co in zino:
    for i, lev in enumerate(pino.lev):
            if np.all(co[i:]):
                amoc_wid.append(lev.values)
                break
amoc_wid = np.stack(amoc_wid)
amoc_pi[('pi', 'amoc_wid')] = amoc_wid

zuki = pino.sel(basin = 1, rlat = slice(30, 50), lev = slice(500., 2000.)).argmax(['rlat', 'lev'])
amoc_max_lev = pino.sel(lev = slice(500., 2000.)).lev[zuki['lev']].values
amoc_pi[('pi', 'amoc_maxlev')] = amoc_max_lev

amax = pino.sel(basin = 0, lev = slice(500., 3000.), rlat = slice(-90, -60)).min(['rlat', 'lev']).values
amoc_pi[('pi', 'aabw_max')] = amax
zuki = pino.sel(basin = 0, lev = slice(500., 3000.), rlat = slice(-90, -60)).argmin(['rlat', 'lev'])
amoc_pi[('pi', 'aabw_maxlev')] = pino.sel(lev = slice(500., 3000.)).lev[zuki['lev']].values

amax = pino.sel(basin = 0, lev = slice(2500., 4500.)).min(['rlat', 'lev']).values
amoc_pi[('pi', 'aby_max')] = amax
zuki = pino.sel(basin = 0, lev = slice(2500., 4500.)).argmin(['rlat', 'lev'])
amoc_pi[('pi', 'aby_maxlev')] = pino.sel(lev = slice(2500., 4500.)).lev[zuki['lev']].values

#amax = pino.sel(basin = 0, lev = slice(500., 3000.), rlat = slice(-90, -40)).max(['rlat', 'lev']).values
#amax = pino.sel(basin = 0, lev = slice(2000., None), rlat = slice(-90, -35)).min(['rlat', 'lev']).values
#amax = pino.sel(basin = 0, lev = slice(2000., None), rlat = slice(-90, -35)).mean('lev').min('rlat').values
amax = pino.sel(basin = 0, lev = slice(3000., 4000.), rlat = slice(-50, -30)).mean('lev').min('rlat').values
amoc_pi[('pi', 'smoc_max')] = amax
zuki = pino.sel(basin = 0, lev = slice(3000., 4000.), rlat = slice(-50, -30)).argmin(['rlat', 'lev'])
maxlev = pino.sel(lev = slice(3000., 4000.)).lev[zuki['lev']].values
amoc_pi[('pi', 'smoc_maxlev')] = maxlev

pickle.dump(amoc_pi, open(cart_out + 'amoc_pi.p', 'wb'))

nyea = 50

amoc_pi = pickle.load(open(cart_out + 'amoc_pi.p', 'rb'))
amoc_all.update(amoc_pi)

# for ke in amoc_all:
#     if len(amoc_all[ke]) > 500:
#         print(ke, len(amoc_all[ke]))
#         amoc_all[ke] = amoc_all[ke][:500]

for cos, lab in zip(['amoc_max', 'aabw_max', 'aby_max', 'smoc_max'], ['AMOC (Sv)', 'AABW cell (Sv)', 'Abyssal cell (Sv)', 'Southern Ocean MOC (Sv)']):
    fac = 1/1.e9
    cose = [amoc_all[('pi', cos)]*fac] + [amoc_all[(ru, cos)][i1:i2]*fac for ru in allru for (i1, i2) in [(0, nyea), (-nyea, None)]]
    names = ['pi'] + allru + allru

    edgecol = np.append(['black'], np.concatenate([('lightslategray', col) for col in colors]))
    fullcol = np.append(['black'], np.concatenate([(col, col) for col in colors]))

    positions = [0.]
    posticks = [0.]
    for i in range(len(allru)):
        positions.append(positions[-1]+0.9+0.5)
        positions.append(positions[-1]+0.9)
        posticks.append(np.mean(positions[-2:]))

    fig, ax = ctl.boxplot(cose, names, fullcol, positions = positions, edge_colors = edgecol, plot_mean = False, plot_minmax = False, plot_ensmeans = False)
    ax.set_xticks(posticks)
    ax.set_xticklabels(['pi'] + allru)
    ax.grid(axis = 'y')

    ax.set_ylabel(lab)

    fig.savefig(cart_out + '{}_boxplot.pdf'.format(cos))

    ##### Adding I experiments
    continue

    cose = [amoc_all[('pi', cos)]*fac] + [amoc_all[(ru, cos)]*fac for ru in allru_wI]
    names = ['pi'] + allru_wI

    edgecol = ['black'] + colors_wI
    fullcol = ['black'] + colors_wI

    positions = [0.]
    posticks = [0.]
    for i in range(len(allru)):
        positions.append(positions[-1]+0.7)
        posticks.append(positions[-1])

    fig, ax = ctl.boxplot(cose, names, fullcol, positions = positions, edge_colors = edgecol, plot_mean = False, plot_minmax = False, plot_ensmeans = False)
    ax.set_xticks(posticks)
    ax.set_xticklabels(['pi'] + allru_wI)
    ax.grid(axis = 'y')

    ax.set_ylabel(lab)
    fig.savefig(cart_out + '{}_boxplot_wI.pdf'.format(cos))

### now the amoc depth
for cos, lab in zip(['amoc_wid', 'amoc_maxlev', 'aabw_maxlev', 'aby_maxlev', 'smoc_maxlev'], ['Lower level at half maximum (m)']+4*['Depth of max (m)']):
    cose = [amoc_all[('pi', cos)]] + [amoc_all[(ru, cos)][i1:i2] for ru in allru for (i1, i2) in [(0, nyea), (-nyea, None)]]
    names = ['pi'] + allru

    fullcol = np.append(['black'], np.concatenate([(col, col) for col in colors]))
    markers = ['o'] + len(colors)*['<','o']

    positions = [0.]
    posticks = [0.]
    for i in range(len(allru)):
        positions.append(positions[-1]+0.9+0.5)
        positions.append(positions[-1]+0.9)
        posticks.append(np.mean(positions[-2:]))

    fig, ax = plt.subplots(figsize = (12,8))

    for val, col, pos, mar in zip(cose, fullcol, positions, markers):
        ax.scatter(pos, np.mean(val), color = col, marker = mar, s = 100)

    ax.set_xticks(posticks)
    ax.set_xticklabels(['pi'] + allru)
    ax.set_ylabel(lab)
    ax.invert_yaxis()
    ax.grid(axis = 'y')

    fig.savefig(cart_out + '{}_scatter.pdf'.format(cos))

    ### now for the I exps
    continue

    cose = [amoc_all[('pi', cos)]] + [amoc_all[(ru, cos)] for ru in allru_wI]
    names = ['pi'] + allru_wI

    fullcol = ['black'] + colors_wI
    markers = ['o'] + len(colors)*['o']

    positions = [0.]
    posticks = [0.]
    for i in range(len(allru)):
        positions.append(positions[-1]+0.7)
        posticks.append(positions[-1])

    fig, ax = plt.subplots(figsize = (12,8))

    for val, col, pos, mar in zip(cose, fullcol, positions, markers):
        ax.scatter(pos, np.mean(val), color = col, marker = mar, s = 100)

    ax.set_xticks(posticks)
    ax.set_xticklabels(['pi'] + allru_wI)
    ax.set_ylabel(lab)
    ax.invert_yaxis()
    ax.grid(axis = 'y')

    fig.savefig(cart_out + '{}_scatter_wI.pdf'.format(cos))


# timeseries
fig = plt.figure(figsize = (16,9))

for ru, col in zip(allru, colors):
    amind = amoc_all[(ru, 'amoc_max')]/1.e9
    plt.plot(amind, ls = ":", lw = 0.1, color = col)
    amoc_low = ctl.running_mean(amind, 20)
    plt.plot(amoc_low, lw = 2, color = col, label = ru)

plt.legend()
plt.grid()
fig.savefig(cart_out + 'amoc_timeseries_1000.pdf')
