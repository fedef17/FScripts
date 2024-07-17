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
    cart = '/home/fabiano/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'xaru':
    cart = '/home/fedef/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'tintin':
    cart = '/home/fabiano/work/lavori/BOTTINO/'

cart_in = cart + 'seasmean/'
cart_out = cart + 'trends/'
ctl.mkdir(cart_out)

allru = ['pi', 'b025', 'b050', 'b100']
allnams = ['piControl', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']

colors = ['black', 'forestgreen', 'orange', 'violet']

allnams2 = allnams + ['ssp585', 'historical']
allru2 = allru + ['ssp585', 'hist']
colors2 = colors + ['indianred', 'steelblue']

allnams3 = ['stabilization-hist-1990'] + allnams[1:]
allru3 = ['b990'] + allru[1:]
colors3 = ['teal'] + colors[1:]

##########################################################

cart_in = '/home/fabiano/Research/lavori/BOTTINO/tediato/output_standalone_2022-04-26T16:29:29.678755/'
cart_out = cart_in + '../tedplots/'

fils = [fi for fi in os.listdir(cart_in) if '.nc' in fi]

res = dict()
for fi in fils:
    realm, _, _, exp, time  = fi.removesuffix('.nc').split('_')
    gigi = xr.load_dataset(cart_in + fi)
    var_name = list(gigi.keys())[0]
    latname = list(gigi.coords)[0]
    lat = gigi[latname].values

    res[(realm, exp, time)] = gigi[var_name].values

allrealms = ['total', 'atmos', 'ocean', 'latent', 'wmb']
alllabs = ['Total heat transport (W)', 'Atmos. enthalpy transport (W)', 'Ocean heat transport (W)', 'Latent heat transport (W)', 'Water mass transport (Kg/s)']

for realm, tit in zip(allrealms, alllabs):
    fig, ax = plt.subplots(figsize = (16,9))
    #plt.title(realm + ' transport')
    for ru, col in zip(allru3, colors3):
        plt.plot(lat, res[(realm, ru, 'end')], color = col, label = ru, lw = 2)
        plt.plot(lat, res[(realm, ru, 'start')], color = col, ls = ':', lw = 2)

    ax.set_xlabel('Latitude')
    ax.set_ylabel(tit)

    plt.grid()
    plt.legend()
    fig.savefig(cart_out + realm + '_transport.pdf')

    fig, ax = plt.subplots()
    plt.title(realm + ' transport evolution')
    for ru, col in zip(allru3, colors3):
        plt.plot(lat, res[(realm, ru, 'end')]-res[(realm, ru, 'start')], color = col, label = ru)
    plt.grid()
    plt.legend()
    fig.savefig(cart_out + realm + '_end-start.pdf')


fig, axs = plt.subplots(2, 2, figsize = (16,9))
#plt.title(realm + ' transport')
for realm, tit, ax in zip(allrealms[:-1], alllabs[:-1], axs.flatten()):
    for ru, col in zip(allru3, colors3):
        ax.plot(lat, res[(realm, ru, 'end')], color = col, label = ru, lw = 2)
        ax.plot(lat, res[(realm, ru, 'start')], color = col, ls = ':', lw = 2)

    if ax in axs[1, :]: ax.set_xlabel('Latitude')
    ax.set_ylabel(tit)

    ax.grid()
# plt.legend()
ctl.custom_legend(fig, colors3, allru3, ncol = 4)
plt.subplots_adjust(bottom = 0.15)
fig.savefig(cart_out + 'all_transports.pdf')

fig, axs = plt.subplots(1, 2, figsize = (16,6))
#plt.title(realm + ' transport')
for realm, tit, ax in zip(allrealms[1:2], alllabs[1:2], axs.flatten()):
    for ru, col in zip(allru3, colors3):
        ax.plot(lat, res[(realm, ru, 'end')], color = col, label = ru, lw = 2)
        ax.plot(lat, res[(realm, ru, 'start')], color = col, ls = ':', lw = 2)

    if ax in axs[1, :]: ax.set_xlabel('Latitude')
    ax.set_ylabel(tit)

    ax.grid()
# plt.legend()
ctl.custom_legend(fig, colors3, allru3, ncol = 4)
plt.subplots_adjust(bottom = 0.15)
fig.savefig(cart_out + 'atm_oce_transport.pdf')


glomeans, pimean = pickle.load(open(cart_in + '../../seasmean/bottino_glomeans.p', 'rb'))

fig = plt.figure()
for realm, tit, mar in zip(allrealms[:-1], alllabs[:-1], ['o', 's', 'D', '*']):
    for ru, col in zip(allru3, colors3):
        sttemp = glomeans[(ru, 'tas')][1][:10].mean()-pimean['tas']
        entemp = glomeans[(ru, 'tas')][1][-10:].mean()-pimean['tas']
        plt.scatter(sttemp, np.max(res[(realm, ru, 'start')]), color = col, facecolor = 'none', marker = mar)
        plt.scatter(entemp, np.max(res[(realm, ru, 'end')]), color = col, marker = mar)

plt.ylabel('Meridional transport (W)')
plt.xlabel('GTAS (K)')
plt.grid()
fig.savefig(cart_out + 'all_transports_vs_gtas.pdf')
