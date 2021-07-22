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

import cartopy.crs as ccrs
import cartopy.util as cutil

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
    cart_out = '/home/fabiano/Research/lavori/BOTTINO/yearmean/'
elif os.uname()[1] == 'xaru':
    cart_out = '/home/fedef/Research/lavori/BOTTINO/yearmean/'
elif os.uname()[1] == 'tintin':
    cart_out = '/home/fabiano/work/lavori/BOTTINO/yearmean/'

ctl.mkdir(cart_out)

allru = ['pi', 'b025', 'b050', 'b100']
allnams = ['piControl', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']

colors = ['black', 'forestgreen', 'orange', 'violet']

####################################################################

fi = 'bottino_seamean_wind.p'
koso = pickle.load(open(cart_out+fi,'rb'))

# mapa = koso[('pi', 'ua')]['seamean'].sel(season = 'DJFM')
# ctl.plot_map_contour(mapa)

# fig, ax = ctl.get_cartopy_fig_ax()
#
# for ru, col in zip(allru, colors):
#     lats_N = []
#     lats_S = []
#     spd_N = []
#     spd_S = []
#
#     if ru == 'pi':
#         mapa = koso[('pi', 'ua', 'DJFM')].mean('year')
#         mapa2 = koso[('pi', 'ua', 'JJAS')].mean('year')
#     else:
#         mapa = koso[(ru, 'ua', 'DJFM')]
#         lye = mapa.year.values[-1]
#         mapa = mapa.sel(year = slice(lye-200, lye)).mean('year')
#         mapa2 = koso[(ru, 'ua', 'JJAS')]
#         lye = mapa2.year.values[-1]
#         mapa2 = mapa2.sel(year = slice(lye-200, lye)).mean('year')
#
#     for lo in mapa['lon'].values:
#         gip = mapa.sel(lon = slice(lo-0.5, lo+0.5)).sel(lat = slice(20,90)).values
#         latok = mapa['lat'].sel(lat = slice(20,90)).values
#         maxla = latok[np.argmax(gip)]
#         spd_N.append(np.max(gip))
#         lats_N.append(maxla)
#
#         gip = mapa2.sel(lon = slice(lo-0.5, lo+0.5)).sel(lat = slice(-90,-20)).values
#         latok = mapa2['lat'].sel(lat = slice(-90,-20)).values
#         maxla = latok[np.argmax(gip)]
#         lats_S.append(maxla)
#         spd_S.append(np.max(gip))
#
#     lats_N = np.array(lats_N)
#     lats_S = np.array(lats_S)
#     lats_N_smo = ctl.running_mean(lats_N, 5, cyclic = True)
#     lats_S_smo = ctl.running_mean(lats_S, 5, cyclic = True)
#
#     spd_N = np.array(spd_N)
#     spd_S = np.array(spd_S)
#
#     # lats_N_smo[spd_N < 5] = np.nan
#     # lats_S_smo[spd_S < 5] = np.nan
#
#     lons = mapa['lon'].values
#     lats_N_smo, lonu = cutil.add_cyclic_point(lats_N_smo, coord = lons)
#     lats_S_smo, lonu = cutil.add_cyclic_point(lats_S_smo, coord = lons)
#
#     ax.plot(lonu, lats_N_smo, color = col, transform = ccrs.PlateCarree())
#     ax.plot(lonu, lats_S_smo, color = col, transform = ccrs.PlateCarree())
#
# fig.savefig(cart_out + 'mean_jetstream_HR.pdf')


# fi = 'bottino_seamean_wind.p'
# koso = pickle.load(open(cart_out+fi,'rb'))

fig, ax = ctl.get_cartopy_fig_ax()

for ru, col in zip(allru, colors):
    lats_N = []
    lats_S = []
    spd_N = []
    spd_S = []

    if ru == 'pi':
        mapa = koso[('pi', 'ua', 'DJFM')]
        mapa2 = koso[('pi', 'ua', 'JJAS')]
    else:
        mapa = koso[(ru, 'ua', 'DJFM')]
        mapa2 = koso[(ru, 'ua', 'JJAS')]

        lye = mapa.year.values[-1]
        mapa = mapa.sel(year = slice(lye-200, lye))
        lye = mapa2.year.values[-1]
        mapa2 = mapa2.sel(year = slice(lye-200, lye))

    for lo in mapa['lon'].values:
        gip = mapa.sel(lon = slice(lo-0.5, lo+0.5)).sel(lat = slice(20,90)).values
        latok = mapa['lat'].sel(lat = slice(20,90)).values
        maxla = latok[np.argmax(gip, axis = 1)]
        spd_N.append(np.max(gip, axis = 1))
        lats_N.append(maxla)

        gip = mapa2.sel(lon = slice(lo-0.5, lo+0.5)).sel(lat = slice(-90,-20)).values
        latok = mapa2['lat'].sel(lat = slice(-90,-20)).values
        maxla = latok[np.argmax(gip, axis = 1)]
        lats_S.append(maxla)
        spd_S.append(np.max(gip, axis = 1))

    lats_N = np.stack(lats_N)
    lats_S = np.stack(lats_S)

    percs = [10, 25, 50, 75, 90]

    lindict = dict()
    for prc in percs:
        lindict[('N', prc)] = np.percentile(lats_N, prc, axis = 1).squeeze()
        lindict[('S', prc)] = np.percentile(lats_S, prc, axis = 1).squeeze()

    lons = mapa['lon'].values
    for pol in ['N', 'S']:
        for prc in percs:
            lindict[(pol, prc)] = ctl.running_mean(lindict[(pol, prc)], 5, cyclic = True)
            lindict[(pol, prc)], lonu = cutil.add_cyclic_point(lindict[(pol, prc)], coord = lons)

    # spd_N = np.array(spd_N)
    # spd_S = np.array(spd_S)

    # lats_N_smo[spd_N < 5] = np.nan
    # lats_S_smo[spd_S < 5] = np.nan

    for pol in ['N', 'S']:
        ax.fill_between(lonu, lindict[(pol, 10)], lindict[(pol, 90)], color = col, transform = ccrs.PlateCarree(), alpha = 0.2)
        ax.plot(lonu, lindict[(pol, 50)], color = col, transform = ccrs.PlateCarree())
        ax.plot(lonu, lindict[(pol, 25)], color = col, transform = ccrs.PlateCarree(), linewidth = 0.5, linestyle = ':')
        ax.plot(lonu, lindict[(pol, 75)], color = col, transform = ccrs.PlateCarree(), linewidth = 0.5, linestyle = ':')

fig.savefig(cart_out + 'seasvar_jetstream.pdf')

latcose = dict()
spdcose = dict()

for sea in ['DJFM', 'JJAS']:
    for ru, col in zip(allru, colors):
        if ru == 'pi':
            mapa = koso[('pi', 'ua', sea)]
        else:
            mapa = koso[(ru, 'ua', sea)]

            lye = mapa.year.values[-1]
            mapa1 = mapa.sel(year = slice(lye-200, lye))
            sye = mapa.year.values[0]
            mapa2 = mapa.sel(year = slice(sye, sye+50))

        #for lo in mapa['lon'].values:
        # area
        for area in ('NATL', 'NEPAC', 'SH'):
            lo1, lo2, la1, la2 = ctl.sel_area_translate(area)
            if lo1 < 0:
                lo1 = 300.
                lo2 = 360.

            if ru == 'pi':
                gip = mapa.sel(lon = slice(lo1, lo2)).sel(lat = slice(la1, la2)).mean('lon').values
                latok = mapa['lat'].sel(lat = slice(la1, la2)).values
                maxla = latok[np.argmax(gip, axis = 1)]
                maxsp = np.max(gip, axis = 1)

                latcose[(ru, sea, area)] = maxla
                spdcose[(ru, sea, area)] = maxsp
            else:
                gip = mapa1.sel(lon = slice(lo1, lo2)).sel(lat = slice(la1, la2)).mean('lon').values
                latok = mapa1['lat'].sel(lat = slice(la1, la2)).values
                maxla = latok[np.argmax(gip, axis = 1)]
                maxsp = np.max(gip, axis = 1)

                latcose[(ru + '_st', sea, area)] = maxla
                spdcose[(ru + '_st', sea, area)] = maxsp

                gip = mapa2.sel(lon = slice(lo1, lo2)).sel(lat = slice(la1, la2)).mean('lon').values
                latok = mapa2['lat'].sel(lat = slice(la1, la2)).values
                maxla = latok[np.argmax(gip, axis = 1)]
                maxsp = np.max(gip, axis = 1)

                latcose[(ru + '_tr', sea, area)] = maxla
                spdcose[(ru + '_tr', sea, area)] = maxsp


colors_vtr = ['black', 'lightgreen', 'forestgreen', 'moccasin', 'orange', 'thistle', 'violet']
areas = ('NATL', 'NEPAC', 'SH')

for co, na, lona in zip([latcose, spdcose], ['jetlat', 'jetspeed'], ['Jet latitude', 'Jet speed (m/s)']):
    fig, axes = plt.subplots(1, 3, figsize = (24,8))

    for area, ax in zip(areas, axes):
        if 'N' in area:
            sea = 'DJFM'
        else:
            sea = 'JJAS'

        allpercs = dict()
        for nu in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(nu)] = [np.percentile(co[('pi', sea, area)], nu)] + [np.percentile(co[(ru+'_'+vers, sea, area)], nu) for ru in allru[1:] for vers in ['tr', 'st']]
        allpercs['mean'] = [np.mean(co[('pi', sea, area)])] + [np.mean(co[(ru+'_'+vers, sea, area)]) for ru in allru[1:] for vers in ['tr', 'st']]
        allpercs['min'] = [np.min(co[('pi', sea, area)])] + [np.min(co[(ru+'_'+vers, sea, area)]) for ru in allru[1:] for vers in ['tr', 'st']]
        allpercs['max'] = [np.max(co[('pi', sea, area)])] + [np.max(co[(ru+'_'+vers, sea, area)]) for ru in allru[1:] for vers in ['tr', 'st']]

        nams = ['pi'] + [ru+'_'+vers for ru in allru[1:] for vers in ['tr', 'st']]
        edgecol = np.append(['black'], np.concatenate([(col, col) for col in colors[1:]]))

        positions = [0.]
        posticks = [0.]
        for i in range(len(allru[1:])):
            positions.append(positions[-1]+0.7+0.4)
            positions.append(positions[-1]+0.7)
            posticks.append(np.mean(positions[-2:]))

        ctl.boxplot_on_ax(ax, allpercs, nams, colors_vtr, positions = positions, edge_colors = edgecol, plot_mean = True, plot_minmax = False, plot_ensmeans = False)#, obsperc = obsperc, obs_color = 'black', obs_name = 'pi')
        # ax.axhline(0, color = 'gray', linewidth = 0.5)
        ax.set_xticks(posticks)
        ax.set_xticklabels(allru)

        ax.set_ylabel(lona)
        ax.set_title(area)

    fig.savefig(cart_out + '{}_boxplot.pdf'.format(na))
