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

mapa = koso[('pi', 'ua')]['seamean'].sel(season = 'DJFM')
ctl.plot_map_contour(mapa)

fig, ax = ctl.get_cartopy_fig_ax()

for ru, col in zip(allru, colors):
    lats_N = []
    lats_S = []
    spd_N = []
    spd_S = []

    if ru == 'pi':
        mapa = koso[('pi', 'ua')]['seamean'].sel(season = 'DJFM')
        mapa2 = koso[('pi', 'ua')]['seamean'].sel(season = 'JJAS')
    else:
        mapa = koso[(ru, 'ua', 'stab')]['seamean'].sel(season = 'DJFM')
        mapa2 = koso[(ru, 'ua', 'stab')]['seamean'].sel(season = 'JJAS')

    for lo in mapa['lon'].values:
        gip = mapa.sel(lon = slice(lo-0.5, lo+0.5)).sel(lat = slice(20,90)).values
        latok = mapa['lat'].sel(lat = slice(20,90)).values
        maxla = latok[np.argmax(gip)]
        spd_N.append(np.max(gip))
        lats_N.append(maxla)

        gip = mapa2.sel(lon = slice(lo-0.5, lo+0.5)).sel(lat = slice(-90,-20)).values
        latok = mapa2['lat'].sel(lat = slice(-90,-20)).values
        maxla = latok[np.argmax(gip)]
        lats_S.append(maxla)
        spd_S.append(np.max(gip))

    lats_N = np.array(lats_N)
    lats_S = np.array(lats_S)
    lats_N_smo = ctl.running_mean(lats_N, 5, cyclic = True)
    lats_S_smo = ctl.running_mean(lats_S, 5, cyclic = True)

    spd_N = np.array(spd_N)
    spd_S = np.array(spd_S)

    # lats_N_smo[spd_N < 5] = np.nan
    # lats_S_smo[spd_S < 5] = np.nan

    lons = mapa['lon'].values
    lats_N_smo, lonu = cutil.add_cyclic_point(lats_N_smo, coord = lons)
    lats_S_smo, lonu = cutil.add_cyclic_point(lats_S_smo, coord = lons)

    ax.plot(lonu, lats_N_smo, color = col, transform = ccrs.PlateCarree())
    ax.plot(lonu, lats_S_smo, color = col, transform = ccrs.PlateCarree())

fig.savefig(cart_out + 'mean_jetstream.pdf')
