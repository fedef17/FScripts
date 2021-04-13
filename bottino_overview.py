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
import xarray as xr
import xclim as xcl
import glob

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 28
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.axisbelow'] = True

#############################################################################
import glob
import xarray as xr
import xclim

cart_out = '/home/fabiano/Research/lavori/BOTTINO/analisi/'

filna = '/nas/BOTTINO/{}/cmorized/cmor_*/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/stabilization-ssp585-2{}/r1i1p1f1/{}/{}/*/v20210315/*nc'

# filist = glob.glob('/nas/BOTTINO/b100/cmorized/cmor_2*/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/stabilization-ssp585-2100/r1i1p1f1/SImon/siconc/gn/v20210315/*nc')

allru = ['b025', 'b050', 'b100']
colors = ['forestgreen', 'orange', 'violet']

#############################################################################
## SEA ICE
areacelfi = '/nas/BOTTINO/areas.nc'
acel = xr.open_dataset(areacelfi)
areaT = np.array(acel['O1t0.srf'].data)

miptab = 'SImon'
var = 'siconc'

# fig = plt.figure(figsize = (24,12))
# ax1 = plt.subplot(1, 2, 1)
# ax2 = plt.subplot(1, 2, 2)
#
# for ru, col in zip(allru, colors):
#     filist = glob.glob(filna.format(ru, ru[1:], miptab, var))
#     gigi = xr.open_mfdataset(filist, use_cftime=True)
#
#     seaice = np.array(gigi.siconc.data)
#     lat = np.array(gigi.latitude.data)
#     okslat = lat > 40.
#
#     areaok = areaT[okslat]
#     oksi = seaice[:, okslat]
#     oksi[oksi < 15.] = 0.
#     oksi[oksi > 15.] = 1.
#     oksiarea = oksi*areaok[np.newaxis, :]
#     seaicearea = np.nansum(oksiarea, axis = 1)
#
#     date = np.array(gigi.time.data)
#     okmarch = np.array([da.month == 3 for da in date])
#     oksept = np.array([da.month == 9 for da in date])
#
#     ax1.plot_date(date[okmarch], seaicearea[okmarch], linestyle='solid', marker = 'None', color = col, label = ru)
#     ax2.plot_date(date[oksept], seaicearea[oksept], linestyle='solid', marker = 'None', color = col, label = ru)
#
# ax1.set_title('March')
# ax2.set_title('September')
# ax2.legend()
# fig.savefig(cart_out + 'bottseaice.pdf')

#pinuc = xclim.indices.sea_ice_area(gigi, acel.data)

### Mean state temperature e varianza?
### Mean state precipitazione e varianza
### Mean state wind e varianza
miptab = 'Amon'
var = 'tas'

resdict = dict()
for var in ['tas', 'pr', 'uas']:
    for ru, col in zip(allru, colors):
        filist = glob.glob(filna.format(ru, ru[1:], miptab, var))
        gigi = xr.open_mfdataset(filist, use_cftime=True)

        var = np.array(gigi[var].data)
        lat = np.array(gigi.lat.data)
        lon = np.array(gigi.lon.data)
        dates = np.array(gigi.time.data)

        yearall = np.array([da.year for da in dates])
        years = np.unique(yearall)
        glomean = []
        for ye in years:
            okye = yearall == ye
            glomean.append(np.mean(ctl.global_mean(var[okye], lat)))
        glomean = np.array(glomean)

        ok200 = np.array([da.year > da[-1].year-200 for da in dates])
        varok = var[ok200]
        dateok = dates[ok200]

        resdict[(ru, var, 'mean200')], resdict[(ru, var, 'std200')] = ctl.seasonal_climatology(varok, dateok, 'year')
        resdict[(ru, var, 'glomean')] = glomean


for var in ['tas', 'pr', 'uas']:
    fiout = cart_out + '{}_mean200.pdf'.format(var)
    varall = [resdict[(ru, var, 'mean200')] for ru in allru] + [resdict[(ru, var, 'std200')] for ru in allru]
    ctl.plot_multimap_contour(varall, lat, lon, fiout, plot_anomalies = False, fix_subplots_shape = (2, 3))

    fiout = cart_out + '{}_mean200_rsp025.pdf'.format(var)
    varall = [resdict[(ru, var, 'mean200')]-resdict[('b025', var, 'mean200')] for ru in allru[1:]] + [resdict[(ru, var, 'std200')] for ru in allru[1:]]
    ctl.plot_multimap_contour(varall, lat, lon, fiout, plot_anomalies = False, fix_subplots_shape = (2, 2))

    fig = plt.figure(figsize = (16,12))
    for ru, col in zip(allru, colors):
        plt.plot(years, resdict[(ru, var, 'glomean')], color = col, label = ru)
    plt.legend()
    fig.savefig(cart_out + '{}_glomean.pdf'.format(var))
