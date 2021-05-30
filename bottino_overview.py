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

#filna = '/nas/BOTTINO/{}/cmorized/cmor_*/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/stabilization-ssp585-2{}/r1i1p1f1/{}/{}/*/v20210315/*nc'
filna = '/nas/BOTTINO/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/{}/r1i1p1f1/{}/{}/*nc'

# filist = glob.glob('/nas/BOTTINO/b100/cmorized/cmor_2*/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/stabilization-ssp585-2100/r1i1p1f1/SImon/siconc/gn/v20210315/*nc')

allru = ['pi', 'b025', 'b050', 'b100']
allnams = ['piControl', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']

colors = ['black', 'forestgreen', 'orange', 'violet']

#############################################################################
## SEA ICE
areacelfi = '/nas/BOTTINO/areas.nc'
acel = xr.open_dataset(areacelfi)
areaT = np.array(acel['O1t0.srf'].data)

miptab = 'SImon'
varnam = 'siconc'

resdict = dict()

fig = plt.figure(figsize = (24,12))
ax1 = plt.subplot(1, 2, 1)
ax2 = plt.subplot(1, 2, 2)

for na, ru, col in zip(allnams, allru, colors):
    filist = glob.glob(filna.format(na, miptab, varnam))
    gigi = xr.open_mfdataset(filist, use_cftime=True)

    seaice = np.array(gigi.siconc.data)
    lat = np.array(gigi.latitude.data)
    okslat = lat > 40.

    areaok = areaT[okslat]
    oksi = seaice[:, okslat]
    oksi[oksi < 15.] = 0.
    oksi[oksi > 15.] = 1.
    oksiarea = oksi*areaok[np.newaxis, :]
    seaicearea = np.nansum(oksiarea, axis = 1)

    dates = np.array(gigi.time.data)
    okmarch = np.array([da.month == 3 for da in dates])
    oksept = np.array([da.month == 9 for da in dates])

    resdict[(ru, varnam, 'glomean', 'mar')] = seaicearea[okmarch]
    resdict[(ru, varnam, 'glomean', 'sep')] = seaicearea[oksept]

    if ru != 'pi':
        ok200 = np.array([da.year > dates[-1].year-200 for da in dates])
        varok = var[ok200]
        dateok = dates[ok200]
        yeaok = [da.year for da in dates[ok200]]
    else:
        varok = var
        yeaok = [da.year-2256+2015 for da in dates]
        dateok = dates

    resdict[(ru, varnam, 'mean', 'mar')], resdict[(ru, varnam, 'std', 'mar')] = ctl.seasonal_climatology(varok, dateok, 'Mar')
    resdict[(ru, varnam, 'mean', 'sep')], resdict[(ru, varnam, 'std', 'sep')] = ctl.seasonal_climatology(varok, dateok, 'Sep')

    ax1.plot(yeaok[okmarch], seaicearea[okmarch], linestyle='solid', marker = 'None', color = col, label = ru)
    ax2.plot(yeaok[oksept], seaicearea[oksept], linestyle='solid', marker = 'None', color = col, label = ru)

ax1.set_title('March')
ax2.set_title('September')
ax1.set_ylabel(r'Sea ice extent (m$^2$)')
ax2.set_ylabel(r'Sea ice extent (m$^2$)')
ax2.legend()
fig.savefig(cart_out + 'bottseaice.pdf')

sys.exit()
#pinuc = xclim.indices.sea_ice_area(gigi, acel.data)

#To describe the stratospheric polar vortex (SPV), we follow Wu et al. (2019) and compute the average zonal wind velocity over 60–75 ◦ N but at 20 hPa instead of 10 hPa.


### Mean state temperature e varianza?
### Mean state precipitazione e varianza
### Mean state wind e varianza
miptab = 'Amon'
var = 'tas'

for varnam in ['tas', 'pr', 'uas']:
    print(varnam)
    for ru, col in zip(allru, colors):
        print(ru)
        filist = glob.glob(filna.format(ru, ru[1:], miptab, varnam))
        gigi = xr.open_mfdataset(filist, use_cftime=True)

        var = np.array(gigi[varnam].data)
        lat = np.array(gigi.lat.data)
        lon = np.array(gigi.lon.data)
        dates = np.array(gigi.time.data)

        varye, datye = ctl.yearly_average(var, dates)
        glomean = ctl.global_mean(varye, lat)
        resdict[(ru, varnam, 'glomean')] = glomean

        ok200 = np.array([da.year > dates[-1].year-200 for da in dates])
        varok = var[ok200]
        dateok = dates[ok200]

        resdict[(ru, varnam, 'mean200')], resdict[(ru, varnam, 'std200')] = ctl.seasonal_climatology(varok, dateok, 'year')
        resdict[(ru, varnam, 'mean200', 'DJFM')], resdict[(ru, varnam, 'std200', 'DJFM')] = ctl.seasonal_climatology(varok, dateok, 'DJFM')
        resdict[(ru, varnam, 'mean200', 'JJAS')], resdict[(ru, varnam, 'std200', 'JJAS')] = ctl.seasonal_climatology(varok, dateok, 'JJAS')


# 3D vars
for varnam in ['ta', 'ua', 'zg']:
    print(varnam)
    for ru, col in zip(allru, colors):
        print(ru)
        filist = glob.glob(filna.format(ru, ru[1:], miptab, varnam))
        gigi = xr.open_mfdataset(filist, use_cftime=True)

        var = np.array(gigi[varnam].data)
        lat = np.array(gigi.lat.data)
        lon = np.array(gigi.lon.data)
        dates = np.array(gigi.time.data)

        varye, datye = ctl.yearly_average(var, dates)
        glomean = ctl.global_mean(varye, lat)
        resdict[(ru, varnam, 'glomean')] = glomean

        ok200 = np.array([da.year > dates[-1].year-200 for da in dates])
        varok = var[ok200]
        dateok = dates[ok200]

        resdict[(ru, varnam, 'mean200')], resdict[(ru, varnam, 'std200')] = ctl.seasonal_climatology(varok, dateok, 'year')
        resdict[(ru, varnam, 'mean200', 'DJFM')], resdict[(ru, varnam, 'std200', 'DJFM')] = ctl.seasonal_climatology(varok, dateok, 'DJFM')
        resdict[(ru, varnam, 'mean200', 'JJAS')], resdict[(ru, varnam, 'std200', 'JJAS')] = ctl.seasonal_climatology(varok, dateok, 'JJAS')


pickle.dump(glomean, open(cart_out + 'mean_climate.p', 'wb'))


for varnam in ['tas', 'pr', 'uas']:
    fiout = cart_out + '{}_mean200.pdf'.format(varnam)
    varall = [resdict[(ru, varnam, 'mean200')] for ru in allru] + [resdict[(ru, varnam, 'std200')] for ru in allru]
    ctl.plot_multimap_contour(varall, lat, lon, fiout, plot_anomalies = False, fix_subplots_shape = (2, 3))

    fiout = cart_out + '{}_mean200_rsp025.pdf'.format(varnam)
    varall = [resdict[(ru, varnam, 'mean200')]-resdict[('b025', varnam, 'mean200')] for ru in allru[1:]] + [resdict[(ru, varnam, 'std200')] for ru in allru[1:]]
    ctl.plot_multimap_contour(varall, lat, lon, fiout, plot_anomalies = False, fix_subplots_shape = (2, 2))

    fig = plt.figure(figsize = (16,12))
    for ru, col in zip(allru, colors):
        plt.plot(years, resdict[(ru, varnam, 'glomean')], color = col, label = ru)
    plt.legend()
    fig.savefig(cart_out + '{}_glomean.pdf'.format(varnam))


for season in ['DJFM', 'JJAS']:
    for varnam in ['tas', 'pr', 'uas']:
        fiout = cart_out + '{}_mean200_{}.pdf'.format(varnam, season)
        varall = [resdict[(ru, varnam, 'mean200', season)] for ru in allru] + [resdict[(ru, varnam, 'std200', season)] for ru in allru]
        ctl.plot_multimap_contour(varall, lat, lon, fiout, plot_anomalies = False, fix_subplots_shape = (2, 3))

        fiout = cart_out + '{}_mean200_rsp025_{}.pdf'.format(varnam, season)
        varall = [resdict[(ru, varnam, 'mean200', season)]-resdict[('b025', varnam, 'mean200', season)] for ru in allru[1:]] + [resdict[(ru, varnam, 'std200', season)] for ru in allru[1:]]
        ctl.plot_multimap_contour(varall, lat, lon, fiout, plot_anomalies = False, fix_subplots_shape = (2, 2))
