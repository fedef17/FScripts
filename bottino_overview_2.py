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
import glob
import xclim

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 28
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.axisbelow'] = True

#############################################################################

cart_out = '/home/fabiano/Research/lavori/BOTTINO/analisi/'

#filna = '/nas/BOTTINO/{}/cmorized/cmor_*/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/stabilization-ssp585-2{}/r1i1p1f1/{}/{}/*/v20210315/*nc'

#bcart = '/nas/BOTTINO/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/'
filna = '/nas/BOTTINO/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/{}/{}i1p1f1/{}/{}/*nc'

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

for na, ru, col in zip(allnams + ['ssp585'], allru + ['ssp585'], colors + ['indianred']):
    mem = 'r1'
    if na == 'ssp585': mem = 'r4'
    filist = glob.glob(filna.format(na, mem, miptab, varnam))
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

    ax1.plot_date(dates[okmarch], seaicearea[okmarch], linestyle='solid', marker = 'None', color = col, label = ru)
    ax2.plot_date(dates[oksept], seaicearea[oksept], linestyle='solid', marker = 'None', color = col, label = ru)

    # if ru == 'ssp585':
    #     continue
    #
    # if ru != 'pi':
    #     ok200 = np.array([da.year > dates[-1].year-200 for da in dates])
    #     varok = var[ok200]
    #     dateok = dates[ok200]
    # else:
    #     varok = var
    #     dateok = dates
    #
    # resdict[(ru, varnam, 'mean', 'mar')], resdict[(ru, varnam, 'std', 'mar')] = ctl.seasonal_climatology(varok, dateok, 'Mar')
    # resdict[(ru, varnam, 'mean', 'sep')], resdict[(ru, varnam, 'std', 'sep')] = ctl.seasonal_climatology(varok, dateok, 'Sep')


ax1.set_title('March')
ax2.set_title('September')
ax1.set_ylabel(r'Sea ice extent (m$^2$)')
ax2.set_ylabel(r'Sea ice extent (m$^2$)')
ax2.legend()
fig.savefig(cart_out + 'bottseaice.pdf')

# coso = gigi['siconc'].sel(time = gigi['time.month'] == 3).mean('time')
# fig, ax = ctl.get_cartopy_fig_ax('nearside', (90., 0.), 30.)
# coso.plot.pcolormesh(ax=ax, x='longitude', y='latitude', vmin = 0.1, vmax = 1.0, transform = ccrs.PlateCarree(), extend = None);
# gogo = gigi.rename({'latitude': 'lat', 'longitude': 'lon'})
# import xesmf as xe
# ds_out = xe.util.grid_global(1, 1)
# regridder = xe.Regridder(gogo, ds_out, 'bilinear')


#pinuc = xclim.indices.sea_ice_area(gigi, acel.data)

#To describe the stratospheric polar vortex (SPV), we follow Wu et al. (2019) and compute the average zonal wind velocity over 60–75 ◦ N but at 20 hPa instead of 10 hPa.

####################################################################################################

miptab = 'Amon'
allvars_2D = 'clt  pr  psl  rlut  rsdt  rsut  tas  uas'.split()

exp = 'b025'
'stabilization-ssp585-2025'
for exp in
if exp == 'pi'
fils = np.concatenate([glob.glob(filna.format(exp, mem, miptab, var)) for var in allvars_2D])

kose = xr.open_mfdataset(fils, use_cftime = True)
kose = kose.drop_vars('time_bnds')

sys.exit()

var_glob_mean = 'tas pr clt net_toa'.split()  # plot global timeseries, including ssp585
var_map_200 = 'clt pr tas rlut uas'.split()  # plot last 200 mean map, stddev, low/high var wrt pi


for var in var_glob_mean:
    ### Add global mean timeseries
    fig = plt.figure()

    for na, ru, col in zip(allnams2, allru2, colors2):
        mem = 'r1'
        if na == 'ssp585': mem = 'r4'

        fils = glob.glob(filna.format(na, mem, miptab, var))

        kose = xr.open_mfdataset(fils, use_cftime = True)
        kose = kose.drop_vars('time_bnds')

        coso = kose[var].mean('lon').groupby("time.year").mean()
        glomean = np.average(coso, weights = abs(np.cos(np.deg2rad(coso.lat))), axis = -1)
        if ru == 'pi':
            years = coso.year.data-2256+2015
        else:
            years = coso.year.data
        plt.plot(coso.year.data, glomean, label = ru, color = col)

    plt.grid()
    plt.title(var)
    plt.legend()
    figs_global.append(fig)

sys.exit()

### Mean state temperature e varianza?
### Mean state precipitazione e varianza
### Mean state wind e varianza
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



####################################################################################################

allvars_3D = 'ta ua zg'.split()
lev_map_3D = [85000, 50000]


AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA