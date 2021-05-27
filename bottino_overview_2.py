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
import tunlib as tl

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

cart_out = '/home/fabiano/Research/lavori/BOTTINO/seasmean/'
ctl.mkdir(cart_out)

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

# fig = plt.figure(figsize = (24,12))
# ax1 = plt.subplot(1, 2, 1)
# ax2 = plt.subplot(1, 2, 2)
#
# for na, ru, col in zip(allnams + ['ssp585'], allru + ['ssp585'], colors + ['indianred']):
#     mem = 'r1'
#     if na == 'ssp585': mem = 'r4'
#     filist = glob.glob(filna.format(na, mem, miptab, varnam))
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
#     dates = np.array(gigi.time.data)
#     okmarch = np.array([da.month == 3 for da in dates])
#     oksept = np.array([da.month == 9 for da in dates])
#
#     resdict[(ru, varnam, 'glomean', 'mar')] = seaicearea[okmarch]
#     resdict[(ru, varnam, 'glomean', 'sep')] = seaicearea[oksept]
#
#     ax1.plot_date(dates[okmarch], seaicearea[okmarch], linestyle='solid', marker = 'None', color = col, label = ru)
#     ax2.plot_date(dates[oksept], seaicearea[oksept], linestyle='solid', marker = 'None', color = col, label = ru)
#
#     # if ru == 'ssp585':
#     #     continue
#     #
#     # if ru != 'pi':
#     #     ok200 = np.array([da.year > dates[-1].year-200 for da in dates])
#     #     varok = var[ok200]
#     #     dateok = dates[ok200]
#     # else:
#     #     varok = var
#     #     dateok = dates
#     #
#     # resdict[(ru, varnam, 'mean', 'mar')], resdict[(ru, varnam, 'std', 'mar')] = ctl.seasonal_climatology(varok, dateok, 'Mar')
#     # resdict[(ru, varnam, 'mean', 'sep')], resdict[(ru, varnam, 'std', 'sep')] = ctl.seasonal_climatology(varok, dateok, 'Sep')
#
#
# ax1.set_title('March')
# ax2.set_title('September')
# ax1.set_ylabel(r'Sea ice extent (m$^2$)')
# ax2.set_ylabel(r'Sea ice extent (m$^2$)')
# ax2.legend()
# fig.savefig(cart_out + 'bottseaice.pdf')

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

miptab = 'Amon_r25'
allvars_2D = 'clt  pr  psl  rlut  rsdt  rsut  tas  uas'.split()
allvars_3D = 'ta ua'.split()

var_glob_mean = 'tas pr clt rlut rsut net_toa'.split()  # plot global timeseries, including ssp585
var_map_200 = 'clt pr tas rlut uas'.split()  # plot last 200 mean map, stddev, low/high var wrt pi

allnams2 = allnams + ['ssp585']
allru2 = allru + ['ssp585']
colors2 = colors + ['indianred']

figs_glob = []
axs_glob = []
pimean = dict()
glomeans = dict()
mapmean = dict()

for var in var_glob_mean:
    fig, ax = plt.subplots(figsize = (12,8))
    axs_glob.append(ax)
    figs_glob.append(fig)
    ax.set_title(var)

fig_greg, ax_greg = plt.subplots(figsize = (12,8))

#for na, ru, col in zip(allnams, allru, colors):
for na, ru, col in zip(allnams2, allru2, colors2):
    mem = 'r1'
    if na == 'ssp585': mem = 'r4'

    fils = np.concatenate([glob.glob(filna.format(na, mem, miptab, var)) for var in allvars_2D[:-1]])

    kose = xr.open_mfdataset(fils, use_cftime = True)
    kose = kose.drop_vars('time_bnds')

    try:
        kose = kose.assign(net_toa = kose.rsdt - kose.rlut - kose.rsut) # net downward energy flux at TOA
    except Exception as exc:
        print(exc)
        pass

    # Separate for uas
    fils = glob.glob(filna.format(na, mem, miptab, allvars_2D[-1]))
    kosettt = xr.open_mfdataset(fils, use_cftime = True)
    kosettt = kosettt.drop_vars('time_bnds')
    kosettt = kosettt.drop_vars('height')
    kose = kose.assign(uas = kosettt.uas)

    for var, fig, ax in zip(var_glob_mean, figs_glob, axs_glob):
        if var not in kose:
            if ru == 'pi':
                pimean[var] = 0.
            continue

        coso = kose[var].mean('lon').groupby("time.year").mean()
        glomean = np.average(coso, weights = abs(np.cos(np.deg2rad(coso.lat))), axis = -1)
        if ru == 'pi':
            years = coso.year.data-2256+2015
            pimean[var] = np.mean(glomean)
        else:
            years = coso.year.data

        glomeans[(ru, var)] = (years, glomean)
        ax.plot(years, glomean-pimean[var], label = ru, color = col)

        if ru == 'b100':
            ax.legend()
            ax.grid()

    # gregory
    try:
        #ax_greg.plot(glomeans[(ru, 'tas')][1]-pimean['tas'], glomeans[(ru, 'net_toa')][1], label = ru, color = col)
        tl.gregplot_on_ax(ax_greg, glomeans[(ru, 'tas')][1]-pimean['tas'], glomeans[(ru, 'net_toa')][1], color = col, label = ru, calc_ERF = False, calc_ECS = False)
    except Exception as exc:
        print(exc)
        pass

    if ru == 'b100':
        ax_greg.legend()
        ax_greg.grid()

    if ru == 'ssp585':
        continue

    if ru != 'pi':
        kose = kose.sel(time = kose['time.year'] >= kose['time.year'].data[-1]-200)

    for var in var_map_200:
        kose_sclim = ctl.seasonal_climatology(kose[var])
        mapmean[(ru, var)] = kose_sclim

    # kose_smean = kose.groupby("time.season").mean()
    # kose_sstd = kose.groupby("time.season").std()
    # kose_p90 = kose.groupby("time.season").percentile(90)
    # kose_p10 = kose.groupby("time.season").percentile(10)
    #
    # for var in var_map_200:
    #     mapmean[(ru, var, 'mean')] = kose_smean
    #     mapmean[(ru, var, 'std')] = kose_sstd
    #     mapmean[(ru, var, 'p90')] = kose_p90
    #     mapmean[(ru, var, 'p10')] = kose_p10

    # for var in var_map_200:
    #     print(var)
    #     vmax = np.nanpercentile(flux_season[var], 98)
    #     fig = plt.figure()
    #     guplo = flux_season[var].plot.contourf(col = 'season', col_wrap = 2, levels = 11, vmax = vmax, transform = proj, figsize = (16,12), subplot_kws = {"projection": proj})
    #     guplo.map(lambda: plt.gca().coastlines())
    #     plt.title(var)
    #     plt.savefig(cart + '{}_seas_{}.pdf'.format(mod, var))


ctl.plot_pdfpages(cart_out + 'bottino_glomeans.pdf', figs_glob, True, )
fig_greg.savefig(cart_out + 'bottino_gregory.pdf')

pickle.dump([glomeans, pimean, mapmean], open(cart_out + 'bottino_seasmean_2D.p', 'wb'))

# 3D vars
for na, ru, col in zip(allnams, allru, colors):
    mem = 'r1'
    if na == 'ssp585': mem = 'r4'

    fils = np.concatenate([glob.glob(filna.format(na, mem, miptab, var)) for var in allvars_3D])

    kose = xr.open_mfdataset(fils, use_cftime = True)
    kose = kose.drop_vars('time_bnds')

    if ru != 'pi':
        kose = kose.sel(time = kose['time.year'] >= kose['time.year'].data[-1]-200)

    for var in var_map_200:
        kose_sclim = ctl.seasonal_climatology(kose[var])
        mapmean[(ru, var)] = kose_sclim

    # kose_smean = kose.groupby("time.season").mean()
    # kose_sstd = kose.groupby("time.season").std()
    # kose_p90 = kose.groupby("time.season").percentile(90)
    # kose_p10 = kose.groupby("time.season").percentile(10)
    #
    # for var in var_map_200:
    #     mapmean[(ru, var, 'mean')] = kose_smean
    #     mapmean[(ru, var, 'std')] = kose_sstd
    #     mapmean[(ru, var, 'p90')] = kose_p90
    #     mapmean[(ru, var, 'p10')] = kose_p10

pickle.dump([glomeans, pimean, mapmean], open(cart_out + 'bottino_seasmean.p', 'wb'))

allcopls = ['seamean', 'seastd', 'seap10', 'seap90']
###### Plots 2D
figs_map = []
for var in var_map_200:
    for copl in allcopls:
        #mappe = [mapmean[('pi', var, copl)]] + [mapmean[(ru, var, copl)]-mapmean[('pi', var, copl)] for ru in allru[1:]]
        mappe = [mapmean[('pi', var)][copl]] + [mapmean[(ru, var)][copl]-mapmean[('pi', var)][copl] for ru in allru[1:]]

        #mappeseas = [ma.sel(time = ma['time.season'] == seasok). for seasok in ['DJF', 'MAM', 'JJA', 'SON'] for ma in mappe]
        mappeseas = [ma.sel(season = seasok) for seasok in ['DJF', 'MAM', 'JJA', 'SON'] for ma in mappe]

        fig = ctl.plot_multimap_contour(mappeseas, figsize = (12,12))
        figs_map.append(fig)

fignames = [var+'_'+copl for var in var_map_200 for copl in allcopls]
ctl.plot_pdfpages(cart_out + 'bottino_mapmeans.pdf', figs_map, True, fignames)

figs_map = []
for var in var_map_200:
    for copl in allcopls:
        fig, axs = plt.subplots(4, 4, figsize = (12,12))
        #mappe = [mapmean[('pi', var, copl)]] + [mapmean[(ru, var, copl)]-mapmean[('pi', var, copl)] for ru in allru[1:]]
        mappe = [mapmean[('pi', var)][copl]] + [mapmean[(ru, var)][copl]-mapmean[('pi', var)][copl] for ru in allru[1:]]

        #mappeseas = [ma.sel(time = ma['time.season'] == seasok). for seasok in ['DJF', 'MAM', 'JJA', 'SON'] for ma in mappe]
        mappeseas = [ma.sel(season = seasok) for seasok in ['DJF', 'MAM', 'JJA', 'SON'] for ma in mappe]

        for ma, ax in zip(mappeseas, axs.flatten()):
            guplo = ma.plot.contourf(x = 'lat', y = 'plev', ax = ax, levels = 11, ylim = (1.e5, 1.e3), yscale = 'log')#vmax = )

        figs_map.append(fig)

fignames = [var+'_'+copl for var in var_map_200 for copl in allcopls]
ctl.plot_pdfpages(cart_out + 'bottino_crossmeans.pdf', figs_map, True, fignames)
