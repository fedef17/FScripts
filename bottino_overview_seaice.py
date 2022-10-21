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
# import tunlib as tl

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
# cart_out = '/home/fabiano/Research/lavori/BOTTINO/seasmean/'
# ctl.mkdir(cart_out)
#
# filna = '/nas/BOTTINO/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/{}/{}i1p1f1/{}/{}/*nc'
#
# allru = ['pi', 'b025', 'b050', 'b100']
# allnams = ['piControl', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']
#
# colors = ['black', 'forestgreen', 'orange', 'violet']

user = 'ffabiano'
cart_out = '/g100_work/IscrB_QUECLIM/BOTTINO/bottino_an/seasmean/'
cart_run = '/g100_scratch/userexternal/{}/ece3/b00I/runtime/'.format(user)

allruadd2 = ['b065', 'b080']
allnadd2 = ['stabilization-ssp585-2065', 'stabilization-ssp585-2080']
colorsadd2 = ['chocolate', 'maroon']

allruall = ['pi', 'hist', 'ssp585', 'b990', 'b025', 'b050', 'b065', 'b080', 'b100']
okye = [(2015, 2515), (1850, 2015), (2015, 2101), (1990, 2490), (2025, 2525), (2050, 2550), (2065, 2565), (2080, 2580), (2100, 2600)]
#colall = ['black', 'steelblue', 'indianred', 'teal', 'forestgreen', 'orange', 'chocolate', 'maroon', 'violet']
colall = ['black', 'royalblue', 'crimson', 'lightslategray', 'forestgreen', 'orange', 'chocolate', 'maroon', 'violet']

# allru = ['pi', 'b025', 'b050', 'b100']
# allnams = ['piControl', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']
#
# colors = ['black', 'forestgreen', 'orange', 'violet']

#############################################################################
## SEA ICE
#areacelfi = '/nas/BOTTINO/areas.nc'
areacelfi = cart_run + 'areas.nc'
acel = xr.open_dataset(areacelfi)
areaT = np.array(acel['O1t0.srf'].data)

#miptab = 'SImon_r1'
miptab = 'SImon'
varnam = 'siconc'

resdict = pickle.load(open(cart_out + 'seaicearea.p', 'rb'))

fig, axs = plt.subplots(2,2,figsize = (12,12))
#
# allnams2 = ['ssp585', 'historical'] + allnams
# allru2 = ['ssp585', 'hist'] + allru
# colors2 = ['indianred', 'steelblue'] + colors
#
# allnams3 = allnams2 + ['stabilization-hist-1990']
# allru3 = allru2 + ['b990']
# colors3 = colors2 + ['teal']

#for na, ru, col in zip(allnams3, allru3, colors3):
for ru, col in zip(allruadd2, colorsadd2):
    print(ru)
    mem = 'r1'
    if ru in ['ssp585', 'hist']: mem = 'r4'

    datadir = '/g100_scratch/userexternal/{}/ece3/{}/cmorized/'.format(user, ru)
    filna = datadir+'cmor_*/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/*/r1i1p1f1/{}/{}/g*/v*/{}*nc'
    filist = glob.glob(filna.format(miptab, var, var))

    #filist = glob.glob(filna.format(na, mem, miptab, varnam))
    gigi = xr.open_mfdataset(filist, use_cftime=True)

    try:
        lat = np.array(gigi.lat.data)
    except:
        print('lat name is latitude')
        lat = np.array(gigi.latitude.data)

    seaice = np.array(gigi.siconc.data)
    #okslat = lat > 40.
    for ii, okslat in zip([0,1], [lat > 40, lat < -40]):
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

pickle.dump(resdict, open(cart_out + 'seaicearea.p', 'wb'))

for ru, col, (ye1, ye2) in zip(allruall, colall, okye):
    yeaok = np.arange(ye1, ye2)

    sia_march = resdict[(ru, varnam, 'glomean', 'mar')]
    sia_sept = resdict[(ru, varnam, 'glomean', 'sep')]

    sima10 = ctl.running_mean(sia_march, 10)
    sise10 = ctl.running_mean(sia_sept, 10)

    axs[ii,0].plot(yeaok, sima10, linestyle='solid', marker = 'None', color = col, label = ru, linewidth = 2)
    axs[ii,1].plot(yeaok, sise10, linestyle='solid', marker = 'None', color = col, label = ru, linewidth = 2)

    axs[ii,0].plot(yeaok, sia_march, linestyle='solid', color = col, alpha = 0.3, linewidth = 0.5)
    axs[ii,1].plot(yeaok, sia_sept, linestyle='solid', color = col, alpha = 0.3, linewidth = 0.5)

axs[0,0].set_title('March')
axs[0,1].set_title('September')
axs[0,0].set_ylabel(r'Sea ice extent (m$^2$)')
axs[1,0].set_ylabel(r'Sea ice extent (m$^2$)')
#axs[1,1].legend()
ctl.custom_legend(fig, colall, allruall, ncol = 5)
fig.savefig(cart_out + 'bottseaice.pdf')

sys.exit()

#
# miptab = 'SImon_r1'
# var_map_200 = ['siconc', 'sithick']
#
# mapmean = dict()
# for na, ru, col in zip(allnams, allru, colors):
#     mem = 'r1'
#     filist = np.concatenate([glob.glob(filna.format(na, mem, miptab, varnam)) for varnam in var_map_200])
#     gigi = xr.open_mfdataset(filist, use_cftime=True)
#
#     if ru != 'pi':
#         gigi = gigi.sel(time = gigi['time.year'] >= gigi['time.year'].data[-1]-200)
#
#     for var in var_map_200:
#         print(var)
#         gigi_sclim = ctl.seasonal_climatology(gigi[var], allseasons = ['Mar', 'Sep'])
#         mapmean[(ru, var)] = gigi_sclim
#
# pickle.dump(mapmean, open(cart_out + 'seaicemap.p', 'wb'))

var_map_200 = ['siconc', 'sithick']
mapmean = pickle.load(open(cart_out + 'seaicemap.p', 'rb'))

allcopls = ['seamean', 'seastd', 'seap10', 'seap90']
###### Plots 2D
figs_map = []
for var in var_map_200:
    for copl in allcopls:
        #mappe = [mapmean[('pi', var, copl)]] + [mapmean[(ru, var, copl)]-mapmean[('pi', var, copl)] for ru in allru[1:]]
        #mappe = [mapmean[('pi', var)][copl]] + [mapmean[(ru, var)][copl]-mapmean[('pi', var)][copl] for ru in allru[1:]]
        mappe = [mapmean[(ru, var)][copl] for ru in allru]
        # for ma in mappe:
        #     ma[ma == 0.] = np.nan

        #mappeseas = [ma.sel(time = ma['time.season'] == seasok). for seasok in ['DJF', 'MAM', 'JJA', 'SON'] for ma in mappe]
        mappeseas = [ma.sel(season = seasok) for seasok in ['Mar', 'Sep'] for ma in mappe]
        mapcont = [None if ru == 'pi' else mapmean[('pi', var)][copl].sel(season = seasok) for seasok in ['Mar', 'Sep'] for ru in allru]
        subtitles = ['{} - {}'.format(ru, seasok) for seasok in ['Mar', 'Sep'] for ru in allru]

        cba = None
        cma = 'viridis'
        if var == 'siconc':
            cba = (0,1)
            cma = 'Blues_r'
        #for pol in ['Npolar', 'Spolar']:
        for clat in [90, -90]:
            print(var, copl, clat)
            fig = ctl.plot_multimap_contour(mappeseas, figsize = (24,12), plot_anomalies = False, subtitles = subtitles, title = var+' - '+copl, add_contour_field = mapcont, add_contour_plot_anomalies = False, visualization = 'nearside', cmap = cma, cbar_range = cba, central_lat_lon = (clat, 0), bounding_lat = 40 * np.sign(clat), fix_subplots_shape = (2,4))

            figs_map.append(fig)

figs_map = np.concatenate(figs_map)
fignames = [var+'_'+copl+'_'+pole for var in var_map_200 for copl in allcopls for pole in ['N', 'S']]
ctl.plot_pdfpages(cart_out + 'bottino_seaicemap.pdf', figs_map, True, fignames)
