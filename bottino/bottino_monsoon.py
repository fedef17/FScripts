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

import cartopy.crs as ccrs

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
    cart_out = '/home/fedef/Research/lavori/BOTTINO/monsoon/'
elif os.uname()[1] == 'tintin':
    cart_out = '/home/fabiano/work/lavori/BOTTINO/yearmean/'

ctl.mkdir(cart_out)

filna = '/nas/BOTTINO/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/{}/{}i1p1f1/{}/{}/*nc'

allru = ['pi', 'b025', 'b050', 'b100']
allnams = ['piControl', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']

colors = ['black', 'forestgreen', 'orange', 'violet']

####################################################################################################

miptab = 'Amon'

figs_glob = []
axs_glob = []

seamean = dict()

allseasons = ['NDJFM', 'MJJAS', 'year']

# for na, ru, col in zip(allnams, allru, colors):
#     print(ru)
#     mem = 'r1'
#
#     fils = glob.glob(filna.format(na, mem, miptab, 'pr'))
#
#     kose = xr.open_mfdataset(fils, use_cftime = True)
#     kose = kose.drop_vars('time_bnds')
#
#     for sea in allseasons:
#         seamean[(ru, sea)] = ctl.seasonal_set(kose['pr'], season = sea, seasonal_stat = 'sum')
#
# pickle.dump(seamean, open(cart_out + 'bottino_seamean_pr.p', 'wb'))

# seamean = pickle.load(open(cart_out + 'bottino_seamean_pr.p', 'rb'))
#
# secuse = dict()
# for ru in allru:
#     for sea in allseasons:
#         for iiy in range(5):
#             pino = seamean[(ru, sea)].isel(year = slice(100*iiy, 100*(iiy+1) )).mean('year')
#             secuse[(ru, sea, iiy)] = pino
#
#         pino = seamean[(ru, sea)].isel(year = slice(-200, None)).mean('year')
#         secuse[(ru, sea, 'stab')] = pino
#
#         if ru == 'pi':
#             secuse[(ru, sea, 'stab')] = seamean[(ru, sea)].mean('year')
#
# pickle.dump(secuse, open(cart_out + 'bottino_seamean_pr_secular.p', 'wb'))

secuse = pickle.load(open(cart_out + 'bottino_seamean_pr_secular.p', 'rb'))

### index is MPI = (MJJAS - NDJFM)/annual mean prec
## reversed for southern hemisphere

spy = 86400.*365

mpi = dict()
monreg = dict()

for ru in allru:
    for iiy in list(range(5)) + ['stab']:
        # NORTH
        diff = (secuse[(ru, 'MJJAS', iiy)] - secuse[(ru, 'NDJFM', iiy)])/secuse[(ru, 'year', iiy)]

        coso = xr.merge([diff.sel(lat = slice(0, None)), -1*diff.sel(lat = slice(None, 0))])['sea_sum']

        yemm = secuse[(ru, 'year', iiy)]/12.*spy # from kg m2 s-1 to mm yr-1

        coso.values[yemm < 300.] = np.nan
        mpi[(ru, iiy)] = coso.copy()

        coso.values[coso < 0.5] = np.nan
        coso.values[~np.isnan(coso)] = 1
        coso.values[np.isnan(coso)] = 0
        monreg[(ru, iiy)] = coso

#### Figurez
# only stab
fig, ax = ctl.get_cartopy_fig_ax(coast_lw = 0.1)
proj = ccrs.PlateCarree()

for ru, col in zip(allru, colors):
    coso = monreg[(ru, 'stab')]
    ctl.plot_mapc_on_ax(ax, coso.values, coso.lat.values, coso.lon.values, proj, None, (0,1), plot_type = 'binary_contour', lw_contour=2., line_color = col)

ctl.custom_legend(fig, colors, allru, ncol = 4)
fig.savefig(cart_out + 'monsoon_stab.pdf')

cosi = [mpi[(ru, 'stab')] for ru in allru]
for co in cosi: co.values[co < 0.5] = np.nan

ctl.plot_multimap_contour(cosi, filename = cart_out + 'monsoon_all.pdf',  cbar_range=(0.5, 1.), subtitles = allru, cmap = 'viridis')
