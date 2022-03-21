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

cart = '/home/fedef/Research/lavori/COSP_for_forum/'
gigi = xr.load_dataset(cart + 'ICMGG_surface_rg.nc')

#ctl.plot_multimap_contour(gigi['var169'], gigi.lat, gigi.lon)
mgri = np.meshgrid(gigi.lat, gigi.lon)

loctime6 = (6 + mgri[1]/360.*24) % 24
loctime12 = (12 + mgri[1]/360.*24) % 24
loctime18 = (18 + mgri[1]/360.*24) % 24
loctime24 = (mgri[1]/360.*24) % 24

loctime = np.stack([loctime6, loctime12, loctime18, loctime24])
loctime = loctime.swapaxes(1,2)

for yea in range(2008, 2017):
    for mo in range(1, 13):
        if mo < 12:
            dates = ctl.date_series(datetime(yea, mo, 1, 0, 0), datetime(yea, mo+1, 1, 0, 0), freq = '6H')[:-1]
        else:
            dates = ctl.date_series(datetime(yea, mo, 1, 0, 0), datetime(yea+1, 1, 1, 0, 0), freq = '6H')[:-1]

        loc_ok = np.vstack(len(dates)//4*[loctime])

        pino = (loc_ok <= 12) & (loc_ok >= 6)
        pino = pino.astype('int8')

        coords = dict([('time', dates), ('lat', gigi.lat), ('lon', gigi.lon)])
        data_vars = dict([('daytime', (('time', 'lat', 'lon'), pino))])

        locxr = xr.Dataset(data_vars=data_vars, coords = coords)
        locxr.to_netcdf(cart + 'daytime_{}{:02d}.nc'.format(yea, mo))

        pino = (loc_ok <= 24) & (loc_ok >= 18)
        pino = pino.astype('int8')

        coords = dict([('time', dates), ('lat', gigi.lat), ('lon', gigi.lon)])
        data_vars = dict([('nighttime', (('time', 'lat', 'lon'), pino))])

        locxr = xr.Dataset(data_vars=data_vars, coords = coords)
        locxr.to_netcdf(cart + 'nighttime_{}{:02d}.nc'.format(yea, mo))
