#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os

from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import pandas as pd

from numpy import linalg as LA
from eofs.standard import Eof
from scipy import stats, optimize
import itertools as itt

from sklearn.cluster import KMeans
from datetime import datetime
import pickle

import climtools_lib as ctl
#import climdiags as cd

from importlib import reload

import xarray as xr
import xesmf as xe

##############################################################

filon = '/data-hobbes/fabiano/CMIP6/CMIP6/model-output/MRI/MRI-ESM2-0/ssp585/atmos/Amon/r1i1p1f1/rsut/rsut_Amon_MRI-ESM2-0_ssp585_r1i1p1f1_gn_210101-230012.nc'

pino = xr.load_dataset(filon, engine = 'netcdf4')
gigi = pino.rsut
gigi

matdat = gigi.data
matdat.shape

plt.ion()
fig = plt.figure(figsize=(14,6))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_global()
gigi[0].plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), x='lon', y='lat', add_colorbar=True)
fig.savefig('testfig.pdf')

# regrid to 2.5 grid
ds_out = xr.Dataset({'lat': (['lat'], np.arange(-90, 90.1, 2.5)), 'lon': (['lon'], np.arange(0, 360, 2.5)), })
regrid = xe.Regridder(pino, ds_out, 'bilinear')
regrid  # print basic regridder information.

gigi_25 = regrid(gigi)
gigi_25

fig = plt.figure(figsize=(14,6))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_global()
gigi_25[0].plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), x='lon', y='lat', add_colorbar=True)
fig.savefig('testfig2.pdf')
