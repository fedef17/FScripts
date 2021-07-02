#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os

from matplotlib import pyplot as plt
import matplotlib.cm as cm
import matplotlib.patheffects as PathEffects
import matplotlib.animation as animation
from matplotlib.animation import ImageMagickFileWriter

import netCDF4 as nc
import cartopy.crs as ccrs
import cartopy.util as cutil
import pandas as pd

from numpy import linalg as LA
from eofs.standard import Eof
from scipy import stats, optimize
import itertools as itt

from sklearn.cluster import KMeans

from datetime import datetime
import pickle
import iris

import glob

import climtools_lib as ctl
import climdiags as cd

import xarray as xr
import xesmf as xe
import xclim

from importlib import reload

#######################################################

cart = '/home/federico/work/reforge/'

plevels = np.array([100000, 85000, 70000, 50000, 25000, 10000, 5000, 1000])

##################### daily 3D
#figs_facet = []
#figs_ovmol = []
figs_facet_3d = []
allvars = ['ta', 'ua', 'va', 'wap']

#fir_HR = '/home/paolo/work/data/REFORGE/EC-Earth3-TL799/rfrg-orog255-noparam/r2i1p1f1/day/{}/{}_day_EC-Earth3-TL799_rfrg-orog255-noparam_r2i1p1f1_r144x73_*nc' ### the 799 wap starts from 2000 instead than 1999
fir_HR = '/home/paolo/work/data/REFORGE/EC-Earth3-TL511/rfrg-orog255-noparam/r2i1p1f1/day/{}/{}_day_EC-Earth3-TL511_rfrg-orog255-noparam_r2i1p1f1_r144x73_*nc'
fir_LR = '/home/paolo/work/data/REFORGE/EC-Earth3/rfrg-ctrl-noparam/r1i1p1f1/day/{}/{}_day_EC-Earth3_rfrg-ctrl-noparam_r1i1p1f1_r144x73_*nc'

fils_HR = np.concatenate([glob.glob(fir_HR.format(var, var)) for var in allvars])
fils_LR = np.concatenate([glob.glob(fir_LR.format(var, var)) for var in allvars])

flux_lr = xr.open_mfdataset(fils_LR, use_cftime = True)
flux_lr = flux_lr.drop_vars('time_bnds')
flux_lr = flux_lr.sel(time = slice('1999-01-01', '1999-04-01'))

flux_hr = xr.open_mfdataset(fils_HR, use_cftime = True)
flux_hr = flux_hr.drop_vars('time_bnds')
flux_hr = flux_hr.sel(time = slice('1999-01-01', '1999-04-01'))

#### Now the calculation
coslat = np.cos(np.deg2rad(flux_lr['lat'].values))
Re = 6371.e3 # Mean Earth radius. Correction for altitude level would mean about 30 km difference for the uppermost level (10 hPa)
ka = 0.286
cp = 1000.

dTdt_lr = np.gradient(flux_lr['ta'].values, axis = 0)
ugrad_lr = flux_lr['ua']*1/(Re * coslat[np.newaxis, np.newaxis, :, np.newaxis])*np.gradient(flux_lr['ta'].values, axis = -1)
vgrad_lr = flux_lr['va']*1/Re*np.gradient(flux_lr['ta'].values, axis = -2)
vertgr_lr = flux_lr['wap'] * (ka * flux_lr['ta'].values/plevels[np.newaxis, :, np.newaxis, np.newaxis] - np.gradient(flux_lr['ta'].values, plevels, axis = 1))

Q_lr = dTdt_lr + ugrad_lr + vgrad_lr - vertgr_lr

dTdt_hr = np.gradient(flux_hr['ta'].values, axis = 0)
ugrad_hr = flux_hr['ua']*1/(Re * coslat[np.newaxis, np.newaxis, :, np.newaxis])*np.gradient(flux_hr['ta'].values, axis = -1)
vgrad_hr = flux_hr['va']*1/Re*np.gradient(flux_hr['ta'].values, axis = -2)
vertgr_hr = flux_hr['wap'] * (ka * flux_hr['ta'].values/plevels[np.newaxis, :, np.newaxis, np.newaxis] - np.gradient(flux_hr['ta'].values, plevels, axis = 1))

Q_hr = dTdt_hr + ugrad_hr + vgrad_hr - vertgr_hr

Q_diff = Q_hr-Q_lr

fig = plt.figure()
pino = Q_diff.sel(time = slice('19990101','19990201'), lat = slice(-40, 40)).mean(['time','lon']).plot.contourf(y = 'plev', ylim = (1.e5, 1.e3), yscale = 'log', levels = 11, vmax = 0.25)
pino.colorbar.set_label('Diab. heating (K/day)')
fig.savefig(cart + 'diabheat_2mo_511orog255-ctrl255.pdf')

fig = plt.figure()
pino = Q_hr.sel(time = slice('19990101','19990201'), lat = slice(-40, 40)).mean(['time','lon']).plot.contourf(y = 'plev', ylim = (1.e5, 1.e3), yscale = 'log', levels = 11, vmax = 0.25)
pino.colorbar.set_label('Diab. heating (K/day)')
fig.savefig(cart + 'diabheat_2mo_511orog255.pdf')

fig = plt.figure()
pino = Q_lr.sel(time = slice('19990101','19990201'), lat = slice(-40, 40)).mean(['time','lon']).plot.contourf(y = 'plev', ylim = (1.e5, 1.e3), yscale = 'log', levels = 11, vmax = 0.25)
pino.colorbar.set_label('Diab. heating (K/day)')
fig.savefig(cart + 'diabheat_2mo_ctrl255.pdf')

## Just the temperature trends
fig = plt.figure()
dTdt_diff = dTdt_hr - dTdt_lr
pino = dTdt_diff.sel(time = slice('19990101','19990201'), lat = slice(-40, 40)).mean(['time','lon']).plot.contourf(y = 'plev', ylim = (1.e5, 1.e3), yscale = 'log', levels = 11, vmax = 0.25)
pino.colorbar.set_label('Tot heating (K/day)')
fig.savefig(cart + 'Temptend_2mo_511orog255-ctrl255.pdf')

fig = plt.figure()
pino = dTdt_hr.sel(time = slice('19990101','19990201'), lat = slice(-40, 40)).mean(['time','lon']).plot.contourf(y = 'plev', ylim = (1.e5, 1.e3), yscale = 'log', levels = 11, vmax = 0.25)
pino.colorbar.set_label('Tot heating (K/day)')
fig.savefig(cart + 'Temptend_2mo_511orog255.pdf')

fig = plt.figure()
pino = dTdt_lr.sel(time = slice('19990101','19990201'), lat = slice(-40, 40)).mean(['time','lon']).plot.contourf(y = 'plev', ylim = (1.e5, 1.e3), yscale = 'log', levels = 11, vmax = 0.25)
pino.colorbar.set_label('Tot heating (K/day)')
fig.savefig(cart + 'Temptend_2mo_ctrl255.pdf')
#.plot(y = 'plev', col = 'time', col_wrap = 6, ylim = (1.e5, 1.e3), yscale = 'log')
