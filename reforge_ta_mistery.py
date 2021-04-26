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

# gigi = xr.load_dataset('/home/paolo/work/data/REFORGE/EC-Earth3/rfrg-ctrl-noparam/r1i1p1f1/day/ta/ta_day_EC-Earth3_rfrg-ctrl-noparam_r1i1p1f1_r144x73_19990101-20291231.nc', use_cftime = True)
# gigi2 = xr.load_dataset('/home/paolo/work/data/REFORGE/EC-Earth3-TL799/rfrg-orog255-noparam/r2i1p1f1/day/ta/ta_day_EC-Earth3-TL799_rfrg-orog255-noparam_r2i1p1f1_r144x73_19990101-20231231.nc', use_cftime = True)
# ta_hr = gigi2['ta']
# ta_lr = gigi['ta']
# ta_lr_mean = ta_lr.mean(['time', 'lon'])
# ta_hr_mean = ta_hr.mean(['time', 'lon'])
# clim_diff = ta_hr_mean-ta_lr_mean
#
# ta_diff = ta_hr - ta_lr
# ta_diff_zon = ta_diff.mean('lon')
# #guplo = ta_diff_zon.sel(time = slice('1999-01-01', '1999-01-30')).plot.contourf(x = 'lat', y = 'plev', col = 'time', col_wrap = 6, levels = 11, ylim = (1.e5, 1.e3), yscale = 'log', vmax = 20.)
# #guplo.set_titles(template='{value}', maxchar = 13)
# #plt.savefig(cart + 'ta_diff_799-cntrl_4mo.pdf')
#
# clim_diff_season = ta_diff_zon.groupby("time.season").mean()
# # guplose = clim_diff_season.plot.contourf(x = 'lat', y = 'plev', col = 'season', col_wrap = 2, levels = 11, ylim = (1.e5, 1.e3), yscale = 'log')
# # plt.savefig(cart + 'season_diff_799-cntrl.pdf')
#
# clim_diff_month = ta_diff_zon.groupby("time.month").mean()
# # guplomo = clim_diff_month.plot.contourf(x = 'lat', y = 'plev', col = 'month', col_wrap = 4, levels = 11, ylim = (1.e5, 1.e3), yscale = 'log')
# # plt.savefig(cart + 'month_diff_799-cntrl.pdf')
#
# ta_diff_zon_mon = ta_diff_zon.resample(time = "M").mean()
# jandiffs = ta_diff_zon_mon[ta_diff_zon_mon.groupby("time.month").groups[1]]
# # guplojan = jandiffs.plot.contourf(x = 'lat', y = 'plev', col = 'time', col_wrap = 5, levels = 11, ylim = (1.e5, 1.e3), yscale = 'log', vmax = 20.)
# # guplojan.set_titles(template='{value}', maxchar = 13)
# # plt.savefig(cart + 'Jan_diffs_799-cntrl.pdf')


### RLUT

#srf_net = ssr + str + sshf + slhf
surf_fluxs = ['rsds', 'rlds', 'rsus', 'rlus', 'hfss', 'hfls']
toa_fluxs = ['rlut', 'rsut', 'rsdt']
allvars = surf_fluxs + toa_fluxs

fir_HR = '/home/paolo/work/data/REFORGE/EC-Earth3-TL799/rfrg-orog255-noparam/r2i1p1f1/mon/{}/{}_Amon_EC-Earth3-TL799_rfrg-orog255-noparam_r2i1p1f1_r144x73_*nc'
fir_LR = '/home/paolo/work/data/REFORGE/EC-Earth3/rfrg-ctrl-noparam/r1i1p1f1/mon/{}/{}_Amon_EC-Earth3_rfrg-ctrl-noparam_r1i1p1f1_r144x73_*nc'

fils_HR = np.concatenate([glob.glob(fir_HR.format(var, var)) for var in allvars])
fils_LR = np.concatenate([glob.glob(fir_LR.format(var, var)) for var in allvars])

flux_lr = xr.open_mfdataset(fils_LR, use_cftime = True)
flux_hr = xr.open_mfdataset(fils_HR, use_cftime = True)

flux_hr = flux_hr.assign(net_sfc = flux_hr.rsds + flux_hr.rlds - flux_hr.rsus - flux_hr.rlus - flux_hr.hfss - flux_hr.hfls) # net downward energy flux at surface
flux_hr = flux_hr.assign(net_toa = flux_hr.rsdt - flux_hr.rlut - flux_hr.rsus) # net downward energy flux at TOA
flux_lr = flux_lr.assign(net_sfc = flux_hr.rsds + flux_hr.rlds - flux_hr.rsus - flux_hr.rlus - flux_lr.hfss - flux_lr.hfls) # net downward energy flux at surface
flux_lr = flux_lr.assign(net_toa = flux_lr.rsdt - flux_lr.rlut - flux_lr.rsus) # net downward energy flux at TOA

allvars = allvars + ['net_sfc', 'net_toa']

flux_lr = flux_lr.drop_vars('time_bnds')
flux_hr = flux_hr.drop_vars('time_bnds')

flux_diff = flux_hr-flux_lr
flux_diff_season = flux_diff.groupby("time.season").mean()
for var in allvars:
    vmax = np.nanpercentile(flux_diff_season[var], 98)
    guplo = flux_diff_season[var].plot.contourf(col = 'season', col_wrap = 2, levels = 11, vmax = vmax)
    plt.savefig(cart + '{}_seas_799-cntrl.pdf'.format(var))

    plt.figure()
    guplo2 = flux_diff[var].mean('lon').sel(time = slice('1999-01-01', '1999-12-30')).plot.contourf(x = 'time', y = 'lat', levels = 11, vmax = vmax)
    plt.savefig(cart + '{}_ovmol_1ye_799-cntrl.pdf'.format(var))

    plt.figure()
    guplo2 = flux_diff[var].mean('lon').sel(time = slice('1999-01-01', '1999-04-30')).plot.contourf(x = 'time', y = 'lat', levels = 11, vmax = vmax)
    plt.savefig(cart + '{}_ovmol_4mo_799-cntrl.pdf'.format(var))

    plt.close('all')

# guplo2 = rlut_diff.sel(time = slice('1999-01-01', '1999-01-30')).plot.contourf(col = 'time', col_wrap = 6, levels = 11, vmax = 100)

### TA lat/time
