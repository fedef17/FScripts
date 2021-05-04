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
# guplo = ta_diff_zon.sel(time = slice('1999-01-01', '1999-01-30')).plot.contourf(x = 'lat', y = 'plev', col = 'time', col_wrap = 6, levels = 11, ylim = (1.e5, 1.e3), yscale = 'log', vmax = 20.)
# #guplo.set_titles(template='{value}', maxchar = 13)
# #plt.savefig(cart + 'ta_diff_799-cntrl_4mo.pdf')
#
# clim_diff_season = ta_diff_zon.groupby("time.season").mean()
# guplose = clim_diff_season.plot.contourf(x = 'lat', y = 'plev', col = 'season', col_wrap = 2, levels = 11, ylim = (1.e5, 1.e3), yscale = 'log')
# # plt.savefig(cart + 'season_diff_799-cntrl.pdf')
#
# clim_diff_month = ta_diff_zon.groupby("time.month").mean()
# guplomo = clim_diff_month.plot.contourf(x = 'lat', y = 'plev', col = 'month', col_wrap = 4, levels = 11, ylim = (1.e5, 1.e3), yscale = 'log')
# # plt.savefig(cart + 'month_diff_799-cntrl.pdf')
#
# ta_diff_zon_mon = ta_diff_zon.resample(time = "M").mean()
# jandiffs = ta_diff_zon_mon[ta_diff_zon_mon.groupby("time.month").groups[1]]
# guplojan = jandiffs.plot.contourf(x = 'lat', y = 'plev', col = 'time', col_wrap = 5, levels = 11, ylim = (1.e5, 1.e3), yscale = 'log', vmax = 20.)
# # guplojan.set_titles(template='{value}', maxchar = 13)
# # plt.savefig(cart + 'Jan_diffs_799-cntrl.pdf')


### RLUT

#srf_net = ssr + str + sshf + slhf
surf_fluxs = ['rsds', 'rlds', 'rsus', 'rlus', 'hfss', 'hfls']
toa_fluxs = ['rlut', 'rsut', 'rsdt']
allvars = surf_fluxs + toa_fluxs + ['clt', 'clwvi', 'evspsbl']

fir_HR = '/home/paolo/work/data/REFORGE/EC-Earth3-TL799/rfrg-orog255-noparam/r2i1p1f1/mon/{}/{}_Amon_EC-Earth3-TL799_rfrg-orog255-noparam_r2i1p1f1_r144x73_*nc'
fir_LR = '/home/paolo/work/data/REFORGE/EC-Earth3/rfrg-ctrl-noparam/r1i1p1f1/mon/{}/{}_Amon_EC-Earth3_rfrg-ctrl-noparam_r1i1p1f1_r144x73_*nc'

fils_HR = np.concatenate([glob.glob(fir_HR.format(var, var)) for var in allvars])
fils_LR = np.concatenate([glob.glob(fir_LR.format(var, var)) for var in allvars])

flux_lr = xr.open_mfdataset(fils_LR, use_cftime = True)
flux_hr = xr.open_mfdataset(fils_HR, use_cftime = True)

flux_hr = flux_hr.assign(net_sfc = flux_hr.rsds + flux_hr.rlds - flux_hr.rsus - flux_hr.rlus - flux_hr.hfss - flux_hr.hfls) # net downward energy flux at surface
flux_hr = flux_hr.assign(net_toa = flux_hr.rsdt - flux_hr.rlut - flux_hr.rsut) # net downward energy flux at TOA
flux_hr = flux_hr.assign(in_atm = flux_hr.net_toa - flux_hr.net_sfc)

flux_lr = flux_lr.assign(net_sfc = flux_lr.rsds + flux_lr.rlds - flux_lr.rsus - flux_lr.rlus - flux_lr.hfss - flux_lr.hfls) # net downward energy flux at surface
flux_lr = flux_lr.assign(net_toa = flux_lr.rsdt - flux_lr.rlut - flux_lr.rsut) # net downward energy flux at TOA
flux_lr = flux_lr.assign(in_atm = flux_lr.net_toa - flux_lr.net_sfc)

allvars = allvars + ['net_sfc', 'net_toa', 'in_atm']

flux_lr = flux_lr.drop_vars('time_bnds')
flux_hr = flux_hr.drop_vars('time_bnds')

figs_clim = []
figs_ovmol = []
figs_ovmol_4m = []
figs_facet_mo = []
figs_facet_mo_3d = []

figs_global = []

flux_diff = flux_hr-flux_lr
flux_diff_season = flux_diff.groupby("time.season").mean()
for var in allvars:
    signc = 1
    if var == 'evspsbl':
        flux_diff[var] = flux_hr[var]+flux_lr[var]
        signc = -1

    print(var)
    vmax = np.nanpercentile(flux_diff_season[var], 98)
    if var in surf_fluxs or var in toa_fluxs or var in ['net_sfc', 'net_toa', 'in_atm']: vmax = 10.
    fig = plt.figure()
    guplo = flux_diff_season[var].plot.contourf(col = 'season', col_wrap = 2, levels = 11, vmax = vmax)
    plt.title(var)
    plt.savefig(cart + '{}_seas_799-cntrl.pdf'.format(var))
    figs_clim.append(guplo.fig)

    vmax = np.nanpercentile(flux_diff[var].mean('lon'), 98)
    if var in surf_fluxs or var in toa_fluxs or var in ['net_sfc', 'net_toa', 'in_atm']: vmax = 10.
    fig = plt.figure()
    guplo2 = flux_diff[var].mean('lon').sel(time = slice('1999-01-01', '1999-12-30')).plot.contourf(x = 'time', y = 'lat', levels = 11, vmax = vmax)
    plt.title(var)
    plt.savefig(cart + '{}_ovmol_1ye_799-cntrl.pdf'.format(var))
    figs_ovmol.append(fig)

    fig = plt.figure()
    guplo2 = flux_diff[var].mean('lon').sel(time = slice('1999-01-01', '1999-04-30')).plot.contourf(x = 'time', y = 'lat', levels = 11, vmax = vmax)
    plt.title(var)
    plt.savefig(cart + '{}_ovmol_4mo_799-cntrl.pdf'.format(var))
    figs_ovmol_4m.append(fig)

    fig = plt.figure()
    vmax = np.nanpercentile(flux_diff[var].sel(time = slice('1999-01-01', '1999-12-30')), 98)
    if var in surf_fluxs or var in toa_fluxs or var in ['net_sfc', 'net_toa', 'in_atm']: vmax = 10.
    guplo2 = flux_diff[var].sel(time = slice('1999-01-01', '1999-12-30')).plot.contourf(levels = 11, vmax = vmax, col = 'time', col_wrap = 4)
    guplo2.set_titles(template='{value}', maxchar = 13)
    plt.title(var)
    plt.savefig(cart + '{}_facet_1ye_799-cntrl.pdf'.format(var))
    figs_facet_mo.append(guplo2.fig)

    ### Add global mean timeseries
    fig = plt.figure()
    coso = flux_lr[var].mean('lon').sel(time = slice('1999-01-01', '1999-12-30'))
    glomean = np.average(coso, weights = abs(np.cos(np.deg2rad(flux_lr.lat))), axis = -1)
    plt.plot_date(coso.time.data, signc*glomean, label = 'LR', color = 'forestgreen')
    coso = flux_hr[var].mean('lon').sel(time = slice('1999-01-01', '1999-12-30'))
    glomean = np.average(coso, weights = abs(np.cos(np.deg2rad(flux_hr.lat))), axis = -1)
    plt.plot_date(coso.time.data, glomean, label = 'HR', color = 'indianred')
    plt.grid()
    plt.title(var)
    plt.legend()
    figs_global.append(fig)

    ### Add global mean timeseries
    fig = plt.figure()
    coso = flux_lr[var].mean('lon').groupby("time.year").mean()
    glomean = np.average(coso, weights = abs(np.cos(np.deg2rad(flux_lr.lat))), axis = -1)
    plt.plot(coso.year.data, signc*glomean, label = 'LR', color = 'forestgreen')
    coso = flux_hr[var].mean('lon').groupby("time.year").mean()
    glomean = np.average(coso, weights = abs(np.cos(np.deg2rad(flux_hr.lat))), axis = -1)
    plt.plot(coso.year.data, glomean, label = 'HR', color = 'indianred')
    plt.grid()
    plt.title(var)
    plt.legend()
    figs_global.append(fig)


ctl.plot_pdfpages(cart + 'month_clim.pdf', figs_clim)
ctl.plot_pdfpages(cart + 'month_ovmol_1e.pdf', figs_ovmol)
ctl.plot_pdfpages(cart + 'month_ovmol_4m.pdf', figs_ovmol_4m)
#plt.close('all')


#sys.exit()

################# 3d month

allvars = ['ta', 'hus', 'wap', 'stream']

fir_HR = '/home/paolo/work/data/REFORGE/EC-Earth3-TL799/rfrg-orog255-noparam/r2i1p1f1/mon/{}/{}_Amon_EC-Earth3-TL799_rfrg-orog255-noparam_r2i1p1f1_r144x73_*nc'
fir_LR = '/home/paolo/work/data/REFORGE/EC-Earth3/rfrg-ctrl-noparam/r1i1p1f1/mon/{}/{}_Amon_EC-Earth3_rfrg-ctrl-noparam_r1i1p1f1_r144x73_*nc'

fils_HR = np.concatenate([glob.glob(fir_HR.format(var, var)) for var in allvars])
fils_LR = np.concatenate([glob.glob(fir_LR.format(var, var)) for var in allvars])

flux_lr = xr.open_mfdataset(fils_LR, use_cftime = True)
flux_hr = xr.open_mfdataset(fils_HR, use_cftime = True)

flux_lr = flux_lr.drop_vars('time_bnds')
flux_hr = flux_hr.drop_vars('time_bnds')

flux_diff = flux_hr-flux_lr
#flux_diff_season = flux_diff.groupby("time.season").mean()
for var in allvars:
    print(var)
    # fig = plt.figure()
    # vmax = np.nanpercentile(flux_diff[var].mean('lon').sel(plev = 5000., time = slice('1999-01-01', '1999-03-01')), 98)
    # guplo2 = flux_diff[var].mean('lon').sel(time = slice('1999-01-01', '1999-03-01')).plot.contourf(x = 'time', y = 'lat', levels = 11, vmax = vmax)
    # plt.title(var)
    # plt.savefig(cart + '{}_ovmol_2mo_799-cntrl.pdf'.format(var))
    # figs_ovmol.append(fig)

    # fig = plt.figure()
    # vmax = np.nanpercentile(np.abs(flux_diff[var].sel(plev = 5000., time = slice('1999-01-01', '1999-12-30'))), 98)
    # guplo2 = flux_diff[var].sel(plev = 5000., time = slice('1999-01-01', '1999-12-30')).plot.contourf(levels = 11, vmax = vmax, col = 'time', col_wrap = 4)
    # guplo2.set_titles(template='{value}', maxchar = 13)
    # plt.title(var)
    # plt.savefig(cart + '{}_facet_1ye_799-cntrl.pdf'.format(var))
    # figs_facet_mo.append(guplo2.fig)
    #
    # fig = plt.figure()
    # vmax = np.nanpercentile(np.abs(flux_diff[var].mean('lon').sel(time = slice('1999-01-01', '1999-12-30'))), 98)
    # guplo2 = flux_diff[var].mean('lon').sel(time = slice('1999-01-01', '1999-12-30')).plot.contourf(x = 'lat', y = 'plev', col = 'time', col_wrap = 4, levels = 11, ylim = (1.e5, 1.e3), yscale = 'log', vmax = vmax)
    # guplo2.set_titles(template='{value}', maxchar = 13)
    # plt.title(var)
    # plt.savefig(cart + '{}_facet_1ye_3d_799-cntrl.pdf'.format(var))
    # figs_facet_mo_3d.append(guplo2.fig)

    ### Add global mean timeseries
    fig = plt.figure()
    coso = flux_lr[var].mean('lon').sel(plev = 5000., time = slice('1999-01-01', '1999-12-30'))
    glomean = np.average(coso, weights = abs(np.cos(np.deg2rad(flux_lr.lat))), axis = -1)
    plt.plot_date(coso.time.data, signc*glomean, label = 'LR', color = 'forestgreen')
    coso = flux_hr[var].mean('lon').sel(time = slice('1999-01-01', '1999-12-30'))
    glomean = np.average(coso, weights = abs(np.cos(np.deg2rad(flux_hr.lat))), axis = -1)
    plt.plot_date(coso.time.data, glomean, label = 'HR', color = 'indianred')
    plt.grid()
    plt.title(var + '- 50 hPa')
    plt.legend()
    figs_global.append(fig)

    ### Add global mean timeseries
    fig = plt.figure()
    coso = flux_lr[var].sel(plev = 5000.).mean('lon').groupby("time.year").mean()
    glomean = np.average(coso, weights = abs(np.cos(np.deg2rad(flux_lr.lat))), axis = -1)
    plt.plot(coso.year.data, signc*glomean, label = 'LR', color = 'forestgreen')
    coso = flux_hr[var].mean('lon').groupby("time.year").mean()
    glomean = np.average(coso, weights = abs(np.cos(np.deg2rad(flux_hr.lat))), axis = -1)
    plt.plot(coso.year.data, glomean, label = 'HR', color = 'indianred')
    plt.grid()
    plt.title(var + '- 50 hPa')
    plt.legend()
    figs_global.append(fig)

#ctl.plot_pdfpages(cart + 'month_facet_1y.pdf', figs_facet_mo+figs_facet_mo_3d)
#ctl.plot_pdfpages(cart + 'day_ovmol_2m.pdf', figs_ovmol)

ctl.plot_pdfpages(cart + 'global_timeseries.pdf', figs_global)
plt.close('all')


# guplo2 = rlut_diff.sel(time = slice('1999-01-01', '1999-01-30')).plot.contourf(col = 'time', col_wrap = 6, levels = 11, vmax = 100)

### TA lat/time

##################### daily
figs_global = []
figs_facet = []
figs_ovmol = []
allvars = ['rsds', 'rlds', 'rsus', 'rlus', 'rlut', 'clt', 'sfcWindmax']

fir_HR = '/home/paolo/work/data/REFORGE/EC-Earth3-TL799/rfrg-orog255-noparam/r2i1p1f1/day/{}/{}_day_EC-Earth3-TL799_rfrg-orog255-noparam_r2i1p1f1_r144x73_*nc'
fir_LR = '/home/paolo/work/data/REFORGE/EC-Earth3/rfrg-ctrl-noparam/r1i1p1f1/day/{}/{}_day_EC-Earth3_rfrg-ctrl-noparam_r1i1p1f1_r144x73_*nc'

fils_HR = np.concatenate([glob.glob(fir_HR.format(var, var)) for var in allvars])
fils_LR = np.concatenate([glob.glob(fir_LR.format(var, var)) for var in allvars])

flux_lr = xr.open_mfdataset(fils_LR, use_cftime = True)
flux_hr = xr.open_mfdataset(fils_HR, use_cftime = True)

flux_hr = flux_hr.assign(in_sfc = flux_hr.rsds + flux_hr.rlds - flux_hr.rsus - flux_hr.rlus) # downward radiative flux at surface
flux_hr = flux_hr.assign(in_toa = -flux_hr.rlut - flux_hr.rsus + flux_hr.rsds) # approx. net downward energy flux at TOA
flux_lr = flux_lr.assign(in_sfc = flux_lr.rsds + flux_lr.rlds - flux_lr.rsus - flux_lr.rlus)
flux_lr = flux_lr.assign(in_toa = -flux_lr.rlut - flux_lr.rsus + flux_lr.rsds)

allvars = allvars + ['in_sfc', 'in_toa']

flux_lr = flux_lr.drop_vars('time_bnds')
flux_hr = flux_hr.drop_vars('time_bnds')

figs_global = []
flux_diff = flux_hr-flux_lr
#flux_diff_season = flux_diff.groupby("time.season").mean()
for var in allvars:
    print(var)
    # fig = plt.figure()
    # vmax = np.nanpercentile(np.abs(flux_diff[var].mean('lon').sel(time = slice('1999-01-01', '1999-03-01'))), 98)
    # guplo2 = flux_diff[var].mean('lon').sel(time = slice('1999-01-01', '1999-03-01')).plot.contourf(x = 'time', y = 'lat', levels = 11, vmax = vmax)
    # plt.title(var)
    # plt.savefig(cart + '{}_ovmol_2mo_799-cntrl.pdf'.format(var))
    # figs_ovmol.append(fig)
    #
    # fig = plt.figure()
    # vmax = np.nanpercentile(np.abs(flux_diff[var].sel(time = slice('1999-01-01', '1999-01-30'))), 98)
    # guplo2 = flux_diff[var].sel(time = slice('1999-01-01', '1999-01-30')).plot.contourf(levels = 11, vmax = vmax, col = 'time', col_wrap = 6)
    # guplo2.set_titles(template='{value}', maxchar = 13)
    # plt.title(var)
    # plt.savefig(cart + '{}_facet_1mo_799-cntrl.pdf'.format(var))
    # figs_facet.append(guplo2.fig)

    ### Add global mean timeseries
    fig = plt.figure()
    coso = flux_lr[var].mean('lon').sel(time = slice('1999-01-01', '1999-04-01'))
    glomean = np.average(coso, weights = abs(np.cos(np.deg2rad(flux_lr.lat))), axis = -1)
    plt.plot_date(coso.time.data, signc*glomean, label = 'LR', color = 'forestgreen')
    coso = flux_hr[var].mean('lon').sel(time = slice('1999-01-01', '1999-04-01'))
    glomean = np.average(coso, weights = abs(np.cos(np.deg2rad(flux_hr.lat))), axis = -1)
    plt.plot_date(coso.time.data, glomean, label = 'HR', color = 'indianred')
    plt.grid()
    plt.title(var)
    plt.legend()
    figs_global.append(fig)

ctl.plot_pdfpages(cart + 'global_timeseries_daily.pdf', figs_global)

sys.exit()
##################### daily 3D
#figs_facet = []
#figs_ovmol = []
figs_facet_3d = []
allvars = ['ta', 'zg']

fir_HR = '/home/paolo/work/data/REFORGE/EC-Earth3-TL799/rfrg-orog255-noparam/r2i1p1f1/day/{}/{}_day_EC-Earth3-TL799_rfrg-orog255-noparam_r2i1p1f1_r144x73_*nc'
fir_LR = '/home/paolo/work/data/REFORGE/EC-Earth3/rfrg-ctrl-noparam/r1i1p1f1/day/{}/{}_day_EC-Earth3_rfrg-ctrl-noparam_r1i1p1f1_r144x73_*nc'

fils_HR = np.concatenate([glob.glob(fir_HR.format(var, var)) for var in allvars])
fils_LR = np.concatenate([glob.glob(fir_LR.format(var, var)) for var in allvars])

flux_lr = xr.open_mfdataset(fils_LR, use_cftime = True)
flux_hr = xr.open_mfdataset(fils_HR, use_cftime = True)

flux_lr = flux_lr.drop_vars('time_bnds')
flux_hr = flux_hr.drop_vars('time_bnds')

flux_diff = flux_hr-flux_lr
#flux_diff_season = flux_diff.groupby("time.season").mean()
for var in allvars:
    print(var)
    if var == 'zg':
        flux_diff_season = flux_diff.groupby("time.season").mean()
        vmax = np.nanpercentile(np.abs(flux_diff_season[var]), 98)
        fig = plt.figure()
        guplo = flux_diff_season[var].mean('lon').plot.contourf(x = 'lat', y ='plev', col = 'season', col_wrap = 2, levels = 11, vmax = vmax, ylim = (1.e5, 1.e3), yscale = 'log')
        plt.title(var)
        plt.savefig(cart + '{}_seas_799-cntrl.pdf'.format(var))

    fig = plt.figure()
    vmax = np.nanpercentile(np.abs(flux_diff[var].mean('lon').sel(plev = 5000., time = slice('1999-01-01', '1999-03-01'))), 98)
    guplo2 = flux_diff[var].mean('lon').sel(plev = 5000., time = slice('1999-01-01', '1999-03-01')).plot.contourf(x = 'time', y = 'lat', levels = 11, vmax = vmax)
    plt.title(var)
    plt.savefig(cart + '{}_ovmol_2mo_799-cntrl.pdf'.format(var))
    figs_ovmol.append(fig)

    fig = plt.figure()
    vmax = np.nanpercentile(np.abs(flux_diff[var].sel(plev = 5000., time = slice('1999-01-01', '1999-01-30'))), 98)
    guplo2 = flux_diff[var].sel(plev = 5000., time = slice('1999-01-01', '1999-01-30')).plot.contourf(levels = 11, vmax = vmax, col = 'time', col_wrap = 6)
    guplo2.set_titles(template='{value}', maxchar = 13)
    plt.title(var)
    plt.savefig(cart + '{}_facet_1mo_799-cntrl.pdf'.format(var))
    figs_facet.append(guplo2.fig)

    fig = plt.figure()
    vmax = np.nanpercentile(np.abs(flux_diff[var].mean('lon').sel(time = slice('1999-01-01', '1999-01-30'))), 98)
    guplo2 = flux_diff[var].mean('lon').sel(time = slice('1999-01-01', '1999-01-30')).plot.contourf(x = 'lat', y = 'plev', col = 'time', col_wrap = 6, levels = 11, ylim = (1.e5, 1.e3), yscale = 'log', vmax = vmax)
    guplo2.set_titles(template='{value}', maxchar = 13)
    plt.title(var)
    plt.savefig(cart + '{}_facet_1mo_3d_799-cntrl.pdf'.format(var))
    figs_facet_3d.append(guplo2.fig)

ctl.plot_pdfpages(cart + 'day_facet_1m.pdf', figs_facet+figs_facet_3d)
ctl.plot_pdfpages(cart + 'day_ovmol_2m.pdf', figs_ovmol)

plt.close('all')
