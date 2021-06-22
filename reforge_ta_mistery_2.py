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

figs_prof = []

flux_diff = flux_hr-flux_lr
#flux_diff_season = flux_diff.groupby("time.season").mean()
for var in allvars:
    print(var)
    # if var == 'zg':
    #     flux_diff_season = flux_diff.groupby("time.season").mean()
    #     vmax = np.nanpercentile(np.abs(flux_diff_season[var]), 98)
    #     fig = plt.figure()
    #     guplo = flux_diff_season[var].mean('lon').plot.contourf(x = 'lat', y ='plev', col = 'season', col_wrap = 2, levels = 11, vmax = vmax, ylim = (1.e5, 1.e3), yscale = 'log')
    #     plt.title(var)
    #     plt.savefig(cart + '{}_seas_799-cntrl.pdf'.format(var))
    #
    # fig = plt.figure()
    # vmax = np.nanpercentile(np.abs(flux_diff[var].sel(plev = 5000., time = slice('1999-01-01', '1999-03-01')).mean('lon')), 98)
    # guplo2 = flux_diff[var].sel(plev = 5000., time = slice('1999-01-01', '1999-03-01')).mean('lon').plot.contourf(x = 'time', y = 'lat', levels = 11, vmax = vmax)
    # plt.title(var)
    # plt.savefig(cart + '{}_ovmol_2mo_799-cntrl.pdf'.format(var))
    # figs_ovmol.append(fig)
    #
    # fig = plt.figure()
    # vmax = np.nanpercentile(np.abs(flux_diff[var].sel(plev = 5000., time = slice('1999-01-01', '1999-01-30'))), 98)
    # guplo2 = flux_diff[var].sel(plev = 5000., time = slice('1999-01-01', '1999-01-30')).plot.contourf(levels = 11, vmax = vmax, col = 'time', col_wrap = 6, transform = proj, figsize = (16,12), subplot_kws = {"projection": proj})
    # guplo2.set_titles(template='{value}', maxchar = 13)
    # guplo2.map(lambda: plt.gca().coastlines())
    # plt.title(var)
    # plt.savefig(cart + '{}_facet_1mo_799-cntrl.pdf'.format(var))
    # figs_facet.append(guplo2.fig)
    #
    # fig = plt.figure()
    # vmax = np.nanpercentile(np.abs(flux_diff[var].mean('lon').sel(time = slice('1999-01-01', '1999-01-30'))), 98)
    # guplo2 = flux_diff[var].sel(time = slice('1999-01-01', '1999-01-30')).mean('lon').plot.contourf(x = 'lat', y = 'plev', col = 'time', col_wrap = 6, levels = 11, ylim = (1.e5, 1.e3), yscale = 'log', vmax = vmax)
    # guplo2.set_titles(template='{value}', maxchar = 13)
    # plt.title(var)
    # plt.savefig(cart + '{}_facet_1mo_3d_799-cntrl.pdf'.format(var))
    # figs_facet_3d.append(guplo2.fig)

    #fig = plt.figure()
    coso_lr = flux_lr[var].sel(time = slice('1999-01-01', '1999-01-30'), lat = slice(-20., 20.)).mean(['lon', 'lat'])
    coso_hr = flux_hr[var].sel(time = slice('1999-01-01', '1999-01-30'), lat = slice(-20., 20.)).mean(['lon', 'lat'])
    # guplo3 = coso_lr.plot(y = 'plev', col = 'time', col_wrap = 6, ylim = (1.e5, 1.e3), yscale = 'log', color = 'forestgreen')
    # guplo2 = coso_hr.plot(ylim = (1.e5, 1.e3), yscale = 'log', color = 'indianred')
    # plt.title(var)
    # plt.savefig(cart + '{}_profs_1mo_799-cntrl.pdf'.format(var))
    # figs_prof.append(fig)

    lrko = coso_lr.data.compute()
    lrko_diff = np.diff(lrko, axis = 0)
    hrko = coso_hr.data.compute()
    hrko_diff = np.diff(hrko, axis = 0)
    levs = coso_hr.plev.data

    cols = ctl.color_set(29, sns_palette = 'crest')
    # plt.figure()
    # for ii, co in zip(range(29), cols):
    #     plt.plot(hrko_diff[ii, :], levs, linestyle = '-', color = co)
    #     plt.plot(lrko_diff[ii, :], levs, linestyle = '--', color = co)
    # plt.yscale('log')
    # plt.ylim(1.e5, 1.e3)
    hrko1 = np.mean(hrko_diff[:10, :], axis = 0)
    hrko2 = np.mean(hrko_diff[10:20, :], axis = 0)
    hrko3 = np.mean(hrko_diff[20:, :], axis = 0)
    lrko1 = np.mean(lrko_diff[:10, :], axis = 0)
    lrko2 = np.mean(lrko_diff[10:20, :], axis = 0)
    lrko3 = np.mean(lrko_diff[20:, :], axis = 0)
    colok = [cols[0], cols[15], cols[-1]]
    fig = plt.figure()
    plt.plot(hrko1, levs, linestyle = '-', color = colok[0], label = 'days 1-10')
    plt.plot(hrko2, levs, linestyle = '-', color = colok[1], label = 'days 10-20')
    plt.plot(hrko3, levs, linestyle = '-', color = colok[2], label = 'days 20-30')
    plt.plot(lrko1, levs, linestyle = '--', color = colok[0])
    plt.plot(lrko2, levs, linestyle = '--', color = colok[1])
    plt.plot(lrko3, levs, linestyle = '--', color = colok[2])
    plt.yscale('log')
    plt.ylim(1.e5, 1.e3)
    plt.grid()
    plt.legend()
    plt.xlabel('daily tendency')
    plt.title(var)
    plt.savefig(cart + '{}_profs_1mo_799-cntrl.pdf'.format(var))
    figs_prof.append(fig)

# ctl.plot_pdfpages(cart + 'day_facet_1m.pdf', figs_facet+figs_facet_3d)
# ctl.plot_pdfpages(cart + 'day_ovmol_2m.pdf', figs_ovmol)
ctl.plot_pdfpages(cart + 'day_prof_1m.pdf', figs_prof)

plt.close('all')
