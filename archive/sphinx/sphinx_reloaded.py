#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os

from matplotlib import pyplot as plt
import netCDF4 as nc
import pandas as pd
from numpy import linalg as LA
from scipy import stats
from scipy import interpolate as itrp
import itertools as itt

from datetime import datetime
import pickle

import climtools_lib as ctl
import climdiags as cd

from copy import deepcopy as cp
###############################################

# Leggo un po' di variabili e per tutte faccio:
# - zonal mean (prima e dopo il salto + diff): 2085 - 2095, 2115 - 2125
# - lo stesso con la weighted zonal mean

cart = '/data-woodstock/SPHINX/extended/'
cart = '/data-woodstock/fede/'

cart_out = '/home/fabiano/Research/lavori/SPHINX_for_lisboa/sphinx_reloaded/'

varlist = ['hcc', 'mcc', 'lcc', 'tas', 'rlus', 'rsus', 'rlds', 'rsds', 'hfls', 'hfss', 'rlut', 'rsut', 'rsdt']
#filename = '1850-2160-{}-monthly.nc'
filename = '1850-2160-{}-monthly_remap25_rechunk.nc'

ens_mem = ['lcb0', 'lcs0', 'lcb1', 'lcs1', 'lcb2', 'lcs2']

figures = []
lat, lon = ctl.genlatlon(73, 144)

fields = dict()

for varnam in varlist:
    print(varnam)
    ifile = cart + filename.format(varnam)
    if varnam == 'tas':
        ifile = cart + '1950-2160-tas-monthly_remap25_rechunk.nc'
    var, lat, lon, dates, time_units, var_units, time_cal = ctl.readxDncfield(ifile)
    print(var.keys())
    weights = abs(np.cos(np.deg2rad(lat)))
    # var, coords, aux_info = ctl.read_iris_nc(ifile)
    # lat = coords['lat']
    # lon = coords['lon']
    # dates = coords['dates']

    for ens in ens_mem:
        # zonal mean
        for season in ['year', 'DJF', 'JJA']:
            var_pre, dates_pre = ctl.sel_time_range(var[ens], dates, ctl.range_years(2085, 2095))

            var_pre_seas, _ = ctl.seasonal_climatology(var_pre, dates_pre, season)
            var_zon = ctl.zonal_mean(var_pre_seas)

            var_zon_wcos = weights * var_zon

            fields[(varnam, ens, 'pre', season, 'zonal')] = var_zon
            fields[(varnam, ens, 'pre', season, 'zonal_wcos')] = var_zon_wcos

            var_pre, dates_pre = ctl.sel_time_range(var[ens], dates, ctl.range_years(2110, 2120))
            var_pre_seas, _ = ctl.seasonal_climatology(var_pre, dates_pre, season)
            var_zon = ctl.zonal_mean(var_pre_seas)

            var_zon_wcos = weights * var_zon

            fields[(varnam,  ens,'post', season, 'zonal')] = var_zon
            fields[(varnam,  ens,'post', season, 'zonal_wcos')] = var_zon_wcos

            fields[(varnam,  ens,'diff', season, 'zonal')] = fields[(varnam,  ens,'post', season, 'zonal')] - fields[(varnam,  ens,'pre', season, 'zonal')]
            fields[(varnam,  ens,'diff', season, 'zonal_wcos')] = fields[(varnam,  ens,'post', season, 'zonal_wcos')] - fields[(varnam,  ens,'pre', season, 'zonal_wcos')]

#### figures
# for varnam in varlist:
#     for season in ['year', 'DJF', 'JJA']:
#         for cos in ['pre', 'post', 'diff']:
#             for ens in ens_mem:
#                 fields[(varnam, ens, cos, season, 'zonal_wcos')] = weights * fields[(varnam, ens, cos, season, 'zonal')]

for ens in ens_mem:
    for season in ['year', 'DJF', 'JJA']:
        for cos in ['pre', 'post', 'diff']:
            for zop in ['zonal', 'zonal_wcos']:
                fields[('surf_rad_bal', ens, cos, season, zop)] = fields[('rsds', ens, cos, season, zop)] + fields[('rlds', ens, cos, season, zop)] - fields[('rsus', ens, cos, season, zop)] - fields[('rlus', ens, cos, season, zop)]
                fields[('surf_tot_bal', ens, cos, season, zop)] = fields[('surf_rad_bal', ens, cos, season, zop)] - fields[('hfss', ens, cos, season, zop)] - fields[('hfls', ens, cos, season, zop)]

###
for varnam in varlist + ['surf_rad_bal', 'surf_tot_bal']:
    for season in ['year', 'DJF', 'JJA']:
        for cos in ['pre', 'post', 'diff']:
            for zop in ['zonal', 'zonal_wcos']:
                fields[(varnam, 'base_mean', cos, season, zop)] = np.mean([fields[(varnam, ens, cos, season, zop)] for ens in ens_mem if 'lcb' in ens], axis = 0)
                fields[(varnam, 'stoc_mean', cos, season, zop)] = np.mean([fields[(varnam, ens, cos, season, zop)] for ens in ens_mem if 'lcs' in ens], axis = 0)
                fields[(varnam, 'base_std', cos, season, zop)] = np.std([fields[(varnam, ens, cos, season, zop)] for ens in ens_mem if 'lcb' in ens], axis = 0)
                fields[(varnam, 'stoc_std', cos, season, zop)] = np.std([fields[(varnam, ens, cos, season, zop)] for ens in ens_mem if 'lcs' in ens], axis = 0)

pickle.dump(fields, open(cart_out + 'sphinx_reloaded_2.p', 'wb'))

fields = pickle.load(open(cart_out + 'sphinx_reloaded_2.p'))

figures_dict = dict()
colors = ctl.color_set(3)
colors = ['indianred', 'steelblue', 'forestgreen']
for season in ['year', 'DJF', 'JJA']:
    figures = []
    for varnam in varlist + ['surf_rad_bal', 'surf_tot_bal']:
        fig = plt.figure()
        plt.title('{} - {}'.format(varnam, season))
        for cos, col in zip(['base', 'stoc'], colors):
            plt.plot(lat, fields[(varnam, cos + '_mean', 'pre', season, 'zonal')], label = cos+' (pre)', color = col, linestyle = ':')
            plt.plot(lat, fields[(varnam, cos + '_mean', 'post', season, 'zonal')], label = cos+' (post)', color = col, linestyle = '-')
        mea = 10*(fields[(varnam, 'stoc_mean', 'pre', season, 'zonal')]-fields[(varnam, 'base_mean', 'pre', season, 'zonal')])
        plt.plot(lat, mea, label = '(stoc-base) x 10 (pre)', color = colors[2], linestyle = ':')
        mea = 10*(fields[(varnam, 'stoc_mean', 'post', season, 'zonal')]-fields[(varnam, 'base_mean', 'post', season, 'zonal')])
        plt.plot(lat, mea, label = '(stoc-base) x 10 (post)', color = colors[2], linestyle = '-')
        plt.axhline(0., color = 'grey')
        plt.legend(fontsize = 'small')
        figures.append(fig)
        figures_dict[(varnam, season, cos)] = fig

        fig = plt.figure()
        plt.title('{} - base vs stoc wdiff - {}'.format(varnam, season))
        for i, cos in enumerate(['base', 'stoc']):
            mea = fields[(varnam, cos + '_mean', 'diff', season, 'zonal_wcos')]
            std = fields[(varnam, cos + '_std', 'diff', season, 'zonal_wcos')]
            plt.fill_between(lat, mea - std, mea + std, color = colors[i], alpha = 0.2)
            plt.plot(lat, mea, label = 'w. diff post-pre - {}'.format(cos), color = colors[i])
        plt.axhline(0., color = 'grey')
        plt.legend()
        figures.append(fig)
        figures_dict[(varnam, season, 'wdiff')] = fig

    filefig = cart_out + 'figures_SPHINX_reloaded_{}.pdf'.format(season)
    ctl.plot_pdfpages(filefig, figures, save_single_figs = True)

# for season in ['year', 'DJF', 'JJA']:
#     figures = []
#     for varnam in varlist + ['surf_rad_bal', 'surf_tot_bal']:
#         fig = plt.figure()
#         plt.title('{} - {}'.format(varnam, season))
#         for cos, col in zip(['base', 'stoc'], colors):
#             plt.plot(lat, fields[(varnam, cos + '_mean', 'pre', season, 'zonal')], label = cos+' (pre)', color = col, linestyle = ':')
#             plt.plot(lat, fields[(varnam, cos + '_mean', 'post', season, 'zonal')], label = cos+' (post)', color = col, linestyle = '-')
