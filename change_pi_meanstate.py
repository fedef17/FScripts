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

from matplotlib.colors import LogNorm
from datetime import datetime

from scipy import stats
import pandas as pd
import glob

import tunlib as tl

from scipy.optimize import Bounds, minimize, least_squares
import itertools as itt
from multiprocessing import Process, Queue, Pool
import random
#from multiprocessing import set_start_method
#from multiprocessing import get_context
#set_start_method("spawn")

import climtools_lib as ctl


plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 18

################################################################################

cart_in = '/home/fabiano/Research/ecearth/TunECS/'
cart_out = '/home/fabiano/Research/lavori/TunECS/results/'

#exps = ['pic0', 'pic5', 'pic6', 'pic8', 'pic9', 'c4c5', 'c4c9']
#exps = ['pic0', 'pic5', 'pic9', 'c4c5', 'c4c9']
exps = ['pic5', 'pic9', 'c4c5', 'c4c9']
scens = ['piControl', 'piControl', 'abrupt-4xCO2', 'abrupt-4xCO2']
colors = ctl.color_set(4)
coldic = dict()
coldic['pic5'] = 'steelblue'
coldic['pic9'] = 'indianred'
coldic['c4c5'] = 'steelblue'
coldic['c4c9'] = 'indianred'

# map of change in pi toa_net, tas and pr
allvars = ['tas', 'pr', 'rlut', 'rsut']

#cart_in = '/scratch/ms/it/ccff/ece3/'
cart_in = '/data-hobbes/fabiano/TunECS/coupled/'
finam = cart_in + '{expname}/cmorized/cmor_{year}/CMIP6/CMIP/EC-Earth-Consortium/EC-Earth3/{scenario}/r1i1p{expname[3]}f1/Amon/{varnam}/gr/v*/{varnam}_Amon_EC-Earth3_{scenario}_r1i1p{expname[3]}f1_gr_{year}01-{year}12.nc'

# cose = dict()
#
# for exp, scen in zip(exps, scens):
#     if scen == 'piControl':
#         year_range = (1881, 1900)
#     else:
#         year_range = (1901, 1920)
#     print(exp, scen, year_range)
#     for varnam in allvars:
#         listafil = []
#         for year in np.arange(year_range[0], year_range[1]+1):
#             finamok = finam.format(varnam = varnam, expname = exp, year = year, scenario = scen)
#             listafilye = glob.glob(finamok)
#             print(finamok, len(listafilye))
#             if len(listafilye) == 1:
#                 listafil.append(listafilye[0])
#             else:
#                 #if not (exp == 'pic0' & year > 1870
#                 raise ValueError('MISSING FILE or too many!')
#         var, coords, aux_info = ctl.read_ensemble_iris(listafil, select_var = varnam)
#
#         mean_field = np.mean(var, axis = 0)
#         zon_mean, zon_std = ctl.zonal_seas_climatology(var, coords['dates'], 'year')
#         glob_mean, glob_std = ctl.global_seas_climatology(var, coords['lat'], coords['dates'], 'year')
#
#         cose[(exp, varnam, 'mean_field')] = mean_field
#         cose[(exp, varnam, 'zon_mean')] = zon_mean
#         cose[(exp, varnam, 'zon_std')] = zon_std
#         cose[(exp, varnam, 'glob_mean')] = glob_mean
#         cose[(exp, varnam, 'glob_std')] = glob_std
#
#
# lat = coords['lat']
# lon = coords['lon']
#
# pickle.dump([cose, lat, lon], open(cart_out + 'mean_state_coupled.p', 'wb'))
cose, lat, lon = pickle.load(open(cart_out + 'mean_state_coupled.p', 'rb'))

cart_fig = cart_out + 'mean_state/'
ctl.mkdir(cart_fig)

cblabels = ['Temp (K)', 'pr', 'rlut (W/m2)', 'rsut (W/m2)']

#### Figure mean field
couples = [('pic9', 'pic5'), ('c4c5', 'pic5'), ('c4c9', 'pic9'), ('c4c9', 'c4c5')]
for varnam, clab in zip(allvars, cblabels):
    if varnam == 'pr':
        cbar_range = (-1e-5, 1e-5)
    else:
        cbar_range = None

    for co in couples:
        field = cose[(co[0], varnam, 'mean_field')]-cose[(co[1], varnam, 'mean_field')]
        filename = cart_fig + 'mean_state_{}_vs_{}_{}.pdf'.format(co[0], co[1], varnam)
        ctl.plot_map_contour(field, lat, lon, filename = filename, visualization = 'standard', cmap = 'RdBu_r', title = None, plot_anomalies = True, draw_grid = True, plot_type = 'filled_contour', add_hatching = None, cb_label = clab, color_percentiles = (2, 98), cbar_range = cbar_range)

        # filename = cart_fig + 'trimean_state_{}_vs_{}_{}.pdf'.format(co[0], co[1], varnam)
        # ctl.plot_triple_sidebyside(cose[(co[0], varnam, 'mean_field')], cose[(co[1], varnam, 'mean_field')], lat, lon, filename = filename, visualization = 'standard', cmap = 'RdBu_r', title = None, plot_anomalies = False, draw_grid = True, plot_type = 'filled_contour', add_hatching = None, cb_label = clab, color_percentiles = (2, 98))

    field = (cose[('c4c9', varnam, 'mean_field')]-cose[('pic9', varnam, 'mean_field')])-(cose[('c4c5', varnam, 'mean_field')]-cose[('pic5', varnam, 'mean_field')])
    filename = cart_fig + 'mean_state_change_9vs5_{}.pdf'.format(varnam)
    ctl.plot_map_contour(field, lat, lon, filename = filename, visualization = 'standard', cmap = 'RdBu_r', title = None, plot_anomalies = True, draw_grid = True, plot_type = 'filled_contour', add_hatching = None, cb_label = clab, color_percentiles = (2, 98), cbar_range = cbar_range)


#### Figure zonal
for varnam in allvars:
    fig = plt.figure(figsize=(24,12))
    ax = plt.subplot(1, 2, 1)
    for exp in ['pic5', 'pic9']:
        field = cose[(exp, varnam, 'zon_mean')]
        err = cose[(exp, varnam, 'zon_std')]

        ax.fill_between(lat, field-err, field+err, color = coldic[exp], alpha = 0.3)
        ax.plot(lat, field, color = coldic[exp], label = exp, linewidth = 2)

    ax.set_xlabel('Latitude')
    ax.set_ylabel(varnam)
    ax.set_title('pre-industrial mean state')
    ax.legend()

    ax = plt.subplot(1, 2, 2)
    for exp in ['pic5', 'pic9']:
        field = cose[('c4'+exp[2:], varnam, 'zon_mean')]-cose[(exp, varnam, 'zon_mean')]
        err = np.mean([cose[('c4'+exp[2:], varnam, 'zon_std')], cose[(exp, varnam, 'zon_std')]], axis = 0)

        ax.fill_between(lat, field-err, field+err, color = coldic[exp], alpha = 0.3)
        ax.plot(lat, field, color = coldic[exp], label = exp[2:], linewidth = 2)

    ax.set_xlabel('Latitude')
    ax.set_ylabel(varnam)
    ax.set_title('change 4xCO2-PI')
    ax.legend()

    filename = cart_fig + 'zon_mean_9vs5_{}.pdf'.format(varnam)
    fig.savefig(filename)

    for exp in exps:
        print(varnam, exp, cose[(exp, varnam, 'glob_mean')])
