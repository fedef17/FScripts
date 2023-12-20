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

#import tunlib as tl

from scipy.optimize import Bounds, minimize, least_squares
import itertools as itt
from multiprocessing import Process, Queue, Pool
import random
#from multiprocessing import set_start_method
#from multiprocessing import get_context
#set_start_method("spawn")

import climtools_lib as ctl
import xarray as xr


plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 18

################################################################################

cart_out = '/home/fabiano/Research/lavori/TunECS/results/amip_vs_ceres/'

exps = ['fml1', 'amc5', 'amc9']
allvars = ['tcc', 'ttr', 'ttrc', 'tsr', 'tsrc']
longnames = ['Cloud cover', 'OLR (W/m2)', 'clear-sky OLR (W/m2)', 'OSR (W/m2)', 'clear-sky OSR (W/m2)']

cart_in = '/data-hobbes/fabiano/TunECS/AMIP_worlds/{}/post/mon/'
finam = '{}_clim00-14_{}_r1.nc'

cose = dict()

for exp in exps:
    for varnam in allvars:
        var, coords, auxv = ctl.read_xr(cart_in.format(exp) + finam.format(exp, varnam))

        mean_field = np.mean(var, axis = 0)
        # zon_mean, zon_std = ctl.zonal_seas_climatology(var, coords['dates'], 'year')
        # glob_mean, glob_std = ctl.global_seas_climatology(var, coords['lat'], coords['dates'], 'year')

        cose[(exp, varnam)] = mean_field
        # cose[(exp, varnam, 'zon_mean')] = zon_mean
        # cose[(exp, varnam, 'zon_std')] = zon_std
        # cose[(exp, varnam, 'glob_mean')] = glob_mean
        # cose[(exp, varnam, 'glob_std')] = glob_std

lat = coords['lat']
lon = coords['lon']

pickle.dump([cose, lat, lon], open(cart_out + 'yearmean_rad_clouds.p', 'wb'))
cose, lat, lon = pickle.load(open(cart_out + 'yearmean_rad_clouds.p', 'rb'))

### Read CERES
cart_cer = '/nas/reference/CERES/'
pino = xr.load_dataset(cart_cer + 'CERES_2000-2015_monclim.nc')
pinostd = xr.load_dataset(cart_cer + 'CERES_2000-2015_monstd.nc')

for exp in exps:
    cose[(exp, 'ttr')] = -cose[(exp, 'ttr')] # cambio segno
    cose[(exp, 'tsr')] = pino['toa_solar_all_mon'].mean('time') - cose[(exp, 'tsr')] # voglio outgoing

    cose[(exp, 'ttrc')] = -cose[(exp, 'ttrc')] # cambio segno
    cose[(exp, 'tsrc')] = pino['toa_solar_all_mon'].mean('time') - cose[(exp, 'tsrc')] # voglio outgoing

    cose[(exp, 'swcre')] = cose[(exp, 'tsr')] - cose[(exp, 'tsrc')]
    cose[(exp, 'lwcre')] = cose[(exp, 'ttr')] - cose[(exp, 'ttrc')]

    cose[(exp, 'tcc')] = 100*cose[(exp, 'tcc')] # %


cernam = dict()
cernam['tcc'] = 'cldarea_total_mon'
cernam['ttr'] = 'toa_lw_all_mon'
cernam['ttrc'] = 'toa_lw_clr_mon'
cernam['tsr'] = 'toa_sw_all_mon'
cernam['tsrc'] = 'toa_sw_clr_mon'

#pino['toa_sw_net_mon'] = pino['toa_sw_all_mon']

biases = dict()
for exp in exps:
    for var in allvars:
        biases[(exp, var)] = cose[(exp, var)] - pino[cernam[var]].mean('time')

    biases[(exp, 'swcre')] = biases[(exp, 'tsr')] - biases[(exp, 'tsrc')]
    biases[(exp, 'lwcre')] = biases[(exp, 'ttr')] - biases[(exp, 'ttrc')]


allvars = ['tcc', 'ttr', 'ttrc', 'tsr', 'tsrc']
lims = 5*[(-15, 15)]
cont_lims = [(10, 90), (130, 290), (140, 320), (50, 170), (20, 170)]
cont_step = [15, 30, 30, 20, 20]
subtitles = ['ctrl bias (vs CERES)', 'amc5-ctrl', 'amc9-ctrl']

#### Figure mean field
for var, clab, lim, clim, cste in zip(allvars, longnames, lims, cont_lims, cont_step):
    # CERES mean state as contour
    # ctrl bias, amc5 change, amc9 change
    # (later) stipling

    filename = cart_out + 'ceres_{}.pdf'.format(var)
    contfield = pino[cernam[var]].mean('time')
    ctl.plot_map_contour(contfield, filename = filename, visualization = 'standard', cmap = 'viridis', title = None, plot_anomalies = False, draw_grid = True, plot_type = 'filled_contour', add_hatching = None, cb_label = clab, figsize = (9,5), cbar_range = clim)#, color_percentiles = (2, 98), cbar_range = cbar_range)

    filename = cart_out + 'tri_{}.pdf'.format(var)

    fields = [biases[('fml1', var)].values, (biases[('amc5', var)]-biases[('fml1', var)]).values, (biases[('amc9', var)]-biases[('fml1', var)]).values]

    ctl.plot_multimap_contour(fields, lat, lon, filename = filename, fix_subplots_shape = (1,3), visualization = 'standard', cmap = 'RdBu_r', title = None, cbar_range = lim, plot_anomalies = True, draw_grid = True, plot_type = 'filled_contour', add_hatching = None, cb_label = clab, subtitles = subtitles, figsize = (18,7))#, add_contour_range = clim, add_contour_lines_step = cste, add_contour_field = 3*[contfield.values], lw_contour = 0.2)

##############################
var = 'swcre'
clab = 'SW CRE (W/m2)'
clim = (-120, 30)
cste = 30
lim = (-10, 10)
filename = cart_out + 'ceres_{}.pdf'.format(var)
contfield = pino[cernam['tsrc']].mean('time')-pino[cernam['tsr']].mean('time')
ctl.plot_map_contour(contfield, filename = filename, visualization = 'standard', cmap = 'viridis', title = None, plot_anomalies = False, draw_grid = True, plot_type = 'filled_contour', add_hatching = None, cb_label = clab, figsize = (9,5))#, color_percentiles = (2, 98), cbar_range = cbar_range)

filename = cart_out + 'tri_{}.pdf'.format(var)

fields = [biases[('fml1', var)].values, (biases[('amc5', var)]-biases[('fml1', var)]).values, (biases[('amc9', var)]-biases[('fml1', var)]).values]

ctl.plot_multimap_contour(fields, lat, lon, filename = filename, fix_subplots_shape = (1,3), visualization = 'standard', cmap = 'RdBu_r', title = None, plot_anomalies = True, draw_grid = True, plot_type = 'filled_contour', add_hatching = None, cb_label = clab, subtitles = subtitles, figsize = (18,7))#, add_contour_range = clim, add_contour_lines_step = cste, add_contour_field = 3*[contfield.values], lw_contour = 0.2)


var = 'lwcre'
clab = 'LW CRE (W/m2)'
clim = (0, 80)
cste = 20
lim = (-10, 10)
filename = cart_out + 'ceres_{}.pdf'.format(var)
contfield = pino[cernam['ttrc']].mean('time')-pino[cernam['ttr']].mean('time')
ctl.plot_map_contour(contfield, filename = filename, visualization = 'standard', cmap = 'viridis', title = None, plot_anomalies = False, draw_grid = True, plot_type = 'filled_contour', add_hatching = None, cb_label = clab, figsize = (9,5))#, color_percentiles = (2, 98), cbar_range = cbar_range)

filename = cart_out + 'tri_{}.pdf'.format(var)

fields = [biases[('fml1', var)].values, (biases[('amc5', var)]-biases[('fml1', var)]).values, (biases[('amc9', var)]-biases[('fml1', var)]).values]

ctl.plot_multimap_contour(fields, lat, lon, filename = filename, fix_subplots_shape = (1,3), visualization = 'standard', cmap = 'RdBu_r', title = None, plot_anomalies = True, draw_grid = True, plot_type = 'filled_contour', add_hatching = None, cb_label = clab, subtitles = subtitles, figsize = (18,7))#, add_contour_range = clim, add_contour_lines_step = cste, add_contour_field = 3*[contfield.values], lw_contour = 0.2)
