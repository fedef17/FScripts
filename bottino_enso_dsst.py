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

if os.uname()[1] == 'hobbes':
    cart = '/home/fabiano/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'xaru':
    cart = '/home/fedef/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'tintin':
    cart = '/home/fabiano/work/lavori/BOTTINO/'

cart_in = cart + 'seasmean/'
cart_out = cart + 'indices/enso/'
ctl.mkdir(cart_out)

allru = ['pi', 'b025', 'b050', 'b100']
allnams = ['piControl', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']

colors = ['black', 'forestgreen', 'orange', 'violet']

allnams2 = allnams + ['ssp585', 'historical']
allru2 = allru + ['ssp585', 'hist']
colors2 = colors + ['indianred', 'steelblue']
####################################################################################################

glomeans, pimean, yeamean, mapmean = pickle.load(open(cart_in + 'bottino_seasmean_2D.p', 'rb'))

cartind = '/nas/BOTTINO/indices/enso500_xr/'
enso = dict()

for ru in allru:
    if ru == 'pi':
        enso[ru] = xr.load_dataset(cartind + '{}_enso_360day.nc'.format(ru), use_cftime = True)
    else:
        enso[ru] = xr.load_dataset(cartind + '{}_enso_360day.nc'.format(ru), use_cftime = True)

        firstye = enso[ru].time.values[0].year
        lastye = enso[ru].time.values[-1].year
        enso[ru+'_st'] = enso[ru].sel(time = slice('{}-01-01'.format(lastye-200), '{}-12-30'.format(lastye)))
        enso[ru+'_tr'] = enso[ru].sel(time = slice('{}-01-01'.format(firstye), '{}-12-30'.format(firstye+50)))

#########################################
#
# fig_o, axs_o = plt.subplots(2, 2, figsize = (16,12))
# fig_a, axs_a = plt.subplots(2, 2, figsize = (16,12))
#
# figli, axsli = plt.subplots(2, 2, figsize = (16,12))

dsst_all = dict()

for ru, col in zip(allru, colors):
    print(ru)
    piuz = enso[ru]['tos'].groupby('time.year').mean()
    years, gtas = glomeans[(ru, 'tas')]
    # coef3, covmat3 = np.polyfit(years, gtas, deg = 3, cov = True)
    # fitco3 = np.polyval(coef3, years)
    #g10 = ctl.butter_filter(gtas, 10)
    yeme = yeamean[(ru, 'tas')]

    if ru == 'pi':
        g10 = g10[:-1]
        yeme = yeme[:-1]

    dsst = yeme.sel(lat = slice(-5, 5), lon = slice(100, 180)).mean(['lat', 'lon']) - yeme.sel(lat = slice(-5, 5), lon = slice(200, 280)).mean(['lat', 'lon'])

    dsst_all[ru] = dsst

pickle.dump(dsst_all, open(cart_out + 'dsst_all.p', 'wb'))
