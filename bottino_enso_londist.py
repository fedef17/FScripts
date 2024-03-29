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

fig_o, axs_o = plt.subplots(2, 2, figsize = (16,12))
fig_a, axs_a = plt.subplots(2, 2, figsize = (16,12))

figli, axsli = plt.subplots(2, 2, figsize = (16,12))

for ru, axo, axa, axl, col in zip(allru, axs_o.flatten(), axs_a.flatten(), axsli.flatten(), colors):
    print(ru)
    piuz = enso[ru]['tos'].groupby('time.year').mean()
    years, gtas = glomeans[(ru, 'tas')]
    # coef3, covmat3 = np.polyfit(years, gtas, deg = 3, cov = True)
    # fitco3 = np.polyval(coef3, years)
    g10 = ctl.butter_filter(gtas, 10)
    yeme = yeamean[(ru, 'tas')]

    if ru == 'pi':
        g10 = g10[:-1]
        yeme = yeme[:-1]

    #tasdetr = yeamean[(ru, 'tas')] - fitco3[:, np.newaxis, np.newaxis]
    tasdetr = yeme - g10[:, np.newaxis, np.newaxis]
    tasanom = tasdetr - tasdetr.mean('year')

    #ctl.plot_map_contour(tasanom[0], central_lat_lon=(0, 205), draw_grid = True, add_rectangles=[160, 270, -5, 5])
    tasanom_pac = tasanom.sel(lat = slice(-5, 5), lon = slice(160, 270)).mean('lat') # 270 for nino3, 240 for nino34

    nino_yrs = piuz > np.percentile(piuz, 90)#0.5
    nina_yrs = piuz < np.percentile(piuz, 10)#-0.5
    oknino = nino_yrs.values.flatten()
    oknina = nina_yrs.values.flatten()

    paclons = tasanom_pac.lon.values

    #tasanom_pac[oknino].plot(ax = axl, x = 'lon', hue = 'year', linewidth = 0.1, color = 'indianred')
    tas90 = np.percentile(tasanom_pac[oknino], 90, axis = 0)
    tas10 = np.percentile(tasanom_pac[oknino], 10, axis = 0)
    axl.fill_between(paclons, tas10, tas90, color = 'indianred', alpha = 0.3)
    tasanom_pac[oknino].mean('year').plot(ax = axl, linewidth = 2, color = 'indianred')

    # tasanom_pac[oknina].plot(ax = axl, x = 'lon', hue = 'year', linewidth = 0.1, color = 'steelblue')
    tas90 = np.percentile(tasanom_pac[oknina], 90, axis = 0)
    tas10 = np.percentile(tasanom_pac[oknina], 10, axis = 0)
    axl.fill_between(paclons, tas10, tas90, color = 'steelblue', alpha = 0.3)
    tasanom_pac[oknina].mean('year').plot(ax = axl, linewidth = 2, color = 'steelblue')

    axl.set_title(ru)
    axl.grid()

    ninodist = np.argmax(tasanom_pac[oknino].values, axis = 1)
    ninadist = np.argmin(tasanom_pac[oknina].values, axis = 1)

    axo.hist(paclons[ninodist], color = col, bins = np.arange(190, 271, 5))
    axa.hist(paclons[ninadist], color = col, bins = np.arange(190, 271, 5))

ctl.adjust_ax_scale(axs_o.flatten())
ctl.adjust_ax_scale(axs_a.flatten())
ctl.adjust_ax_scale(axsli.flatten())

fig_o.savefig(cart_out + 'nino_lon_hist.pdf')
fig_a.savefig(cart_out + 'nina_lon_hist.pdf')
figli.savefig(cart_out + 'nino_lon_anom.pdf')
