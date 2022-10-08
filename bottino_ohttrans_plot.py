#!/usr/bin/python
# -*- coding: utf-8 -*-


import numpy as np
import sys
import os
from matplotlib import pyplot as plt
from matplotlib import cm

import pickle
#import netCDF4 as nc

import climtools_lib as ctl
#import climdiags as cd
#from tunlib import gregplot_on_ax

#from matplotlib.colors import LogNorm
#from datetime import datetime

#from scipy import stats
import xarray as xr
import glob
#import xclim

import multiprocessing as mp
import psutil

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 22
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['axes.axisbelow'] = True

#############################################################################

if os.uname()[1] == 'hobbes':
    cart_out = '/home/fabiano/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'xaru':
    cart_out = '/home/fedef/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'tintin':
    cart_out = '/home/fabiano/work/lavori/BOTTINO/'

#cart_out = '/g100_work/IscrB_QUECLIM/BOTTINO/bottino_an/'
cart_out = cart_out + 'ocean3d/'
ctl.mkdir(cart_out)

cart_in = cart_out + '../seasmean/'
gogo = pickle.load(open(cart_in + 'bottino_seasmean_2D.p', 'rb'))
glomeans, pimean, yeamean, _ = gogo

# allru = ['b990', 'b025', 'b050', 'b100']
# allsyear = [1990, 2025, 2050, 2100]
# allnams = ['stabilization-hist-1990', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']
# colors = ['teal', 'forestgreen', 'orange', 'violet']

allru = ['b990', 'b025', 'b050', 'b065', 'b080', 'b100']
colors = ['lightslategray', 'forestgreen', 'orange', 'chocolate', 'maroon', 'violet']

cp0 = 3989.245 # J/kg/K
rho = 1025

oht_all = dict()

read_ts = False

fig, axs = plt.subplots(1, 3, figsize = (18,6))
for ru, col in zip(allru, colors):
    if not read_ts:
        oht_lev = []
        filo = open(cart_out + 'ohtrans_{}.p'.format(ru), 'rb')
        for i in range(500):
            try:
                gigi = pickle.load(filo)
            except:
                break
            oht_lev.append(gigi[0])

        filo.close()

        if ru == 'b065':
            print('AAAAA MISSING 1 YEAR FOR b065!! copying last year twice for now')
            oht_lev.append(gigi[0])

        oht_lev = xr.concat(oht_lev, dim = 'year')*cp0*rho

        oht100 = oht_lev.sel(lev = slice(96., 98.)).mean('lev')
        oht700 = oht_lev.sel(lev = slice(696., 698.)).mean('lev')
        oht2000 = oht_lev.sel(lev = slice(1944., 1946.)).mean('lev')

        oht_all[(ru, 100)] = oht100
        oht_all[(ru, 700)] = oht700
        oht_all[(ru, 2000)] = oht2000

    gtas = glomeans[(ru, 'tas')][1] - pimean['tas']
    yeas = np.arange(500)

    grun = ctl.running_mean(gtas, 20, remove_nans = True)
    for cosu, ax in zip([oht100, oht700, oht2000], axs.flatten()):
        if ru != 'b100':
            ax.scatter(gtas, cosu, s = 5, color = col, label = ru)
            larun = ctl.running_mean(cosu, 20, remove_nans = True)
            ax.plot(grun, larun, color = col, label = ru, lw = 2)
        else:
            ax.scatter(gtas[106:], cosu[:394], s = 5, color = col, label = ru)
            larun = ctl.running_mean(cosu[:394], 20, remove_nans = True)
            ax.plot(grun, larun, color = col, label = ru, lw = 2)

for ax, tit in zip(axs.flatten(), ['100 m', '700 m', '2000 m']):
    ax.grid()
    ax.set_title(tit)
    ax.set_xlabel('GTAS (K)')

axs[0,2].legend()
axs[0,0].set_ylabel('Downward OHT (J/s)')

fig.savefig(cart_out + 'ohtrans_vs_gtas.pdf')

pickle.dump(oht_all, open(cart_out + 'ohtrans_mean.p', 'wb'))

#############################################################

#############################################################
### maps of OHT trends
lats = np.linspace(-89.5, 89.5, 180)
lons = np.linspace(0, 359, 360)

oht_patt = dict()

for ru, col in zip(allru, colors):
    print(ru)
    filo = open(cart_out + 'ohtrans_{}.p'.format(ru), 'rb')

    oht100 = []
    oht700 = []
    oht2000 = []
    for i in range(500):
        try:
            oht_lev, oht100_i, oht700_i, oht2000_i = pickle.load(filo)
        except:
            break
        oht100.append(oht100_i)
        oht700.append(oht700_i)
        oht2000.append(oht2000_i)

    if ru == 'b065':
        print('AAAAA MISSING 1 YEAR FOR b065!! copying last year twice for now')
        oht100.append(oht100_i)
        oht700.append(oht700_i)
        oht2000.append(oht2000_i)

    filo.close()

    oht100 = np.stack(oht100)
    oht700 = np.stack(oht700)
    oht2000 = np.stack(oht2000)

    for var, lab in zip([oht100, oht700, oht2000], [100, 700, 2000]):
        var[var == 0.0] = np.nan

        oht_patt[(ru, lab, 'ini')] = var[:50].mean(axis = 0)*cp0*rho
        oht_patt[(ru, lab, 'fin')] = var[-50:].mean(axis = 0)*cp0*rho
        oht_patt[(ru, lab, 'change')] = oht_patt[(ru, lab, 'fin')] - oht_patt[(ru, lab, 'ini')]

pickle.dump(oht_patt, open(cart_out + 'ohtrans_patt.p', 'wb'))
oht_patt = pickle.load(open(cart_out + 'ohtrans_patt.p', 'rb'))

plpa = []
subt = []
hatch = []
#for lev, tit in zip([700, 2000, 'deep'], ['0-700m', '700-2000m', '> 2000 m']):
for lev, tit in zip([100, 700, 2000], ['100 m', '700 m', '2000 m']):
    for ru in allru:
        plpa.append(oht_patt[(ru, lev, 'fin')])
        subt.append(ru + ': ' + tit)

[fig] = ctl.plot_multimap_contour(plpa, lats, lons, visualization = 'Robinson', central_lat_lon = (0., -120.), filename = cart_out + 'ohtrans_patt_fin.pdf', subtitles = subt, plot_anomalies = True, cmap = ctl.heatmap(), figsize = (16,9), fix_subplots_shape = (3,4), cb_label = 'Downward heat trasport (J/s)')#, add_hatching = hatch, hatch_styles = ['///', '', ''])

for ax in fig.axes[:-1]:
    ax.set_facecolor('gainsboro')
fig.savefig(cart_out + 'ohtrans_patt_fin.pdf')
