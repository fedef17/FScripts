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
import xarray as xr
import glob
import xclim

import cartopy.crs as ccrs

plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
titlefont = 22
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.axisbelow'] = True
plt.rcParams['legend.fontsize'] = 18

#############################################################################

if os.uname()[1] == 'hobbes':
    cart_out = '/home/fabiano/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'xaru':
    cart_out = '/home/fedef/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'tintin':
    cart_out = '/home/fabiano/work/lavori/BOTTINO/'

colors = ['black', 'royalblue', 'lightslategray', 'forestgreen', 'orange', 'chocolate', 'maroon', 'violet', 'crimson']
allru = ['pi', 'hist', 'b990', 'b025', 'b050', 'b065', 'b080', 'b100', 'ssp585']

####################################################################################################

cart_out = cart_out + 'simple_maps/'

yeamean = pickle.load(open(cart_out + '../seasmean/bottino_seasmean_2D.p', 'rb'))
coso = yeamean[2][('b025', 'tas')]
del yeamean

anom_maps, patt_maps = pickle.load(open(cart_out + 'all_maps.p', 'rb'))

# map of temp trends
ziup = patt_maps[('tas', 'b100', 'stab')]
piuk = patt_maps[('tas', 'ssp585', 'fin')]

ctl.plot_map_contour(ziup/piuk, coso.lat, coso.lon, visualization = 'Robinson', cbar_range = (0, 2), cmap = ctl.heatmap(), filename = cart_out + 'taspatt_ssp_vs_b100_ratio.pdf')


### map of prec trends
ziup = patt_maps[('pr_rel', 'b100', 'stab')]
piuk = patt_maps[('pr_rel', 'ssp585', 'fin')]

fig = plt.figure()
plt.scatter(piuk.flatten(), ziup.flatten(), color = 'black', s = 1)
plt.grid()
piuk[np.abs(piuk) < 1] = np.nan

diago = lambda x : x
plt.plot(np.linspace(-20,150,20), diago(np.linspace(-20,150,20)), color = 'grey')

diago2 = lambda x : x + 1
diago3 = lambda x : x - 1
#plt.plot(np.linspace(-20,150,20), diago3(np.linspace(-20,150,20)), color = 'red')
#plt.plot(np.linspace(-20,150,20), diago2(np.linspace(-20,150,20)), color = 'green')

dw = (piuk < 0) & (ziup > 0)
wd = (ziup < diago3(piuk)) & (piuk > 0)
wdd = (piuk > 0) & (ziup < 0)
ww = (ziup > diago2(piuk)) & (piuk > 0)
colors = ['forestgreen', 'orange', 'indianred', 'steelblue']
labels = ['dry-wet', 'wet slower', 'wet-dry', 'wet faster']

#for cos, col, lab in zip([dw, wd, wdd, ww], colors, labels):
for cos, col, lab in zip([dw, wdd], ['forestgreen', 'indianred'], ['dry-wet', 'wet-dry']):
    plt.scatter(piuk[cos], ziup[cos], color = col, label = lab, s = 1)

plt.xlabel('SSP5-8.5 trend (%/K)')
plt.ylabel('b100 trend (%/K)')
fig.savefig(cart_out + 'prssp_prb100_scatter.png')

long, latg = np.meshgrid(coso.lon, coso.lat)

fig, ax = ctl.get_cartopy_fig_ax(visualization='Robinson')
#for cos, col, lab in zip([dw, wd, wdd, ww], colors, labels):
for cos, col, lab in zip([dw, wdd], ['forestgreen', 'indianred'], ['dry-wet', 'wet-dry']):
    ax.scatter(long[cos], latg[cos], s = 2, color = col, transform = ccrs.PlateCarree(), label = lab)
#ctl.custom_legend(fig, colors, labels, ncol = 4)
ctl.custom_legend(fig, ['forestgreen', 'indianred'], ['dry-wet', 'wet-dry'], ncol = 2)
fig.savefig(cart_out + 'map_prssp_prb100_scatter.png')
