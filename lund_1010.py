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
rat1 = patt_maps[('tas', 'b100', 'stab')]/patt_maps[('tas', 'ssp585', 'fin')]
rat2 = patt_maps[('tas', 'b025', 'stab')]/patt_maps[('tas', 'ssp585', 'fin')]

# ctl.plot_multimap_contour([rat2, rat1], coso.lat, coso.lon, visualization = 'Robinson', subtitles = ['b025/ssp585', 'b100/ssp585'], cbar_range = (0, 2), cmap = ctl.heatmap(), filename = cart_out + 'taspatt_b025_b100_sspratio.pdf', figsize = (16,7), cb_label = 'Ratio of warming patterns', cbar_bottomspace = 0.11)

#ctl.plot_map_contour(rat1, coso.lat, coso.lon, visualization = 'Robinson', cbar_range = (0, 2), cmap = ctl.heatmap(), filename = cart_out + 'taspatt_ssp_vs_b100_ratio.pdf')


fig, axs = ctl.get_cartopy_fig_ax(visualization='Robinson', fix_subplots_shape = (1, 2), figsize = (16,6))

fig2, axs2 = plt.subplots(1, 2, figsize = (16,9))
### map of prec trends

for ru, ax, ax2 in zip(['b025', 'b100'], axs, axs2):
    ziup = patt_maps[('pr_rel', ru, 'stab')]
    piuk = patt_maps[('pr_rel', 'ssp585', 'fin')]

    ax2.scatter(piuk.flatten(), ziup.flatten(), color = 'black', s = 1)
    ax2.grid()
    piuk[np.abs(piuk) < 1] = np.nan

    diago = lambda x : x
    ax2.plot(np.linspace(-20,150,20), diago(np.linspace(-20,150,20)), color = 'grey')

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
        ax2.scatter(piuk[cos], ziup[cos], color = col, label = lab, s = 1)

    ax2.set_xlabel('SSP5-8.5 trend (%/K)')
    ax2.set_ylabel(ru+' trend (%/K)')

    long, latg = np.meshgrid(coso.lon, coso.lat)

    #for cos, col, lab in zip([dw, wd, wdd, ww], colors, labels):
    # for cos, col, lab in zip([dw, wdd], ['forestgreen', 'indianred'], ['dry-wet', 'wet-dry']):
    #     ax.scatter(long[cos], latg[cos], s = 0.5, color = col, transform = ccrs.PlateCarree(), label = lab)

    cmap = ctl.cmap_shading(['forestgreen', 'white', 'indianred'])

    cose = np.zeros(long.shape)
    cose[dw] = -1
    cose[wdd] = 1
    ax.pcolormesh(coso.lon, coso.lat, cose, vmin = -1, vmax = 1, cmap = cmap, transform = ccrs.PlateCarree())

    #ctl.custom_legend(fig, colors, labels, ncol = 4)
    ax.set_title(ru)

ctl.custom_legend(fig, ['forestgreen', 'indianred'], ['dry-wet', 'wet-dry'], ncol = 2)
#ctl.custom_legend(fig, [cm.get_cmap('PiYG_r')(0), cm.get_cmap('PiYG_r')(1)], ['dry-wet', 'wet-dry'], ncol = 2)

fig2.savefig(cart_out + 'prssp_b025_b100_scatter.png')
fig.savefig(cart_out + 'map_prssp_b025_b100_scatter.png')
fig.savefig(cart_out + 'map_prssp_b025_b100_scatter.pdf')
