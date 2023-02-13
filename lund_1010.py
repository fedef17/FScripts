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

from matplotlib import colors as mcolors
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

exp1 = 'b025'
exp2 = 'b100'

yeamean = pickle.load(open(cart_out + '../seasmean/bottino_seasmean_2D.p', 'rb'))
coso = yeamean[2][(exp1, 'tas')]
del yeamean

anom_maps, patt_maps = pickle.load(open(cart_out + 'all_maps_1000.p', 'rb'))

# map of temp trends
rat1 = patt_maps[('tas', exp2, 'stab')]/patt_maps[('tas', 'ssp585', 'fin')]
rat2 = patt_maps[('tas', exp1, 'stab')]/patt_maps[('tas', 'ssp585', 'fin')]

divnorm = mcolors.TwoSlopeNorm(vmin=0.05, vcenter=1., vmax=2.75)
ctl.plot_multimap_contour([rat2, rat1], coso.lat, coso.lon, visualization = 'Robinson', subtitles = ['{}/ssp585'.format(exp1), '{}/ssp585'.format(exp2)], cbar_range = (0.05, 2.75), color_norm = divnorm, cmap = ctl.heatmap(), filename = cart_out + 'tas_patt_NEW_p2.pdf', figsize = (16,7), cb_label = 'Ratio of warming patterns', cbar_bottomspace = 0.11, n_color_levels = 28, draw_grid = True)

#### Patt map with only ssp585, b025, b100
cmap = ctl.heatmap()
cbran = (0.25, 2.95)
divnorm = mcolors.TwoSlopeNorm(vmin=0.25, vcenter=1., vmax=2.95)
plotanom = False
cblab = 'Temperature change per degree of global warming'
okfigs = [patt_maps[('tas', 'ssp585', 'fin')], patt_maps[('tas', exp1, 'stab')], patt_maps[('tas', exp2, 'stab')]]
ctl.plot_multimap_contour(okfigs, coso.lat, coso.lon, filename = cart_out + 'tas_patt_NEW_p1.pdf', cmap = cmap, plot_anomalies = plotanom, cbar_range = cbran, figsize = (18, 7), subtitles = ['ssp585', exp1, exp2], color_norm = divnorm, cb_label = cblab, visualization = 'Robinson', central_lat_lon = (0.,0.), n_color_levels = 28, fix_subplots_shape = (1, 3), draw_grid = True)

ctl.plot_map_contour(rat1, coso.lat, coso.lon, visualization = 'Robinson', cbar_range = (0.05, 2.75), color_norm = divnorm, cmap = ctl.heatmap(), filename = cart_out + 'taspatt_ssp_vs_{}_ratio.pdf'.format(exp2), n_color_levels = 28, draw_grid = True)

#### Patt map with only ssp585, b025, b100
c5 = -11
c95 = 21
divnorm = mcolors.TwoSlopeNorm(vmin=c5, vcenter=0., vmax=c95)
cmap = 'BrBG'
cbran = (c5, c95)
cblab = 'Relative precipitation change per degree of global warming (%/K)'
okfigs = [patt_maps[('pr_rel', 'ssp585', 'fin')], patt_maps[('pr_rel', exp1, 'stab')], patt_maps[('pr_rel', exp2, 'stab')]]
ctl.plot_multimap_contour(okfigs, coso.lat, coso.lon, filename = cart_out + 'pr_patt_NEW_p1.pdf', cmap = cmap, plot_anomalies = plotanom, cbar_range = cbran, figsize = (18, 7), subtitles = ['ssp585', exp1, exp2], color_norm = divnorm, cb_label = cblab, visualization = 'Robinson', central_lat_lon = (0.,0.), n_color_levels = 17, fix_subplots_shape = (1, 3), draw_grid = True)


fig, axs = ctl.get_cartopy_fig_ax(visualization='Robinson', fix_subplots_shape = (1, 2), figsize = (16,6))

fig2, axs2 = plt.subplots(1, 2, figsize = (16,9))
### map of prec trends

for ru, ax, ax2 in zip([exp1, exp2], axs, axs2):
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

fig2.savefig(cart_out + 'prssp_{}_{}_scatter.png'.format(exp1, exp2))
fig.savefig(cart_out + 'map_prssp_{}_{}_scatter.png'.format(exp1, exp2))
fig.savefig(cart_out + 'map_prssp_{}_{}_scatter.pdf'.format(exp1, exp2))
