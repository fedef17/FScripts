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
from matplotlib import colors as mcolors
import matplotlib.gridspec as gridspec

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
plt.rcParams['hatch.linewidth'] = 0.1
plt.rcParams['hatch.color'] = 'grey'

#############################################################################

if os.uname()[1] == 'hobbes':
    cart = '/home/fabiano/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'xaru':
    cart = '/home/fedef/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'tintin':
    cart = '/home/fabiano/work/lavori/BOTTINO/'

cart_in = cart + 'seasmean/'
cart_out = cart + 'simple_maps/'
ctl.mkdir(cart_out)

allru = ['pi', 'b025', 'b050', 'b100']
allnams = ['piControl', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']

colors = ['black', 'forestgreen', 'orange', 'violet']

allnams2 = allnams + ['ssp585', 'historical']
allru2 = allru + ['ssp585', 'hist']
colors2 = colors + ['indianred', 'steelblue']

allnams3 = ['stabilization-hist-1990'] + allnams[1:]
allru3 = ['b990'] + allru[1:]
colors3 = ['teal'] + colors[1:]

colors = ['black', 'steelblue', 'indianred', 'teal', 'forestgreen', 'orange', 'chocolate', 'maroon', 'violet']
allru = ['pi', 'hist', 'ssp585', 'b990', 'b025', 'b050', 'b065', 'b080', 'b100']
allru3 = allru[3:]
colors3 = colors[3:]

####################################################################################################

anom_maps = dict() # maps of anomalies
patt_maps = dict() # maps of anomalies divided by change in global tas

ypre_trans = 30#50
ypre_stab = 50

cb_labels = ['Temp. anomaly (K)', 'Prec. anomaly (mm/year)', 'Relative change']

glomeans, pimean = pickle.load(open(cart_in + 'bottino_glomeans_1000.p', 'rb'))

hist_sd = dict()
delta_stab = dict()

for var in ['tas', 'pr']:
    ### Estimate historical variability over 30 years from historical ensemble spread
    exp = 'historical'
    yeamean_hs, seamean_hs = pickle.load(open(cart_in + '../yearmean/bottino_yeamean_3_{}_{}.p'.format(exp, var), 'rb'))
    mean_clim = []
    for mem in yeamean_hs[('historical', 'members')]:
        print(mem)
        mean_clim.append(np.mean(yeamean_hs[('historical', mem, var)][-30:], axis = 0))
    mean_clim
    hist_sd[var] = np.std(mean_clim, axis = 0)

    yeamean = pickle.load(open(cart_in + 'bottino_seasmean_2D_{}_1000.p'.format(var), 'rb'))
    for ru in ['pi', 'hist', 'ssp585']:
        yeamean_add = pickle.load(open(cart_in + 'bottino_seasmean_2D_{}.p'.format(ru), 'rb'))

        for ke in list(yeamean_add):
            if var not in ke: del yeamean_add[ke]
        yeamean.update(yeamean_add)

    coso = yeamean[('b025', var)]

    # This is to calculate the ensemble mean of historical and ssp5
    for exp in ['ssp585', 'historical']:
        yeamean_hs, seamean_hs = pickle.load(open(cart_in + '../yearmean/bottino_yeamean_3_{}_{}.p'.format(exp, var), 'rb'))
        if exp == 'ssp585':
            memok = yeamean_hs[(exp, 'members')]

        yeamean[(exp+'_mean', var)] = yeamean_hs[(exp, 'ensmean', var)]
        yeamean[(exp+'_std', var)] = yeamean_hs[(exp, 'ensstd', var)]

        continue
        print(yeamean[(exp+'_mean', var)].shape)

        print(memok)
        lenmin = np.min([len(yeamean_hs[(exp, mem, var)]) for mem in memok])
        cosoye = np.mean([yeamean_hs[(exp, mem, var)][-lenmin:] for mem in memok], axis = 0)

        # cosok = xr.DataArray(data = cosoye, dims = ['year', 'lat', 'lon'], coords = [sssssssssssss], name = 'tos')

        yeamean[(exp+'_mean', var)] = cosoye
        cosoye = np.std([yeamean_hs[(exp, mem, var)][-lenmin:] for mem in memok], axis = 0)
        yeamean[(exp+'_std', var)] = cosoye


    if var == 'tas':
        glotas_hissp = np.concatenate([ctl.global_mean(yeamean[('historical_mean', 'tas')], coso.lat), ctl.global_mean(yeamean[('ssp585_mean', 'tas')], coso.lat)])

        for ru, y0 in zip(allru3, [1990, 2025, 2050, 2065, 2080, 2100]):
            deltatas = glomeans[(ru, 'tas')][1][-yrang:].mean()-pimean['tas']
            deltatas_ini = glotas_hissp[(yeahissp >= y0-ypre_trans) & (yeahissp < y0)].mean()-pimean['tas']
            delta_stab[ru] = deltatas-deltatas_ini

        for ru in ['hist', 'ssp585']:
            deltatas = glomeans[(ru, 'tas')][1][-yrang:].mean()-pimean['tas']
            delta_stab[ru] = deltatas


#     #for varnam, var, cblab in zip(['tas', 'pr', 'pr_rel'], ['tas', 'pr', 'pr'], cb_labels):
#     for tip in ['', '_rel']:
#         if tip == '_rel' and var == 'tas': continue
#         varnam = var + tip

#         cosopi = yeamean[('pi', var)]
#         cosohist = yeamean[('historical_mean', var)]
#         cosossp = yeamean[('ssp585_mean', var)]

#         cosohistssp = np.concatenate([cosohist, cosossp], axis = 0)
#         yeahissp = np.arange(1970, 2101)

#         fact = 1
#         if var == 'pr':
#             fact = 60*60*24*365

#         pimap = fact*cosopi.mean('year').values

#         pirun = ctl.running_mean(cosopi, ypre_stab)
#         pistd = (cosopi-pirun).std('year').values # this can be used for the significance

#         if var == 'pr':
#             thres = pimap < 50.

#         #for ru in ['hist', 'ssp585'] + allru3:
#         for ru in allru:
#             coso = yeamean[(ru, var)]
#             y0 = coso.year[0].values

#             if ru in ['hist', 'ssp585']:
#                 yrang = ypre_trans # using shorter period for transient
#             else:
#                 yrang = ypre_stab
            
#             equi = fact*coso[-yrang:].mean('year').values
#             deltatas = glomeans[(ru, 'tas')][1][-yrang:].mean()-pimean['tas']

#             if varnam == 'pr_rel':
#                 equi[thres] = np.nan
#                 mappa = 100*((equi-pimap)/pimap)
#             else:
#                 mappa = equi-pimap

#             anom_maps[(varnam, ru, 'fin')] = mappa
#             patt_maps[(varnam, ru, 'fin')] = mappa/deltatas

#             if ru in allru3:
#                 deltatas_ini = glotas_hissp[(yeahissp >= y0-ypre_trans) & (yeahissp < y0)].mean()-pimean['tas']
#                 transient = fact*np.mean(cosohistssp[(yeahissp >= y0-ypre_trans) & (yeahissp < y0), ...], axis = 0)

#                 if varnam == 'pr_rel':
#                     transient[thres] = np.nan
#                     mappa_trans = 100*(transient-pimap)/pimap
#                     mappa_stab = 100*(equi-transient)/pimap
#                 else:
#                     mappa_trans = transient-pimap
#                     mappa_stab = equi-transient

#                 anom_maps[(varnam, ru, 'trans')] = mappa_trans
#                 patt_maps[(varnam, ru, 'trans')] = mappa_trans/deltatas_ini

#                 anom_maps[(varnam, ru, 'stab')] = mappa_stab
#                 patt_maps[(varnam, ru, 'stab')] = mappa_stab/(deltatas-deltatas_ini)

#     print('Finished ', var)
#     del yeamean

# pickle.dump([anom_maps, patt_maps], open(cart_out + 'all_maps_1000.p', 'wb'))

yeamean = pickle.load(open(cart_in + 'bottino_seasmean_2D_b025.p', 'rb'))
coso = yeamean[('b025', 'tas')]

anom_maps, patt_maps = pickle.load(open(cart_out + 'all_maps_1000.p', 'rb'))

var = 'pr'
yeamean_pi = pickle.load(open(cart_in + 'bottino_seasmean_2D_pi.p', 'rb'))
cosopi = yeamean_pi[('pi', var)]
fact = 60*60*24*365
pimap = fact*cosopi.mean('year').values

hist_sd['pr'] = fact*hist_sd['pr']
hist_sd['pr_rel'] = 100*hist_sd['pr']/pimap


def make_fig_axes():
    proj = ctl.def_projection(visualization = 'Robinson', central_lat_lon = (0.,0.))
    fig = plt.figure(figsize=(16,12))
    gs = gridspec.GridSpec(3, 6)
    axes = []
    axes.append(fig.add_subplot(gs[0, 1:3], projection = proj))
    axes.append(fig.add_subplot(gs[0, 3:5], projection = proj))
    for i in range(1,3):
        axes.append(fig.add_subplot(gs[i, 0:2], projection = proj))
        axes.append(fig.add_subplot(gs[i, 2:4], projection = proj))
        axes.append(fig.add_subplot(gs[i, 4:], projection = proj))

    return fig, axes


## now the figures
plotanom = False
for varnam in ['tas', 'pr', 'pr_rel']:
#for varnam in ['pr', 'pr_rel']:
    # maps of total change at the end of the sims
    if varnam == 'tas':
        cmap = ctl.heatmap_mono()
        cbran = (0., 20)
        n_color_levels = 21
        divnorm = None
        plotanom = False
        cblab = 'Temperature anomaly (K)'
    elif varnam == 'pr_rel':
        c5 = -55
        c95 = 105
        n_color_levels = 17
        divnorm = mcolors.TwoSlopeNorm(vmin=c5, vcenter=0., vmax=c95)
        cmap = 'BrBG'
        cbran = (c5, c95)
        # cmap = ctl.wetmap()
        # cbran = (-80, 80)
        # plotanom = True
        cblab = 'Relative precipitation change (%)'
    elif varnam == 'pr':
        n_color_levels = 17
        c5 = -550
        c95 = 1050
        divnorm = mcolors.TwoSlopeNorm(vmin=c5, vcenter=0., vmax=c95)
        cmap = 'BrBG'
        cbran = (c5, c95)
        # divnorm = None
        # plotanom = True
        # cmap = ctl.wetmap()
        # cbran = (-1000., 1000.)
        cblab = 'Precipitation anomaly (mm/year)'

    okru = ['hist', 'ssp585', 'b990', 'b025', 'b050', 'b065', 'b080', 'b100']
    okfi = [anom_maps[(varnam, ru, 'fin')] for ru in okru]
    

    # fig, axes = make_fig_axes()
    # ctl.plot_multimap_contour(okfi, coso.lat, coso.lon, filename = cart_out + 'All_anom_fin_{}_newproj_1000.pdf'.format(varnam), cmap = cmap, plot_anomalies = plotanom, cbar_range = cbran, figsize = (16, 9), subtitles = okru, color_norm = divnorm, cb_label = cblab, visualization = 'Robinson', central_lat_lon = (0.,0.), axs_external = axes, fig_external = fig, n_color_levels = n_color_levels)

    okha = [np.abs(anom_maps[(varnam, ru, 'fin')]) > hist_sd[varnam] for ru in okru]
    for ii, ma in enumerate(okfi):
        okha[ii][np.isnan(ma)] = np.nan

    fig, axes = make_fig_axes()
    ctl.plot_multimap_contour(okfi, coso.lat, coso.lon, filename = cart_out + 'All_anom_fin_{}_newproj_1000_vs_histvar.pdf'.format(varnam), cmap = cmap, plot_anomalies = plotanom, cbar_range = cbran, figsize = (16, 9), subtitles = okru, color_norm = divnorm, cb_label = cblab, visualization = 'Robinson', central_lat_lon = (0.,0.), axs_external = axes, fig_external = fig, n_color_levels = n_color_levels, add_hatching = okha, hatch_styles = ['////', '', ''])
    #fig.savefig(cart_out + 'All_anom_fin_{}_newproj.pdf'.format(varnam))
    ###########################

    # maps of change per degree of global warming during the simulation (-> stab for the bottins)
    if varnam == 'tas':
        n_color_levels = 29
        cmap = ctl.heatmap()
        c5 = 0.25
        c95 = 3.05
        cbran = (c5, c95)
        divnorm = mcolors.TwoSlopeNorm(vmin=c5, vcenter=1., vmax=c95)
        plotanom = False
        cblab = 'Temperature change per degree of global warming'
    elif varnam == 'pr_rel':
        n_color_levels = 17
        c5 = -11
        c95 = 21
        divnorm = mcolors.TwoSlopeNorm(vmin=c5, vcenter=0., vmax=c95)
        cmap = 'BrBG'
        cbran = (c5, c95)
        #cmap = ctl.wetmap()
        #cbran = (-30, 30)
        #plotanom = True
        cblab = 'Relative precipitation change per degree of global warming (%/K)'
    elif varnam == 'pr':
        n_color_levels = 27
        c5 = -150
        c95 = 300
        cbran = (c5, c95)
        divnorm = mcolors.TwoSlopeNorm(vmin=c5, vcenter=0., vmax=c95)
        cmap = 'BrBG'
        # divnorm = None
        # plotanom = True
        #cmap = ctl.wetmap()
        # cbran = (-300., 300.)
        cblab = 'Precipitation change per degree of global warming (mm/year/K)'

    okfi = [patt_maps[(varnam, ru, 'fin')] if ru in ['hist', 'ssp585'] else patt_maps[(varnam, ru, 'stab')] for ru in okru]

    #fig, axes = make_fig_axes()
    #ctl.plot_multimap_contour(okfi, coso.lat, coso.lon, filename = cart_out + 'All_patt_stab_{}_newproj_1000.pdf'.format(varnam), cmap = cmap, plot_anomalies = plotanom, cbar_range = cbran, figsize = (16, 9), subtitles = okru, color_norm = divnorm, cb_label = cblab, visualization = 'Robinson', central_lat_lon = (0.,0.), n_color_levels = n_color_levels, axs_external = axes, fig_external = fig)

    scaled_hist_sds = [hist_sd[varnam]/delta_stab[ru] for ru in okru]

    #okha = [np.abs(co) > hist_sd[varnam] for co in okfi]
    okha = [np.abs(co) > shisd for co, shisd in zip(okfi, scaled_hist_sds)]
    for ii, ma in enumerate(okfi):
        okha[ii][np.isnan(ma)] = np.nan

    fig, axes = make_fig_axes()
    ctl.plot_multimap_contour(okfi, coso.lat, coso.lon, filename = cart_out + 'All_patt_stab_{}_newproj_1000_vs_histvar.pdf'.format(varnam), cmap = cmap, plot_anomalies = plotanom, cbar_range = cbran, figsize = (16, 9), subtitles = okru, color_norm = divnorm, cb_label = cblab, visualization = 'Robinson', central_lat_lon = (0.,0.), axs_external = axes, fig_external = fig, n_color_levels = n_color_levels, add_hatching = okha, hatch_styles = ['////', '', ''])

    #################
    # trans vs stab pattern for each bottin

    continue

    okru = ['hist', 'ssp585', 'b990', 'b025', 'b050', 'b065', 'b080', 'b100']
    okru = allru3 + 4*['']
    okfi = [patt_maps[(varnam, ru, 'trans')] for ru in allru3] + [patt_maps[(varnam, ru, 'stab')] for ru in allru3]

    fig, axes = make_fig_axes()
    ctl.plot_multimap_contour(okfi, coso.lat, coso.lon, filename = cart_out + 'All_patt_transvsstab_{}_newproj_1000.pdf'.format(varnam), cmap = cmap, plot_anomalies = plotanom, cbar_range = cbran, figsize = (16, 9), fix_subplots_shape = (2,4), subtitles = okru, color_norm = divnorm, cb_label = cblab, visualization = 'Robinson', central_lat_lon = (0.,0.), n_color_levels = n_color_levels)
    fig.savefig(cart_out + 'All_patt_transvsstab_{}_newproj_1000.pdf'.format(varnam))
