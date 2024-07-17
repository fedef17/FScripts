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
import matplotlib.gridspec as gs

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
####################################################################################################

glomeans, pimean, yeamean, mapmean = pickle.load(open(cart_in + 'bottino_seasmean_2D.p', 'rb'))

plot_all = False

for var in ['tas', 'pr']:
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


ext_pr = (-500., 500.)
ext_pr_rel = (-0.5, 0.5)
#ext_tas = (0., 15.)
ext_tas = (-10., 10.)
cmappa_tas = ctl.heatmap()

ext_relcha = (0.25, 0.75)
ext_relcha_pr = (-0.5, 0.5)

#cmappa_tas = cm.get_cmap('viridis').copy()#ctl.heatmap()
cmappa_tas2 = mcolors.LinearSegmentedColormap.from_list('bau2', ['white', 'indianred', 'brown', 'violet', 'white'])
cmappa_tas2.set_under('steelblue')
cmappa_tas2.set_over('white')
cmappa_tas2.set_bad('lightslategray', alpha = 0.5)

cmappa_pr = cm.get_cmap('BrBG').copy()
cmappa_pr.set_bad('lightslategray', alpha = 0.5)

ypre_trans = 30#50
ypre_stab = 30

cb_labels = ['Temp. anomaly (K)', 'Prec. anomaly (mm/year)', 'Relative change']
coso = yeamean[('b025', 'tas')]

glotas_hissp = np.concatenate([ctl.global_mean(yeamean[('historical_mean', 'tas')], coso.lat), ctl.global_mean(yeamean[('ssp585_mean', 'tas')], coso.lat)])

#fig, axs = plt.subplots(1, 3, figsize = (18,5))
for varnam, var, cmappa, ext, ext_rel, cblab in zip(['tas', 'pr', 'pr_rel'], ['tas', 'pr', 'pr'], [cmappa_tas, cmappa_pr, cmappa_pr], [ext_tas, ext_pr, ext_pr_rel], [ext_relcha, ext_relcha_pr, ext_relcha_pr], cb_labels):
    cosopi = yeamean[('pi', var)]
    # cosohist = yeamean[('hist', var)]
    # cosossp = yeamean[('ssp585', var)]
    cosohist = yeamean[('historical_mean', var)]
    cosossp = yeamean[('ssp585_mean', var)]
    #cosohistssp = xr.concat([cosohist, cosossp], dim = 'year')
    cosohistssp = np.concatenate([cosohist, cosossp], axis = 0)
    yeahissp = np.arange(1970, 2101)

    fact = 1
    if varnam == 'pr':
        fact = 60*60*24*365

    pimap = fact*cosopi.mean('year').values

    pirun = ctl.running_mean(cosopi, 50)
    pistd = (cosopi-pirun).std('year').values

    if varnam == 'pr':
        thres = pimap < 50.

    fig_patt_trans = []
    fig_patt_stab = []
    fig_patt_stab_relcha = []
    fig_abs_start = []
    fig_abs_end = []

    fig_ratio = []
    fig_std = []
    fig_tot = []

    for ru in ['hist', 'ssp585'] + allru3:
        coso = yeamean[(ru, var)]
        y0 = coso.year[0].values

        #transient = fact*xr.concat([cosohistssp.sel(year = slice(None, y0)), coso[:ypre_trans//2+1]], dim = 'year').sel(year = slice(y0-ypre_trans//2, y0+ypre_trans//2)).mean('year').values

        # finestra a cavallo
        #transient = fact*np.concatenate([cosohistssp.sel(year = slice(y0-ypre_trans//2, y0)).values, coso.sel(year = slice(y0, y0+ypre_trans//2)).values], axis = 0).mean(axis = 0)
        # finestra tutta prima
        #transient = fact*cosohistssp.sel(year = slice(y0-ypre_trans, y0)).mean('year').values
        transient = fact*np.mean(cosohistssp[(yeahissp >= y0-ypre_trans) & (yeahissp < y0), ...], axis = 0)

        equi = fact*coso[-ypre_stab:].mean('year').values

        deltatas = glomeans[(ru, 'tas')][1][-ypre_stab:].mean()-pimean['tas']
        deltatas_ini = glotas_hissp[(yeahissp >= y0-ypre_trans) & (yeahissp < y0)].mean()-pimean['tas']

        #mappe = [transient-pimap, equi-pimap, equi-transient]
        if ru in allru3:
            if plot_all:
                if varnam == 'pr_rel':
                    transient[thres] = np.nan
                    equi[thres] = np.nan
                    mappe = [(transient-pimap)/pimap, (equi-transient)/pimap]
                elif varnam == 'tas':
                    mappe = [transient-pimap, equi-transient]
                    #mappe = [ma - ctl.global_mean(ma, coso.lat) for ma in mappe]
                else:
                    mappe = [transient-pimap, equi-transient]

                ctl.plot_multimap_contour(mappe, coso.lat, coso.lon, filename = cart_out + 'trans_vs_stab_{}_{}.pdf'.format(varnam, ru), subtitles = ['transient', 'stabilization'], cmap = cmappa, plot_anomalies = True, cbar_range = ext, figsize = (16,6), cb_label = cblab)

                fig_patt_trans.append(mappe[0])
                fig_patt_stab.append(mappe[1])
                fig_patt_stab_relcha.append(mappe[1]/(deltatas-deltatas_ini))

                mappe += [mappe[0]+mappe[1]]
                ctl.plot_multimap_contour(mappe, coso.lat, coso.lon, filename = cart_out + 'triple_trans_vs_stab_{}_{}.pdf'.format(varnam, ru), subtitles = ['transient', 'stabilization', 'final'], cmap = cmappa, plot_anomalies = True, cbar_range = ext, figsize = (18,5), fix_subplots_shape = (1,3), cb_label = cblab)

                ctl.plot_multimap_contour(mappe, coso.lat, coso.lon, filename = cart_out + 'triple_trans_vs_stab_{}_{}_medzoom.pdf'.format(varnam, ru), subtitles = ['transient', 'stabilization', 'final'], cmap = cmappa, plot_anomalies = True, cbar_range = ext, figsize = (18,5), plot_margins = (-20, 45., 27., 50), fix_subplots_shape = (1,3), cb_label = cblab)

                if varnam == 'tas':
                    mappe = [transient-pimap, equi-pimap]
                    cmap = cmappa_tas2
                    cbran = (0., 20.)
                elif varnam == 'pr_rel':
                    mappe = [(transient-pimap)/pimap, (equi-pimap)/pimap]
                    cmap = cmappa
                    cbran = ext
                elif varnam == 'pr':
                    mappe = [transient-pimap, equi-pimap]
                    cmap = cmappa
                    cbran = ext

                ctl.plot_multimap_contour(mappe, coso.lat, coso.lon, filename = cart_out + 'abs_trans_vs_stab_{}_{}.pdf'.format(varnam, ru), subtitles = ['start', 'end'], cmap = cmap, plot_anomalies = True, cbar_range = ext, figsize = (16, 6), cb_label = cblab)

                fig_abs_start.append(mappe[0])
                fig_abs_end.append(mappe[1])

        if varnam == 'tas':
            #mappe = (equi-pimap)/deltatas
            mappe = (equi-pimap)
            cmap = ctl.heatmap_mono()
            #cmap = 'gist_ncar'
            #cbran = (0.5, 3.5)
            cbran = (0, 12)
            divnorm = None
            plotanom = False
        elif varnam == 'pr_rel':
            c5 = -10
            c95 = 20
            divnorm = mcolors.TwoSlopeNorm(vmin=c5, vcenter=0., vmax=c95)
            mappe = 100*((equi-pimap)/pimap)/deltatas
            mappe[thres] = np.nan
            cmap = 'BrBG'
            cbran = (c5, c95)
            plotanom = False
        elif varnam == 'pr':
            # c5 = -0.1
            # c95 = 0.3
            # divnorm = mcolors.TwoSlopeNorm(vmin=c5, vcenter=0., vmax=c95)
            divnorm = None
            plotanom = True
            mappe = (equi-pimap)/deltatas
            cmap = 'BrBG'
            cbran = ext

        ctl.plot_map_contour(mappe, coso.lat, coso.lon, filename = cart_out + 'tot_stab_{}_{}.pdf'.format(varnam, ru), cmap = cmap, plot_anomalies = plotanom, cbar_range = cbran, figsize = (16, 6), color_norm = divnorm, cb_label = cblab)

        fig_tot.append(mappe)

        # mappeanom = [ma-ctl.global_mean(ma, coso.lat) for ma in mappe]
        # ctl.plot_multimap_contour(mappeanom, coso.lat, coso.lon, cmap = cmappa, plot_anomalies = True, cbar_range = [-10., 10.])

        cmap = mcolors.LinearSegmentedColormap.from_list('bau', ['violet', 'white', 'lightgreen'])
        cmap.set_under('violet')
        cmap.set_over('lightgreen')

        if plot_all and ru in allru3:
            if varnam == 'tas':
                print('entro tas')
                mappa = (equi-transient)/(equi-pimap)
                ctl.plot_map_contour(mappa, coso.lat, coso.lon, filename = cart_out + 'stabtransratio_{}_{}.pdf'.format(varnam, ru), plot_anomalies = True, cbar_range = ext_rel, cmap = cmap, cb_label = 'Fraction of residual warming')
                fig_ratio.append(mappa)
                print(len(fig_ratio))
            # elif varnam == 'pr_rel':
            else:
                # cmap = cm.get_cmap('viridis').copy()
                # cmap.set_under('violet')
                # cmap.set_bad('lightslategray', alpha = 0.5)

                #mappa = (equi-transient)/(transient-pimap)
                mappa = (equi-transient)/(equi-pimap)

                if varnam == 'pr_rel':
                    #mappa[np.abs(transient-pimap)/pimap < 0.02] = np.nan
                    mappa[np.abs(equi-pimap)/pimap < 0.02] = np.nan
                    tramap = (transient-pimap)/pimap
                    tramap[thres] = np.nan
                    lst = 0.25
                elif varnam == 'pr':
                    #mappa[np.abs(transient-pimap) < 20] = np.nan
                    mappa[np.abs(equi-pimap) < 20] = np.nan
                    tramap = transient-pimap
                    lst = 200

                #cblab = 'stabilization/transient ratio'
                cblabU = 'Ratio of change during stabilization to total change'
                ctl.plot_map_contour(mappa, coso.lat, coso.lon, filename = cart_out + 'stabtransratio_{}_{}.pdf'.format(varnam, ru), plot_anomalies = True, cbar_range = ext_rel, cmap = cmap, cb_label = cblabU)#, add_contour_field = tramap, add_contour_lines_step = lst)
                fig_ratio.append(mappa)

                ctl.plot_map_contour(mappa, coso.lat, coso.lon, filename = cart_out + 'stabtransratio_{}_{}_medzoom.pdf'.format(varnam, ru), plot_anomalies = True, cbar_range = ext_rel, cmap = cmap, plot_margins = (-20, 45., 27., 50), cb_label = cblabU)#, add_contour_field = tramap, add_contour_lines_step = lst)

            cosorun = ctl.running_mean(coso, 50)
            equistd = (coso-cosorun)[-ypre_stab:].std('year').values

            mappa = equistd/pistd
            if varnam == 'pr_rel':
                mappa[thres] = np.nan

            ctl.plot_map_contour(mappa, coso.lat, coso.lon, cmap = cmappa, plot_anomalies = True, cbar_range = [0., 2.], filename = cart_out + 'equistdratio_{}_{}.pdf'.format(varnam, ru), cb_label = 'Interannual variability ratio')
            fig_std.append(mappa)


    #######

    if varnam == 'tas':
        cmap = ctl.heatmap_mono()
        cbran = (0., 12)
        divnorm = None
        plotanom = False
        cblab = 'Temperature anomaly (K)'
    elif varnam == 'pr_rel':
        # c5 = -0.1
        # c95 = 0.2
        # divnorm = mcolors.TwoSlopeNorm(vmin=c5, vcenter=0., vmax=c95)
        # cmap = 'BrBG'
        # cbran = (c5, c95)
        cmap = ctl.wetmap()
        cbran = (-30, 30)
        plotanom = True
        cblab = 'Relative precipitation change per degree of global warming (%/K)'
    elif varnam == 'pr':
        # c5 = -0.1
        # c95 = 0.3
        # divnorm = mcolors.TwoSlopeNorm(vmin=c5, vcenter=0., vmax=c95)
        divnorm = None
        plotanom = True
        cmap = 'BrBG'
        cbran = ext
        cblab = 'Precipitation anomaly (mm/year)'

    okru = ['hist', 'ssp585'] + allru3
    okfi = fig_tot
    ctl.plot_multimap_contour(okfi, coso.lat, coso.lon, filename = cart_out + 'All_tot_stab_{}.pdf'.format(varnam), cmap = cmap, plot_anomalies = plotanom, cbar_range = cbran, figsize = (16, 9), fix_subplots_shape = (2,3), subtitles = okru, color_norm = divnorm, cb_label = cblab)

    ctl.plot_multimap_contour(okfi, coso.lat, coso.lon, filename = cart_out + 'All_tot_stab_{}_newproj.pdf'.format(varnam), cmap = cmap, plot_anomalies = plotanom, cbar_range = cbran, figsize = (16, 9), fix_subplots_shape = (2,3), subtitles = okru, color_norm = divnorm, cb_label = cblab, visualization = 'Robinson', central_lat_lon = (0.,0.))
    ###########################3

    if plot_all:
        ctl.plot_multimap_contour(fig_patt_trans+fig_patt_stab, coso.lat, coso.lon, filename = cart_out + 'All_trans_vs_stab_{}.pdf'.format(varnam), subtitles = allru3+4*[None], cmap = cmappa, plot_anomalies = True, cbar_range = ext, figsize = (16,6), fix_subplots_shape = (2,4), cb_label = cblab)

        if varnam == 'tas':
            cmap = cmappa_tas2
            cbran = (0., 20.)
        else:
            cmap = cmappa
            cbran = ext
        ctl.plot_multimap_contour(fig_abs_start+fig_abs_end, coso.lat, coso.lon, filename = cart_out + 'All_abs_trans_vs_stab_{}.pdf'.format(varnam), subtitles = allru3+4*[None], cmap = cmap, plot_anomalies = True, cbar_range = cbran, figsize = (16,6), fix_subplots_shape = (2,4), cb_label = cblab)

        cmap = mcolors.LinearSegmentedColormap.from_list('bau', ['violet', 'white', 'lightgreen'])
        cmap.set_under('violet')
        cmap.set_over('lightgreen')

        if varnam == 'tas':
            addcont = None
            lst = 1.
            cla2 = 'Fraction of residual warming'
        else:
            addcont = fig_patt_trans
            #cla2 = 'stabilization/transient ratio'
            cla2 = 'Ratio of change during stabilization to total change'

        # fig = plt.figure(figsize = (16, 12))
        # gsco = gs.GridSpec(2, 4)
        # ax1 = plt.subplot(gsco[0, :2], projection = proj)
        # ax2 = plt.subplot(gsco[0, 2:], projection = proj)
        # ax3 = plt.subplot(gsco[1, 1:3], projection = proj)
        #
        # axs = [ax1, ax2, ax3]
        # for ax in axs:
        #     ax.set_global()
        #     ax.coastlines(linewidth = 1.)
        #
        # print(len(fig_ratio))
        # ctl.plot_multimap_contour(fig_ratio, coso.lat, coso.lon, filename = cart_out + 'All_stabtransratio_{}.pdf'.format(varnam), plot_anomalies = True, cbar_range = ext_rel, cmap = cmap, fig_external = fig, axs_external = axs, subtitles = allru3, cb_label = cla2)#figsize = (18,5), fix_subplots_shape = (1,3))
        ctl.plot_multimap_contour(fig_ratio, coso.lat, coso.lon, filename = cart_out + 'All_stabtransratio_{}.pdf'.format(varnam), plot_anomalies = True, cbar_range = ext_rel, cmap = cmap, subtitles = allru3, cb_label = cla2)#figsize = (18,5), fix_subplots_shape = (1,3))

        ctl.plot_multimap_contour(fig_std, coso.lat, coso.lon, cmap = cmappa, plot_anomalies = True, cbar_range = [0., 2.], filename = cart_out + 'All_equistdratio_{}.pdf'.format(varnam), figsize = (18,5), subtitles = allru3, cb_label = 'Interannual variability ratio')
        fig_std.append(mappa)
