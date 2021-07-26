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
    cart_out = '/home/fabiano/Research/lavori/BOTTINO/yearmean/'
elif os.uname()[1] == 'xaru':
    cart_out = '/home/fedef/Research/lavori/BOTTINO/yearmean/'
elif os.uname()[1] == 'tintin':
    cart_out = '/home/fabiano/work/lavori/BOTTINO/yearmean/'

ctl.mkdir(cart_out)

filna = '/nas/BOTTINO/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/{}/{}i1p1f1/{}/{}/*nc'

allru = ['pi', 'b025', 'b050', 'b100']
allnams = ['piControl', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']

colors = ['black', 'forestgreen', 'orange', 'violet']

####################################################################################################

miptab = 'Amon'
allvars_2D = 'clt pr psl rlut rsut tas uas'.split()
allvars_3D = 'ta ua'.split()

var_map = 'clt pr tas rlut uas'.split()  # plot last 200 mean map, stddev, low/high var wrt pi
allnams2 = allnams + ['ssp585']
allru2 = allru + ['ssp585']
colors2 = colors + ['indianred']

yef = 50
yel = 200


figs_glob = []
axs_glob = []
# pimean = dict()
# glomeans = dict()
# yeamean = dict()
# mapmean = dict()
#
# for na, ru, col in zip(allnams, allru, colors):
# #for na, ru, col in zip(allnams2, allru2, colors2):
#     print(ru)
#     mem = 'r1'
#     if na == 'ssp585': mem = 'r4'
#
#     fils = np.concatenate([glob.glob(filna.format(na, mem, miptab, var)) for var in allvars_2D[:-1]])
#
#     kose = xr.open_mfdataset(fils, use_cftime = True)
#     kose = kose.drop_vars('time_bnds')
#
#     # Separate for uas
#     fils = glob.glob(filna.format(na, mem, miptab, allvars_2D[-1]))
#     if len(fils) > 0:
#         kosettt = xr.open_mfdataset(fils, use_cftime = True)
#         kosettt = kosettt.drop_vars('time_bnds')
#         kosettt = kosettt.drop_vars('height')
#         kose = kose.assign(uas = kosettt.uas)
#
#     for var in allvars_2D:
#         print(var)
#         if var not in kose:
#             continue
#
#         cosoye = kose[var].groupby("time.year").mean().compute()
#         yeamean[(ru, var)] = cosoye
#
#
# # 3D vars
# for na, ru, col in zip(allnams, allru, colors):
#     mem = 'r1'
#     if na == 'ssp585': mem = 'r4'
#
#     fils = np.concatenate([glob.glob(filna.format(na, mem, miptab, var)) for var in allvars_3D])
#
#     kose = xr.open_mfdataset(fils, use_cftime = True)
#     kose = kose.drop_vars('time_bnds')
#
#     for var in allvars_3D:
#         print(var)
#         cosoye = kose[var].groupby("time.year").mean().compute()
#         yeamean[(ru, var)] = cosoye
#
#
# pickle.dump(yeamean, open(cart_out + 'bottino_yeamean_fullres.p', 'wb'))
#
# sys.exit()

yeamean, seamean = pickle.load(open(cart_out + 'bottino_yeamean.p', 'rb'))

#cart_out_seas = cart_out + '../seasmean/'
#glomeans, pimean, yeamean, mapmean = pickle.load(open(cart_out_seas + 'bottino_seasmean.p', 'rb'))

#var_map = 'tas pr rlut clt'.split()  # plot last 200 mean map, stddev, low/high var wrt pi

###### Plots 2D
figs_map = []
for var in var_map + allvars_3D:
    # Mean pattern
    mappe = []
    pimean = yeamean[('pi', var)].mean('year')

    subtitles = []
    for ru in allru[1:]:
        yefirst = yeamean[(ru, var)].year.values[0]
        pinko = yeamean[(ru, var)].sel(year = slice(yefirst, yefirst + yef)).mean('year') - pimean
        mappe.append(pinko)
        subtitles.append(ru + ' - tr')

    for ru in allru[1:]:
        yefirst = yeamean[(ru, var)].year.values[0]
        yelast = yeamean[(ru, var)].year.values[-1]
        pinko = yeamean[(ru, var)].sel(year = slice(yelast - yel, yelast)).mean('year') - yeamean[(ru, var)].sel(year = slice(yefirst, yefirst + yef)).mean('year')
        mappe.append(pinko)
        subtitles.append(ru + ' - eq')

    cbr0 = ctl.get_cbar_range([ma.values for ma in mappe[:3]], True, (2,98))
    cbr1 = ctl.get_cbar_range([ma.values for ma in mappe[3:]], True, (2,98))
    cbar_range = 3*[cbr0] + 3*[cbr1]

    # subtitles = ['{} - {}'.format(ru, seasok) for seasok in ['DJF', 'MAM', 'JJA', 'SON'] for ru in allru]

    if var in var_map:
        fig = ctl.plot_multimap_contour(mappe, figsize = (16,8), title = var, plot_anomalies = True, add_contour_field = 6*[pimean], add_contour_plot_anomalies = False, add_contour_same_levels = False, fix_subplots_shape = (2,3), subtitles = subtitles, use_different_cbars = True, cbar_range = cbar_range) #, cmap = cmaps, cbar_range = cbar_range, use_different_cbars = True, use_different_cmaps = True)
        figs_map.append(fig[0])
    else:
        # 3D vars
        fig, axs = plt.subplots(2, 3, figsize = (12,8))
        for ma, ax, subt in zip(mappe, axs.flatten(), subtitles):
            if 'tr' in subt:
                vma = 10
            else:
                vma = 5
            guplo = ma.mean('lon').plot.contourf(x = 'lat', y = 'plev', ax = ax, levels = 11, ylim = (1.e5, 1.e3), yscale = 'log', vmax = vma)
            try:
                guplo.colorbar.set_label('')
            except Exception as exc:
                print(exc)
                pass

            guplo2 = pimean.mean('lon').plot.contour(x = 'lat', y = 'plev', ax = ax, levels = 11, ylim = (1.e5, 1.e3), yscale = 'log', colors="k", add_colorbar=False)

        for i in range(2):
            for j in range(3):
                if j > 0:
                    axs[i,j].set_ylabel('')
                if i == 0:
                    axs[i,j].set_xlabel('')

        figs_map.append(fig)

#figs_map = np.concatenate(figs_map)
fignames = [var for var in var_map+allvars_3D]
ctl.plot_pdfpages(cart_out + 'bott3_map_f{}vsl{}.pdf'.format(yef, yel), figs_map, True, fignames)


#### Now, this is the ratio between the two warmings. Equilibration warming/transient (which is larger for most regions). This is always positive for temp, but maybe not for prec/rlut/clt.
# Other possibility is to plot the two contributions separated but in terms of ratio to the global warming. Ok doing this. AAA ok.

###### Plots 2D
figs_map = []
for var in var_map + allvars_3D:
    # Mean pattern
    mappe = []
    mappediv = []
    pimean = yeamean[('pi', var)].mean('year')

    totchan = []
    for ru in allru[1:]:
        yelast = yeamean[(ru, var)].year.values[-1]
        pinko = yeamean[(ru, var)].sel(year = slice(yelast - yel, yelast)).mean('year') - pimean
        totchan.append(pinko)

    subtitles = []
    for ru, totc in zip(allru[1:], totchan):
        yefirst = yeamean[(ru, var)].year.values[0]
        pinko = yeamean[(ru, var)].sel(year = slice(yefirst, yefirst + yef)).mean('year') - pimean
        mappe.append(pinko)

        gigino = pinko/totc
        if var != 'tas': ## where the final change is less than 1%, set to nan
            sputu = np.abs(totc/pimean) < 0.01
            gigino.values[sputu.values] = np.nan
            # gigino[np.abs(totc/pimean) < 0.01] = np.nan # does not work

        mappediv.append(gigino)
        subtitles.append(ru + ' - tr')

    for ru, totc in zip(allru[1:], totchan):
        yefirst = yeamean[(ru, var)].year.values[0]
        yelast = yeamean[(ru, var)].year.values[-1]
        pinko = yeamean[(ru, var)].sel(year = slice(yelast - yel, yelast)).mean('year') - yeamean[(ru, var)].sel(year = slice(yefirst, yefirst + yef)).mean('year')
        mappe.append(pinko)

        gigino = pinko/totc
        if var != 'tas': ## where the final change is less than 1%, set to nan
            sputu = np.abs(totc/pimean) < 0.01
            gigino.values[sputu.values] = np.nan
            # gigino[np.abs(totc/pimean) < 0.01] = np.nan # does not work
        mappediv.append(gigino)
        subtitles.append(ru + ' - eq')

    cbr0 = ctl.get_cbar_range([ma.values for ma in mappe[:3]], True, (2,98))
    cbr1 = ctl.get_cbar_range([ma.values for ma in mappe[3:]], True, (2,98))
    cbar_range = 3*[cbr0] + 3*[cbr1]

    # subtitles = ['{} - {}'.format(ru, seasok) for seasok in ['DJF', 'MAM', 'JJA', 'SON'] for ru in allru]

    #if var == 'tas':
    if var in var_map:
        if var == 'tas':
            cbran = [0., 1.0]
        else:
            cbran = [-2., 2.]
        fig = ctl.plot_multimap_contour(mappediv[:3], figsize = (24,8), title = None, plot_anomalies = True, add_contour_field = 3*[pimean], add_contour_plot_anomalies = False, add_contour_same_levels = False, fix_subplots_shape = (1,3), subtitles = allru[1:], use_different_cbars = False, cbar_range = cbran, cmap = 'viridis', cb_label = var)
        figs_map.append(fig[0])
    #elif var in var_map:
        # zonal mean
        fig, axs = plt.subplots(1, 3, figsize = (24,8))
        for iii, (ax, subt, totc) in enumerate(zip(axs.flatten(), allru[1:], totchan)):
            ma = mappe[iii]
            guplo = ma.mean('lon').plot(x = 'lat', ax = ax, label = 'transient', color = 'orange')
            ma2 = mappe[iii+3]
            guplo = ma2.mean('lon').plot(x = 'lat', ax = ax, label = 'stabilization', color = 'blue')
            guplo = totc.mean('lon').plot(x = 'lat', ax = ax, label = 'total', color = 'black')
            ax.legend()
            ax.grid()
            ax.set_title(subt.split()[0])

        figs_map.append(fig)
    else:
        # 3D vars
        fig, axs = plt.subplots(2, 3, figsize = (12,8))
        for ma, ax, subt, totc in zip(mappe, axs.flatten(), subtitles, 2*totchan):
            vma = 1
            # if 'tr' in subt:
            #     vma = 10
            # else:
            #     vma = 5
            guplo = (ma.mean('lon')/totc.mean('lon')).plot.contourf(x = 'lat', y = 'plev', ax = ax, levels = 11, ylim = (1.e5, 1.e3), yscale = 'log', vmax = vma)
            try:
                guplo.colorbar.set_label('')
            except Exception as exc:
                print(exc)
                pass

            guplo2 = pimean.mean('lon').plot.contour(x = 'lat', y = 'plev', ax = ax, levels = 11, ylim = (1.e5, 1.e3), yscale = 'log', colors="k", add_colorbar=False)

        for i in range(2):
            for j in range(3):
                if j > 0:
                    axs[i,j].set_ylabel('')
                if i == 0:
                    axs[i,j].set_xlabel('')

        figs_map.append(fig)

#figs_map = np.concatenate(figs_map)
fignames = [var for var in var_map+allvars_3D]
ctl.plot_pdfpages(cart_out + 'bott3_map_f{}vsl{}_rel.pdf'.format(yef, yel), figs_map, True, fignames)


#### for JJAS and DJFM
###### Plots 2D

for sea in ['DJFM', 'JJAS']:
    figs_map = []
    for var in var_map + allvars_3D:
        # Mean pattern
        mappe = []
        pimean = seamean[('pi', var)]['seamean'].sel(season = sea)

        subtitles = []
        for ru in allru[1:]:
            pinko = seamean[(ru, var, 'trans')]['seamean'].sel(season = sea) - pimean
            mappe.append(pinko)
            subtitles.append(ru + ' - tr')

        for ru in allru[1:]:
            pinko = seamean[(ru, var, 'stab')]['seamean'].sel(season = sea) - seamean[(ru, var, 'trans')]['seamean'].sel(season = sea)
            mappe.append(pinko)
            subtitles.append(ru + ' - eq')

        cbr0 = ctl.get_cbar_range([ma.values for ma in mappe[:3]], True, (2,98))
        cbr1 = ctl.get_cbar_range([ma.values for ma in mappe[3:]], True, (2,98))
        cbar_range = 3*[cbr0] + 3*[cbr1]

        # subtitles = ['{} - {}'.format(ru, seasok) for seasok in ['DJF', 'MAM', 'JJA', 'SON'] for ru in allru]

        if var in var_map:
            fig = ctl.plot_multimap_contour(mappe, figsize = (16,8), title = var, plot_anomalies = True, add_contour_field = 6*[pimean], add_contour_plot_anomalies = False, add_contour_same_levels = False, fix_subplots_shape = (2,3), subtitles = subtitles, use_different_cbars = True, cbar_range = cbar_range) #, cmap = cmaps, cbar_range = cbar_range, use_different_cbars = True, use_different_cmaps = True)
            figs_map.append(fig[0])
        else:
            # 3D vars
            fig, axs = plt.subplots(2, 3, figsize = (12,8))
            for ma, ax, subt in zip(mappe, axs.flatten(), subtitles):
                if 'tr' in subt:
                    vma = 10
                else:
                    vma = 5
                guplo = ma.mean('lon').plot.contourf(x = 'lat', y = 'plev', ax = ax, levels = 11, ylim = (1.e5, 1.e3), yscale = 'log', vmax = vma)
                try:
                    guplo.colorbar.set_label('')
                except Exception as exc:
                    print(exc)
                    pass

                guplo2 = pimean.mean('lon').plot.contour(x = 'lat', y = 'plev', ax = ax, levels = 11, ylim = (1.e5, 1.e3), yscale = 'log', colors="k", add_colorbar=False)

            for i in range(2):
                for j in range(3):
                    if j > 0:
                        axs[i,j].set_ylabel('')
                    if i == 0:
                        axs[i,j].set_xlabel('')

            figs_map.append(fig)

    #figs_map = np.concatenate(figs_map)
    fignames = [var for var in var_map+allvars_3D]
    ctl.plot_pdfpages(cart_out + 'bott3_map_f{}vsl{}_{}.pdf'.format(yef, yel, sea), figs_map, True, fignames)


    #### Now, this is the ratio between the two warmings. Equilibration warming/transient (which is larger for most regions). This is always positive for temp, but maybe not for prec/rlut/clt.
    # Other possibility is to plot the two contributions separated but in terms of ratio to the global warming. Ok doing this. AAA ok.

    ###### Plots 2D
    figs_map = []
    for var in var_map + allvars_3D:
        # Mean pattern
        mappe = []
        mappediv = []
        pimean = seamean[('pi', var)]['seamean'].sel(season = sea)

        totchan = []
        for ru in allru[1:]:
            pinko = seamean[(ru, var, 'stab')]['seamean'].sel(season = sea) - pimean
            totchan.append(pinko)

        subtitles = []
        for ru, totc in zip(allru[1:], totchan):
            pinko = seamean[(ru, var, 'trans')]['seamean'].sel(season = sea) - pimean
            mappe.append(pinko)

            gigino = pinko/totc
            if var != 'tas': ## where the final change is less than 1%, set to nan
                sputu = np.abs(totc/pimean) < 0.01
                gigino.values[sputu.values] = np.nan
                # gigino[np.abs(totc/pimean) < 0.01] = np.nan # does not work
            mappediv.append(gigino)
            subtitles.append(ru + ' - tr')

        for ru, totc in zip(allru[1:], totchan):
            pinko = seamean[(ru, var, 'stab')]['seamean'].sel(season = sea) - seamean[(ru, var, 'trans')]['seamean'].sel(season = sea)
            mappe.append(pinko)

            gigino = pinko/totc
            if var != 'tas': ## where the final change is less than 1%, set to nan
                sputu = np.abs(totc/pimean) < 0.01
                gigino.values[sputu.values] = np.nan
                # gigino[np.abs(totc/pimean) < 0.01] = np.nan # does not work
            mappediv.append(gigino)
            subtitles.append(ru + ' - eq')

        cbr0 = ctl.get_cbar_range([ma.values for ma in mappe[:3]], True, (2,98))
        cbr1 = ctl.get_cbar_range([ma.values for ma in mappe[3:]], True, (2,98))
        cbar_range = 3*[cbr0] + 3*[cbr1]

        # subtitles = ['{} - {}'.format(ru, seasok) for seasok in ['DJF', 'MAM', 'JJA', 'SON'] for ru in allru]

        #if var == 'tas':
        if var in var_map:
            if var == 'tas':
                cbran = [0., 1.0]
            else:
                cbran = [-2., 2.]
            fig = ctl.plot_multimap_contour(mappediv[:3], figsize = (24,8), title = None, plot_anomalies = True, add_contour_field = 3*[pimean], add_contour_plot_anomalies = False, add_contour_same_levels = False, fix_subplots_shape = (1,3), subtitles = allru[1:], use_different_cbars = False, cbar_range = cbran, cmap = 'viridis', cb_label = var)
            figs_map.append(fig[0])
        #elif var in var_map:
            # zonal mean
            fig, axs = plt.subplots(1, 3, figsize = (24,8))
            for iii, (ax, subt, totc) in enumerate(zip(axs.flatten(), allru[1:], totchan)):
                ma = mappe[iii]
                guplo = ma.mean('lon').plot(x = 'lat', ax = ax, label = 'transient', color = 'orange')
                ma2 = mappe[iii+3]
                guplo = ma2.mean('lon').plot(x = 'lat', ax = ax, label = 'stabilization', color = 'blue')
                guplo = totc.mean('lon').plot(x = 'lat', ax = ax, label = 'total', color = 'black')
                ax.legend()
                ax.grid()
                ax.set_title(subt.split()[0])
                ax.set_ylabel('')
            fig.suptitle(var)

            figs_map.append(fig)
        else:
            # 3D vars
            fig, axs = plt.subplots(2, 3, figsize = (12,8))
            for ma, ax, subt, totc in zip(mappe, axs.flatten(), subtitles, 2*totchan):
                vma = 1
                # if 'tr' in subt:
                #     vma = 10
                # else:
                #     vma = 5
                guplo = (ma.mean('lon')/totc.mean('lon')).plot.contourf(x = 'lat', y = 'plev', ax = ax, levels = 11, ylim = (1.e5, 1.e3), yscale = 'log', vmax = vma)
                try:
                    guplo.colorbar.set_label('')
                except Exception as exc:
                    print(exc)
                    pass

                guplo2 = pimean.mean('lon').plot.contour(x = 'lat', y = 'plev', ax = ax, levels = 11, ylim = (1.e5, 1.e3), yscale = 'log', colors="k", add_colorbar=False)

            for i in range(2):
                for j in range(3):
                    if j > 0:
                        axs[i,j].set_ylabel('')
                    if i == 0:
                        axs[i,j].set_xlabel('')

            figs_map.append(fig)

    #figs_map = np.concatenate(figs_map)
    fignames = [var for var in var_map+allvars_3D]
    ctl.plot_pdfpages(cart_out + 'bott3_map_f{}vsl{}_rel_{}.pdf'.format(yef, yel, sea), figs_map, True, fignames)
