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
#import tunlib as tl

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
    cart_out = '/home/fabiano/Research/lavori/BOTTINO/seasmean/'
elif os.uname()[1] == 'xaru':
    cart_out = '/home/fedef/Research/lavori/BOTTINO/seasmean/'

ctl.mkdir(cart_out)

filna = '/nas/BOTTINO/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/{}/{}i1p1f1/{}/{}/*nc'

allru = ['pi', 'b025', 'b050', 'b100']
allnams = ['piControl', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']

colors = ['black', 'forestgreen', 'orange', 'violet']

####################################################################################################

miptab = 'Amon_r25'
allvars_2D = 'clt  pr  psl  rlut  rsdt  rsut  tas  uas'.split()
allvars_3D = 'ta ua'.split()

var_glob_mean = 'tas pr clt rlut rsut net_toa'.split()  # plot global timeseries, including ssp585

var_map_200 = 'clt pr tas rlut uas'.split()  # plot last 200 mean map, stddev, low/high var wrt pi
allnams2 = allnams + ['ssp585']
allru2 = allru + ['ssp585']
colors2 = colors + ['indianred']

# figs_glob = []
# axs_glob = []
# pimean = dict()
# glomeans = dict()
# yeamean = dict()
# mapmean = dict()
#
# for var in var_glob_mean:
#     fig, ax = plt.subplots(figsize = (12,8))
#     axs_glob.append(ax)
#     figs_glob.append(fig)
#     ax.set_title(var)
#
# fig_greg, ax_greg = plt.subplots(figsize = (12,8))
#
# #for na, ru, col in zip(allnams, allru, colors):
# for na, ru, col in zip(allnams2, allru2, colors2):
#     print(ru)
#     mem = 'r1'
#     if na == 'ssp585': mem = 'r4'
#
#     fils = np.concatenate([glob.glob(filna.format(na, mem, miptab, var)) for var in allvars_2D[:-1]])
#
#     kose = xr.open_mfdataset(fils, use_cftime = True)
#     kose = kose.drop_vars('time_bnds')
#
#     try:
#         kose = kose.assign(net_toa = kose.rsdt - kose.rlut - kose.rsut) # net downward energy flux at TOA
#     except Exception as exc:
#         print(exc)
#         pass
#
#     # Separate for uas
#     fils = glob.glob(filna.format(na, mem, miptab, allvars_2D[-1]))
#     if len(fils) > 0:
#         kosettt = xr.open_mfdataset(fils, use_cftime = True)
#         kosettt = kosettt.drop_vars('time_bnds')
#         kosettt = kosettt.drop_vars('height')
#         kose = kose.assign(uas = kosettt.uas)
#
#     for var, fig, ax in zip(var_glob_mean, figs_glob, axs_glob):
#         print(var)
#         if var not in kose:
#             if ru == 'pi':
#                 pimean[var] = 0.
#             continue
#
#         cosoye = kose[var].groupby("time.year").mean().compute()
#         yeamean[(ru, var)] = cosoye
#
#         coso = cosoye.mean('lon')
#         glomean = np.average(coso, weights = abs(np.cos(np.deg2rad(coso.lat))), axis = -1)
#         if ru == 'pi':
#             years = coso.year.data-2256+2015
#             pimean[var] = np.mean(glomean)
#         else:
#             years = coso.year.data
#
#         glomeans[(ru, var)] = (years, glomean)
#         ax.plot(years, glomean-pimean[var], label = ru, color = col)
#
#         if ru == 'b100':
#             ax.legend()
#             ax.grid()
#
#     # gregory
#     try:
#         #ax_greg.plot(glomeans[(ru, 'tas')][1]-pimean['tas'], glomeans[(ru, 'net_toa')][1], label = ru, color = col)
#         tl.gregplot_on_ax(ax_greg, glomeans[(ru, 'tas')][1]-pimean['tas'], glomeans[(ru, 'net_toa')][1], color = col, label = ru, calc_ERF = False, calc_ECS = False)
#     except Exception as exc:
#         print(exc)
#         pass
#
#     if ru == 'b100':
#         ax_greg.legend()
#         ax_greg.grid()
#
#     if ru == 'ssp585':
#         continue
#
#     if ru != 'pi':
#         kose = kose.sel(time = kose['time.year'] >= kose['time.year'].data[-1]-200)
#
#     for var in var_map_200:
#         print(var)
#         kose_sclim = ctl.seasonal_climatology(kose[var])
#         mapmean[(ru, var)] = kose_sclim
#
#     # kose_smean = kose.groupby("time.season").mean()
#     # kose_sstd = kose.groupby("time.season").std()
#     # kose_p90 = kose.groupby("time.season").percentile(90)
#     # kose_p10 = kose.groupby("time.season").percentile(10)
#     #
#     # for var in var_map_200:
#     #     mapmean[(ru, var, 'mean')] = kose_smean
#     #     mapmean[(ru, var, 'std')] = kose_sstd
#     #     mapmean[(ru, var, 'p90')] = kose_p90
#     #     mapmean[(ru, var, 'p10')] = kose_p10
#
#     # for var in var_map_200:
#     #     print(var)
#     #     vmax = np.nanpercentile(flux_season[var], 98)
#     #     fig = plt.figure()
#     #     guplo = flux_season[var].plot.contourf(col = 'season', col_wrap = 2, levels = 11, vmax = vmax, transform = proj, figsize = (16,12), subplot_kws = {"projection": proj})
#     #     guplo.map(lambda: plt.gca().coastlines())
#     #     plt.title(var)
#     #     plt.savefig(cart + '{}_seas_{}.pdf'.format(mod, var))
#
#
# ctl.plot_pdfpages(cart_out + 'bottino_glomeans.pdf', figs_glob, True, )
# fig_greg.savefig(cart_out + 'bottino_gregory.pdf')
#
# pickle.dump([glomeans, pimean, yeamean, mapmean], open(cart_out + 'bottino_seasmean_2D.p', 'wb'))
#
# # 3D vars
# for na, ru, col in zip(allnams2, allru2, colors2):
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
#     if ru == 'ssp585':
#         continue
#
#     if ru != 'pi':
#         kose = kose.sel(time = kose['time.year'] >= kose['time.year'].data[-1]-200)
#
#     for var in allvars_3D:
#         print(var)
#         kose_sclim = ctl.seasonal_climatology(kose[var])
#         mapmean[(ru, var)] = kose_sclim
#
#     # kose_smean = kose.groupby("time.season").mean()
#     # kose_sstd = kose.groupby("time.season").std()
#     # kose_p90 = kose.groupby("time.season").percentile(90)
#     # kose_p10 = kose.groupby("time.season").percentile(10)
#     #
#     # for var in var_map_200:
#     #     mapmean[(ru, var, 'mean')] = kose_smean
#     #     mapmean[(ru, var, 'std')] = kose_sstd
#     #     mapmean[(ru, var, 'p90')] = kose_p90
#     #     mapmean[(ru, var, 'p10')] = kose_p10
#
# pickle.dump([glomeans, pimean, yeamean, mapmean], open(cart_out + 'bottino_seasmean.p', 'wb'))

glomeans, pimean, yeamean, mapmean = pickle.load(open(cart_out + 'bottino_seasmean.p', 'rb'))

for var in var_glob_mean:
    fig, ax = plt.subplots(figsize = (12,8))
    axs_glob.append(ax)
    figs_glob.append(fig)
    ax.set_title(var)

fig_greg, ax_greg = plt.subplots(figsize = (12,8))

#for na, ru, col in zip(allnams, allru, colors):
for na, ru, col in zip(allnams2, allru2, colors2):
    print(ru)

    for var, fig, ax in zip(var_glob_mean, figs_glob, axs_glob):
        print(var)
        if (ru, var) not in glomeans.keys():
            print('NOT found')
            continue

        cosoye = yeamean[(ru, var)]
        years, glomean = glomeans[(ru, var)]

        ax.plot(years, glomean-pimean[var], label = ru, color = col)

        if ru == 'b100':
            ax.legend()
            ax.grid()

    # gregory
    try:
        #ax_greg.plot(glomeans[(ru, 'tas')][1]-pimean['tas'], glomeans[(ru, 'net_toa')][1], label = ru, color = col)
        tl.gregplot_on_ax(ax_greg, glomeans[(ru, 'tas')][1]-pimean['tas'], glomeans[(ru, 'net_toa')][1], color = col, label = ru, calc_ERF = False, calc_ECS = False)
    except Exception as exc:
        print(exc)
        pass

    if ru == 'b100':
        ax_greg.legend()
        ax_greg.grid()

ctl.plot_pdfpages(cart_out + 'bottino_glomeans.pdf', figs_glob, True, )

ax_greg.set_xlabel('Global mean tas (K)')
ax_greg.set_ylabel('Global net incoming TOA flux (W/m2)')
fig_greg.savefig(cart_out + 'bottino_gregory.pdf')

ax_greg.set_xlim((-0.5, 0.5))
ax_greg.set_ylim((-0.5, 0.5))
fig_greg.savefig(cart_out + 'bottino_gregory_pizoom.pdf')


var_map_200 = 'clt pr tas rlut uas'.split()  # plot last 200 mean map, stddev, low/high var wrt pi

allcopls = ['seamean', 'seastd', 'seap10', 'seap90']
for ru in allru:
    zup = mapmean[(ru, 'tas')]
    mapmean[(ru, 'tas_patt')] = zup
    zupme = ctl.global_mean(zup['seamean'])
    for copl in allcopls:
        mapmean[(ru, 'tas_patt')][copl] = zup[copl]-zupme[:, np.newaxis, np.newaxis]

    # zup = mapmean[(ru, 'pr')]
    # mapmean[(ru, 'pr_perc')] = zup
    # zupme = mapmean[('pi', 'pr')]['seamean']
    # for copl in allcopls:
    #     mapmean[(ru, 'pr_perc')][copl] = (zup[copl]-zupme)/zupme

print('CHECK! -> ', mapmean[(ru, 'tas')] is mapmean[(ru, 'tas_patt')])
sys.exit()

var_map_200 += ['tas_patt', 'pr_perc']
###### Plots 2D
figs_map = []
for var in var_map_200:
    for copl in allcopls:
        #mappe = [mapmean[('pi', var, copl)]] + [mapmean[(ru, var, copl)]-mapmean[('pi', var, copl)] for ru in allru[1:]]
        if var == 'pr_perc':
            zupme = mapmean[('pi', 'pr')]['seamean']
            mappe = [mapmean[('pi', 'pr')][copl]] + [(mapmean[(ru, 'pr')][copl]-zupme)/zupme for ru in allru[1:]]
            mapcont = [None if ru == 'pi' else mapmean[('pi', 'pr')][copl].sel(season = seasok) for seasok in ['DJF', 'MAM', 'JJA', 'SON'] for ru in allru]
        else:
            mappe = [mapmean[('pi', var)][copl]] + [mapmean[(ru, var)][copl]-mapmean[('pi', var)][copl] for ru in allru[1:]]
            mapcont = [None if ru == 'pi' else mapmean[('pi', var)][copl].sel(season = seasok) for seasok in ['DJF', 'MAM', 'JJA', 'SON'] for ru in allru]

        plot_anoms = [False if ru == 'pi' else True for se in range(4) for ru in allru]
        cmaps = ['viridis' if ru == 'pi' else 'RdBu_r' for se in range(4) for ru in allru]

        cbr0 = ctl.get_cbar_range(mappe[0].values, False, (2,98))
        cbr1 = ctl.get_cbar_range([ma.values for ma in mappe[1:]], True, (2,98))
        cbar_range = [cbr0 if ru == 'pi' else cbr1 for se in range(4) for ru in allru]

        #mappeseas = [ma.sel(time = ma['time.season'] == seasok). for seasok in ['DJF', 'MAM', 'JJA', 'SON'] for ma in mappe]
        mappeseas = [ma.sel(season = seasok) for seasok in ['DJF', 'MAM', 'JJA', 'SON'] for ma in mappe]

        subtitles = ['{} - {}'.format(ru, seasok) for seasok in ['DJF', 'MAM', 'JJA', 'SON'] for ru in allru]

        fig = ctl.plot_multimap_contour(mappeseas, figsize = (20,12), cmap = cmaps, cbar_range = cbar_range, use_different_cbars = True, use_different_cmaps = True, subtitles = subtitles, title = var+' - '+copl, add_contour_field = mapcont, add_contour_plot_anomalies = False)
        figs_map.append(fig)

figs_map = np.concatenate(figs_map)
fignames = [var+'_'+copl for var in var_map_200 for copl in allcopls]
ctl.plot_pdfpages(cart_out + 'bottino_mapmeans.pdf', figs_map, True, fignames)

figs_map = []
for var in allvars_3D:
    for copl in allcopls:
        fig, axs = plt.subplots(4, 4, figsize = (20,12))
        #mappe = [mapmean[('pi', var, copl)]] + [mapmean[(ru, var, copl)]-mapmean[('pi', var, copl)] for ru in allru[1:]]
        mappe = [mapmean[('pi', var)][copl]] + [mapmean[(ru, var)][copl]-mapmean[('pi', var)][copl] for ru in allru[1:]]

        #mappeseas = [ma.sel(time = ma['time.season'] == seasok). for seasok in ['DJF', 'MAM', 'JJA', 'SON'] for ma in mappe]
        mappeseas = [ma.sel(season = seasok) for seasok in ['DJF', 'MAM', 'JJA', 'SON'] for ma in mappe]
        subtitles = ['{} - {}'.format(ru, seasok) for seasok in ['DJF', 'MAM', 'JJA', 'SON'] for ru in allru]

        for ma, ax, subt in zip(mappeseas, axs.flatten(), subtitles):
            guplo = ma.mean('lon').plot.contourf(x = 'lat', y = 'plev', ax = ax, levels = 11, ylim = (1.e5, 1.e3), yscale = 'log')#vmax = )
            try:
                guplo.colorbar.set_label('')
            except Exception as exc:
                print(exc)
                pass

            if 'pi' not in subt:
                seasok = subt.split('-')[1].strip()
                guplo2 = mapmean[('pi', var)][copl].sel(season = seasok).mean('lon').plot.contour(x = 'lat', y = 'plev', ax = ax, levels = 11, ylim = (1.e5, 1.e3), yscale = 'log', colors="k", add_colorbar=False)#vmax = )
            #guplo.set_titles(template='{value}', maxchar = 13, fontsize = 12)
            ax.set_title(subt)

        for i in range(4):
            for j in range(4):
                if j > 0:
                    ax.set_ylabel('')

        fig.suptitle(var+' - '+copl)
        plt.tight_layout()

        figs_map.append(fig)

fignames = [var+'_'+copl for var in var_map_200 for copl in allcopls]
ctl.plot_pdfpages(cart_out + 'bottino_crossmeans.pdf', figs_map, True, fignames)
