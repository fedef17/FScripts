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
cart_out = cart + 'trends/'
ctl.mkdir(cart_out)

allru = ['pi', 'b025', 'b050', 'b100']
allnams = ['piControl', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']

colors = ['black', 'forestgreen', 'orange', 'violet']

allnams2 = allnams + ['ssp585', 'historical']
allru2 = allru + ['ssp585', 'hist']
colors2 = colors + ['indianred', 'steelblue']

allnams3 = allnams2 + ['stabilization-hist-1990']
allru3 = allru2 + ['b990']
colors3 = colors2 + ['teal']
####################################################################################################

glomeans, pimean, yeamean, mapmean = pickle.load(open(cart_in + 'bottino_seasmean_2D.p', 'rb'))

var = 'tas'

latbins = np.arange(-90, 91, 20)
lacol = ctl.color_set(9, use_seaborn=False)

for var in ['tas', 'pr']:
    fig, axs = plt.subplots(2, 4, figsize = (16,9))
    for ru, ax in zip(allru3, axs.flatten()):
        coso = yeamean[(ru, var)]
        coso1 = coso[:10].mean(axis = 0)
        cosoanom = coso-coso1
        ax.set_title(ru)
        smut = 50
        if ru in ['ssp585', 'hist']:
            smut = 20

        for la1, la2, col in zip(latbins[:-1], latbins[1:], lacol):
            print(la1,la2)
            cosolat = cosoanom.sel(lat = slice(la1, la2)).mean(['lat','lon'])
            cosmu = ctl.butter_filter(cosolat, smut)
            ax.plot(coso.year, cosmu, color = col)

    fig.savefig(cart_out + '{}_latbin_evol.pdf'.format(var))


## Climate stripes plot (lat/time)
import matplotlib as mpl
#cmappa = cm.get_cmap('viridis', 13)

# colo = '#a50026 #d73027 #f46d43 #fdae61 #fee090 #ffffff #e0f3f8 #abd9e9 #74add1 #4575b4 #313695'
# colo = colo.split()
# colo = colo[::-1]
# cmappa = mpl.colors.ListedColormap(colo)
cmappa = cm.get_cmap('RdBu_r', 17)
#cmappa.set_over('#800026') #662506
#cmappa.set_under('#023858') #542788

cmappa2 = cm.get_cmap('BrBG', 17)
# cmappa.set_over('#800026') #662506
# cmappa.set_under('#023858') #542788

latbins = np.arange(-90, 91, 10)
lacol = ctl.color_set(18, use_seaborn=False)


def difun(years, exp = 3, idiv = 100, imax = 550):
    outye = np.empty_like(years)
    okye = years < idiv
    outye[okye] = np.log(years[okye])
    okye2 = years >= idiv
    outye[okye2] = np.log(idiv) + exp * (np.log(years[okye2])-np.log(idiv))/(np.log(imax)-np.log(idiv))

    return outye

def invfun(outye, exp = 3, idiv = 100, imax = 550):
    years = np.empty_like(outye)
    yexp = np.exp(outye)
    okye = yexp < idiv
    years[okye] = np.exp(outye[okye])
    okye2 = ~okye

    years[okye2] = np.exp(np.log(idiv) + (np.log(imax)-np.log(idiv)) * (outye[okye2] - np.log(idiv))/exp)

    return years


funcsca = (difun, invfun)

masfi = '/nas/BOTTINO/masks.nc'
cose = xr.load_dataset(masfi)
oce_mask = cose['RnfA.msk'].values.astype('bool') # ocean mask: 1 over ocean, 0 over land

add_ssp = True
plot_b990 = True

cases = [True, True, False, False]
cases2 = [True, False, True, False]

for add_ssp, plot_b990 in zip(cases, cases2):
    assp = ''
    ab90 = ''
    if add_ssp: assp = '_wssp50'
    if plot_b990: ab90 = '_wb990'

    for tip in ['all', 'land', 'oce']:
        var = 'tas'
        cosopi = yeamean[('pi', var)]
        cosohist = yeamean[('hist', var)]
        cosossp = yeamean[('ssp585', var)]
        cosohistssp = xr.concat([cosohist, cosossp], dim = 'year')

        ypre = 30

        if plot_b990:
            fig, axs = plt.subplots(1, 4, figsize = (20,7))
            okru = ['b990'] + allru[1:]
        else:
            fig, axs = plt.subplots(1, 3, figsize = (16,7))
            okru = allru[1:]

        if tip == 'oce':
            cosohistssp = cosohistssp.where(oce_mask)
            cosopi = cosopi.where(oce_mask)
        elif tip == 'land':
            cosohistssp = cosohistssp.where(~oce_mask)
            cosopi = cosopi.where(~oce_mask)

        for ru, ax in zip(okru, axs.flatten()):
            print(ru)
            coso = yeamean[(ru, var)]

            if tip == 'oce':
                coso = coso.where(oce_mask)
            elif tip == 'land':
                coso = coso.where(~oce_mask)

            if add_ssp:
                # if ru == 'b990':
                #     coso = coso.assign_coords({'lat': cosohistssp.lat})
                #     coso = xr.concat([cosohistssp[-75:-25], coso], dim = 'year')
                # else:
                y0 = coso.year.values[0]
                cosolu = cosohistssp.sel(year = slice(y0-50, y0))
                coso = coso.assign_coords({'lat': cosohistssp.lat})
                coso = xr.concat([cosolu, coso], dim = 'year')

            # pio = cosossp.sel(year = slice(2015, int('2'+ru[1:])))
            # coso = xr.concat([pio, yeamean[(ru, var)]], dim = 'year')

            ax.set_title(ru)
            smut = 50

            matrix = []

            for la1, la2, col in zip(latbins[:-1], latbins[1:], lacol):
                print(la1,la2)
                cosolat = coso.sel(lat = slice(la1, la2)).mean(['lat','lon'])
                cosmu = ctl.butter_filter(cosolat, smut)

                # if add_ssp:
                #     base = float(cosobase.sel(year = slice(1984, 2014), lat = slice(la1, la2)).mean())
                # else:
                #     base = float(cosobase.sel(year = slice(coso.year[0]-ypre, coso.year[0]), lat = slice(la1, la2)).mean())
                base = float(cosopi.sel(lat = slice(la1, la2)).mean())

                relcha = (cosmu-base)/(np.mean(cosmu[-ypre:]) - base)
                matrix.append(relcha)

            matrix = np.stack(matrix)

            pizz = ax.imshow(matrix, vmin = 0, vmax = 1, aspect = 0.02, origin = 'lower', extent = [10, len(coso), -90, 90], cmap = cmappa)
            #ax.set_xscale('logit')

            ax.set_xscale('function', functions = funcsca)
            ax.set_xticks([20., 50, 100, 200, 300, 500])
            if add_ssp:
                ax.set_xticklabels([-30, 0, 50, 150, 250, 450])
                ax.axvline(0., color = 'grey', ls = '--')

            if ru == okru[0]:
                ax.set_yticks(np.arange(-90, 91, 30))
            else:
                ax.set_yticks([])

            #ax.imshow(matrix, vmin = 0, vmax = 1, aspect = 3, origin = 'lower', extent = [0, 500, -90, 90], cmap = cmappa)

        cax = plt.axes([0.1, 0.1, 0.8, 0.05])
        cb = plt.colorbar(pizz, cax=cax, orientation='horizontal', extend = 'both')
        cb.ax.tick_params(labelsize=18)
        cb.set_label('Realized change', fontsize=20)
        plt.subplots_adjust(left=0.1, bottom=0.17, right=0.98, top=0.86, wspace=0.05, hspace=0.20)

        fig.savefig(cart_out + '{}_ovmol_lattime_{}{}{}.pdf'.format(var, tip, assp, ab90))

        var = 'pr'

        cosopi = yeamean[('pi', var)]
        cosohist = yeamean[('hist', var)]
        cosossp = yeamean[('ssp585', var)]
        cosohistssp = xr.concat([cosohist, cosossp], dim = 'year')

        ypre = 30

        if plot_b990:
            #fig, axs = plt.subplots(2, 2, figsize = (16,9))
            fig, axs = plt.subplots(1, 4, figsize = (20,7))
            okru = ['b990'] + allru[1:]
        else:
            fig, axs = plt.subplots(1, 3, figsize = (16,7))
            okru = allru[1:]

        if tip == 'oce':
            cosohistssp = cosohistssp.where(oce_mask)
            cosopi = cosopi.where(oce_mask)
        elif tip == 'land':
            cosohistssp = cosohistssp.where(~oce_mask)
            cosopi = cosopi.where(~oce_mask)

        for ru, ax in zip(okru, axs.flatten()):
        #fig, axs = plt.subplots(1, 3, figsize = (16,7))
        #for ru, ax in zip(allru[1:], axs.flatten()):
        # for ru in allru[1:]:
        #     fig, ax = plt.subplots(figsize = (16,7))
            print(ru)
            coso = yeamean[(ru, var)]

            if tip == 'oce':
                coso = coso.where(oce_mask)
            elif tip == 'land':
                coso = coso.where(~oce_mask)

            if add_ssp:
                y0 = coso.year.values[0]
                cosolu = cosohistssp.sel(year = slice(y0-50, y0))
                coso = coso.assign_coords({'lat': cosohistssp.lat})
                coso = xr.concat([cosolu, coso], dim = 'year')

            ax.set_title(ru)
            smut = 50

            matrix = []

            for la1, la2, col in zip(latbins[:-1], latbins[1:], lacol):
                print(la1,la2)
                cosolat = coso.sel(lat = slice(la1, la2)).mean(['lat','lon'])
                cosmu = ctl.butter_filter(cosolat, smut)
                #ax.plot(coso.year, cosmu, color = col)

                #base = cosmu[0]
                #base = float(cosopi.sel(lat = slice(la1, la2)).mean())
                #base = float(cosohist[-20:].sel(lat = slice(la1, la2)).mean())

                # if add_ssp:
                #     base = float(cosobase.sel(year = slice(1984, 2014), lat = slice(la1, la2)).mean())
                # else:
                #     base = float(cosobase.sel(year = slice(coso.year[0]-ypre, coso.year[0]), lat = slice(la1, la2)).mean())
                base = float(cosopi.sel(lat = slice(la1, la2)).mean())

                # relcha = (cosmu-base)/(np.mean(cosmu[-ypre:]) - base)
                # if np.mean(cosmu[-ypre:]) - base < 0.001*base:
                #     relcha = np.nan*relcha
                relcha = (cosmu-base)/base
                matrix.append(relcha)

            matrix = np.stack(matrix)

            #divnorm = mpl.colors.TwoSlopeNorm(vmin=-0.2, vcenter=0., vmax=0.5)

            pizz = ax.imshow(matrix, aspect = 0.02, origin = 'lower', extent = [10, len(coso), -90, 90], cmap = cmappa2, vmin = -0.4, vmax = 0.4)#, norm = divnorm)
            #ax.set_xscale('logit')

            ax.set_xscale('function', functions = funcsca)
            ax.set_xticks([20., 50, 100, 200, 300, 500])
            if add_ssp:
                ax.set_xticklabels([-30, 0, 50, 150, 250, 450])
                ax.axvline(0., color = 'grey', ls = '--')

            if ru == okru[0]:
                ax.set_yticks(np.arange(-90, 91, 30))
            else:
                ax.set_yticks([])

            #ax.imshow(matrix, vmin = 0, vmax = 1, aspect = 3, origin = 'lower', extent = [0, 500, -90, 90], cmap = cmappa)
        cax = plt.axes([0.1, 0.1, 0.8, 0.05])
        cb = plt.colorbar(pizz, cax=cax, orientation='horizontal', extend = 'both')
        cb.ax.tick_params(labelsize=18)
        cb.set_label('Relative change', fontsize=20)
        plt.subplots_adjust(left=0.1, bottom=0.17, right=0.98, top=0.86, wspace=0.05, hspace=0.20)

        fig.savefig(cart_out + '{}_ovmol_lattime_{}{}{}.pdf'.format(var, tip, assp, ab90))


# ##### NOW for each sector
#
# var = 'tas'
# cosopi = yeamean[('pi', var)]
# cosohist = yeamean[('hist', var)]
# cosossp = yeamean[('ssp585', var)]
# cosobase = xr.concat([cosohist, cosossp], dim = 'year')
#
# ypre = 30
#
# fig, axs = plt.subplots(3, 3, figsize = (16,16))
# for iiru, ru in enumerate(allru[1:]):
#     coso = yeamean[(ru, var)]
#
#     # pio = cosossp.sel(year = slice(2015, int('2'+ru[1:])))
#     # coso = xr.concat([pio, yeamean[(ru, var)]], dim = 'year')
#
#     for ise, sect in enumerate([(-80, 40), (40, 160), (160, 280)]):
#         ax = axs[iiru, ise]
#
#         if iiru == 0:
#             ax.set_title('lon: {} -> {}'.format(sect[0], sect[1]))
#
#         #ax.set_title(ru)
#         smut = 50
#
#         matrix = []
#
#         print(coso.year[0]-ypre, coso.year[0])
#         for la1, la2, col in zip(latbins[:-1], latbins[1:], lacol):
#             #print(la1,la2)
#             if sect[0] > 0:
#                 cosolat = coso.sel(lat = slice(la1, la2), lon = slice(sect[0], sect[1])).mean(['lat','lon'])
#             else:
#                 cosolat1 = coso.sel(lat = slice(la1, la2), lon = slice(sect[0], 360))
#                 cosolat2 = coso.sel(lat = slice(la1, la2), lon = slice(0, sect[1]))
#                 cosolat = xr.concat([cosolat1, cosolat2], dim = 'lon').mean(['lat','lon'])
#
#             cosmu = ctl.butter_filter(cosolat, smut)
#
#             base = float(cosossp.sel(year = slice(coso.year[0]-ypre, coso.year[0]), lat = slice(la1, la2)).mean())
#             tutlat = coso.sel(lat = slice(la1, la2))[-ypre:].mean(['lat','lon','year']).values
#
#             # if sect[0] > 0:
#             #     base = float(cosossp.sel(year = slice(coso.year[0]-ypre, coso.year[0]), lat = slice(la1, la2), lon = slice(sect[0], sect[1])).mean())
#             # else:
#             #     co1 = cosossp.sel(year = slice(coso.year[0]-ypre, coso.year[0]), lat = slice(la1, la2), lon = slice(sect[0], 360)).mean()
#             #     co2 = cosossp.sel(year = slice(coso.year[0]-ypre, coso.year[0]), lat = slice(la1, la2), lon = slice(0, sect[1])).mean()
#             #     base = float(np.mean([co1, co2]))
#
#             relcha = (cosmu-base)/(tutlat - base)
#             #relcha = (cosmu-base)
#             matrix.append(relcha)
#
#         matrix = np.stack(matrix)
#
#         pizz = ax.imshow(matrix, vmin = 0, vmax = 1, aspect = 0.02, origin = 'lower', extent = [10, len(coso), -90, 90], cmap = cmappa)
#         #ax.set_xscale('logit')
#
#         ax.set_xscale('function', functions = funcsca)
#         ax.set_xticks([20., 50, 100, 200, 300, 500])
#         ax.set_yticks(np.arange(-90, 91, 30))
#
#         #ax.imshow(matrix, vmin = 0, vmax = 1, aspect = 3, origin = 'lower', extent = [0, 500, -90, 90], cmap = cmappa)
#
# cax = plt.axes([0.1, 0.1, 0.8, 0.05])
# cb = plt.colorbar(pizz, cax=cax, orientation='horizontal', extend = 'both')
# cb.ax.tick_params(labelsize=18)
# cb.set_label('Realized change', fontsize=20)
# plt.subplots_adjust(left=0.1, bottom=0.17, right=0.98, top=0.86, wspace=0.05, hspace=0.20)
#
# fig.savefig(cart_out + '{}_ovmol_lattime_sectorial.pdf'.format(var))
#
# var = 'pr'
#
# cosopi = yeamean[('pi', var)]
# cosohist = yeamean[('hist', var)]
# cosossp = yeamean[('ssp585', var)]
# cosobase = xr.concat([cosohist, cosossp], dim = 'year')
#
# ypre = 30
#
# fig, axs = plt.subplots(3, 3, figsize = (16,16))
# for iiru, ru in enumerate(allru[1:]):
#     coso = yeamean[(ru, var)]
#
#     # pio = cosossp.sel(year = slice(2015, int('2'+ru[1:])))
#     # coso = xr.concat([pio, yeamean[(ru, var)]], dim = 'year')
#
#     for ise, sect in enumerate([(-80, 40), (40, 160), (160, 280)]):
#         ax = axs[iiru, ise]
#
#         if iiru == 0:
#             ax.set_title('lon: {} -> {}'.format(sect[0], sect[1]))
#
#         #ax.set_title(ru)
#         smut = 50
#
#         matrix = []
#
#         for la1, la2, col in zip(latbins[:-1], latbins[1:], lacol):
#             print(la1,la2)
#             if sect[0] > 0:
#                 cosolat = coso.sel(lat = slice(la1, la2), lon = slice(sect[0], sect[1])).mean(['lat','lon'])
#             else:
#                 cosolat1 = coso.sel(lat = slice(la1, la2), lon = slice(sect[0], 360))
#                 cosolat2 = coso.sel(lat = slice(la1, la2), lon = slice(0, sect[1]))
#                 cosolat = xr.concat([cosolat1, cosolat2], dim = 'lon').mean(['lat','lon'])
#
#             cosmu = ctl.butter_filter(cosolat, smut)
#
#             if sect[0] > 0:
#                 base = float(cosossp.sel(year = slice(coso.year[0]-ypre, coso.year[0]), lat = slice(la1, la2)).mean())
#             else:
#                 co1 = cosossp.sel(year = slice(coso.year[0]-ypre, coso.year[0]), lat = slice(la1, la2), lon = slice(sect[0], 360)).mean()
#                 co2 = cosossp.sel(year = slice(coso.year[0]-ypre, coso.year[0]), lat = slice(la1, la2), lon = slice(0, sect[1])).mean()
#                 base = float(np.mean([co1, co2]))
#
#             relcha = (cosmu-base)/base
#             matrix.append(relcha)
#
#         matrix = np.stack(matrix)
#
#         divnorm = mpl.colors.TwoSlopeNorm(vmin=-0.05, vcenter=0., vmax=0.2)
#
#         pizz = ax.imshow(matrix, aspect = 0.02, origin = 'lower', extent = [10, len(coso), -90, 90], cmap = cmappa, norm = divnorm)
#
#         ax.set_xscale('function', functions = funcsca)
#         ax.set_xticks([20., 50, 100, 200, 300, 500])
#         ax.set_yticks(np.arange(-90, 91, 30))
#
#         #ax.imshow(matrix, vmin = 0, vmax = 1, aspect = 3, origin = 'lower', extent = [0, 500, -90, 90], cmap = cmappa)
#
# cax = plt.axes([0.1, 0.1, 0.8, 0.05])
# cb = plt.colorbar(pizz, cax=cax, orientation='horizontal', extend = 'both')
# cb.ax.tick_params(labelsize=18)
# cb.set_label('Relative change', fontsize=20)
# plt.subplots_adjust(left=0.1, bottom=0.17, right=0.98, top=0.86, wspace=0.05, hspace=0.20)
#
# fig.savefig(cart_out + '{}_ovmol_lattime_sectorial.pdf'.format(var))
