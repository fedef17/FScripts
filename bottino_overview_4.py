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
####################################################################################################

glomeans, pimean, yeamean, mapmean = pickle.load(open(cart_in + 'bottino_seasmean_2D.p', 'rb'))

# century trends. or bicentury?
from matplotlib import colors

# trendz = dict()
# for var in ['tas', 'pr', 'clt', 'rlut']:
#     print(var)
#     for ru in allru2:
#         print(ru)
#         gtas = glomeans[(ru, var)][1]
#         g50 = ctl.butter_filter(gtas, 50)
#
#         if ru[0] == 'b':
#             for cent, start, end in zip(['1st', '2nd', '5th'], [0, 100, 200], [100, 200, 500]):
#                 print(cent)
#                 coso = yeamean[(ru, var)][start:end]
#                 # print(coso.shape, gtas.shape)
#
#                 var_trend, var_intercept, var_trend_err, var_intercept_err, var_pval = ctl.calc_trend_climatevar(g50[start:end], coso)
#
#                 trendz[(var, ru, cent)] = var_trend
#                 trendz[(var, ru, cent, 'err')] = var_trend_err
#                 trendz[(var, ru, cent, 'pval')] = var_pval
#
#         coso = yeamean[(ru, var)]
#
#         var_trend, var_intercept, var_trend_err, var_intercept_err, var_pval = ctl.calc_trend_climatevar(g50, coso)
#
#         trendz[(var, ru)] = var_trend
#         trendz[(var, ru, 'err')] = var_trend_err
#         trendz[(var, ru, 'pval')] = var_pval
#
# pickle.dump(trendz, open(cart_out + 'trendz.p', 'wb'))

trendz2 = dict()
for var in ['tas', 'pr', 'clt', 'rlut']:
    print(var)
    for ru in allru2:
        print(ru)
        if ru[0] != 'b': continue

        gtas = glomeans[(ru, var)][1]
        g50 = ctl.butter_filter(gtas, 50)

        coso = yeamean[(ru, var)]

        var_trend, var_intercept, var_trend_err, var_intercept_err, var_pval = ctl.calc_trend_climatevar(g50, coso)

        trendz2[(var, ru)] = var_trend
        trendz2[(var, ru, 'err')] = var_trend_err
        trendz2[(var, ru, 'pval')] = var_pval

trendz = pickle.load(open(cart_out + 'trendz.p', 'rb'))
trendz.update(trendz2)

pickle.dump(trendz, open(cart_out + 'trendz.p', 'wb'))

# trendz = pickle.load(open(cart_out + 'trendz.p', 'rb'))

coso = yeamean[('b025', 'tas')]

#Figures
vcenall = [1., 0., 0., 0.]
cmapall = ['RdBu_r', 'BrBG', 'BrBG', 'RdBu_r']
for var, vcen, cmap in zip(['tas', 'pr', 'clt', 'rlut'], vcenall, cmapall):
    figtrend = []
    fighatch = []
    tits = []

    for ru in allru[1:]:
        for cent, start, end in zip(['1st', '2nd', '5th'], [0, 100, 200], [100, 200, 500]):
            var_trend = trendz[(var, ru, cent)]
            figtrend.append(var_trend)
            fighatch.append(trendz[(var, ru, cent, 'pval')] < 0.05)
            tits.append(ru + ' - ' + cent)

    c5 = np.percentile(figtrend, 5)
    c95 = np.percentile(figtrend, 95)
    divnorm = colors.TwoSlopeNorm(vmin=c5, vcenter=vcen, vmax=c95)

    ctl.plot_multimap_contour(figtrend, coso.lat, coso.lon, filename = cart_out+var+'_trendz.pdf', fix_subplots_shape = (3,3), figsize = (16, 9), cbar_range = (c5, c95), color_norm = divnorm, subtitles = tits, add_hatching = fighatch, cmap = cmap, n_color_levels = 13, hatch_styles = ['////', '', ''])

    okfi = [fi for fi, tit in zip(figtrend, tits) if '1st' in tit]+[fi for fi, tit in zip(figtrend, tits) if '5th' in tit]
    oktit = [tit for tit in tits if '1st' in tit]+[tit for tit in tits if '5th' in tit]
    okha = [fi for fi, tit in zip(fighatch, tits) if '1st' in tit]+[fi for fi, tit in zip(fighatch, tits) if '5th' in tit]

    ctl.plot_multimap_contour(okfi, coso.lat, coso.lon, filename = cart_out+var+'_trendz_fastslow.pdf', fix_subplots_shape = (2,3), figsize = (16, 9), cbar_range = (c5, c95), color_norm = divnorm, subtitles = oktit, add_hatching = okha, cmap = cmap, n_color_levels = 13, hatch_styles = ['////', '', ''])

    for ru in ['pi']+allru2[-2:]:
        var_trend = trendz[(var, ru)]
        var_hatch = trendz[(var, ru, 'pval')] < 0.05

        ctl.plot_map_contour(var_trend, coso.lat, coso.lon, filename = cart_out+var+'_trend_{}.pdf'.format(ru), figsize = (16, 9), cbar_range = (c5, c95), color_norm = divnorm, add_hatching = var_hatch, cmap = cmap, n_color_levels = 13, hatch_styles = ['///', '', ''])

    #### Now fig with ssp585 and all botts
    okfi = [trendz[(var, ru)] for ru in ['ssp585'] + allru[1:]]
    okha = [trendz[(var, ru, 'pval')] < 0.05 for ru in ['ssp585'] + allru[1:]]
    oktit = ['ssp585'] + allru[1:]

    ctl.plot_multimap_contour(okfi, coso.lat, coso.lon, filename = cart_out+var+'_trendzfull_vs_ssp.pdf', fix_subplots_shape = (2,3), figsize = (16, 9), cbar_range = (c5, c95), color_norm = divnorm, subtitles = oktit, add_hatching = okha, cmap = cmap, n_color_levels = 13, hatch_styles = ['////', '', ''])
