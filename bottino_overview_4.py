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
#         else:
#             coso = yeamean[(ru, var)]
#
#             var_trend, var_intercept, var_trend_err, var_intercept_err, var_pval = ctl.calc_trend_climatevar(g50, coso)
#
#             trendz[(var, ru)] = var_trend
#             trendz[(var, ru, 'err')] = var_trend_err
#             trendz[(var, ru, 'pval')] = var_pval
#
# pickle.dump(trendz, open(cart_out + 'trendz.p', 'wb'))

trendz = pickle.load(open(cart_out + 'trendz.p', 'rb'))

coso = yeamean[('b025', 'tas')]

#Figures
vcenall = [1., 0., 0., 0.]
for var, vcen in zip(['tas', 'pr', 'clt', 'rlut'], vcenall):
    figtrend = []
    fighatch = []
    tits = []

    for ru in allru[1:]:
        for cent, start, end in zip(['1st', '2nd', '5th'], [0, 100, 200], [100, 200, 500]):
            var_trend = trendz[(var, ru, cent)]
            figtrend.append(var_trend)
            fighatch.append(trendz[(var, ru, cent, 'pval')] < 0.01)
            tits.append(ru + ' - ' + cent)

    c5 = np.percentile(figtrend, 2)
    c95 = np.percentile(figtrend, 98)
    divnorm = colors.TwoSlopeNorm(vmin=c5, vcenter=vcen, vmax=c95)

    ctl.plot_multimap_contour(figtrend, coso.lat, coso.lon, filename = cart_out+var+'_trendz.pdf', fix_subplots_shape = (3,3), figsize = (16, 9), cbar_range = (c5, c95), color_norm = divnorm, subtitles = tits, add_hatching = fighatch)#, n_color_levels = 11)

    okfi = [fi for fi, tit in zip(figtrend, tits) if '1st' in tit]+[fi for fi, tit in zip(figtrend, tits) if '5th' in tit]
    oktit = [tit for tit in tits if '1st' in tit]+[tit for tit in tits if '5th' in tit]
    ctl.plot_multimap_contour(okfi, coso.lat, coso.lon, filename = cart_out+var+'_trendz_fastslow.pdf', fix_subplots_shape = (2,3), figsize = (16, 9), cbar_range = (c5, c95), color_norm = divnorm, subtitles = oktit, add_hatching = fighatch)#, n_color_levels = 11)

    for ru in allru2[-2:]:
        var_trend = trendz[(var, ru)]
        var_hatch = trendz[(var, ru, 'pval')] < 0.01

        ctl.plot_map_contour(var_trend, coso.lat, coso.lon, filename = cart_out+var+'_trend_{}.pdf'.format(ru), figsize = (16, 9), cbar_range = (c5, c95), color_norm = divnorm, add_hatching = var_hatch)#, n_color_levels = 11)