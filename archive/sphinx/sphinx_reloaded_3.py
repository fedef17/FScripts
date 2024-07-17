#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os

from matplotlib import pyplot as plt
import netCDF4 as nc
import pandas as pd
from numpy import linalg as LA
from scipy import stats
from scipy import interpolate as itrp
import itertools as itt

from datetime import datetime
import pickle

import climtools_lib as ctl
import climdiags as cd

from copy import deepcopy as cp
###############################################

# cart = '/data-woodstock/SPHINX/extended/'
# filename = '1850-2160-{}-annual.nc'

cart = '/data-woodstock/fede/'
filename = '1850-2160-{}-annual_remap25_rechunk.nc'

cart_out = '/home/fabiano/Research/lavori/SPHINX_for_lisboa/sphinx_reloaded/'

ens_mem = ['lcb0', 'lcs0', 'lcb1', 'lcs1', 'lcb2', 'lcs2']

# leggo: hcc, mcc, lcc, tas

figures = []
#lat, lon = ctl.genlatlon(73, 144)
lat, lon = ctl.genlatlon(256, 512, lon_limits = (0., 359.3), lat_limits = (-89.46, 89.46))

periods = [(1900, 2005), (2005, 2050), (2050, 2090), (2090, 2110), (2110, 2160)]
allvars = ['tas', 'hcc', 'mcc', 'lcc', 'rlut', 'rsus']

#deg_periods = [(0., 1.), (1.,2.), (2.,3.), (3.,4.), (4., 5.)]

# da fare:
# diff plots hcc, mcc, lcc vs tas con colore dato dall'anno. Qui posso usare gli annual. Separare tropici?
# ovmoller lat/temp ?
# trova warming pattern per grado di temp

varnam = 'tas'
tas, lat, lon, dates, time_units, var_units, time_cal = ctl.readxDncfield(cart + filename.format(varnam))

# periodi

# globalT = dict()
# for ens in ens_mem:
#     globalT[ens] = ctl.global_mean(tas[ens], lat)
#
# globalT['base'] = np.mean([globalT[ens] for ens in ens_mem if 'lcb' in ens], axis = 0)
# globalT['stoc'] = np.mean([globalT[ens] for ens in ens_mem if 'lcs' in ens], axis = 0)
# globalT['base_std'] = np.std([globalT[ens] for ens in ens_mem if 'lcb' in ens], axis = 0)
# globalT['stoc_std'] = np.std([globalT[ens] for ens in ens_mem if 'lcs' in ens], axis = 0)
#
# trenddict = dict()
#
# for varnam in allvars:
#     var, lat, lon, dates, time_units, var_units, time_cal = ctl.readxDncfield(cart + filename.format(varnam))
#
#     for per in periods:
#         for ens in ens_mem:
#             var_ok, dates_ok = ctl.sel_time_range(var[ens], dates, ctl.range_years(per[0], per[1]))
#             glot, _ = ctl.sel_time_range(globalT[ens], dates, ctl.range_years(per[0], per[1]))
#
#             var_trend, var_trend_err = ctl.calc_trend_climatevar(glot, var_ok)
#             trenddict[(varnam, per, ens)] = var_trend
#             trenddict[(varnam, per, ens, 'err')] = var_trend_err
#
#         trenddict[(varnam, per, 'base')] = np.mean([trenddict[(varnam, per, ens)] for ens in ens_mem if 'lcb' in ens], axis = 0)
#         trenddict[(varnam, per, 'stoc')] = np.mean([trenddict[(varnam, per, ens)] for ens in ens_mem if 'lcs' in ens], axis = 0)
#         trenddict[(varnam, per, 'diff')] = trenddict[(varnam, per, 'stoc')] - trenddict[(varnam, per, 'base')]
#         trenddict[(varnam, per, 'base_std')] = np.std([trenddict[(varnam, per, ens)] for ens in ens_mem if 'lcb' in ens], axis = 0)
#         trenddict[(varnam, per, 'stoc_std')] = np.std([trenddict[(varnam, per, ens)] for ens in ens_mem if 'lcs' in ens], axis = 0)
#         trenddict[(varnam, per, 'diff_std')] = trenddict[(varnam, per, 'stoc_std')] + trenddict[(varnam, per, 'base_std')]
#
# pickle.dump([globalT, trenddict], open(cart_out + 'warming_patt_v2.p', 'wb'))
globalT, trenddict = pickle.load(open(cart_out + 'warming_patt_v2.p'))

for per in periods:
    for cos in ['base', 'stoc', 'diff']:
        trenddict[('toa_net', per, cos)] = -(trenddict[('rsus', per, cos)]+trenddict[('rlut', per, cos)])

for varnam in allvars + ['toa_net']:
    figs = []
    for per in periods:
        for cos in ['base', 'stoc', 'diff']:
            if varnam == 'tas':
                cbar_range = (-2., 2.)
            elif 'cc' in varnam:
                cbar_range = (-0.09, 0.09)
            else:
                cbar_range = (-20., 20.)
            if cos == 'diff' and cbar_range is not None: cbar_range = (cbar_range[0]/4., cbar_range[1]/4.)
            fig = ctl.plot_map_contour(trenddict[(varnam, per, cos)], lat, lon, visualization= 'Robinson', central_lat_lon=(0., 180.), title = '{} - {} - {}'.format(varnam, per, cos), cbar_range = cbar_range)
            figs.append(fig)

    ctl.plot_pdfpages(cart_out + 'trend_patterns_v2_{}.pdf'.format(varnam), figs)
