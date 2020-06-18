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
import pandas as pd
#import pymannkendall as mk

#############################################################################
if os.uname()[1] == 'hobbes':
    cart_in = '/home/fabiano/Research/lavori/WeatherRegimes/'
    cart_out = '/home/fabiano/Research/lavori/prima_D45/'
elif os.uname()[1] == 'ff-clevo':
    cart_out = '/home/fedefab/Scrivania/Research/Post-doc/lavori/prima_D45/'

fil_pres = 'prima_D45_pres/out_prima_D45_pres_DJF_EAT_4clus_4pcs_1979-2014_refEOF_dtr.p'
fil_fut = 'prima_D45_fut/out_prima_D45_fut_DJF_EAT_4clus_4pcs_2015-2050_refEOF_dtr.p'

results_pres, results_ref = ctl.load_wrtool(cart_in + fil_pres)
results_fut, _ = ctl.load_wrtool(cart_in + fil_fut)

cart_data = '/nas/PRIMAVERA/Stream1/'
filtas = 'highresSST-{}/{}/{}/day/{}/{}_day_{}_highresSST-{}_{}_gr_*_r25_rc.nc'
#filtas.format(temp, mod, mem, varnam, varnam, mod, temp, mem)

composites = dict()
mod = 'ref'
temp = 'present'
file_ref = dict()
file_ref['tas'] = '/data-hobbes/fabiano/OBS/ERA/ERAInterim/ERAInt_daily_1979-2018_167_r25.nc'
file_ref['pr'] = '/data-hobbes/fabiano/OBS/ERA/ERAInterim/ERAInt_daily_1979-2018_228_pr_daysum_ok_r25.nc'
for varnam in ['tas', 'pr']:
    var, coords, aux_info = ctl.read_iris_nc(file_ref[varnam])
    var_season, dates_season = ctl.sel_season(var, coords['dates'], 'DJF')

    for comp_moment in ['mean', 'std', 'percentile']:
        comps = ctl.composites_regimes_daily(coords['lat'], coords['lon'], var_season, dates_season, results_ref['labels'], results_ref['dates'], comp_moment = comp_moment, detrend_global = True, area_dtr = 'NML')
        composites[(temp, mod, varnam, comp_moment)] = comps

# results_all = dict()
# for ke in results_pres:
#     mod = ke.split('_')[0]
#     results_all[mod+'_pres'] = results_pres[ke]
# for ke in results_fut:
#     mod = ke.split('_')[0]
#     results_all[mod+'_fut'] = results_pres[ke]

for ke in results_pres:
    temp = 'present'
    mod, mem = ke.split('_')
    for varnam in ['tas', 'pr']:
        listafilo = ctl.find_file(cart_data, filtas.format(temp, mod, mem, varnam, varnam, mod, temp, mem))
        if len(listafilo) > 0:
            filo = listafilo[0]
        else:
            print('NOT FOUND: ', temp, mod, varnam)
            continue

        var, coords, aux_info = ctl.read_iris_nc(filo)
        var_season, dates_season = ctl.sel_season(var, coords['dates'], 'DJF')

        for comp_moment in ['mean', 'std', 'percentile']:
            comps = ctl.composites_regimes_daily(coords['lat'], coords['lon'], var_season, dates_season, results_pres[ke]['labels'], results_pres[ke]['dates'], comp_moment = comp_moment, detrend_global = True, area_dtr = 'NML')
            composites[(temp, mod, varnam, comp_moment)] = comps


for ke in results_fut:
    temp = 'future'
    mod, mem = ke.split('_')
    for varnam in ['tas', 'pr']:
        listafilo = ctl.find_file(cart_data, filtas.format(temp, mod, mem, varnam, varnam, mod, temp, mem))
        if len(listafilo) > 0:
            filo = listafilo[0]
        else:
            print('NOT FOUND: ', temp, mod, varnam)
            continue

        var, coords, aux_info = ctl.read_iris_nc(filo)
        var_season, dates_season = ctl.sel_season(var, coords['dates'], 'DJF')

        for comp_moment in ['mean', 'std', 'percentile']:
            comps = ctl.composites_regimes_daily(coords['lat'], coords['lon'], var_season, dates_season, results_fut[ke]['labels'], results_fut[ke]['dates'], comp_moment = comp_moment, detrend_global = True, area_dtr = 'NML')
            composites[(temp, mod, varnam, comp_moment)] = comps

with open(cart_out + 'composites_taspr.p', 'wb') as filonz:
    pickle.dump(composites, filonz)
