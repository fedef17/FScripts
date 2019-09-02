#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
from matplotlib import pyplot as plt
from matplotlib import cm
import pickle

import climtools_lib as ctl
import climdiags as cd

from matplotlib.colors import LogNorm

from scipy import stats

#######################################

season = 'DJF'
area = 'EAT'

cart_out = '/home/fabiano/Research/lavori/WeatherRegimes/test_ERA_reclustering/'
if not os.path.exists(cart_out): os.mkdir(cart_out)

ctl.openlog(cart_out, tag = 'test_ERA_stable')

# ctl.openlog('.')
file_in = '/data-hobbes/fabiano/OBS/ERA/ERA40+Int_daily_1957-2018_zg500_remap25_meters.nc'

print(ctl.datestamp())

kwar = dict()
kwar['numclus'] = 4
kwar['run_significance_calc'] = False
kwar['numpcs'] = 4
kwar['detrended_eof_calculation'] = False # detrendo io all'inizio
kwar['detrended_anom_for_clustering'] = False
kwar['nrsamp_sig'] = 500

var, coords, aux_info = ctl.read_iris_nc(file_in, extract_level_hPa = 500)
lat = coords['lat']
lon = coords['lon']
dates = coords['dates']

var, dates = ctl.sel_time_range(var, dates, ctl.range_years(1957, 2014))

var_anoms = ctl.anomalies_daily_detrended(var, dates)
var_season, dates_season = ctl.sel_season(var_anoms, dates, season)
all_years = np.arange(dates[0].year, dates[-1].year+1)

results_ref = cd.WRtool_core(var_season, lat, lon, dates_season, area, heavy_output = True, **kwar)
kwar['ref_solver'] = results_ref['solver']
kwar['ref_patterns_area'] = results_ref['cluspattern_area']

all_results = dict()
for i in range(100):
    print(i)
    result = cd.WRtool_core(var_season, lat, lon, dates_season, area, heavy_output = False, **kwar)
    all_results['try {}'.format(i)] = result

all_names = ['ref'] + ['try {}'.format(i) for i in range(20)]
all_results['ref'] = results_ref

pickle.dump(all_results, open(cart_out + 'test_ERA_stable.p', 'wb'))

# file_res = cart_out + 'results_test_ERA_stable.dat'
# cd.out_WRtool_mainres(file_res, all_results, results_ref, kwar)

ctl.plot_multimodel_regime_pdfs(all_results, model_names = all_names, filename = cart_out + 'regimes_testERA_stable.pdf', eof_proj = [(0,1), (1,2), (2,3), (3,0)], figsize = (20,12), reference = 'ref', eof_axis_lim = (-3000, 3000), check_for_eofs = True)
