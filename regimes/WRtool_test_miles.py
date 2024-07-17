#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
from matplotlib import pyplot as plt
import pickle

import climtools_lib as ctl
import climdiags as cd

#######################################

cart_in = '/home/fabiano/Research/lavori/test_ensclus/miles_block/preproc/miles_diagnostics_preproc1_zg/'
cart_out = '/home/fabiano/Research/lavori/test_ensclus/'
if not os.path.exists(cart_out): os.mkdir(cart_out)

filo = cart_in + 'CMIP5_EC-EARTH_day_historical_r1i1p1_T2Ds_zg_1980-1989.nc'
file_era = cart_in + 'OBS_ERA-Interim_reanaly_1_T2Ds_zg_1980-1989.nc'

seas = 'DJF'
area = 'EAT'

ERA_ref = cd.WRtool_from_file(file_era, seas, area, extract_level_4D = 50000., numclus = 4, heavy_output = True, run_significance_calc = False)
model = cd.WRtool_from_file(filo, seas, area, extract_level_4D = 50000., numclus = 4, heavy_output = True, run_significance_calc = False, ref_solver = ERA_ref['solver'], ref_patterns_area = ERA_ref['cluspattern_area'])

pickle.dump([ERA_ref, model], open(cart_out+'out_WRtool.p', 'wb'))

ERA_ref, model = pickle.load(open(cart_out+'out_WRtool.p', 'rb'))
cd.plot_WRtool_results(cart_out, 'test_miles', 1, model, ERA_ref, model_name = 'ECE', obs_name = 'ERA', patnames = ['Sc.Blocking', 'NAO -', 'NAO +', 'Atl. Ridge'], patnames_short = ['BL', 'NN', 'NP', 'AR'])
