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

#############################################################################
cart_in = '/data-hobbes/fabiano/WR_CMIP6/'
cart_out_orig = '/home/fabiano/Research/lavori/CMIP6/Results_histssp_reground/'
ctl.mkdir(cart_out_orig)

file_hist = cart_in + 'out_NEW_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refCLUS_dtr.p'
file_hist_light = cart_in + 'out_NEW_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refCLUS_dtr_light.p'
gen_file_ssp = cart_in + 'out_NEW_cmip6_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr.p'
file_rebase = cart_in + 'out_NEW_cmip6_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr_histrebase.p'
file_light = cart_in + 'out_NEW_cmip6_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr_light.p'

filref = '/home/fabiano/Research/lavori/WeatherRegimes/ERA_ref_r25_v4/out_ERA_NDJFM_{}_4clus_4pcs_1964-2014_dtr.p'

numclus = 4
reg_names_area = dict()
reg_names_area['EAT'] = ['NAO+', 'SBL', 'AR', 'NAO-']
reg_names_area['PNA'] = ['PT', 'PNA+', 'PNA-', 'AR']

area = 'EAT'
allssps = 'ssp126 ssp245 ssp370 ssp585'.split()
for area in ['EAT']:#, 'PNA']:
    results_hist, results_ref = ctl.load_wrtool(file_hist.format(area))
    results_ref = pickle.load(open(filref.format(area), 'rb'))
    # del results_ref['var_glob']
    # del results_ref['var_area']
    # del results_ref['solver']
    for mod in results_hist.keys():
        del results_hist[mod]['var_glob']
        del results_hist[mod]['var_area']
        del results_hist[mod]['solver']

    restot = dict()
    restot['models'] = results_hist
    restot['reference'] = results_ref
    with open(file_hist_light.format(area), 'wb') as fillo:
        pickle.dump(restot, fillo)

    for ssp in allssps:
        print(ssp)
        results_ssp, _ = ctl.load_wrtool(gen_file_ssp.format(ssp, area))
        res_rebase = dict()
        for mod in results_ssp.keys():
            print(mod)
            if mod not in results_hist:
                print('skipping '+mod)
                continue

            # riaggiungi global mean futura e togli quella storica. Riassegna gli stati e salva tutto in res_rebase
            var_mod = results_ssp[mod]['var_glob']
            climate_rebase = results_hist[mod]['climate_mean']-results_ssp[mod]['climate_mean']
            #var_mod_new = ctl.anomalies_daily(var_mod, results_ssp[mod]['dates'], climate_mean = climate_rebase, dates_climate_mean = results_ssp[mod]['climate_mean_dates'])

            reres = cd.WRtool_core(var_mod, results_ssp[mod]['lat'], results_ssp[mod]['lon'], results_ssp[mod]['dates'], area, wnd_days = 20, numpcs = 4, numclus = 4, ref_solver = results_ref['solver'], ref_patterns_area = results_ref['cluspattern_area'], heavy_output = False, run_significance_calc = False, use_reference_eofs = True, use_reference_clusters = True, ref_clusters_centers = results_ref['centroids'], climate_mean = climate_rebase, dates_climate_mean = results_ssp[mod]['climate_mean_dates'])

            res_rebase[mod] = reres

        restot = dict()
        restot['models'] = res_rebase
        restot['reference'] = results_ref
        with open(file_rebase.format(ssp, area), 'wb') as fillo:
            pickle.dump(restot, fillo)

        for mod in results_ssp.keys():
            del results_ssp[mod]['var_glob']
            del results_ssp[mod]['var_area']
            del results_ssp[mod]['solver']

        restot = dict()
        restot['models'] = results_ssp
        restot['reference'] = results_ref
        with open(file_light.format(ssp, area), 'wb') as fillo:
            pickle.dump(restot, fillo)
