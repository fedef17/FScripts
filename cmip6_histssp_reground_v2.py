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
file_light = cart_in + 'out_NEW_cmip6_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr_light.p'

gen_file_ssp = cart_in + 'out_NEW_cmip6_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr.p'
#gen_file_ssp = cart_in + 'out_NEW_cmip6_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr_histrebase.p'

file_refit = cart_in + 'out_NEW_cmip6_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr_refit.p'
file_refit_2 = cart_in + 'out_NEW_cmip6_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr_refit_rebasetot.p'
fil_ece_ssp = cart_in + 'out_eceens_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr.p'

filref = '/home/fabiano/Research/lavori/WeatherRegimes/ERA_ref_r25_v4/out_ERA_NDJFM_{}_4clus_4pcs_1964-2014_dtr.p'

numclus = 4
reg_names_area = dict()
reg_names_area['EAT'] = ['NAO+', 'SBL', 'AR', 'NAO-']
reg_names_area['PNA'] = ['PT', 'PNA+', 'PNA-', 'AR']

area = 'EAT'
allssps = 'ssp126 ssp245 ssp370 ssp585'.split()
for area in ['EAT', 'PNA']:
    results_hist, results_ref = ctl.load_wrtool(file_hist.format(area))
    results_ref = pickle.load(open(filref.format(area), 'rb'))
    # del results_ref['var_glob']
    # del results_ref['var_area']
    # del results_ref['solver']
    for mod in results_hist.keys():
        #del results_hist[mod]['var_glob']
        del results_hist[mod]['var_area']
        del results_hist[mod]['solver']

    restot = dict()
    restot['models'] = results_hist
    restot['reference'] = results_ref


    for ssp in ['ssp585']:
        print(ssp)
        results_ssp, _ = ctl.load_wrtool(gen_file_ssp.format(ssp, area))
        for mod in results_ssp.keys():
            del results_ssp[mod]['var_area']
            del results_ssp[mod]['solver']

        res_rebase = dict()
        res_rebase_tot = dict()

        # Adding the ece ensemble
        ece_ssp, _ = ctl.load_wrtool(fil_ece_ssp.format(ssp, area))
        for mod in ece_ssp.keys():
            results_hist[mod] = results_hist['EC-Earth3_r1i1p1f1']
        del ece_ssp['EC-Earth3_r1i1p1f1']
        results_ssp.update(ece_ssp)

        for mod in results_ssp.keys():
            print(mod)
            if mod not in results_hist:
                print('skipping '+mod)
                continue

            bau = results_hist[mod]['var_dtr']
            bauda = np.arange(1965, 2015)

            if len(bau) != len(bauda):
                dates_set, _ = ctl.seasonal_set(results_hist[mod]['dates'], results_hist[mod]['dates'], 'NDJFM', seasonal_average = False)
                bauda = np.array([da[0].year for da in dates_set])

            gigi = results_ssp[mod]['var_dtr']
            gigida = np.arange(2015, 2100)
            if len(gigi) != len(gigida):
                dates_set, _ = ctl.seasonal_set(results_ssp[mod]['dates'], results_ssp[mod]['dates'], 'NDJFM', seasonal_average = False)
                gigida = np.array([da[0].year for da in dates_set])

            annette = np.concatenate([bauda, gigida])
            cosette = np.concatenate([bau, gigi])
            coeffs, covmat = np.polyfit(annette, cosette, deg = 3, cov = True)
            new_fit = np.polyval(coeffs, annette)

            hist_var = ctl.restore_fullfield_from_anomalies_daily(results_hist[mod]['var_glob'], results_hist[mod]['dates'], results_hist[mod]['climate_mean'], results_hist[mod]['climate_mean_dates'])
            ssp_var = ctl.restore_fullfield_from_anomalies_daily(results_ssp[mod]['var_glob'], results_ssp[mod]['dates'], results_ssp[mod]['climate_mean'], results_ssp[mod]['climate_mean_dates'])

            var_mod = np.concatenate([hist_var, ssp_var], axis = 0)
            old_fit = np.concatenate([np.polyval(results_hist[mod]['coeffs_dtr'], bauda), np.polyval(results_ssp[mod]['coeffs_dtr'], gigida)])

            diffit = new_fit - old_fit

            dates_mod = np.concatenate([results_hist[mod]['dates'], results_ssp[mod]['dates']])

            var_set, dates_set = ctl.seasonal_set(var_mod, dates_mod, 'NDJFM', seasonal_average = False)
            var_mod_new = []
            for va, ye, cos in zip(var_set, annette, diffit):
                var_mod_new.append(va + cos)

            var_mod = np.concatenate(var_mod_new, axis = 0)
            print(np.nanmean(var_mod))

            # Corretto la discrepanza nel fit. Ora faccio due prove: una con la hist climate mean, l'altra facendo ricalcolare la climate_mean, quindi con una climate_mean media tra hist e ssp.
            climate_rebase = results_hist[mod]['climate_mean']

            reres = cd.WRtool_core(var_mod, results_ssp[mod]['lat'], results_ssp[mod]['lon'], dates_mod, area, wnd_days = 20, numpcs = 4, numclus = 4, ref_solver = results_ref['solver'], ref_patterns_area = results_ref['cluspattern_area'], heavy_output = False, run_significance_calc = False, use_reference_eofs = True, use_reference_clusters = True, ref_clusters_centers = results_ref['centroids'], climate_mean = climate_rebase, dates_climate_mean = results_hist[mod]['climate_mean_dates'])

            res_rebase[mod] = reres

            reres = cd.WRtool_core(var_mod, results_ssp[mod]['lat'], results_ssp[mod]['lon'], dates_mod, area, wnd_days = 20, numpcs = 4, numclus = 4, ref_solver = results_ref['solver'], ref_patterns_area = results_ref['cluspattern_area'], heavy_output = False, run_significance_calc = False, use_reference_eofs = True, use_reference_clusters = True, ref_clusters_centers = results_ref['centroids'], climate_mean = None, dates_climate_mean = None)

            res_rebase_tot[mod] = reres

            if 'EC-Earth3' not in mod: del results_hist[mod]['var_glob']
            del results_ssp[mod]['var_glob']

        restot = dict()
        restot['models'] = res_rebase
        restot['reference'] = results_ref
        with open(file_refit.format(ssp, area), 'wb') as fillo:
            pickle.dump(restot, fillo)

        restot = dict()
        restot['models'] = res_rebase_tot
        restot['reference'] = results_ref
        with open(file_refit_2.format(ssp, area), 'wb') as fillo:
            pickle.dump(restot, fillo)
