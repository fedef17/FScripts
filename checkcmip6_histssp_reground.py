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
    del results_ref['var_glob']
    del results_ref['var_area']
    del results_ref['solver']

    for mod in results_hist.keys():
        if mod == 'MPI-ESM1-2-LR_r1i1p1f1':
            continue
        del results_hist[mod]['var_glob']
        del results_hist[mod]['var_area']
        del results_hist[mod]['solver']

    results_ssp, _ = ctl.load_wrtool(gen_file_ssp.format('ssp585', area))
    for mod in results_ssp.keys():
        if mod == 'MPI-ESM1-2-LR_r1i1p1f1':
            continue
        del results_ssp[mod]['var_glob']
        del results_ssp[mod]['var_area']
        del results_ssp[mod]['solver']

    # riaggiungi global mean futura e togli quella storica. Riassegna gli stati e salva tutto in res_rebase
    mod = 'MPI-ESM1-2-LR_r1i1p1f1'
    var_mod = results_ssp[mod]['var_glob']
    climate_rebase = results_hist[mod]['climate_mean']-results_ssp[mod]['climate_mean']

    fig = plt.figure()
    plt.plot(results_hist[mod]['climate_mean'][:, 50:70, -8])
    plt.plot(results_ssp[mod]['climate_mean'][:, 50:70, -8], linestyle = '--')
    plt.plot(climate_rebase[:, 50:70, -8], linestyle = ':')
    fig.savefig(cart_out_orig + 'baugigicheck.pdf')

    bau = results_hist[mod]['var_dtr']
    bauda = np.arange(1965, 2015)
    gigi = results_ssp[mod]['var_dtr']
    gigida = np.arange(2015, 2100)
    fig = plt.figure()
    plt.plot(bauda, bau)
    plt.plot(bauda, np.polyval(results_hist[mod]['coeffs_dtr'], bauda), linestyle = '--')
    plt.plot(gigida, gigi)
    plt.plot(gigida, np.polyval(results_ssp[mod]['coeffs_dtr'], gigida), linestyle = '--')
    fig.savefig(cart_out_orig + 'vardtr_check.pdf')

    fig = plt.figure()
    cols = ctl.color_set(len(results_ssp.keys()))
    for mod, co in zip(results_ssp.keys(), cols):
        try:
            bau = results_hist[mod]['var_dtr']
            bauda = np.arange(1965, 2015)
            gigi = results_ssp[mod]['var_dtr']
            gigida = np.arange(2015, 2100)
            plt.scatter(bauda, bau, color = co, s = 2)
            plt.plot(bauda, np.polyval(results_hist[mod]['coeffs_dtr'], bauda), color = co, linewidth = 0.5)
            plt.scatter(gigida, gigi, color = co, s = 2)
            plt.plot(gigida, np.polyval(results_ssp[mod]['coeffs_dtr'], gigida), color = co, linewidth = 0.5)
        except:
            pass
    fig.savefig(cart_out_orig + 'vardtr_check_allmods.pdf')