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
import pymannkendall as mk
import glob

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 18
#############################################################################
if os.uname()[1] == 'hobbes':
    cart_in = '/home/fabiano/Research/lavori/CMIP6/'
elif os.uname()[1] == 'ff-clevo':
    cart_in = '/home/fedefab/Scrivania/Research/Post-doc/lavori/CMIP6/'

cart_tas = '/data-woodstock/CMIP6/tasOld/ssp585/remap/'
fil_tas = cart_tas + 'tas_Amon_{}_ssp585_{}_201501-210012_remap.nc'

cart_out = cart_in + 'Results_SST_corrmap/'
ctl.mkdir(cart_out)

cart_data = '/data-hobbes/fabiano/WR_CMIP6/'
file_hist_refEOF = cart_data + 'out_NEW_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refEOF_dtr.p'
file_hist = cart_data + 'out_NEW_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refCLUS_dtr_light.p'
gen_file_ssp = cart_data + 'out_NEW_cmip6_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr_histrebase.p'

yr10 = 10 # length of running mean

numclus = 4
reg_names_area = dict()
reg_names_area['EAT'] = ['NAO+', 'SBL', 'AR', 'NAO-']
reg_names_area['PNA'] = ['PT', 'PNA+', 'PNA-', 'AR']
#okmods = ['ACCESS-CM2_r1i1p1f1', 'BCC-CSM2-MR_r1i1p1f1', 'CanESM5_r1i1p1f1', 'CESM2-WACCM_r1i1p1f1','CNRM-CM6-1-HR_r1i1p1f2', 'CNRM-CM6-1_r1i1p1f2','CNRM-ESM2-1_r1i1p1f2', 'EC-Earth3_r1i1p1f1', 'FGOALS-g3_r1i1p1f1','INM-CM4-8_r1i1p1f1', 'INM-CM5-0_r1i1p1f1', 'IPSL-CM6A-LR_r1i1p1f1','MIROC6_r1i1p1f1', 'MPI-ESM1-2-HR_r1i1p1f1','MPI-ESM1-2-LR_r1i1p1f1', 'MRI-ESM2-0_r1i1p1f1','NorESM2-LM_r1i1p1f1', 'NorESM2-MM_r1i1p1f1', 'UKESM1-0-LL_r1i1p1f2']
# mancano = ['BCC-ESM1_r1i1p1f1', 'CESM2_r1i1p1f1', 'GFDL-CM4_r1i1p1f1', 'HadGEM3-GC31-LL_r1i1p1f3', 'KACE-1-0-G_r1i1p1f1', 'MPI-ESM-1-2-HAM_r1i1p1f1']

allssps = 'ssp126 ssp245 ssp370 ssp585'.split()
allsimcol = ['hist', 'ssp126', 'ssp245', 'bau', 'ssp370', 'ssp585', 'rcp85_cmip5']
coldic = dict(zip(allsimcol, ctl.color_set(7)))
colssp = [coldic[ssp] for ssp in allssps]

# tas_anom = dict()
# tas_trends = dict()
# for ssp in ['ssp585']:
#     print('SSP '+ssp)
#
#     area = 'EAT'
#     cart_lui = cart_in + 'Results_v5_rebase/{}_NDJFM/'.format(area)
#     freqs, residtimes, patterns = pickle.load(open(cart_lui + 'allresults_dicts_{}_v3.p'.format(area), 'rb'))
#
#     okmods = [ke[1] for ke in freqs if ssp in ke and 'tot50' in ke and 'all' not in ke and 'rel' not in ke]
#     if ssp == 'ssp126':
#         okmods = [mod for mod in okmods if mod != 'FGOALS-g3_r1i1p1f1']
#     print(okmods)
#
#     for mod in okmods:
#         # lisf = glob.glob(fil_tas.format(mod, '*'))
#         # if len(lisf) > 1:
#         #     raise ValueError('too many models matching: {}'.format(lisf))
#         # elif len(lisf) == 1:
#         #     okfil = lisf[0]
#         # else:
#         #     print('Model {} not found!'.format(mod))
#         #     continue
#         okfil = fil_tas.format(*(mod.split('_')))
#         # leggo tas, faccio detrend, calcolo anom e ok.
#         var, coords, aux_info = ctl.readxDncfield(okfil)
#         lat = coords['lat']
#         lon = coords['lon']
#         dates = coords['dates']
#
#         tas, coeffs, var_regional, dates_seas = ctl.remove_global_polytrend(lat, lon, var, dates, 'NDJFM', deg = 3, area = 'global', print_trend = True)
#         tas_anom[(ssp, mod)] = tas
#
#         trendmat, errtrendmat, cmat, errcmat = ctl.local_lineartrend_climate(lat, lon, var, dates, 'NDJFM', print_trend = True, remove_global_trend = False, global_deg = 3, global_area = 'global')
#         tas_trends[(ssp, mod)] = trendmat
#
# pickle.dump([tas_anom, tas_trends], open(cart_out + 'tas_anom_ssp585.p', 'wb'))
tas_anom, tas_trends = pickle.load(open(cart_out + 'tas_anom_ssp585.p', 'rb'))

area = 'EAT'
ssp = 'ssp585'
cose = dict()

for area in ['EAT', 'PNA']:
    results_hist, results_ref = ctl.load_wrtool(file_hist.format(area))
    results_hist_refEOF, _ = ctl.load_wrtool(file_hist_refEOF.format(area))

    cart_lui = cart_in + 'Results_v5_rebase/{}_NDJFM/'.format(area)
    freqs, residtimes, patterns = pickle.load(open(cart_lui + 'allresults_dicts_{}_v3.p'.format(area), 'rb'))
    trend_ssp, residtime_ssp = pickle.load(open(cart_lui + 'trends_wrfreq_e_restime_{}.p'.format(area), 'rb'))
    seasfreq, runfreq = pickle.load(open(cart_lui + 'seasfreqs_{}_v4.p'.format(area), 'rb'))

    #for ssp in allssps:
    for ssp in ['ssp585']:
        print('SSP '+ssp)
        results_ssp, _ = ctl.load_wrtool(gen_file_ssp.format(ssp, area))

        okmods = [ke[1] for ke in freqs if ssp in ke and 'tot50' in ke and 'all' not in ke and 'rel' not in ke]
        if ssp == 'ssp126':
            okmods = [mod for mod in okmods if mod != 'FGOALS-g3_r1i1p1f1']
        print(okmods)

        for mod in okmods:
            for reg in range(4):
                cose[(ssp, area, mod, 'freq', reg)] = freqs[(ssp, mod, 'tot50')][reg]
                cose[(ssp, area, mod, 'trend', reg)]= trend_ssp[(ssp, mod, 'trend', 'seafreq', reg)]


#### ok.
### ORA ho tas_anom e cose
# provo
# NO
# alur. io voglio per ogni punto vedere come correlano i trend delle frequenze con i trend delle SST
# e sarebbe bello capire chi viene prima ma non so se funziona
# però ecco vedere se i modelli con WR trend maggiore si scaldano anche di più in qualche zona sarebbe carino
