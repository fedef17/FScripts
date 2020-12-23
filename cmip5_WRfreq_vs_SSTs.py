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

cart_tas = '/data-hobbes/fabiano/cmip5_tas/'
fil_tas = cart_tas + 'tas_Amon_{}_rcp85_{}_200601-210012.nc'

cart_cmip5 = '/home/fabiano/Research/lavori/CMIP6/Results_cmip5/{}_NDJFM/'
cart_out = cart_in + 'Results_SST_corrmap/'
ctl.mkdir(cart_out)

yr10 = 10 # length of running mean

numclus = 4
reg_names_area = dict()
reg_names_area['EAT'] = ['NAO+', 'SBL', 'AR', 'NAO-']
reg_names_area['PNA'] = ['PT', 'PNA+', 'PNA-', 'BR']
#okmods = ['ACCESS-CM2_r1i1p1f1', 'BCC-CSM2-MR_r1i1p1f1', 'CanESM5_r1i1p1f1', 'CESM2-WACCM_r1i1p1f1','CNRM-CM6-1-HR_r1i1p1f2', 'CNRM-CM6-1_r1i1p1f2','CNRM-ESM2-1_r1i1p1f2', 'EC-Earth3_r1i1p1f1', 'FGOALS-g3_r1i1p1f1','INM-CM4-8_r1i1p1f1', 'INM-CM5-0_r1i1p1f1', 'IPSL-CM6A-LR_r1i1p1f1','MIROC6_r1i1p1f1', 'MPI-ESM1-2-HR_r1i1p1f1','MPI-ESM1-2-LR_r1i1p1f1', 'MRI-ESM2-0_r1i1p1f1','NorESM2-LM_r1i1p1f1', 'NorESM2-MM_r1i1p1f1', 'UKESM1-0-LL_r1i1p1f2']
# mancano = ['BCC-ESM1_r1i1p1f1', 'CESM2_r1i1p1f1', 'GFDL-CM4_r1i1p1f1', 'HadGEM3-GC31-LL_r1i1p1f3', 'KACE-1-0-G_r1i1p1f1', 'MPI-ESM-1-2-HAM_r1i1p1f1']

allssps = 'ssp126 ssp245 ssp370 ssp585'.split()
allsimcol = ['hist', 'ssp126', 'ssp245', 'bau', 'ssp370', 'ssp585', 'rcp85_cmip5']
coldic = dict(zip(allsimcol, ctl.color_set(7)))
colssp = [coldic[ssp] for ssp in allssps]

tas_anom = dict()
tas_trends = dict()
cose = dict()
for ssp in ['rcp85_cmip5']:
    print('SSP '+ssp)

    area = 'EAT'
    cart_lui = cart_cmip5.format(area)

    freqs, trend_ssp, residtimes = pickle.load(open(cart_cmip5.format(area) + 'freqs_cmip5_{}.p'.format(area), 'rb'))
    seasfreq, runfreq = pickle.load(open(cart_cmip5.format(area) + 'seasfreqs_cmip5_{}.p'.format(area), 'rb'))

    okmods = [ke[1] for ke in freqs if ssp in ke and 'tot50' in ke and 'all' not in ke and 'rel' not in ke]
    print(okmods)

    for mod in okmods:
        for reg in range(4):
            cose[(ssp, area, mod, 'freq', reg)] = freqs[(ssp, mod, 'tot50')][reg]
            cose[(ssp, area, mod, 'trend', reg)]= trend_ssp[(ssp, mod, 'trend', 'seafreq', reg)]

    for mod in okmods:
        okfil = fil_tas.format(*(mod.split('_')))
        # leggo tas, faccio detrend, calcolo anom e ok.
        var, coords, aux_info = ctl.readxDncfield(okfil)
        lat = coords['lat']
        lon = coords['lon']
        dates = coords['dates']

        tas, coeffs, var_regional, dates_seas = ctl.remove_global_polytrend(lat, lon, var, dates, 'NDJFM', deg = 3, area = 'global', print_trend = True)
        tas_anom[(ssp, mod, 'NDJFM')] = tas

        trendmat, errtrendmat, cmat, errcmat = ctl.local_lineartrend_climate(lat, lon, var, dates, 'NDJFM', print_trend = True, remove_global_trend = False, global_deg = 3, global_area = 'global')
        tas_trends[(ssp, mod, 'NDJFM')] = trendmat

        tas, coeffs, var_regional, dates_seas = ctl.remove_global_polytrend(lat, lon, var, dates, None, deg = 3, area = 'global', print_trend = True)
        tas_anom[(ssp, mod, 'year')] = tas

        trendmat, errtrendmat, cmat, errcmat = ctl.local_lineartrend_climate(lat, lon, var, dates, None, print_trend = True, remove_global_trend = False, global_deg = 3, global_area = 'global')
        tas_trends[(ssp, mod, 'year')] = trendmat

pickle.dump([tas_anom, tas_trends], open(cart_out + 'tas_anom_rcp85.p', 'wb'))
pickle.dump([okmods, cose], open(cart_out + 'cose_rcp85.p', 'wb'))
