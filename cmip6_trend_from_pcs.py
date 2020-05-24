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
if os.uname()[1] == 'hobbes':
    cart_in = '/home/fabiano/Research/lavori/CMIP6/'
elif os.uname()[1] == 'ff-clevo':
    cart_in = '/home/fedefab/Scrivania/Research/Post-doc/lavori/CMIP6/'

cart_out_orig = cart_in + 'Results_trendssp585/'
ctl.mkdir(cart_out_orig)

file_hist = cart_in + 'cmip6_hist/out_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refCLUS_dtr.p'
#file_hist_refEOF = cart_in + 'cmip6_hist/out_cmip6_hist_NDJFM_EAT_4clus_4pcs_1964-2014_refEOF.p'
gen_file_ssp = cart_in + 'cmip6_{}/out_cmip6_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr.p'

numclus = 4
reg_names_area = dict()
reg_names_area['EAT'] = ['NAO+', 'SBL', 'AR', 'NAO-']
reg_names_area['PNA'] = ['PT', 'PNA+', 'PNA-', 'AR']


area = 'EAT'
ssp = 'ssp585'

results_hist, results_ref = pickle.load(open(file_hist.format(area), 'rb'))
results_ssp = pickle.load(open(gen_file_ssp.format(ssp, ssp, area), 'rb'))['models']

cart_lui = cart_in + 'Results_v4/{}_NDJFM/'.format(area)
freqs, residtimes = pickle.load(open(cart_lui + 'allresults_dicts_{}_v3.p'.format(area), 'rb'))

okmods = [ke[1] for ke in freqs if 'ssp585' in ke and 'tot50' in ke and 'all' not in ke and 'rel' not in ke]
#['BCC-CSM2-MR_r1i1p1f1', 'CanESM5_r1i1p1f1', 'CESM2-WACCM_r1i1p1f1\', 'CNRM-CM6-1_r1i1p1f2', 'CNRM-ESM2-1_r1i1p1f2', 'EC-Earth3_r1i1p1f1', 'FGOALS-g3_r1i1p1f1', 'INM-CM4-8_r1i1p1f1', 'INM-CM5-0_r1i1p1f1', 'IPSL-CM6A-LR_r1i1p1f1', 'MIROC6_r1i1p1f1', 'MPI-ESM1-2-HR_r1i1p1f1', 'MRI-ESM2-0_r1i1p1f1', 'UKESM1-0-LL_r1i1p1f2']

for nu in range(4):
    for mod in okmods:
        var, dat = ctl.seasonal_set(results_ssp['EC-Earth3_r1i1p1f1']['pcs'][:, 0], results_ssp['EC-Earth3_r1i1p1f1']['dates'], 'DJF', seasonal_average = True)
        years = np.array([da.year for da in dat])
        m, c, err_m, err_c = ctl.linear_regre_witherr(years, var)

        daynum = np.arange(len(results_ssp['EC-Earth3_r1i1p1f1']['dates']))
        m, c, err_m, err_c = ctl.linear_regre_witherr(daynum, results_ssp['EC-Earth3_r1i1p1f1']['pcs'])

cartmon_hist = '/data-hobbes/fabiano/CMIP6/historical_mon_zg/'
cartmon_ssp = '/data-hobbes/fabiano/CMIP6/ssp585_mon_zg/'

filssp = 'zg_Amon_ssp585_{}_{}_2015-2100.nc'
filhist = 'zg_{}_mon.nc'

season = 'NDJFM'

cose = dict()
for modmem in okmods:
    mod, mem = modmem.split('_')
    zgssp, coords, _ = ctl.readxDncfield(cartmon_ssp + filssp.format(mod, mem))
    lat = coords['lat']
    lon = coords['lon']
    dates = coords['dates']
    trendssp, errtrendssp = ctl.local_lineartrend_climate(lat, lon, zgssp, dates, season)

    zghist, coords, _ = ctl.readxDncfield(cartmon_hist + filhist.format(mod))
    lat = coords['lat']
    lon = coords['lon']
    dates = coords['dates']
    zgmean, zgstd = ctl.seasonal_climatology(zghist, dates, season, dates_range = ctl.range_years(1964,2014))

    cose[('trend', mod)] = trendssp
    cose[('errtrend', mod)] = errtrendssp
    cose[('hist mean', mod)] = zgmean
    cose[('hist std', mod)] = zgstd

pickle.dump(cose, open(cart_out_orig + 'zgtrends_ssp585.p', 'wb'))
