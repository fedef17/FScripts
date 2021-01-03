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
    cart_out = cart_in + 'Climate_mean/'
    ctl.mkdir(cart_out)
elif os.uname()[1] == 'wilma':
    cart_in = 'baugigi'
    cart_out = '/home/federico/work/CMIP6/Climate_mean/'
    ctl.mkdir(cart_out)

print('ciaociao sto partendo')

#cart_data = '/data-hobbes/fabiano/WR_CMIP6/'
cart_data = '/home/federico/work/CMIP6/data_25deg/historical/'
yr10 = 10 # length of running mean

numclus = 4
reg_names_area = dict()
reg_names_area['EAT'] = ['NAO+', 'SBL', 'AR', 'NAO-']
reg_names_area['PNA'] = ['PT', 'PNA+', 'PNA-', 'BR']

okmods = ['BCC-CSM2-MR_r1i1p1f1', 'ACCESS-CM2_r1i1p1f1', 'CanESM5_r1i1p1f1', 'CESM2-WACCM_r1i1p1f1', 'CNRM-CM6-1-HR_r1i1p1f2', 'CNRM-CM6-1_r1i1p1f2', 'CNRM-ESM2-1_r1i1p1f2', 'EC-Earth3_r1i1p1f1', 'FGOALS-g3_r1i1p1f1', 'INM-CM4-8_r1i1p1f1', 'INM-CM5-0_r1i1p1f1', 'IPSL-CM6A-LR_r1i1p1f1', 'MIROC6_r1i1p1f1', 'MPI-ESM1-2-HR_r1i1p1f1', 'MPI-ESM1-2-LR_r1i1p1f1', 'MRI-ESM2-0_r1i1p1f1', 'NorESM2-LM_r1i1p1f1', 'NorESM2-MM_r1i1p1f1', 'UKESM1-0-LL_r1i1p1f2']
okmods_mo = [co.split('_')[0] for co in okmods]

ref_file = '/home/federico/work/CMIP6/ERA/ERA40+Int_daily_1957-2018_zg500_remap25_meters.nc'

fieldnam = 'zg'
climate_mean = dict()
climate_std = dict()
num_members = dict()
climate_mean_dates = dict()
levok = 500

with open(cart_data + 'lista_all_hist.dat', 'r') as fillo:
    filli = [fi.rstrip() for fi in fillo.readlines()]

all_mods = np.unique([fi.split('/')[0] for fi in filli])
all_mems = dict()
for mod in all_mods:
   all_mems[mod] = [fi.split('/')[1] for fi in filli if fi.split('/')[0] == mod]

#all_data = ctl.check_available_cmip6_data(fieldnam, 'day', 'historical')
#all_mods = [co[1] for co in all_data]

for mod in okmods_mo:
    print(mod)
    if mod not in all_mods:
        print('NO data for {}'.format(mod))
        continue

    climmeans = dict()
    climmeans['EAT'] = []
    climmeans['PNA'] = []

    for mem in all_mems[mod]:
        okfil = [fi for fi in filli if mod in fi and mem in fi][0]
        try:
            #var, coords, aux_info = ctl.read_cmip6_data(fieldnam, 'day', 'historical', mod, sel_member = mem, extract_level_hPa = levok, regrid_to_reference_file = ref_file, sel_yr_range = (1964, 2014), select_season_first = True, season = 'ONDJFMA', select_area_first = True, area = 'NML')
            var, coords, aux_info = ctl.readxDncfield(cart_data + okfil)
        except Exception as exp:
            print('Unable to read data for {}, going on with next model..'.format(mod + '_' + mem))
            print(exp)
            continue

        lat = coords['lat']
        lon = coords['lon']
        dates = coords['dates']

        zg_dtr, coeffs, var_regional, dates_seas = ctl.remove_global_polytrend(lat, lon, var, dates, 'NDJFM', deg = 1, area = 'NML', print_trend = True)

        for area in ['EAT', 'PNA']:
            zg_ok, _, _ = ctl.sel_area(lat, lon, zg_dtr, area)
            cmean, dates_cm, _ = ctl.daily_climatology(zg_ok, dates_seas, window = 20)
            climmeans[area].append(cmean)
        #trendmat, errtrendmat, cmat, errcmat = ctl.local_lineartrend_climate(lat, lon, var, dates, 'NDJFM', print_trend = True, remove_global_trend = False, global_deg = 3, global_area = 'global')
        #field_trends[(ssp, mod, 'NDJFM')] = trendmat

    for area in ['EAT', 'PNA']:
        climate_mean[(area, mod)] = np.mean(climmeans[area], axis = 0)
        climate_std[(area, mod)] = np.std(climmeans[area], axis = 0)

    climate_mean_dates[mod] = dates_cm
    num_members[mod] = len(climmeans)

pickle.dump([climate_mean, climate_mean_dates, climate_std, num_members], open(cart_out + 'climate_mean_hist.p', 'wb'))

for area in ['EAT', 'PNA']:
    res = dict()
    for mod in okmods_mo:
        res[mod+'_ensmean'] = dict()
        res[mod+'_ensmean']['climate_mean'] = climate_mean[(area, mod)]
        res[mod+'_ensmean']['dates_climate_mean'] = climate_mean_dates[mod]

    pickle.dump([res, dict()], open(cart_out + 'dict_climate_mean_hist_{}.p'.format(area), 'wb'))
