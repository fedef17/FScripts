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
    cart_out = cart_in + 'Results_uaprec_corrmap/'
    ctl.mkdir(cart_out)
elif os.uname()[1] == 'ff-clevo':
    cart_in = '/home/fedefab/Scrivania/Research/Post-doc/lavori/CMIP6/'
elif os.uname()[1] == 'wilma':
    cart_in = 'baugigi'
    cart_out = '/home/federico/work/CMIP6/Results_ua_corrmap/'
    ctl.mkdir(cart_out)

cart_data = '/data-hobbes/fabiano/WR_CMIP6/'
file_hist_refEOF = cart_data + 'out_NEW_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refEOF_dtr.p'
file_hist = cart_data + 'out_NEW_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refCLUS_dtr_light.p'
gen_file_ssp = cart_data + 'out_NEW_cmip6_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr_histrebase.p'

yr10 = 10 # length of running mean

numclus = 4
reg_names_area = dict()
reg_names_area['EAT'] = ['NAO+', 'SBL', 'AR', 'NAO-']
reg_names_area['PNA'] = ['PT', 'PNA+', 'PNA-', 'BR']

okmods = ['ACCESS-CM2_r1i1p1f1', 'BCC-CSM2-MR_r1i1p1f1', 'CanESM5_r1i1p1f1', 'CESM2-WACCM_r1i1p1f1', 'CNRM-CM6-1-HR_r1i1p1f2', 'CNRM-CM6-1_r1i1p1f2', 'CNRM-ESM2-1_r1i1p1f2', 'EC-Earth3_r1i1p1f1', 'FGOALS-g3_r1i1p1f1', 'INM-CM4-8_r1i1p1f1', 'INM-CM5-0_r1i1p1f1', 'IPSL-CM6A-LR_r1i1p1f1', 'MIROC6_r1i1p1f1', 'MPI-ESM1-2-HR_r1i1p1f1', 'MPI-ESM1-2-LR_r1i1p1f1', 'MRI-ESM2-0_r1i1p1f1', 'NorESM2-LM_r1i1p1f1', 'NorESM2-MM_r1i1p1f1', 'UKESM1-0-LL_r1i1p1f2']
okmods_mo = [co.split('_')[0] for co in okmods]

allssps = 'ssp126 ssp245 ssp370 ssp585'.split()
allsimcol = ['hist', 'ssp126', 'ssp245', 'bau', 'ssp370', 'ssp585', 'rcp85_cmip5']
coldic = dict(zip(allsimcol, ctl.color_set(7)))
colssp = [coldic[ssp] for ssp in allssps]

# ref_file = '/home/federico/work/CMIP6/ERA/ERA40+Int_daily_1957-2018_zg500_remap25_meters.nc'
#
# for fieldnam in ['ua', 'pr']:
#     print(fieldnam)
#
#     field_anom = dict()
#     field_trends = dict()
#     if fieldnam == 'ua':
#         levok = 300
#     else:
#         levok = None
#
#     all_data = ctl.check_available_cmip6_data(fieldnam, 'Amon', 'ssp585')
#     all_modmem = [co[1]+'_'+co[2] for co in all_data]
#
#     for ssp in ['ssp585']:
#         print('SSP '+ssp)
#
#         for modmem in okmods:
#             print(modmem)
#             if modmem not in all_modmem:
#                 print('NO data for {}'.format(modmem))
#                 continue
#
#             mod, mem = modmem.split('_')
#             try:
#                 var, coords, aux_info = ctl.read_cmip6_data(fieldnam, 'Amon', 'ssp585', mod, sel_member = mem, extract_level_hPa = levok, regrid_to_reference_file = ref_file, sel_yr_range = (2015, 2100))
#             except Exception as exp:
#                 print('Unable to read data for {}, going on with next model..'.format(modmem))
#                 print(exp)
#                 continue
#
#             lat = coords['lat']
#             lon = coords['lon']
#             dates = coords['dates']
#
#             tas, coeffs, var_regional, dates_seas = ctl.remove_global_polytrend(lat, lon, var, dates, 'NDJFM', deg = 3, area = 'global', print_trend = True)
#             field_anom[(ssp, mod, 'NDJFM')] = tas
#
#             trendmat, errtrendmat, cmat, errcmat = ctl.local_lineartrend_climate(lat, lon, var, dates, 'NDJFM', print_trend = True, remove_global_trend = False, global_deg = 3, global_area = 'global')
#             field_trends[(ssp, mod, 'NDJFM')] = trendmat
#
#             tas, coeffs, var_regional, dates_seas = ctl.remove_global_polytrend(lat, lon, var, dates, None, deg = 3, area = 'global', print_trend = True)
#             field_anom[(ssp, mod, 'year')] = tas
#
#             trendmat, errtrendmat, cmat, errcmat = ctl.local_lineartrend_climate(lat, lon, var, dates, None, print_trend = True, remove_global_trend = False, global_deg = 3, global_area = 'global')
#             field_trends[(ssp, mod, 'year')] = trendmat
#
#     pickle.dump([field_anom, field_trends], open(cart_out + '{}_anom_ssp585.p'.format(fieldnam), 'wb'))



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


tas_anom, tas_trends = pickle.load(open(cart_out + '../Results_SST_corrmap/tas_anom_ssp585.p', 'rb'))
#### ok.
### ORA ho tas_anom e cose
for fieldnam in ['ua', 'pr']:
    field_anom, field_trends = pickle.load(open(cart_out + '{}_anom_ssp585.p'.format(fieldnam), 'rb'))
    print(field_trends.keys())
    okmok = []
    for mod in okmods:
        if ('ssp585', mod.split('_')[0], 'year') in field_trends:
            okmok.append(mod)

    print('Available models:', okmok)

    # provo
    # per ogni punto devo fare corr tra freq trend e sst trend
    # però.. entrambi detrendati per il global warming?
    # sennò non ha molto senso, cioè beccherò le zone che si scaldano di più

    # devo dividere sst_trend e freq_trend per il global tas trend di ogni modello
    corrmaps = dict()
    ssp = 'ssp585'
    for seas in ['NDJFM', 'year']:
        for area in ['EAT', 'PNA']:
            for reg in range(4):
                # trendmat = tas_trends[(ssp, okmods[0])]
                trendmat = field_trends[(ssp, okmok[0].split('_')[0], seas)]
                corr_map = np.empty_like(trendmat)
                pval_map = np.empty_like(trendmat)
                nlat, nlon = trendmat.shape
                lat, lon = ctl.genlatlon(nlat, nlon)

                gw = np.array([ctl.global_mean(tas_trends[(ssp, mod, seas)], lat) for mod in okmok])
                #gw = np.array([ctl.global_mean(tas_trends[(ssp, mod.split('_')[0])], lat) for mod in okmok])
                frok = np.array([cose[(ssp, area, mod, 'trend', reg)] for mod in okmok])

                for la in range(nlat):
                    for lo in range(nlon):
                        tastr = np.array([field_trends[(ssp, mod.split('_')[0], seas)][la, lo] for mod in okmok])
                        pears, pval = stats.pearsonr(frok/gw, tastr/gw)

                        corr_map[la, lo] = pears
                        pval_map[la, lo] = pval

                corrmaps[('corr', area, reg)] = corr_map
                corrmaps[('pval', area, reg)] = pval_map

                fnam = cart_out + '{}_corrmap_{}_{}_{}.pdf'.format(fieldnam, area, reg, seas)
                ctl.plot_map_contour(corr_map, lat, lon, filename = fnam, add_hatching = pval_map <= 0.05, cbar_range = (-1, 1), cb_label = 'Correlation', title = 'area: {}, regime: {}, seas: {}'.format(area, reg_names_area[area][reg], seas), draw_grid = True)

                fnam = cart_out + '{}_pvalmap_{}_{}_{}.pdf'.format(fieldnam, area, reg, seas)
                ctl.plot_map_contour(pval_map, lat, lon, filename = fnam, cbar_range = (0., 0.1), plot_anomalies = False, extend_opt = 'neither', draw_grid = True, cb_label = 'P-value', title = 'area: {}, regime: {}, seas: {}'.format(area, reg_names_area[area][reg], seas))

            cmape = [corrmaps[('corr', area, reg)] for reg in range(4)]
            pvlep = [corrmaps[('pval', area, reg)] <= 0.05 for reg in range(4)]

            fnam = cart_out + '{}_corrmap_{}_allregs_{}.pdf'.format(fieldnam, area, seas)
            ctl.plot_multimap_contour(cmape, lat, lon, filename = fnam, add_hatching = pvlep, cbar_range = (-1, 1), cb_label = 'Correlation', title = 'area: {}, seas: {}'.format(area, seas), subtitles = reg_names_area[area], draw_grid = True, figsize = (18,12))
