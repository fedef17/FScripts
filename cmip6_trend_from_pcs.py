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
    cart = '/home/fabiano/Research/lavori/CMIP6/'
elif os.uname()[1] == 'ff-clevo':
    cart = '/home/fedefab/Scrivania/Research/Post-doc/lavori/CMIP6/'

cart_out_orig = cart + 'Results_trendssp585/'
ctl.mkdir(cart_out_orig)

cart_in = '/data-hobbes/fabiano/WR_CMIP6/'

file_hist_refEOF = cart_in + 'out_NEW_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refEOF_dtr.p'
file_hist = cart_in + 'out_NEW_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refCLUS_dtr_light.p'
gen_file_ssp = cart_in + 'out_NEW_cmip6_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr_histrebase.p'

numclus = 4
reg_names_area = dict()
reg_names_area['EAT'] = ['NAO+', 'SBL', 'AR', 'NAO-']
reg_names_area['PNA'] = ['PT', 'PNA+', 'PNA-', 'AR']

area = 'EAT'
ssp = 'ssp585'

results_hist, results_ref = ctl.load_wrtool(file_hist.format(area))
results_ssp, _ = ctl.load_wrtool(gen_file_ssp.format(ssp, area))
print(results_hist.keys())
print(len(results_hist.keys()))
print(results_ssp.keys())
print(len(results_ssp.keys()))

cart_lui = cart + 'Results_v5_rebase/{}_NDJFM/'.format(area)
freqs, residtimes, patterns = pickle.load(open(cart_lui + 'allresults_dicts_{}_v3.p'.format(area), 'rb'))

okmods = [ke[1] for ke in freqs if 'ssp585' in ke and 'tot50' in ke and 'all' not in ke and 'rel' not in ke]
print(okmods)
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
    #trendssp, errtrendssp = ctl.local_lineartrend_climate(lat, lon, zgssp, dates, season)
    zg_noglob, coeffs, var_regional, dates_seas = ctl.remove_global_polytrend(lat, lon, zgssp, dates, season, deg = 3, area = 'NML')

    trend, errtrend, _, _ = ctl.local_lineartrend_climate(lat, lon, zg_noglob, dates, None)
    zontrend = ctl.zonal_mean(trend)
    zonerrtrend = ctl.zonal_mean(errtrend)

    zg_noglob_zonme = ctl.zonal_mean(zg_noglob)
    zg_se = zg_noglob - zg_noglob_zonme[..., np.newaxis]
    se_trend, se_errtrend, _, _ = ctl.local_lineartrend_climate(lat, lon, zg_se, dates, None)

    try:
        zghist, coords, _ = ctl.readxDncfield(cartmon_hist + filhist.format(mod))
        lat = coords['lat']
        lon = coords['lon']
        dates = coords['dates']
        zgmean, zgstd = ctl.seasonal_climatology(zghist, dates, season, dates_range = ctl.range_years(1964,2014))
        cose[('hist mean', mod)] = zgmean
        cose[('hist std', mod)] = zgstd
    except:
        print('passsssss ' + mod)
        pass

    cose[('trend', mod)] = trend
    cose[('errtrend', mod)] = errtrend
    cose[('se_trend', mod)] = se_trend
    cose[('se_errtrend', mod)] = se_errtrend
    cose[('zontrend', mod)] = zontrend
    cose[('zonerrtrend', mod)] = zonerrtrend

pickle.dump(cose, open(cart_out_orig + 'zgtrends_ssp585.p', 'wb'))


area = (-180, 180, 20, 90)

trendsanom = []
trendsstateddy = []
zontrend = []
zontrend_EAT = []
zontrend_PNA = []

hatchs = []
hatchs_se = []
for modmem in okmods:
    mod, mem = modmem.split('_')
    trend = cose[('trend', mod)].squeeze()
    errtrend = cose[('errtrend', mod)].squeeze()

    # zonme = ctl.zonal_mean(trend)
    # stat_eddy_trend = np.empty_like(trend)
    # for i in range(trend.shape[0]):
    #     stat_eddy_trend[i,:] = trend[i,:]-zonme[i]

    stat_eddy_trend = cose[('se_trend', mod)].squeeze()
    se_errtrend = cose[('se_errtrend', mod)].squeeze()

    # trend_area, lat_area, lon_area = ctl.sel_area(lat, lon, trend, area)
    # trend_anom = trend-np.mean(trend_area)
    trendsanom.append(trend)
    trendsstateddy.append(stat_eddy_trend)
    zontrend.append(cose[('zontrend', mod)].squeeze())

    lata, lona, var = ctl.sel_area(lat, lon, trend, 'EAT')
    zeat = ctl.zonal_mean(var)
    zontrend_EAT.append(zeat)

    lata, lona, var = ctl.sel_area(lat, lon, trend, 'PNA')
    zpna = ctl.zonal_mean(var)
    zontrend_PNA.append(zpna)

    hatchs.append(np.abs(trend) > 2*errtrend)
    hatchs_se.append(np.abs(stat_eddy_trend) > 2*se_errtrend)

NSIG = int(0.8*len(okmods)) # number of models to consider response significant

trendsanom.append(np.mean(trendsanom, axis = 0))
hatchs.append(np.sum([np.sign(tre) == np.sign(trendsanom[-1]) for tre in trendsanom[:-1]], axis = 0) >= NSIG)
trendsstateddy.append(np.mean(trendsstateddy, axis = 0))
hatchs_se.append(np.sum([np.sign(tre) == np.sign(trendsstateddy[-1]) for tre in trendsstateddy[:-1]], axis = 0) >= NSIG)

allmods = [modmem.split('_')[0] for modmem in okmods]
#meanfields = [cose[('hist mean', mod)] for mod in allmods]

meanfields = []
stateddies = []
for mod in allmods:
    mf = cose[('hist mean', mod)].squeeze()
    zonme = ctl.zonal_mean(mf)
    stat_eddy = mf - zonme[:, np.newaxis]

    meanfields.append(mf)
    stateddies.append(stat_eddy)

stateddies.append(np.mean(stateddies, axis = 0))
allmods_MM = allmods + ['Multi-model mean']

filename = cart_out_orig + 'zontrend_ssp585.pdf'
oklats = lat >= 30
fig = plt.figure(figsize = (16,12))
ax = fig.add_subplot(111)
for modmem, co in zip(okmods, zontrend):
    ax.plot(co[oklats], lat[oklats], label = modmem)
ax.plot(np.mean(zontrend, axis = 0)[oklats], lat[oklats], label = 'MMM', color = 'black', linewidth = 3)
ax.axvline(0., color = 'lightslategray', linewidth = 0.2)
plt.legend()
fig.savefig(filename)

cols = ctl.color_set(3)

filename = cart_out_orig + 'zontrend_ssp585_MMM.pdf'
fig = plt.figure(figsize = (16,12))
ax = fig.add_subplot(111)
coso = np.mean(zontrend, axis = 0)[oklats]
coserr = np.std(zontrend, axis = 0)[oklats]

ax.fill_betweenx(lat[oklats], coso-coserr, coso+coserr, color = cols[0], alpha = 0.2)
ax.plot(coso, lat[oklats], color = cols[0], linewidth = 2, label = 'NML')

coso_ = np.mean(zontrend_PNA, axis = 0)
coserr = np.std(zontrend_PNA, axis = 0)
ax.fill_betweenx(lata, coso-coserr, coso+coserr, color = cols[1], alpha = 0.2)
ax.plot(coso, lata, color = cols[1], linewidth = 2, label = 'EAT')

coso = np.mean(zontrend_EAT, axis = 0)
coserr = np.std(zontrend_EAT, axis = 0)
ax.fill_betweenx(lata, coso-coserr, coso+coserr, color = cols[2], alpha = 0.2)
ax.plot(coso, lata, color = cols[2], linewidth = 2, label = 'PNA')

ax.axvline(0., color = 'lightslategray', linewidth = 0.2)
ax.set_xlabel('Zonal trend anomaly (m/yr)')
ax.set_ylabel('Latitude')
plt.legend()
fig.savefig(filename)

filename = cart_out_orig + 'trend_anom_ssp585.pdf'
ctl.plot_multimap_contour(trendsanom, lat, lon, filename, plot_anomalies=True, plot_margins=(-180, 180, 20, 90), cbar_range=(-2,2), add_hatching = hatchs, fix_subplots_shape = (7, 2), figsize = (15,20), subtitles = allmods_MM, cb_label = 'm/year', verbose = True, draw_grid = True)

filename = cart_out_orig + 'trend_anom_ssp585_EAT.pdf'
ctl.plot_multimap_contour(trendsanom, lat, lon, filename, plot_anomalies=True, visualization = 'nearside', central_lat_lon = (65, -30), cbar_range=(-2,2), add_hatching = hatchs, fix_subplots_shape = (3, 5), figsize = (18,12), subtitles = allmods_MM, cb_label = 'm/year', verbose = True, draw_grid = True)

filename = cart_out_orig + 'trend_anom_mmm_vs_hist.pdf'
ctl.plot_map_contour(trendsanom[-1], lat, lon, filename, plot_anomalies=True, plot_margins=(-180, 180, 20, 90), cbar_range=(-1,1), add_contour_field = meanfields[-1], figsize = (24,12), cb_label = 'm/year', draw_grid = True, add_hatching = hatchs[-1], n_lines = 8, add_contour_same_levels = False, add_contour_plot_anomalies = False)

filename = cart_out_orig + 'trend_anom_mmm_vs_hist_EAT.pdf'
ctl.plot_map_contour(trendsanom[-1], lat, lon, filename, plot_anomalies=True, visualization = 'nearside', central_lat_lon = (65, -30), cbar_range=(-1,1), add_contour_field = meanfields[-1], cb_label = 'm/year', draw_grid = True, add_hatching = hatchs[-1], n_lines = 8, add_contour_same_levels = False, add_contour_plot_anomalies = False)

filename = cart_out_orig + 'trend_anom_mmm_vs_hist_polar.pdf'
ctl.plot_map_contour(trendsanom[-1], lat, lon, filename, plot_anomalies=True, visualization = 'nearside', bounding_lat = 10, central_lat_lon = (90, 0), cbar_range=(-1,1), add_contour_field = meanfields[-1], cb_label = 'm/year', draw_grid = True, add_hatching = hatchs[-1], n_lines = 8, add_contour_same_levels = False, add_contour_plot_anomalies = False)

filename = cart_out_orig + 'trend_stateddy_ssp585.pdf'
ctl.plot_multimap_contour(trendsstateddy, lat, lon, filename, plot_anomalies=True, plot_margins=(-180, 180, 20, 90), cbar_range=(-2,2), add_hatching = hatchs_se, fix_subplots_shape = (7, 2), figsize = (15,20), subtitles = allmods_MM, cb_label = 'm/year', draw_grid = True)

filename = cart_out_orig + 'trend_stateddy_ssp585_EAT.pdf'
ctl.plot_multimap_contour(trendsstateddy, lat, lon, filename, plot_anomalies=True, visualization = 'nearside', central_lat_lon = (65, -30), cbar_range=(-2,2), add_hatching = hatchs_se, fix_subplots_shape = (3, 5), figsize = (18,12), subtitles = allmods_MM, cb_label = 'm/year', verbose = True, draw_grid = True)

filename = cart_out_orig + 'trend_stateddy_ssp585_polar.pdf'
ctl.plot_multimap_contour(trendsstateddy, lat, lon, filename, plot_anomalies=True, visualization = 'nearside', bounding_lat = 10, central_lat_lon = (90, 0), cbar_range=(-2,2), add_hatching = hatchs_se, fix_subplots_shape = (3, 5), figsize = (18,12), subtitles = allmods_MM, cb_label = 'm/year', verbose = True, draw_grid = True)

filename = cart_out_orig + 'stateddies_hist.pdf'
ctl.plot_multimap_contour(stateddies, lat, lon, filename, plot_anomalies=True, plot_margins=(-180, 180, 20, 90), cbar_range=(-200,200), fix_subplots_shape = (7, 2), figsize = (15,20), subtitles = allmods_MM, cb_label = 'm', draw_grid = True)

filename = cart_out_orig + 'trend_stateddy_ssp585_whist.pdf'
ctl.plot_multimap_contour(trendsstateddy, lat, lon, filename, plot_anomalies=True, plot_margins=(-180, 180, 20, 90), cbar_range=(-2,2), add_contour_field = stateddies, fix_subplots_shape = (7, 2), figsize = (15,20), subtitles = allmods_MM, cb_label = 'm/year', draw_grid = True, n_lines = 8, add_contour_same_levels = False, add_contour_plot_anomalies = True)

filename = cart_out_orig + 'stateddies_trend_mmm_vs_hist.pdf'
ctl.plot_map_contour(trendsstateddy[-1], lat, lon, filename, plot_anomalies=True, plot_margins=(-180, 180, 20, 90), cbar_range=(-1,1), add_contour_field = stateddies[-1], figsize = (24,12), cb_label = 'm/year', draw_grid = True, add_hatching = hatchs_se[-1], n_lines = 8, add_contour_same_levels = False, add_contour_plot_anomalies = True)

filename = cart_out_orig + 'stateddies_trend_mmm_vs_hist_EAT.pdf'
ctl.plot_map_contour(trendsstateddy[-1], lat, lon, filename, plot_anomalies=True, visualization = 'nearside', central_lat_lon = (65, -30), cbar_range=(-1,1), add_contour_field = stateddies[-1], cb_label = 'm/year', draw_grid = True, add_hatching = hatchs_se[-1], n_lines = 8, add_contour_same_levels = False, add_contour_plot_anomalies = True)

filename = cart_out_orig + 'stateddies_trend_mmm_vs_hist_polar.pdf'
ctl.plot_map_contour(trendsstateddy[-1], lat, lon, filename, plot_anomalies=True, visualization = 'nearside', bounding_lat = 10, central_lat_lon = (90, 0), cbar_range=(-1,1), add_contour_field = stateddies[-1], cb_label = 'm/year', draw_grid = True, add_hatching = hatchs_se[-1], n_lines = 8, add_contour_same_levels = False, add_contour_plot_anomalies = True)
