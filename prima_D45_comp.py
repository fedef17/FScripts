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
#import pymannkendall as mk

#############################################################################
if os.uname()[1] == 'hobbes':
    cart_in = '/home/fabiano/Research/lavori/WeatherRegimes/'
    cart_out = '/home/fabiano/Research/lavori/prima_D45/'
elif os.uname()[1] == 'ff-clevo':
    cart_out = '/home/fedefab/Scrivania/Research/Post-doc/lavori/prima_D45/'

fil_pres = 'prima_D45_pres/out_prima_D45_pres_DJF_EAT_4clus_4pcs_1979-2014_refEOF_dtr.p'
fil_fut = 'prima_D45_fut/out_prima_D45_fut_DJF_EAT_4clus_4pcs_2015-2050_refEOF_dtr.p'

results_pres, results_ref = ctl.load_wrtool(cart_in + fil_pres)
results_fut, _ = ctl.load_wrtool(cart_in + fil_fut)

cart_data = '/nas/PRIMAVERA/Stream1/'
filtas = 'highresSST-{}/{}/{}/day/{}/{}_day_{}_highresSST-{}_{}_*_r25_rc.nc'
#filtas.format(temp, mod, mem, varnam, varnam, mod, temp, mem)

# composites = dict()
# mod = 'ref'
# temp = 'present'
# file_ref = dict()
# file_ref['tas'] = '/data-hobbes/fabiano/OBS/ERA/ERAInterim/ERAInt_daily_1979-2018_167_r25.nc'
# file_ref['pr'] = '/data-hobbes/fabiano/OBS/ERA/ERAInterim/ERAInt_daily_1979-2018_228_pr_daysum_ok_r25.nc'
# for varnam in ['tas', 'pr']:
#     var, coords, aux_info = ctl.read_iris_nc(file_ref[varnam])
#     var_season, dates_season = ctl.sel_season(var, coords['dates'], 'DJF')
#
#     for comp_moment in ['mean', 'std', 'percentile']:
#         comps = ctl.composites_regimes_daily(coords['lat'], coords['lon'], var_season, dates_season, results_ref['labels'], results_ref['dates'], comp_moment = comp_moment, detrend_global = True, area_dtr = 'NML')
#         composites[(temp, mod, varnam, comp_moment)] = comps


# for ke in results_pres:
#     temp = 'present'
#     mod, mem = ke.split('_')
#     for varnam in ['tas', 'pr']:
#         listafilo = ctl.find_file(cart_data, filtas.format(temp, mod, mem, varnam, varnam, mod, temp, mem))
#         if len(listafilo) > 0:
#             filo = listafilo[0]
#         else:
#             print('NOT FOUND: ', temp, mod, varnam)
#             continue
#
#         var, coords, aux_info = ctl.read_iris_nc(filo)
#         var_season, dates_season = ctl.sel_season(var, coords['dates'], 'DJF')
#
#         for comp_moment in ['mean', 'std', 'percentile']:
#             comps = ctl.composites_regimes_daily(coords['lat'], coords['lon'], var_season, dates_season, results_pres[ke]['labels'], results_pres[ke]['dates'], comp_moment = comp_moment, detrend_global = True, area_dtr = 'NML')
#             composites[(temp, mod, varnam, comp_moment)] = comps
#
#
# for ke in results_fut:
#     temp = 'future'
#     mod, mem = ke.split('_')
#     for varnam in ['tas', 'pr']:
#         listafilo = ctl.find_file(cart_data, filtas.format(temp, mod, mem, varnam, varnam, mod, temp, mem))
#         if len(listafilo) > 0:
#             filo = listafilo[0]
#         else:
#             print('NOT FOUND: ', temp, mod, varnam)
#             continue
#
#         var, coords, aux_info = ctl.read_iris_nc(filo)
#         var_season, dates_season = ctl.sel_season(var, coords['dates'], 'DJF')
#
#         for comp_moment in ['mean', 'std', 'percentile']:
#             comps = ctl.composites_regimes_daily(coords['lat'], coords['lon'], var_season, dates_season, results_fut[ke]['labels'], results_fut[ke]['dates'], comp_moment = comp_moment, detrend_global = True, area_dtr = 'NML')
#             composites[(temp, mod, varnam, comp_moment)] = comps
#
# with open(cart_out + 'composites_taspr.p', 'wb') as filonz:
#     pickle.dump(composites, filonz)
# with open(cart_out + 'composites_taspr.p', 'rb') as filonz:
#     composites = pickle.load(filonz)
#
# # Riporto tutto alle anomalie
# for ke in composites:
#     composites[ke] = composites[ke]-np.mean(composites[ke], axis = 0)
# with open(cart_out + 'composites_taspr_anom.p', 'wb') as filonz:
#     pickle.dump(composites, filonz)
with open(cart_out + 'composites_taspr_anom.p', 'rb') as filonz:
    composites = pickle.load(filonz)

# plots
okmods = ['CMCC-CM2-HR4', 'CMCC-CM2-VHR4', 'EC-Earth3P', 'EC-Earth3P-HR', 'HadGEM3-GC31-LM', 'HadGEM3-GC31-MM', 'MPI-ESM1-2-HR', 'MPI-ESM1-2-XR']
okmods_L = ['CMCC-CM2-HR4', 'EC-Earth3P', 'HadGEM3-GC31-LM', 'MPI-ESM1-2-HR']
okmods_H = ['CMCC-CM2-VHR4', 'EC-Earth3P-HR', 'HadGEM3-GC31-MM', 'MPI-ESM1-2-XR']

okmods_hist_L = okmods_L[1:] + ['ECMWF-IFS-LR']
okmods_hist_H = okmods_H[1:] + ['ECMWF-IFS-HR']
# Freq difference
freqs = dict()

cose_H = []
cose_L = []
for ke in results_pres:
    for mod in okmods_H:
        if mod == ke.split('_')[0]:
            cose_H.append(results_pres[ke]['freq_clus'])
    for mod in okmods_L:
        if mod == ke.split('_')[0]:
            cose_L.append(results_pres[ke]['freq_clus'])

freqs[('pres', 'HR')] = np.mean(cose_H, axis = 0)
freqs[('pres', 'LR')] = np.mean(cose_L, axis = 0)
freqs[('pres', 'HRstd')] = np.std(cose_H, axis = 0)
freqs[('pres', 'LRstd')] = np.std(cose_L, axis = 0)

cose_H = []
cose_L = []
for ke in results_fut:
    for mod in okmods_H:
        if mod == ke.split('_')[0]:
            cose_H.append(results_fut[ke]['freq_clus'])
    for mod in okmods_L:
        if mod == ke.split('_')[0]:
            cose_L.append(results_fut[ke]['freq_clus'])

freqs[('fut', 'HR')] = np.mean(cose_H, axis = 0)
freqs[('fut', 'LR')] = np.mean(cose_L, axis = 0)
freqs[('fut', 'HRstd')] = np.std(cose_H, axis = 0)
freqs[('fut', 'LRstd')] = np.std(cose_L, axis = 0)

keall = [('pres', 'LR'), ('pres', 'HR'), ('fut', 'LR'), ('fut', 'HR')]
laball = ['_'.join(ke) for ke in keall]
colors = ctl.color_set(4)
regnames = ['NAO+', 'SBL', 'NAO-', 'AR']

freq_ref = results_ref['freq_clus']
fig = plt.figure(figsize=(16,12))
axes = []
for reg in range(4):
    ax = plt.subplot(2,2,reg+1)
    frqall = np.array([freqs[cos][reg] for cos in keall])
    ax.bar(np.arange(4), frqall-freq_ref[reg], color = colors)

    ax.set_title(regnames[reg])
    ax.set_xticks([])
    axes.append(ax)

ctl.adjust_ax_scale(axes)
ctl.custom_legend(fig, colors, laball, ncol = 2)
fig.savefig(cart_out+'Regime_frequency.pdf')


# Hist composite difference (HR, LR)
compdiffs = dict()
for varnam in ['tas', 'pr']:
    comps = []
    for mod in okmods_hist_L:
        comps.append(composites[('present', mod, varnam, 'mean')]-composites[('present', 'ref', varnam, 'mean')])
        print('CHECK', mod, np.mean(composites[('present', mod, varnam, 'mean')]))
    compdiffs[(varnam, 'LR')] = np.mean(comps, axis = 0)
    comps = []
    for mod in okmods_hist_H:
        comps.append(composites[('present', mod, varnam, 'mean')]-composites[('present', 'ref', varnam, 'mean')])
        print('CHECK', mod, np.mean(composites[('present', mod, varnam, 'mean')]))
    compdiffs[(varnam, 'HR')] = np.mean(comps, axis = 0)

for varnam in ['tas', 'pr']:
    compdiffs[(varnam, 'diff')] = compdiffs[(varnam, 'HR')]-compdiffs[(varnam, 'LR')]


compcomp = dict()
for cos in ['present', 'future']:
    for varnam in ['tas', 'pr']:
        comps = []
        for mod in okmods_hist_L[:-1]:
            comps.append(composites[(cos, mod, varnam, 'mean')])
        compcomp[(cos, varnam, 'LR')] = np.mean(comps, axis = 0)
        comps = []
        for mod in okmods_hist_H[:-1]:
            comps.append(composites[(cos, mod, varnam, 'mean')])
        compcomp[(cos, varnam, 'HR')] = np.mean(comps, axis = 0)


cmaps = dict()
cmaps['tas'] = 'RdBu_r'
cmaps['pr'] = 'BrBG'
cbar_range = dict()
cbar_range['tas'] = (-5, 5)
cbar_range['pr'] = (-3, 3)
lat = results_ref['lat']
lon = results_ref['lon']
cblab = dict()
cblab['tas'] = 'Temperature (K)'
cblab['pr'] = 'Daily prec (mm)'

margs = [-30,70, 20,80]

for varnam in ['tas', 'pr']:
    if varnam == 'tas':
        fields = composites[('present', 'ref', varnam, 'mean')]
    else:
        fields = 1000*composites[('present', 'ref', varnam, 'mean')]
    filnam = cart_out + 'refcomp_{}.pdf'.format(varnam)
    ctl.plot_multimap_contour(fields, lat, lon, filnam, visualization = 'standard', central_lat_lon = (70, -20), plot_margins = margs, cmap = cmaps[varnam], title = '', subtitles = regnames, cb_label = cblab[varnam], bounding_lat = 0., draw_grid = True, n_color_levels = 10, draw_contour_lines = False, lw_contour = 0.7, cbar_range = cbar_range[varnam])#, plot_type = 'pcolormesh')

for cos in ['present', 'future']:
    for resol in ['LR', 'HR']:
        for varnam in ['tas', 'pr']:
            if varnam == 'tas':
                fields = compcomp[(cos, varnam, resol)]
            else:
                fields = 1000*compcomp[(cos, varnam, resol)]
            filnam = cart_out + 'modcomp_{}_{}_{}.pdf'.format(varnam, cos, resol)
            ctl.plot_multimap_contour(fields, lat, lon, filnam, visualization = 'standard', central_lat_lon = (70, -20), plot_margins = margs, cmap = cmaps[varnam], title = '', subtitles = regnames, cb_label = cblab[varnam], bounding_lat = 0., draw_grid = True, n_color_levels = 10, draw_contour_lines = False, lw_contour = 0.7, cbar_range = cbar_range[varnam])#, plot_type = 'pcolormesh')

for varnam in ['tas', 'pr']:
    if varnam == 'tas':
        fields = composites[('present', 'ref', varnam, 'mean')]
    else:
        fields = 1000*composites[('present', 'ref', varnam, 'mean')]
    filnam = cart_out + 'refcomp_{}.pdf'.format(varnam)
    ctl.plot_multimap_contour(fields, lat, lon, filnam, visualization = 'standard', central_lat_lon = (70, -20), plot_margins = margs, cmap = cmaps[varnam], title = '', subtitles = regnames, cb_label = cblab[varnam], bounding_lat = 0., draw_grid = True, n_color_levels = 10, draw_contour_lines = False, lw_contour = 0.7, cbar_range = cbar_range[varnam])#, plot_type = 'pcolormesh')

cbar_range['tas'] = (-2, 2)
cbar_range['pr'] = (-2, 2)

for varnam in ['tas', 'pr']:
    for cos in ['LR', 'HR', 'diff']:
        if varnam == 'tas':
            fields = compdiffs[(varnam, cos)]
        else:
            fields = 1000*compdiffs[(varnam, cos)]
        filnam = cart_out + 'comp_diff_{}_{}.pdf'.format(varnam, cos)
        ctl.plot_multimap_contour(fields, lat, lon, filnam, visualization = 'standard', central_lat_lon = (70, -20), plot_margins = margs, cmap = cmaps[varnam], title = '', subtitles = regnames, cb_label = cblab[varnam], bounding_lat = 0., draw_grid = True, n_color_levels = 10, draw_contour_lines = False, lw_contour = 0.7, cbar_range = cbar_range[varnam])#, plot_type = 'pcolormesh')

# Fut composite - Hist composite (HR, LR)
compfut = dict()
for varnam in ['tas', 'pr']:
    comps = []
    for mod in okmods_hist_L[:-1]:
        comps.append(composites[('future', mod, varnam, 'mean')]-composites[('present', mod, varnam, 'mean')])
        print('CHECK', mod, np.mean(composites[('present', mod, varnam, 'mean')]))
    compfut[(varnam, 'LR')] = np.mean(comps, axis = 0)
    comps = []
    for mod in okmods_hist_H[:-1]:
        comps.append(composites[('future', mod, varnam, 'mean')]-composites[('present', mod, varnam, 'mean')])
        print('CHECK', mod, np.mean(composites[('present', mod, varnam, 'mean')]))
    compfut[(varnam, 'HR')] = np.mean(comps, axis = 0)

for varnam in ['tas', 'pr']:
    compfut[(varnam, 'diff')] = compfut[(varnam, 'HR')]-compfut[(varnam, 'LR')]

for varnam in ['tas', 'pr']:
    for cos in ['LR', 'HR', 'diff']:
        if varnam == 'tas':
            fields = compfut[(varnam, cos)]
        else:
            fields = 1000*compfut[(varnam, cos)]
        filnam = cart_out + 'comp_futchange_{}_{}.pdf'.format(varnam, cos)
        ctl.plot_multimap_contour(fields, lat, lon, filnam, visualization = 'standard', central_lat_lon = (70, -20), plot_margins = margs, cmap = cmaps[varnam], title = '', subtitles = regnames, cb_label = cblab[varnam], bounding_lat = 0., draw_grid = True, n_color_levels = 10, draw_contour_lines = False, lw_contour = 0.7, cbar_range = cbar_range[varnam])#, plot_type = 'pcolormesh')
