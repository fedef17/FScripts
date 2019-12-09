#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
from matplotlib import pyplot as plt
from matplotlib import cm

import pickle

import climtools_lib as ctl
import climdiags as cd

from matplotlib.colors import LogNorm
from datetime import datetime

from scipy import stats
import pandas as pd

#######################################
cart_out = '/home/fabiano/Research/articoli/Papers/primavera_regimes/figures/SST_biases/'
if not os.path.exists(cart_out): os.mkdir(cart_out)

cart_out_maps = cart_out + 'maps/'
if not os.path.exists(cart_out_maps): os.mkdir(cart_out_maps)

cart = '/home/fabiano/Research/lavori/WeatherRegimes/prima_coup_v6/'
filogen = cart + 'out_prima_coup_v6_DJF_EAT_4clus_4pcs_1957-2014_refEOF.p'

model_names = ['AWI-CM-1-0-LR', 'AWI-CM-1-0-HR', 'CMCC-CM2-HR4', 'CMCC-CM2-VHR4', 'CNRM-CM6-1', 'CNRM-CM6-1-HR', 'EC-Earth3P', 'EC-Earth3P-HR', 'ECMWF-IFS-LR', 'ECMWF-IFS-MR', 'ECMWF-IFS-HR', 'MPI-ESM1-2-HR', 'MPI-ESM1-2-XR', 'HadGEM3-GC31-LL', 'HadGEM3-GC31-MM', 'HadGEM3-GC31-HM', 'HadGEM3-GC31-HH']
vers = ['LR', 'HR', 'LR', 'HR', 'LR', 'HR', 'LR', 'HR', 'LR', 'MR', 'HR', 'LR', 'HR', 'LR', 'MR', 'MR', 'HR'] + ['OBS']
#model_names = ['CMCC-CM2-HR4', 'CMCC-CM2-VHR4', 'CNRM-CM6-1', 'CNRM-CM6-1-HR', 'EC-Earth3P', 'EC-Earth3P-HR', 'ECMWF-IFS-LR', 'ECMWF-IFS-HR', 'MPI-ESM1-2-HR', 'MPI-ESM1-2-XR', 'HadGEM3-GC31-LL', 'HadGEM3-GC31-MM']#, 'HadGEM3-GC31-HM']
#ens_mems = ['r1i1p1f002', 'r1i1p1f002', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f2', 'r1i1p1f2', 'r1i1p2f1', 'r1i1p2f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1']
#ens_mems = ['r1i1p1f1', 'r1i1p1f1', 'r1i1p1f2', 'r1i1p1f2', 'r3i1p2f1', 'r1i1p2f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1']

#vers = 7*['LR', 'HR'] + ['OBS']
model_coups = ['AWI-CM-1-0', 'CMCC-CM2', 'CNRM-CM6-1', 'EC-Earth3P', 'ECMWF-IFS', 'MPI-ESM1-2', 'HadGEM3-GC31']
#model_coups = ['CMCC-CM2', 'CNRM-CM6-1', 'EC-Earth3P', 'ECMWF-IFS', 'MPI-ESM1-2', 'HadGEM3-GC31']
model_names_all = model_names + ['ERA']

colors = ctl.color_set(len(model_coups)*2, sns_palette = 'Paired')
mr_cos = np.where(np.array(vers) == 'MR')[0]
for gi in mr_cos:
    colors.insert(gi, np.mean(colors[gi-1:gi+1], axis = 0))

colors_wERA = colors + ['black']
# colors_wERA2 = np.array(colors_wERA)
# colors_wERA2[0:-1:2] = colors_wERA2[1::2]

regnam = ['NAO +', 'Sc. BL', 'AR', 'NAO -']

################################################################################

# results_refEOF, results_ref = pickle.load(open(filogen, 'rb'))
# results_refEOF['ERA'] = results_ref

cart = '/home/fabiano/Research/lavori/WeatherRegimes/panosjet/'
filogen = cart + 'out_panosjet_DJF_EAT_4clus_4pcs_1957-2014_refCLUS.p'
results_refCLUS, results_ref = pickle.load(open(filogen, 'rb'))
results_refCLUS['ERA'] = results_ref

cart = '/home/fabiano/Research/lavori/WeatherRegimes/panosjet/'
filogen = cart + 'out_panosjet_DJF_EAT_4clus_4pcs_1957-2014_refEOF.p'
results_refEOF, results_ref = pickle.load(open(filogen, 'rb'))
results_refEOF['ERA'] = results_ref

##################################################################
cart_sst = '/data-woodstock/PRIMAVERA/stream1/merged/tos/'
filsst = 'tos-{}-1950-2014-{}-remap.nc'

#cart_sst = '/nas/PRIMAVERA/Stream1/hist-1950/{}/{}/tos/'
filera = '/nas/reference/ERA40+Int/sst_Amon_ERA40_195701-201812_1deg.nc'
sstera, datacoords, _ = ctl.readxDncfield(filera)
datesera = datacoords['dates']
lat = datacoords['lat']
lon = datacoords['lon']

sstera = sstera - 273.15
sstera_mean, sstera_sd = ctl.seasonal_climatology(sstera, datesera, 'DJF', dates_range = ctl.range_years(1957, 2014))

area_box = (-80, 10, 20, 80)

sstera_mean_area, latsel, lonsel = ctl.sel_area(lat, lon, sstera_mean, area_box)
okpoera = (sstera_mean_area < -100) | (sstera_mean_area > 500)

# allrms = dict()
# allpatcor = dict()
# #for mod, mem in zip(model_names, ens_mems):
# for mod in model_names:
#     print(mod)
#     filmod_part = 'tos-{}-1950-2014-'.format(mod)
#     if 'AWI' in mod:
#         filmod_part = 'tos-{}-1950-2010-'.format(mod)
#     allfimod = [fi for fi in os.listdir(cart_sst) if filmod_part in fi and 'remap.nc' in fi]
#     print(allfimod)
#
#     for fi in allfimod:
#         filmod = cart_sst + fi
#         i0 = len(filmod_part)
#         mem = fi[i0:i0+8]
#         print(mem)
#         if 'EC-Earth' in mod and mem == 'r1i1p1f1':
#             continue
#         if 'HadGEM' in mod and mem == 'r1i1p2f1':
#             continue
#         sstmod, datacoords, _ = ctl.readxDncfield(filmod)
#         datesmod = datacoords['dates']
#         lat = datacoords['lat']
#         lon = datacoords['lon']
#
#         if sstmod.min() > 200:
#             sstmod = sstmod - 273.15
#
#         sstmod_mean, sstmod_sd = ctl.seasonal_climatology(sstmod, datesmod, 'DJF', dates_range = ctl.range_years(1957, 2014))
#
#         # compare
#         sstmod_mean_area, latsel, lonsel = ctl.sel_area(lat, lon, sstmod_mean, area_box)
#
#         okpomod = (sstmod_mean_area < -100) | (sstmod_mean_area > 500)
#         oktot = (okpomod) | (okpoera)
#
#         sstmod_x = np.ma.masked_array(sstmod_mean_area, mask = oktot)
#         sstera_x = np.ma.masked_array(sstera_mean_area, mask = oktot)
#
#         ctl.plot_triple_sidebyside(sstmod_x, sstera_x, latsel, lonsel, plot_type = 'pcolormesh', filename = cart_out_maps + 'map_EAT_{}_{}.pdf'.format(mod, mem), plot_margins = area_box, title = 'DJF SST bias - {} - {}'.format(mod, mem), stitle_1 = mod, stitle_2 = 'ERA', cb_label = 'SST bias (K)')
#
#         rms = ctl.E_rms(sstmod_x, sstera_x, latitude = latsel, masked = True)
#         patcor = ctl.Rcorr(sstmod_x, sstera_x, latitude = latsel, masked = True)
#         print(rms, patcor)
#         allrms[(mod, mem)] = rms
#         allpatcor[(mod, mem)] = patcor
#
#
# pickle.dump([allrms, allpatcor], open(cart_out + 'sst_bias_rms_djf_eat.p', 'wb'))
allrms, allpatcor = pickle.load(open(cart_out + 'sst_bias_rms_djf_eat.p', 'rb'))

fig = plt.figure(figsize = (16, 12))
ax = plt.subplot()
ax.set_ylabel('RMS SST bias (K)')
ax.set_title('RMS SST bias in North Atlantic')
ax.set_xticks([])
i = 0
wi = 0.6
#for mod, col, vv, mem in zip(model_names, colors, vers, ens_mems):
    #ax.bar(i, allrms[(mod, mem)], width = wi, color = col)
for mod, col, vv in zip(model_names, colors, vers):
    rmsall = []
    for ke in allrms.keys():
        if mod in ke:
            rmsall.append(allrms[ke])

    allrms[mod] = np.mean(rmsall)

    if len(rmsall) == 1:
        ax.bar(i, rmsall[0], width = wi, color = col)
    else:
        print('{} has {} ens members'.format(mod, len(rmsall)))
        okrms = np.mean(rmsall)
        rmserr_lo = np.min(rmsall)
        rmserr_hi = np.max(rmsall)
        print(okrms, rmserr_lo, rmserr_hi)
        yerr = np.array([[okrms-rmserr_lo, rmserr_hi-okrms]]).T
        # yerr.shape
        # # ax.errorbar(1,1, yerr = yerr)
        # # ax.errorbar(1,1, yerr = yerr, capsize=10)
        ax.bar(i, okrms, width = wi, color = col, yerr = yerr, capsize = 5)

    i += 0.7
    if vv == 'HR': i += 0.2

all_LR = np.mean([allrms[mod] for mod, vv in zip(model_names, vers) if vv == 'LR'])
all_HR = np.mean([allrms[mod] for mod, vv in zip(model_names, vers) if vv == 'HR'])

i+=0.2
col_LR = 'grey'
col_HR = 'black'
ax.bar(i, all_LR, width = wi, color = col_LR)
i += 0.7
ax.bar(i, all_HR, width = wi, color = col_HR)

ctl.custom_legend(fig, colors, model_names, ncol = 6)

fig.savefig(cart_out + 'sst_bias_RMS_DJF_EAT.pdf')


fig = plt.figure(figsize = (16, 12))
ax = plt.subplot()
ax.set_ylabel('(1-patcor)')
ax.set_title('Pattern correlation with observed SSTs in the North Atlantic')
ax.set_xticks([])
i = 0
wi = 0.6
#for mod, col, vv, mem in zip(model_names, colors, vers, ens_mems):
    #ax.bar(i, allrms[(mod, mem)], width = wi, color = col)
for mod, col, vv in zip(model_names, colors, vers):
    rmsall = []
    for ke in allpatcor.keys():
        if mod in ke:
            rmsall.append(allpatcor[ke])

    allpatcor[mod] = np.mean(rmsall)

    if len(rmsall) == 1:
        ax.bar(i, 1-rmsall[0], width = wi, color = col)
    else:
        print('{} has {} ens members'.format(mod, len(rmsall)))
        okrms = np.mean(rmsall)
        rmserr_lo = np.min(rmsall)
        rmserr_hi = np.max(rmsall)
        print(okrms, rmserr_lo, rmserr_hi)
        yerr = np.array([[okrms-rmserr_lo, rmserr_hi-okrms]]).T
        # yerr.shape
        # # ax.errorbar(1,1, yerr = yerr)
        # # ax.errorbar(1,1, yerr = yerr, capsize=10)
        ax.bar(i, 1-okrms, width = wi, color = col, yerr = yerr[::-1], capsize = 5)

    i += 0.7
    if vv == 'HR': i += 0.2

all_LR = np.mean([allpatcor[mod] for mod, vv in zip(model_names, vers) if vv == 'LR'])
all_HR = np.mean([allpatcor[mod] for mod, vv in zip(model_names, vers) if vv == 'HR'])

i+=0.2
col_LR = 'grey'
col_HR = 'black'
ax.bar(i, 1-all_LR, width = wi, color = col_LR)
i += 0.7
ax.bar(i, 1-all_HR, width = wi, color = col_HR)

ctl.custom_legend(fig, colors, model_names, ncol = 6)

fig.savefig(cart_out + 'sst_bias_PATCOR_DJF_EAT.pdf')
