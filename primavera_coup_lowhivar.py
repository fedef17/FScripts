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

#######################################

season = 'DJF'
area = 'EAT'
projtype = 'robinson'
projtypemf = 'Nstereo'

cart_out = '/home/fabiano/Research/lavori/PRIMAVERA_bias/hist_1950/'
if not os.path.exists(cart_out): os.mkdir(cart_out)

vtag = 'v7'
cart = '/home/fabiano/Research/lavori/WeatherRegimes/prima_coup_{}/'.format(vtag)
filogen = cart + 'out_prima_coup_{}_DJF_EAT_4clus_4pcs_1957-2014_refEOF.p'.format(vtag)
results, results_ref = pickle.load(open(filogen, 'rb'))

results.pop('HadGEM3-GC31-LL_r1i1p2f1')
results.pop('EC-Earth3P_r1i1p1f1')
results.pop('EC-Earth3P-HR_r1i1p1f1')
allresmembers = results.keys()
del results

##############################################################

model_names = ['AWI-CM-1-1-LR', 'AWI-CM-1-1-HR', 'CMCC-CM2-HR4', 'CMCC-CM2-VHR4', 'CNRM-CM6-1', 'CNRM-CM6-1-HR', 'EC-Earth3P', 'EC-Earth3P-HR', 'ECMWF-IFS-LR', 'ECMWF-IFS-MR', 'ECMWF-IFS-HR', 'MPI-ESM1-2-HR', 'MPI-ESM1-2-XR', 'HadGEM3-GC31-LL', 'HadGEM3-GC31-MM', 'HadGEM3-GC31-HM', 'HadGEM3-GC31-HH']

vers = ['LR', 'HR', 'LR', 'HR', 'LR', 'HR', 'LR', 'HR', 'LR', 'MR', 'HR', 'LR', 'HR', 'LR', 'MR', 'MR', 'HR'] + ['OBS']
model_coups = ['AWI-CM-1-1', 'CMCC-CM2', 'CNRM-CM6-1', 'EC-Earth3P', 'ECMWF-IFS', 'MPI-ESM1-2', 'HadGEM3-GC31']

model_names_all = model_names + ['ERA']

colors = ctl.color_set(len(model_coups)*2, sns_palette = 'Paired')
colorscoup = [np.mean([col1, col2], axis = 0) for col1, col2 in zip(colors[:-1:2], colors[1::2])]
color_main = []
for col2 in colors[1::2]:
    color_main += [col2, col2]

mr_cos = np.where(np.array(vers) == 'MR')[0]
for gi in mr_cos:
    colors.insert(gi, np.mean(colors[gi-1:gi+1], axis = 0))
    color_main.insert(gi, colors[gi+1])

colors_wERA = colors + ['black']
color_main.append('black')

#############################################################################################
#############################################################################################

mean_field_all = dict()
lowfrvar = dict()
highfrvar = dict()
stat_eddy_all = dict()

file_in = '/data-hobbes/fabiano/OBS/ERA/ERA40+Int_daily_1957-2018_zg500_remap25_meters.nc'

#var, coords, aux_info = ctl.read_iris_nc(file_in, extract_level_hPa = 500)
var, coords, aux_info = ctl.readxDncfield(file_in, extract_level = 500)
lat = coords['lat']
lon = coords['lon']
dates = coords['dates']

mean_field_all['ERA'], lowfrvar['ERA'], highfrvar['ERA'], stat_eddy_all['ERA'] = ctl.variability_lowhi(lat, lon, var, dates, season, area = area, dates_range = ctl.range_years(1957, 2014))

cart_in = '/data-hobbes/fabiano/PRIMAVERA/incoming/hist-1950/'
flist = cart_in + 'lista_filez_v7.dat'
fi = open(flist, 'r')
lista_files = [lin.rstrip() for lin in fi.readlines()]
fi.close()


#mean_RMSbias_all = dict()
for filo in lista_files:
    metadata = ctl.cmip6_naming(filo)
    mod_name = metadata['model'] + '_' + metadata['member']

    #var, coords, aux_info = ctl.read_iris_nc(file_in, extract_level_hPa = 500)
    var, coords, aux_info = ctl.readxDncfield(cart_in + filo, extract_level = 500)
    lat = coords['lat']
    lon = coords['lon']
    dates = coords['dates']

    mean_field_all[mod_name], lowfrvar[mod_name], highfrvar[mod_name], stat_eddy_all[mod_name] = ctl.variability_lowhi(lat, lon, var, dates, season, area = area, dates_range = ctl.range_years(1957, 2014))

    #mean_RMSbias_all[mod_name] = ctl.E_rms(mean_field_all[mod_name], mean_field_all['ERA'], latitude = lat)

pickle.dump([mean_field_all, lowfrvar, highfrvar, stat_eddy_all], open(cart_out + 'out_lowhighstat_var.p', 'wb'))
