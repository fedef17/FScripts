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

from scipy import stats
import pandas as pd

#######################################

cart_out = '/home/fabiano/Research/lavori/WeatherRegimes/check_mean_climatology/'

cart_in = '/data-hobbes/fabiano/PRIMAVERA/hist_1950/Stream1/'
filelista = '/data-hobbes/fabiano/PRIMAVERA/hist_1950/Stream1/filelist.in'
filo = open(filelista, 'r')
listafils = [cart_in + fi.strip() for fi in filo.readlines()]
filo.close()

model_names = ['CMCC-CM2-HR4', 'CMCC-CM2-VHR4', 'CNRM-CM6-1', 'CNRM-CM6-1-HR', 'EC-Earth-3-LR', 'EC-Earth-3-HR', 'ECMWF-IFS-LR', 'ECMWF-IFS-HR', 'MPI-ESM1-2-HR', 'MPI-ESM1-2-XR', 'HadGEM3-GC31-MM', 'HadGEM3-GC31-HM', 'HadGEM3-GC31-LL-det', 'HadGEM3-GC31-LL-stoc', 'EC-Earth-3P-LR-det', 'EC-Earth-3P-LR-stoc']

listafils.append('/data-hobbes/fabiano/OBS/ERA/ERA40+Int_daily_1957-2018_zg500_remap25_meters.nc')
model_names.append('ERA')

nya = 30 # number of years for running mean
clm_fullperiod = dict()
datesall = dict()
climat_mean_dict = dict()

allnums = [5, 10, 20]

for ifile, name in zip(listafils, model_names):
    print(name)
    var, lat, lon, dates, time_units, var_units, time_cal = ctl.readxDncfield(ifile, extract_level = 500)
    # var, coords, aux_info = ctl.read_iris_nc(ifile)
    # lat = coords['lat']
    # lon = coords['lon']
    # dates = coords['dates']

    var, dates = ctl.sel_time_range(var, dates, ctl.range_years(1957,2014))

    for num in allnums:
        clm, dtclm, _ = ctl.daily_climatology(var, dates, num)
        clmarea, lat_area, lon_area = ctl.sel_area(lat, lon, clm, 'EAT')
        clm_fullperiod[(name, num)] = clmarea
    dtclmpd = pd.to_datetime(dtclm)
    datesall[name] = dtclmpd

    climat_mean, dates_climate_mean = ctl.trend_daily_climat(var, dates, window_days = 20, window_years = nya)
    difftot, lat_area, lon_area = ctl.sel_area(lat, lon, climat_mean[-1]-climat_mean[0], 'EAT')
    climat_mean_dict[name] = difftot

lat_sect = [(16, None), (8, 16), (0, 8)]
lon_sect = [(0, 16), (16, 32), (32, None)]

colors = ctl.color_set(len(model_names)-1, sns_palette = 'Paired')
colors.append('black')

for num in allnums:
    fig = plt.figure(figsize = (16,12))
    i = 0
    for lla1, lla2 in lat_sect:
        axes = []
        for llo1, llo2 in lon_sect:
            lw = 0.7
            i+=1
            ax = plt.subplot(3,3,i)
            for name, col in zip(model_names, colors):
                if name == 'ERA': lw = 2.0
                ax.plot(datesall[name].dayofyear, np.mean(clm_fullperiod[(name, num)][:, lla1:lla2, llo1:llo2], axis = (1,2)), color = col, linewidth = lw)

            if lla2 is None: lla2 = -1
            if llo2 is None: llo2 = -1
            ax.set_title('lat {:2.0f}/{:2.0f}, lon {:3.0f}/{:3.0f}'.format(lat_area[lla1], lat_area[lla2], lon_area[llo1], lon_area[llo2]))
            axes.append(ax)
        ctl.adjust_ax_scale(axes)

    fig.tight_layout()
    fig = ctl.custom_legend(fig, colors, model_names)
    fig.savefig(cart_out + 'mean_climatology_zg_models_{}days.pdf'.format(num))


fig = plt.figure(figsize = (16,12))
i = 0
for lla1, lla2 in lat_sect:
    axes = []
    for llo1, llo2 in lon_sect:
        lw = 0.7
        i+=1
        ax = plt.subplot(3,3,i)
        for name, col in zip(model_names, colors):
            if name == 'ERA': lw = 2.0
            ax.plot(datesall[name].dayofyear, np.mean(climat_mean_dict[name][:, lla1:lla2, llo1:llo2], axis = (1,2)), color = col, linewidth = lw)

        if lla2 is None: lla2 = -1
        if llo2 is None: llo2 = -1
        ax.set_title('lat {:2.0f}/{:2.0f}, lon {:3.0f}/{:3.0f}'.format(lat_area[lla1], lat_area[lla2], lon_area[llo1], lon_area[llo2]))
        axes.append(ax)
    ctl.adjust_ax_scale(axes)

fig.tight_layout()
fig = ctl.custom_legend(fig, colors, model_names)
fig.savefig(cart_out + 'trend_climatology_zg_models_last30-first30.pdf')
