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

era_ref_file = '/data-hobbes/fabiano/OBS/ERA/ERA40+Int_daily_1957-2018_zg500_remap25_meters.nc'

var, coords, aux_info = ctl.read_iris_nc(era_ref_file, extract_level_hPa = 500)
lat = coords['lat']
lon = coords['lon']
dates = coords['dates']

nya = 30 # number of years for running mean

climat_mean_dict = dict()

for num in [1, 5, 10, 15, 20]:
    climat_mean, dates_climate_mean = ctl.trend_daily_climat(var, dates, window_days = num, window_years = nya)
    climat_mean_area = []
    for cos in climat_mean:
        nucos, lat_area, lon_area = ctl.sel_area(lat, lon, cos, 'EAT')
        climat_mean_area.append(nucos)
    climat_mean_dict[num] = climat_mean_area

all_years_ref = np.array([dat[0].year for dat in dates_climate_mean])

allnums = [1,5,10,15,20]

clm_fullperiod = dict()
for num in [1,5,10,15,20]:
    clm, dtclm, _ = ctl.daily_climatology(var, dates, num)
    clmarea, _, _ = ctl.sel_area(lat, lon, clm, 'EAT')
    clm_fullperiod[num] = clmarea

dtclmpd = pd.to_datetime(dtclm)

lat_sect = [(16, None), (8, 16), (0, 8)]
lon_sect = [(0, 16), (16, 32), (32, None)]

colors = ctl.color_set(len(allnums))

fig = plt.figure(figsize = (16,12))
i = 0
for lla1, lla2 in lat_sect:
    axes = []
    for llo1, llo2 in lon_sect:
        i+=1
        ax = plt.subplot(3,3,i)
        for num, col in zip(allnums, colors):
            ax.plot(dtclmpd.dayofyear, np.mean(clm_fullperiod[num][:, lla1:lla2, llo1:llo2], axis = (1,2)), color = col)

        if lla2 is None: lla2 = -1
        if llo2 is None: llo2 = -1
        ax.set_title('lat {:2.0f}/{:2.0f}, lon {:3.0f}/{:3.0f}'.format(lat_area[lla1], lat_area[lla2], lon_area[llo1], lon_area[llo2]))
        axes.append(ax)
    ctl.adjust_ax_scale(axes)

fig.tight_layout()
fig = ctl.custom_legend(fig, colors, ['{}'.format(i) for i in allnums])
fig.savefig(cart_out + 'mean_climatology_zg.pdf')

figs = []
axaxes = [[], [], []]
colors = ctl.color_set(len(climat_mean_area))

for num in allnums:
    cosone = climat_mean_dict[num]
    ref = clm_fullperiod[num]

    fig = plt.figure(figsize = (16,12))
    i = 0
    for uu, (lla1, lla2) in enumerate(lat_sect):
        for llo1, llo2 in lon_sect:
            i+=1
            ax = plt.subplot(3,3,i)
            for cos, col in zip(cosone, colors):
                ax.plot(dtclmpd.dayofyear, np.mean(cos[:, lla1:lla2, llo1:llo2], axis = (1,2)), color = col)
            ax.plot(dtclmpd.dayofyear, np.mean(ref[:, lla1:lla2, llo1:llo2], axis = (1,2)), color = 'black', linewidth = 2)

            if lla2 is None: lla2 = -1
            if llo2 is None: llo2 = -1
            ax.set_title('lat {:2.0f}/{:2.0f}, lon {:3.0f}/{:3.0f}'.format(lat_area[lla1], lat_area[lla2], lon_area[llo1], lon_area[llo2]))
            axaxes[uu].append(ax)

    fig.tight_layout()
    fig.subplots_adjust(top = 0.9)
    fig.suptitle('{} days running mean'.format(num))

    labels = ['{}'.format(ye) for ye in all_years_ref]
    fig = ctl.custom_legend(fig, colors, labels)
    figs.append(fig)

for uu in range(3):
    ctl.adjust_ax_scale(axaxes[uu])

for fig, nam in zip(figs, allnums):
    fig.savefig(cart_out + 'mean_climatology_zg_trend{}days_{}years.pdf'.format(num, nya))
