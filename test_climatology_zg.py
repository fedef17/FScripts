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

# cart_out = '/home/fabiano/Research/lavori/WeatherRegimes/check_mean_climatology/obs_trends/'
# ref_file = '/data-hobbes/fabiano/OBS/ERA/ERA40+Int_daily_1957-2018_zg500_remap25_meters.nc'
#
# modnam = 'ERA'
# yearange = (1957, 2018)
# cbar_range = [-1.2,1.2]
# cbar_range_notr = [-1.2,1.2]


cart_out = '/home/fabiano/Research/lavori/WeatherRegimes/check_mean_climatology/future_trends/'
ctl.mkdir(cart_out)
ref_file = '/data-hobbes/fabiano/CMIP6/zg_day_EC-Earth3_ssp585_r4i1p1f1_gr_201501-210012_r25_rc.nc'

modnam = 'EC-Earth3P'
yearange = (2015, 2100)
cbar_range = [0.,3.]
cbar_range_notr = [-1,1]
#############################################################################

var, coords, aux_info = ctl.read_iris_nc(ref_file, extract_level_hPa = 500)
lat = coords['lat']
lon = coords['lon']
dates = coords['dates']

var, dates = ctl.sel_time_range(var, dates, ctl.range_years(yearange[0], yearange[1]))


var_set, dates_set = ctl.seasonal_set(var, dates, 'DJF', seasonal_average = True)
years = np.array([da.year for da in dates_set])

############## PLOT GLOBAL TRENDS ######################

fig, ax = plt.subplots()
glob_mea = ctl.global_mean(var_set, lat)
g0 = glob_mea[0]
m,c = ctl.linear_regre(years, glob_mea)
ax.scatter(years, glob_mea-g0, label = 'Global', color = 'blue')
ax.plot(years, c+m*years-g0, color = 'blue')

var_area, lat_area, lon_area = ctl.sel_area(lat, lon, var_set, 'EAT')
eat_mea = ctl.global_mean(var_area, lat_area)
g0 = eat_mea[0]
m,c = ctl.linear_regre(years, eat_mea)
ax.scatter(years, eat_mea-g0, label = 'EAT', color = 'green')
ax.plot(years, c+m*years-g0, color = 'green')

var_area, lat_area, lon_area = ctl.sel_area(lat, lon, var_set, 'NH')
eat_mea = ctl.global_mean(var_area, lat_area)
g0 = eat_mea[0]
m,c = ctl.linear_regre(years, eat_mea)
ax.scatter(years, eat_mea-g0, label = 'NH', color = 'orange')
ax.plot(years, c+m*years-g0, color = 'orange')
ax.legend()

ax.set_title('Trend in the average zg500 during DJF')
ax.set_ylabel('m')
fig.savefig(cart_out + 'global_{}_trend_{}-{}_DJF.pdf'.format(modnam, yearange[0], yearange[1]))


############## PLOT TREND MAPS ######################

m, c, merr, cerr = ctl.calc_trend_climatevar(years, var_set)

ctl.plot_map_contour(m, lat, lon, plot_anomalies = True, cbar_range = cbar_range, n_color_levels = 11, n_lines = 11, draw_contour_lines=True, draw_grid=True, cb_label = 'Trend (m/year)', filename = cart_out + '{}_trend_{}-{}_DJF.pdf'.format(modnam, yearange[0], yearange[1]))

var_set_notr = []
for ye, va, glo in zip(years, var_set, glob_mea):
    var_set_notr.append(va - glo)
var_set_notr = np.stack(var_set_notr)

m, c, merr, cerr = ctl.calc_trend_climatevar(years, var_set_notr)

ctl.plot_map_contour(m, lat, lon, plot_anomalies = True, cbar_range = cbar_range_notr, n_color_levels = 11, n_lines = 11, draw_contour_lines=True, draw_grid=True, cb_label = 'Trend (m/year)', filename = cart_out + '{}_trend_{}-{}_DJF_noglobal.pdf'.format(modnam, yearange[0], yearange[1]))

sys.exit()

#var_20d = ctl.running_mean(var, 20)
var_set, dates_set = ctl.seasonal_set(var, dates, 'NDJFM', seasonal_average = False)
ndays = var_set.shape[1]

m, c, merr, cerr = ctl.calc_trend_climatevar(years, np.mean(var_set, axis = 1))
ctl.plot_map_contour(m, lat, lon, plot_anomalies = True, cbar_range = cbar_range, n_color_levels = 11, n_lines = 11, draw_contour_lines=True, draw_grid=True, cb_label = 'Trend (m/year)', filename = cart_out + '{}_trend_{}-{}_NDJFM.pdf'.format(modnam, yearange[0], yearange[1]))


window_days = 40
daok = np.arange(0, ndays - window_days, 10)

years = np.array([da[0].year for da in dates_set])

allfigs = []
trend_yea = []
intrc_yea = []
for da in daok:
    print(da)
    allpo = np.stack([np.mean(gigi[da:da+window_days], axis = 0) for gigi in var_set])
    m, c, merr, cerr = ctl.calc_trend_climatevar(years, allpo)

    first_day = dates_set[0][da]
    last_day = dates_set[0][da+window_days]
    tit = '{}/{} to {}/{}'.format(first_day.day, first_day.month, last_day.day, last_day.month)
    fig = ctl.plot_map_contour(m, lat, lon, plot_anomalies = True, cbar_range = cbar_range, n_color_levels = 11, n_lines = 11, draw_contour_lines=True, draw_grid=True, cb_label = 'Trend (m/year)', title = tit)
    allfigs.append(fig)
    trend_yea.append(m)
    intrc_yea.append(c)

trend_yea = np.stack(trend_yea)
intrc_yea = np.stack(intrc_yea)

ctl.plot_pdfpages(cart_out + '{}_trend_{}-{}_NDJFM_alongseason.pdf'.format(modnam, yearange[0], yearange[1]), allfigs)

##############################################################################

sys.exit()

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
