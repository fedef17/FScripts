#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
from matplotlib import pyplot as plt
import netCDF4 as nc
#import cartopy.crs as ccrs
from numpy import linalg as LA
import pickle
import pandas as pd

sys.path.insert(0,'/home/fabiano/Research/git/WRtool/CLUS_tool/WRtool/')
from readsavencfield import read3Dncfield, read4Dncfield, save3Dncfield

sys.path.insert(0,'/home/fabiano/Research/git/EnsClus/clus/')
from sel_season_area import sel_season, sel_area

################################################3

parname = 't2m'
season = 'DJF'
year_ok = 2012
area = 'Med'

cart = '/data-hobbes/fabiano/Medscope/ERAInterim_1d5/'
filename = 'ERAInterim_MSMM_167_grid150.nc'

##################################################

if season == 'JJA':
    months = [6,7,8]
elif season == 'DJF':
    months = [12,1,2]
else:
    raise ValueError('Season not included')


print(filename)
filename = cart + filename

var, lat, lon, dates, time_units, var_units = read3Dncfield(filename)

fh = nc.Dataset(filename, mode='r')
time_cal    = fh.variables['time'].calendar

dates_pdh = pd.to_datetime(dates)

init = pd.to_datetime('19921231', format='%Y%m%d')
fin = pd.to_datetime('20161231', format='%Y%m%d')

mask = (dates_pdh > init) & (dates_pdh < fin)
dates_ok = dates_pdh[mask]
var_ok = var[mask,:,:]

climat_mean = []
climat_std = []
dates_clim = []

for mon in range(1,13):
    mask = (dates_ok.month == mon)
    climat_mean.append(np.mean(var_ok[mask,:,:], axis = 0))
    climat_std.append(np.std(var_ok[mask,:,:], axis = 0))
    dates_clim.append(pd.to_datetime('2000{:02d}01'.format(mon), format='%Y%m%d'))

climat_mean = np.array(climat_mean)
climat_std = np.array(climat_std)

var_anom = []
for el, dat in zip(var, dates_pdh):
    anom = el - climat_mean[dat.month-1,:,:]
    var_anom.append(anom)

var_anom = np.array(var_anom)

var_season, dates_season = sel_season(var_anom, dates, season)
dates_seas_pdh = pd.to_datetime(dates_season)
years = np.unique(dates_seas_pdh.year)

var_season_tot = []

if months[0] == 12:
    years == years[1:]

for yea in years:
    if months[0] == 12:
        init = pd.to_datetime('{:04d}{:02d}01'.format(yea-1, months[0]), format='%Y%m%d')
    else:
        init = pd.to_datetime('{:04d}{:02d}01'.format(yea, months[0]), format='%Y%m%d')
    fin = pd.to_datetime('{:04d}{:02d}01'.format(yea, months[-1]), format='%Y%m%d')
    oks = (dates_seas_pdh >= init) & (dates_seas_pdh <= fin)
    var_season_tot.append(np.mean(var_season[oks], axis = 0))

var_season_tot = np.array(var_season_tot)

var_area, _, _ = sel_area(lat, lon, var_season_tot, area)
print(var_area.shape)
var_anom_mean_area = np.array([var.mean() for var in var_area])

fig = plt.figure()
plt.grid()
plt.scatter(years, var_anom_mean_area)
plt.scatter(years[years == year_ok], var_anom_mean_area[years == year_ok], s = 80)
plt.xlabel('Year')
plt.ylabel('Anomaly (K)')
plt.title('{} mean anomalies on area {} (wrt 1993-2016 average)'.format(season, area))
filename = 'Mean_anom_{}_{}_{}.pdf'.format(area, season, year_ok)
fig.savefig(cart + filename)
