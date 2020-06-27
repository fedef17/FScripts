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

parname = 't2m'

cart = '/data-hobbes/fabiano/Medscope/ERAInterim_1d5/'
filename = 'ERAInterim_MSMM_167_grid150.nc'
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

climate_mean = []
climat_std = []
dates_clim = []

for mon in range(1,13):
    mask = (dates_ok.month == mon)
    climate_mean.append(np.mean(var_ok[mask,:,:], axis = 0))
    climat_std.append(np.std(var_ok[mask,:,:], axis = 0))
    dates_clim.append(pd.to_datetime('2000{:02d}01'.format(mon), format='%Y%m%d'))

climate_mean = np.array(climate_mean)
climat_std = np.array(climat_std)

filename = 'ERAInterim_167_grid150_mean_1993-2016.nc'
print(filename)
filename = cart + filename
save3Dncfield(lat,lon,climate_mean,parname,var_units,dates_clim,time_units,time_cal,filename)

filename = 'ERAInterim_167_grid150_std_1993-2016.nc'
print(filename)
filename = cart + filename
save3Dncfield(lat,lon,climat_std,parname,var_units,dates_clim,time_units,time_cal,filename)

var_anom = []
for el, dat in zip(var, dates_pdh):
    anom = el - climate_mean[dat.month-1,:,:]
    var_anom.append(anom)

var_anom = np.array(var_anom)

filename = 'ERAInterim_anomalies_167_grid150.nc'
print(filename)
filename = cart + filename
save3Dncfield(lat,lon,var_anom,parname,var_units,dates,time_units,time_cal,filename)
