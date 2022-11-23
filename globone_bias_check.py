#!/usr/bin/env python
# coding: utf-8

import numpy as np
import sys
import os
from matplotlib import pyplot as plt
import pickle

import climtools_lib as ctl
import climdiags as cd

from scipy import stats
import xarray as xr
import glob

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 22
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['axes.axisbelow'] = True

#############################################################################

### script for producing some maps of clouds, radiation and precipitation
#### inputs:
if len(sys.argv) < 2:
    print('Specify path of experiment runtime and output path in the call (e.g. python globone_bias_check.py cart_in cart_out)')
    sys.exit()
else:
    cart = sys.argv[1] # Name of input file (relative path)
    cart_out = sys.argv[2]

globon = xr.open_mfdataset(cart + 'GLOBONE_atm_6hrs_*.nc')

#### CLOUDS
cloud = globon['clt']
cloud_clim = cloud.mean('time')

ceres_fi = '/nas/reference/CERES/CERES_2000-2015_monclim.nc'
ceres = xr.load_dataset(ceres_fi)
ceres_clim = ceres.mean('time')

cloud_clim = cloud_clim.rename({'longitude': 'lon','latitude': 'lat'})
cloud_clim_rg = ctl.regrid_dataset(cloud_clim, regrid_to_reference=ceres_clim)

era_cloud = xr.open_mfdataset('/nas/reference/ERAInterim/monthly/tcc/tcc*nc')
era_cloud_clim = era_cloud.mean('time')
era_cloud_clim_rg = ctl.regrid_dataset(era_cloud_clim, regrid_to_reference=ceres_clim)

ctl.plot_multimap_contour([cloud_clim_rg, ceres_clim['cldarea_total_mon'], 100*era_cloud_clim_rg['tcc']], figsize = (18,7), subtitles = ['GLOBO KM 312', 'CERES', 'ERAInt'], cb_label = 'Cloud cover', fix_subplots_shape=(1,3), filename = cart_out + 'cloud_check.pdf')

# ## Shortwave flux at TOA
varname = 'rsnt'
var = globon[varname].mean('time').rename({'longitude': 'lon','latitude': 'lat'})
var_ceres = ceres_clim['toa_solar_all_mon']-ceres_clim['toa_sw_all_mon']

var_era = xr.open_mfdataset('/nas/reference/ERAInterim/monthly/rsnt/rsnt*nc')
var_era = var_era['tsr'].mean('time')
var_rg = ctl.regrid_dataset(var, regrid_to_reference=var_ceres)
var_era_rg = ctl.regrid_dataset(var_era, regrid_to_reference=var_ceres)

#ctl.plot_multimap_contour([var_rg, var_ceres, var_era_rg], figsize = (18,7), subtitles = ['GLOBO KM 312', 'CERES', 'ERA Int'], cb_label = 'Net SW at TOA', fix_subplots_shape=(1,3))

#ctl.plot_multimap_contour([var_rg, var_ceres, var_rg-var_ceres], figsize = (18,7), subtitles = ['GLOBO KM 312', 'CERES', 'DIFF'], cb_label = 'Net SW at TOA', fix_subplots_shape=(1,3))
ctl.plot_map_contour(var_rg-var_ceres, cb_label = 'Net SW at TOA bias', plot_anomalies = True, cbar_range=(-50, 50), filename = cart_out + 'net_SW_bias.pdf')

# ## Precipitation
varname = 'pr'
var = globon[varname].mean('time').rename({'longitude': 'lon','latitude': 'lat'})
var_era = xr.open_mfdataset('/nas/reference/ERAInterim/monthly/pr/pr*nc')
var_era = var_era['tp'].rename({'longitude': 'lon','latitude': 'lat'}).mean('time')
var_rg = ctl.regrid_dataset(var, regrid_to_reference=var_era)


ctl.plot_multimap_contour([1.5*var_rg*86400., var_era*1000], figsize = (16,9), subtitles = ['GLOBO KM 312', 'ERA'], cb_label = 'Precipit', fix_subplots_shape=(1,2), cbar_range = (0, 10), cmap='BrBG', filename = cart_out + 'check_pr.pdf')
ctl.plot_map_contour((var_rg*86400-var_era*1000), figsize = (16,9), title = ['GLOBO KM 312 - ERA'], cb_label = 'Precipit', cbar_range = (-5, 5), cmap='BrBG', filename = cart_out + 'check_pr_bias.pdf')
