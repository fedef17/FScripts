#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os

from matplotlib import pyplot as plt
import matplotlib.cm as cm
import matplotlib.patheffects as PathEffects

import netCDF4 as nc
import cartopy.crs as ccrs
import pandas as pd

from numpy import linalg as LA
from eofs.standard import Eof
from scipy import stats
import itertools as iter

from sklearn.cluster import KMeans


from datetime import datetime
import pickle

import climtools_lib as ctl


cart_ERA = '/home/fabiano/Research/lavori/WeatherRegimes/OUT_WRTOOL/ERA_ref/'
cart_NCEP = '/home/fabiano/Research/lavori/WeatherRegimes/OUT_WRTOOL/NCEP_ref/'
cart = '/home/fabiano/Research/lavori/WeatherRegimes/OUT_WRTOOL/OUT_ECMWF-IFS_allHR_1979-2008/'

# ncep_file = 'area_anomaly_zg500_day_NCEPNCAR_144x73_1ens_DJF_EAT_1979-2008_0.nc'
# era_file = 'area_anomaly_zg500_day_ERAInterim_144x73_1ens_DJF_EAT_1979-2008_0.nc'
# mod_file = ''
#
# var_ERA, lat, lon, dates_ERA, time_units, var_units = ctl.readxDncfield(cart_ERA+era_file)
# var_NCEP, lat, lon, dates_NCEP, time_units, var_units = ctl.readxDncfield(cart_NCEP+ncep_file)
#
# indxs_ERA = np.loadtxt(cart_ERA+'indclORD_4clus_zg500_day_ERAInterim_144x73_1ens_DJF_EAT_1979-2008_4pcs.txt')
# indxs_NCEP = np.loadtxt(cart_NCEP+'indclORD_4clus_zg500_day_NCEPNCAR_144x73_1ens_DJF_EAT_1979-2008_4pcs.txt')
#
# resid_times_ERA = ctl.calc_regime_residtimes(indxs_ERA, dates = dates_ERA)
# resid_times_ERA_noinc = ctl.calc_regime_residtimes(indxs_ERA, dates = dates_ERA, count_incomplete = False)
#
# resid_times_NCEP = ctl.calc_regime_residtimes(indxs_NCEP, dates = dates_NCEP)
# resid_times_NCEP_noinc = ctl.calc_regime_residtimes(indxs_NCEP, dates = dates_NCEP, count_incomplete = False)


data_ERA = nc.Dataset(cart_ERA+'cluspatternORD_4clus_zg500_day_ERAInterim_144x73_1ens_DJF_EAT_1979-2008_4pcs.nc')
data_mod = nc.Dataset(cart+'cluspatternORD_4clus_zg500_day_ECMWF-IFS_144x73_4ens_DJF_EAT_1979-2008_4pcs.nc')

solver_ERA = pickle.load(open(cart_ERA+'solver_zg500_day_ERAInterim_144x73_1ens_DJF_EAT_1979-2008_4pcs.p'))

clus_ERA_tot = data_ERA.variables['cluspattern'][:]
clus_mod_tot = data_mod.variables['cluspattern'][:]
lat_era = data_ERA.variables['lat'][:]
lon_era = data_ERA.variables['lon'][:]
latgr_era, longr_era = np.meshgrid(lat_era, lon_era)
lat_mod = data_mod.variables['lat'][:]
lon_mod = data_mod.variables['lon'][:]
latgr_mod, longr_mod = np.meshgrid(lat_era, lon_era)

clus_mod_tot_ok, lat_mod, lon_mod = ctl.check_increasing_latlon(clus_mod_tot, lat_mod, lon_mod)

clus_ERA_PC = []
clus_mod_PC = []
for i, (co,chi) in enumerate(zip(clus_ERA_tot, clus_mod_tot_ok)):
    clus_ERA_PC.append(solver_ERA.projectField(co, neofs=10,eofscaling=0,weighted=True))
    clus_mod_PC.append(solver_ERA.projectField(chi, neofs=10,eofscaling=0,weighted=True))
    print(i, chi.min(), chi.max())

ok_perm = ctl.match_pc_sets(clus_ERA_PC, clus_mod_PC)
print(ok_perm)

fr_eras = [32.57, 27.09, 20.19, 20.15]
fr_ecms = [25.91, 29.10, 20.27, 24.72]

cbar_range = ctl.get_cbar_range(np.stack([clus_ERA_tot, clus_mod_tot_ok]), symmetrical = True, percentiles = (2,98))
for i, p, fr_era, fr_ecm in zip(list(range(4)), ok_perm, fr_eras, fr_ecms):
    ctl.plot_double_sidebyside(clus_ERA_tot[i], clus_mod_tot_ok[p], lat_era, lon_era, filename = cart+'Cluster_{}_vs_ERA.pdf'.format(i), title = 'Cluster {}'.format(i), cb_label = 'Geopotential height anomaly (m)', cbar_range = cbar_range, stitle_1 = 'ERA - {}%'.format(fr_era), stitle_2 = 'ECMWF-IFS-HR - {}%'.format(fr_ecm))
    ctl.plot_double_sidebyside(clus_ERA_tot[i], clus_mod_tot_ok[p], lat_era, lon_era, filename = cart+'Cluster_{}_vs_ERA_orto.pdf'.format(i), title = 'Cluster {}'.format(i), cb_label = 'Geopotential height anomaly (m)', cbar_range = cbar_range, stitle_1 = 'ERA - {}%'.format(fr_era), stitle_2 = 'ECMWF-IFS-HR - {}%'.format(fr_ecm), visualization = 'polar')

# print('\nRMS error\n')
# for i in range(4):
#     pini = []
#     for j in range(4):
#         pino = ctl.E_rms(clus_ERA[i], clus_mod[j])
#         pini.append(pino)
#     print(i,np.argmin(pini),pini)
#
# print('\nCorrelation\n')
# for i in range(4):
#     pini = []
#     for j in range(4):
#         pino = ctl.Rcorr(clus_ERA[i], clus_mod[j])
#         pini.append(pino)
#     print(i,np.argmax(pini),pini)
#
# print('\nCP RMS error\n')
# for i in range(4):
#     pini = []
#     for j in range(4):
#         pino = ctl.E_rms_cp(clus_ERA[i], clus_mod[j])
#         pini.append(pino)
#     print(i,np.argmin(pini),pini)
