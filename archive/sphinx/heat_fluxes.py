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
import itertools as itt

from sklearn.cluster import KMeans



from datetime import datetime
import pickle

import climtools_lib as ctl
import climdiags as cd

from copy import deepcopy as cp

##############################################
### constants
# cp = 4185.5 # J kg-1 K-1 QUESTA E L'ACQUA COJOOOO
cp = 1005.0 # specific enthalpy dry air - J kg-1 K-1
cpw = 1840.0 # specific enthalpy water vapor - J kg-1 K-1
# per la moist air sarebbe cpm = cp + q*cpw, ma conta al massimo per un 3 %
L = 2501000.0 # J kg-1
Lsub = 2835000.0 # J kg-1
g = 9.81 # m s-2
Rearth = 6371.0e3 # mean radius
#####################################

cart_in = '/data-hobbes/fabiano/SPHINX/heat_flux/'
cart_out = '/home/fabiano/Research/lavori/SPHINX_for_lisboa/heat_flux/'

expnames = ['lcb0','lcb1', 'lcb2','lcs0','lcs1','lcs2']
varnames = ['mshf', 'mpef', 'mlhf']
fluxnames = ['Sensible Heat', 'Potential Energy', 'Latent Heat']
factors = [-cp/g, -1., -L/g]

yearsall = np.arange(1950,2101)
all_fluxes_yrl = dict()
monthly_clim_hist = dict()
monthly_clim_fut = dict()
for exp in expnames:
    for varna, flun, fac in zip(varnames, fluxnames, factors):
        # leggi da file e salva in dict
        fil = cart_in + '{}_{}_1950-2100.nc'.format(exp, varna)
        print(fil)
        var, lat, lon, dates, time_units, var_units = ctl.read3Dncfield(fil)
        yearsok = np.unique(pd.to_datetime(dates).year)
        for y in yearsall:
            if y not in yearsok: print('{} {} - {} missing\n'.format(exp, varna, y))
        varfac = var*fac
        all_fluxes_yrl[(exp, varna)], yrl_dates = ctl.yearly_average(varfac, dates)
        monthly_clim_hist[(exp, varna)], mon_dates, _ = ctl.monthly_climatology(varfac, dates, dates_range = ctl.range_years(1950,2025))
        monthly_clim_fut[(exp, varna)], mon_dates, _ = ctl.monthly_climatology(varfac, dates, dates_range = ctl.range_years(2026,2100))

    all_fluxes_yrl[(exp, 'all')] = np.sum(np.stack([all_fluxes_yrl[cos] for cos in all_fluxes_yrl.keys() if exp in cos]), axis=0)

varnames.append('all')
fluxnames.append('Total')

map_means = dict()
zonal_ints = dict()

zonal_factor = 2*np.pi*Rearth*np.cos(np.deg2rad(lat))

for coso, cosotto in zip(['lcb','lcs'], ['base', 'stoc']):
    for varna in varnames:
        cosiok = [all_fluxes_yrl[(exp, varna)] for exp in expnames if coso in exp]
        all_fluxes_yrl[(cosotto, varna)] = np.mean(np.stack(cosiok), axis = 0)

for key in all_fluxes_yrl:
    map_means[key] = np.mean(all_fluxes_yrl[key], axis = 0)
    zonal_ints[key] = np.mean(all_fluxes_yrl[key], axis = -1)*zonal_factor

# for varna, flun in zip(varnames, fluxnames):
#     ctl.plot_double_sidebyside(map_means[('base',varna)], map_means[('stoc',varna)], lat, lon, title = flun)

plt.figure()
plt.title('present day - base vs stoc')
plt.plot(lat, zonal_ints[('base','all')][50], label = 'base')
plt.plot(lat, zonal_ints[('stoc','all')][50], label = 'stoc')
plt.legend()

plt.figure()
plt.title('present day fluxes')
for varna, flun in zip(varnames, fluxnames):
    plt.plot(lat, zonal_ints[('base',varna)][50], label = flun)
plt.legend()

plt.figure()
plt.title('future - base vs stoc')
plt.plot(lat, zonal_ints[('base','all')][-10], label = 'base')
plt.plot(lat, zonal_ints[('stoc','all')][-10], label = 'stoc')
plt.legend()

# zonal_mean?
# stampa anno per anno + plot per artico integrato e non
