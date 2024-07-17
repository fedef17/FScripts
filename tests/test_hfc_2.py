#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os

from matplotlib import pyplot as plt
import netCDF4 as nc
import pandas as pd
from numpy import linalg as LA
from scipy import stats
from scipy import interpolate as itrp
import itertools as itt

from datetime import datetime
import pickle

import climtools_lib as ctl
import climdiags as cd

from copy import deepcopy as cp

###############################################
### constants
cp = 1005.0 # specific enthalpy dry air - J kg-1 K-1
cpw = 1840.0 # specific enthalpy water vapor - J kg-1 K-1
# per la moist air sarebbe cpm = cp + q*cpw, ma conta al massimo per un 3 %
L = 2501000.0 # J kg-1
Lsub = 2835000.0 # J kg-1
g = 9.81 # m s-2
Rearth = 6371.0e3 # mean radius
#####################################

L = 2501000.0
Rearth = 6371.0e3
##

cart_in = '/data-hobbes/fabiano/OBS/ERA/ERAInterim/'
file_list = cart_in+'ERAInterim_6hrs_1988_vatazgq.nc'
file_ps = cart_in+'ERAInterim_mon_1988_ps.nc'
cart_out = '/home/fabiano/Research/lavori/SPHINX_for_lisboa/heat_flux/ERA_ref_6hrs/'

# Loading ERA reference
era_fi = 'prova_heatflux_1988_MM.nc'

vars, lat, lon, dates, time_units, var_units, time_cal = ctl.readxDncfield(cart_in+era_fi)
print(vars.keys())

fi = 'northward_water_flux.nc'
var2, lat, lon, dates, time_units, var_units, time_cal = ctl.readxDncfield(cart_in+fi)
print(var2.keys())
vars[list(var2.keys())[0]] = L*var2[list(var2.keys())[0]]

era_zonal_factor = 2*np.pi*Rearth*np.cos(np.deg2rad(lat))

era_fluxes_maps = dict()
era_fluxes_zonal = dict()

seasons = ['Feb','DJF', 'JJA']
fluxnames = ['tot', 'SH', 'PE', 'LH']
eraname = {'tot': 'p76.162', 'SH': 'p70.162', 'PE': 'p74.162', 'LH': 'p72.162'}

for flun in fluxnames:
    for seas in seasons:
        era_fluxes_maps[(flun, seas)] = np.mean(ctl.sel_season(vars[eraname[flun]], dates, seas, cut = False)[0], axis = 0)
    era_fluxes_maps[(flun, 'year')] = np.mean(vars[eraname[flun]], axis = 0)

for fu in era_fluxes_maps:
     era_fluxes_zonal[fu] = np.mean(era_fluxes_maps[fu], axis = 1)*era_zonal_factor

######
print('Is it the first level???\n')
tag = 'ERAwith1000'
file_list = cart_in + 'all_vtgq_1988_6hrs.nc'

factors = {'SH': cp/g, 'PE': 1., 'LH': L/g}
#factors['PE'] = 1./g

press0row, latdad, londsad, datespress, time_units, var_units = ctl.read3Dncfield(file_ps)
pressok = press0row[1]

flux_levels = dict()
v, level, lat, lon, dates, time_units, var_units, time_cal = ctl.readxDncfield(file_list, select_var = ['v'])

t, level, lat, lon, dates, time_units, var_units, time_cal = ctl.readxDncfield(file_list, select_var = ['t'])
mshf = factors['SH']*v['v']*t['t']
flux_levels['SH'] = np.mean(mshf, axis=0)
del t, mshf

z, level, lat, lon, dates, time_units, var_units, time_cal = ctl.readxDncfield(file_list, select_var = ['z'])
mpef = factors['PE']*v['v']*z['z']
flux_levels['PE'] = np.mean(mpef, axis=0)
del z, mpef

q, level, lat, lon, dates, time_units, var_units, time_cal = ctl.readxDncfield(file_list, select_var = ['q'])
mlhf = factors['LH']*v['v']*q['q']
flux_levels['LH'] = np.mean(mlhf, axis=0)
del q, mlhf

fluxnames = ['SH', 'PE', 'LH']
level = level *100.

results_recalc = dict()
for flun in fluxnames:
    results_recalc[flun] = cd.quant_flux_calc(np.stack([flux_levels[flun], flux_levels[flun]]), lat, lon, level, None, np.stack([pressok,pressok]), None, quantity = None, seasons = [])['zonal']['year']

out = pickle.load(open('/home/fabiano/Research/lavori/SPHINX_for_lisboa/heat_flux/out_hfc_testERA_meanFeb.p'))

results_testERA = dict()
for cos, flun in zip(out, fluxnames):
    results_testERA[flun] = cd.quant_flux_calc(np.stack([cos, cos]), lat, lon, level, None, np.stack([pressok,pressok]), None, quantity = None, seasons = [])['zonal']['year']

figures = []
figure_file = cart_out + 'recalc_vs_testERA_Feb.pdf'
for flun in fluxnames:
    fig = plt.figure()
    plt.title('{} fluxes - ERAwith1000 vs ERA-sphinxlev'.format(flun))
    #cset = ctl.color_set(len(seasons), bright_thres = 1)
    #for seas, col in zip(seasons, cset):
    seas = 'Feb'
    plt.plot(lat, results_recalc[flun], label = 'recalc '+seas, linewidth = 1.5)
    #plt.plot(lat, results_testERA[flun], label = 'testERA '+seas, linewidth = 1.5)
    plt.plot(lat, era_fluxes_zonal[(flun, seas)], label = 'ref '+seas, linewidth = 0.7, linestyle = '--')
    plt.legend()
    plt.grid()
    #plt.ylim(mrgegegeg)
    plt.xlabel('Latitude')
    plt.ylabel('Integrated Net Heat Flux (W)')
    figures.append(fig)

ctl.plot_pdfpages(figure_file, figures)
