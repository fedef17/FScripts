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
from scipy import interpolate as itrp
import itertools as itt

from sklearn.cluster import KMeans



from datetime import datetime
import pickle

import climtools_lib as ctl
import climdiags as cd

from copy import deepcopy as cp

import logging
logging.basicConfig()

##############################################
KtoC = 273.15
cp_air = 1005.0 # specific enthalpy dry air - J kg-1 K-1

cart_in = '/data-hobbes/fabiano/SPHINX/tas_pr_mon/'
cart_in_3d = '/data-hobbes/fabiano/SPHINX/va_ua_ta_mon/'

cart_out = '/home/fabiano/Research/lavori/SPHINX_for_lisboa/mean_climates/'

varss3d = ['va', 'ta', 'ua']

ensmem = ['lcb0','lcb1','lcb2','lcs0','lcs1','lcs2']

seasons = ['DJF', 'JJA']

#Reading lat, lon and level
filena = '{}_mon_{}_{}.nc'.format('lcb0', 1988, 'ta')
varuna, level, lat, lon, datesuna, time_units, var_units, time_cal = ctl.read4Dncfield(cart_in_3d+filena)

years = np.arange(1950,2101)

cblabels = dict()
cblabels['ta'] = 'Temp (K)'
cblabels['ua'] = 'u (m/s)'
cblabels['va'] = 'v (m/s)'

carttemp = '/data-hobbes/fabiano/SPHINX/tas_mon/'
globalme, zonalme = pickle.load(open(carttemp+'global_tasmean_yearly.p', 'rb'))

for key in globalme:
    globalme[key] = globalme[key][100:]

for varna in ['ta', 'ua', 'va']:
    filos = cart_out+'out_mean3d{}_SPHINX.p'.format(varna)
    if os.path.exists(filos):
        cross3d = pickle.load(open(filos,'rb'))
    else:
        cross3d = dict()
        for ens in ensmem:
            for seas in seasons+['year']:
                cross3d[(ens, seas)] = []

            for yea in years:
                filena = '{}_mon_{}_{}.nc'.format(ens, yea, varna)
                var, level, lat, lon, dates, time_units, var_units, time_cal = ctl.read4Dncfield(cart_in_3d+filena)

                cli, datescli, _ = ctl.monthly_climatology(var, dates)
                for seas in seasons:
                    coso = np.mean(ctl.sel_season(cli, datescli, seas, cut = False)[0], axis = 0)
                    cross3d[(ens, seas)].append(np.mean(coso, axis = -1))
                coso = np.mean(cli, axis = 0)
                cross3d[(ens, 'year')].append(np.mean(coso, axis = -1))

            for seas in seasons+['year']:
                cross3d[(ens, seas)] = np.stack(cross3d[(ens, seas)])

        for seas in seasons+['year']:
            cross3d[('base', seas)] = np.mean([cross3d[(ens, seas)] for ens in ensmem[:3]], axis = 0)
            cross3d[('stoc', seas)] = np.mean([cross3d[(ens, seas)] for ens in ensmem[3:]], axis = 0)

        pickle.dump(cross3d, open(filos,'wb'))

    for key in cross3d:
        cross3d[key] = cross3d[key]-KtoC

    figure_file = cart_out+'{}_cross_bvss_greg.pdf'.format(varna)
    figures = []
    figures_zonal = []

    temp_ranges = np.arange(13.5, 18.1, 0.5)
    temp_base = globalme[('base', 'tas')]
    temp_stoc = globalme[('stoc', 'tas')]

    print('piop')
    mino = np.percentile([cross3d[('stoc', seas)]-cross3d[('base', seas)] for seas in seasons+['year']], 1)
    maxo = np.percentile([cross3d[('stoc', seas)]-cross3d[('base', seas)] for seas in seasons+['year']], 99)
    mimax = np.max([abs(mino), abs(maxo)])
    for seas in ['year']+seasons:
        axes = []
        for t1,t2 in zip(temp_ranges[:-1], temp_ranges[1:]):
            print(t1,t2)

            okann_base = (temp_base >=  t1) & (temp_base < t2)
            okann_stoc = (temp_stoc >=  t1) & (temp_stoc < t2)
            stoc = np.mean(cross3d[('stoc', seas)][okann_stoc, ...], axis = 0)
            base = np.mean(cross3d[('base', seas)][okann_base, ...], axis = 0)

            fig, ax = ctl.plot_lat_crosssection(stoc-base, lat, level, title = '3d {} diff stoc-base {}: at global T {} to {} C'.format(varna, seas, t1, t2), cbar_range = (-mimax, mimax), cb_label = cblabels[varna], set_logscale_levels = True, return_ax = True)

            figures.append(fig)

            stoczon = cp_air*np.trapz(stoc, x=level, axis=0)
            basezon = cp_air*np.trapz(base, x=level, axis=0)

            fig = plt.figure()
            ax = plt.subplot(1,1,1)
            plt.title('Total zonal energy - between {} and {} C'.format(t1, t2))

            ax.plot(lat, basezon, label = 'base', linewidth = 2.0)
            ax.plot(lat, stoczon, label = 'stoc', linewidth = 2.0)
            ax.plot(lat, 10*(stoczon-basezon), label = '(diff s-b) x 10', linewidth = 1.0, linestyle = '--')
            plt.xlabel('Latitude')
            plt.ylabel('Thermal energy (J)')
            plt.legend()
            plt.grid()
            figures_zonal.append(fig)
            axes.append(ax)

        ctl.adjust_ax_scale(axes, sel_axis = 'both')

    ctl.plot_pdfpages(figure_file, figures+figures_zonal)
    plt.close('all')
