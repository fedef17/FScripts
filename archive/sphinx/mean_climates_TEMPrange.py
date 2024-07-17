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
tommday = 86400

cart_in = '/data-hobbes/fabiano/SPHINX/tas_mon/'
cart_out = '/home/fabiano/Research/lavori/SPHINX_for_lisboa/mean_climates_sameT/'
if not os.path.exists(cart_out): os.mkdir(cart_out)
#varss = ['tas', 'pr']
varss = ['tas']

cblabels = dict()
cblabels['tas'] = 'Temp (C)'
cblabels['pr'] = 'pr (mm/day)'

ensmem = ['lcb0','lcb1','lcb2','lcs0','lcs1','lcs2']

# climatzon = dict()
# globalme = dict()
#
# for ens in ensmem:
#     # carico MM di tutti gli ensmems
#     for varna in varss:
#         filena = '{}-1850-2100-{}_mon.nc'.format(ens, varna)
#
#         var, lat, lon, dates, time_units, var_units = ctl.read3Dncfield(cart_in+filena)
#         if varna == 'tas':
#             var = var-KtoC
#         elif varna == 'pr':
#             var = var*tommday
#         dates_pdh = pd.to_datetime(dates)
#
#         varzon = ctl.zonal_mean(var)
#         climatzon[(ens, varna)] = ctl.zonal_mean(var)
#
#         # Global stuff
#         varye, _ = ctl.yearly_average(var, dates)
#         global_mean = ctl.global_mean(varye, lat)
#         globalme[(ens, varna)] = global_mean
#
# for varna in varss:
#     globalme[('base', varna)] = np.mean([globalme[(ens, varna)] for ens in ensmem[:3]], axis = 0)
#     globalme[('stoc', varna)] = np.mean([globalme[(ens, varna)] for ens in ensmem[3:]], axis = 0)
#     climatzon[('base', varna)] = np.mean([climatzon[(ens, varna)] for ens in ensmem[:3]], axis = 0)
#     climatzon[('stoc', varna)] = np.mean([climatzon[(ens, varna)] for ens in ensmem[3:]], axis = 0)
#
temp_horiz = [(13.5, 14.5), (14.5, 15.5), (15.5, 16.5), (16.5, 17.5)]
lat_bands = [(-90,-70), (-60,-40), (-20,20), (40,60), (70,90)]
#
# band_clim = dict()
# band_clim_std = dict()
# for exp in ensmem:
#     for varna in varss:
#         for th in temp_horiz:
#             yok = (globalme[(exp, 'tas')] >= th[0]) & (globalme[(exp, 'tas')] <= th[1])
#             yokmon = []
#             for co in yok:
#                 yokmon += 12*[co]
#             okvars = climatzon[(exp, varna)][yokmon, ...]
#             okdates = dates[yokmon]
#             climzonmean, _, climzonstd = ctl.monthly_climatology(okvars, okdates, refyear = 2001, dates_range = None)
#             for band in lat_bands:
#                 band_clim[(exp, varna, th, band)] = ctl.band_mean_from_zonal(climzonmean, lat, band[0], band[1])
#                 band_clim_std[(exp, varna, th, band)] = ctl.band_mean_from_zonal(climzonstd, lat, band[0], band[1])
#
# for varna in varss:
#     for th in temp_horiz:
#         for band in lat_bands:
#             band_clim[('base', varna, th, band)] = np.mean([band_clim[(ens, varna, th, band)] for ens in ensmem[:3]], axis = 0)
#             band_clim[('stoc', varna, th, band)] = np.mean([band_clim[(ens, varna, th, band)] for ens in ensmem[3:]], axis = 0)
#             band_clim_std[('base', varna, th, band)] = np.mean([band_clim[(ens, varna, th, band)] for ens in ensmem[:3]], axis = 0)
#             band_clim_std[('stoc', varna, th, band)] = np.mean([band_clim[(ens, varna, th, band)] for ens in ensmem[3:]], axis = 0)
#
# pickle.dump([globalme, band_clim], open(cart_out+'out_seascycle_SPHINX.p','wb'))
globalme, band_clim = pickle.load(open(cart_out+'out_seascycle_SPHINX.p','rb'))


# ORA i plots

names = ['SP', 'SML', 'TROP', 'NML', 'NP']
figures = []
for bandnam, band in zip(names, lat_bands):
    axes = []
    for coso in ['base', 'stoc']:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('{} - {}'.format(coso, bandnam))
        for th in temp_horiz:
            ax.plot(np.arange(12), band_clim[(coso, 'tas', th, band)], label = th)
        ax.legend()
        figures.append(fig)
        axes.append(ax)

    ctl.adjust_ax_scale(axes, sel_axis = 'both')

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('Diff stoc-base - {}'.format(bandnam))
    for th in temp_horiz:
        ax.plot(np.arange(12), band_clim[('stoc', 'tas', th, band)] - band_clim[('base', 'tas', th, band)], label = th)
    ax.legend()
    ax.grid()
    figures.append(fig)

ctl.plot_pdfpages(cart_out+'seas_cycle_all.pdf', figures)
