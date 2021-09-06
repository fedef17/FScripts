#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
from matplotlib import pyplot as plt
from matplotlib import cm

import pickle
import netCDF4 as nc

import climtools_lib as ctl
import climdiags as cd
#from tunlib import gregplot_on_ax

from matplotlib.colors import LogNorm
from datetime import datetime

from scipy import stats
import xarray as xr
import glob
import xclim

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 22
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['axes.axisbelow'] = True

#############################################################################

if os.uname()[1] == 'hobbes':
    cart_out = '/home/fabiano/Research/lavori/BOTTINO/extreme_risk/'
elif os.uname()[1] == 'xaru':
    cart_out = '/home/fedef/Research/lavori/BOTTINO/extreme_risk/'
elif os.uname()[1] == 'tintin':
    cart_out = '/home/fabiano/work/lavori/BOTTINO/extreme_risk/'

ctl.mkdir(cart_out)

filna = '/nas/BOTTINO/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/{}/{}i1p1f1/{}/{}/*nc'

allru = ['pi', 'b025', 'b050', 'b100']
allnams = ['piControl', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']

colors = ['black', 'forestgreen', 'orange', 'violet']

####################################################################################################

allvars = ['tas', 'pr']

cart_in = '/home/fedef/Research/lavori/BOTTINO/yearmean/'
yeamean = pickle.load(open(cart_in + 'bottino_yeastat_tas.p', 'rb'))
yeamean_pr = pickle.load(open(cart_in + 'bottino_yeastat_pr.p', 'rb'))
yeamean.update(yeamean_pr)

resmap = dict()

figs = []
for ru in allru[1:]:
    ##### Risk of yearly drought
    piperc = np.percentile(yeamean[('pi', 'pr', 'sum')], 5, axis = 0)
    cos = '_prmin'

    coso = yeamean[(ru, 'pr', 'sum')]
    coso = coso.isel(year = slice(300, 500))
    tyt5 = np.sum(coso < piperc, axis = 0)/10.
    fig = ctl.plot_map_contour(tyt5, cmap='viridis', cbar_range=(1,9), n_color_levels=5, title = 'rr_'+ru+cos)
    figs.append(fig)

    resmap[(ru, cos)] = tyt5

    ##### Risk of monthly exceptional precipitation
    cos = '_prmax'
    piperc = np.percentile(yeamean[('pi', 'pr', 'max')], 95, axis = 0)
    tyt95 = np.sum(yeamean[(ru, 'pr', 'max')].isel(year = slice(300, 500)) > piperc, axis = 0)/10.
    fig = ctl.plot_map_contour(tyt95, cmap='viridis', cbar_range=(1,21), n_color_levels=5, title = 'rr_'+ru+cos)
    figs.append(fig)

    resmap[(ru, cos)] = tyt95

    ##### Risk of monthly exceptional temperature
    pime = np.mean(yeamean[('pi', 'tas', 'mean')], axis = 0)
    rume = np.mean(yeamean[(ru, 'tas', 'mean')].isel(year = slice(300, 500)), axis = 0)

    cos = '_tasmax'
    piperc = np.percentile(yeamean[('pi', 'tas', 'max')], 95, axis = 0)
    kazu = ((yeamean[(ru, 'tas', 'max')].isel(year = slice(300, 500)) - rume).values > (piperc - pime).values)
    tet95 = np.sum(kazu, axis = 0)/10.
    # tet95 = np.sum((yeamean[(ru, 'tas', 'max')].isel(year = slice(300, 500)) - rume > piperc - pime), axis = 0)/10.
    fig = ctl.plot_map_contour(tet95, pime.lat.values, pime.lon.values, cmap='viridis', cbar_range=(1,21), n_color_levels=5, title = 'rr_'+ru+cos)
    figs.append(fig)

    resmap[(ru, cos)] = tet95

    ##### Risk of monthly exceptional temperature
    cos = '_tasmin'
    piperc = np.percentile(yeamean[('pi', 'tas', 'min')], 5, axis = 0)
    kazu = ((yeamean[(ru, 'tas', 'min')].isel(year = slice(300, 500)) - rume).values < (piperc - pime).values)
    tet5 = np.sum(kazu, axis = 0)/10.
    #tet5 = np.sum(yeamean[(ru, 'tas', 'min')].isel(year = slice(300, 500)) - rume < piperc - pime, axis = 0)/10.
    fig = ctl.plot_map_contour(tet5, pime.lat.values, pime.lon.values, cmap='viridis', cbar_range=(0., 2.), n_color_levels=5, title = 'rr_'+ru+cos)
    figs.append(fig)

    resmap[(ru, cos)] = tet5

allcose = ['_prmin', '_prmax', '_tasmax', '_tasmin']
fnams = ['rr_'+ru+cos for ru in allru for cos in allcose]

allcbran = [(1,9), (1,21), (1,21), (0,2)]

ctl.plot_pdfpages(cart_out + 'risk_ratio_all.pdf', figs, save_single_figs = True, fig_names = fnams)

for cos, cbran in zip(allcose, allcbran):
    mape = [resmap[(ru, cos)] for ru in allru[1:]]
    fig = ctl.plot_multimap_contour(mape, pime.lat.values, pime.lon.values, filename = cart_out + 'all_RR{}.pdf'.format(cos), fix_subplots_shape = (1, 3), figsize = (24, 8), cmap='viridis', cbar_range = cbran, n_color_levels = 5, subtitles = allru[1:])
