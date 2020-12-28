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

from matplotlib.colors import LogNorm
from datetime import datetime

from scipy import stats
import pandas as pd

plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 28
plt.rcParams['axes.labelsize'] = 22

#############################################################################

yr10 = 10 # length of running mean
#dtrtyp = 'light'
dtrtyp = 'histrebase'
cart_out_orig = '/home/fabiano/Research/lavori/CMIP6/Results_v6_eceens/'
ctl.mkdir(cart_out_orig)

cart_in = '/data-hobbes/fabiano/WR_CMIP6/'
#file_hist = cart_in + 'out_NEW_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refCLUS_dtr_light.p'
file_hist = cart_in + 'out_eceens_hist_NDJFM_{}_4clus_4pcs_1964-2014_refCLUS_dtr.p'
file_hist_refEOF = cart_in + 'out_eceens_hist_NDJFM_{}_4clus_4pcs_1964-2014_refEOF_dtr.p'
gen_file_ssp = cart_in + 'out_eceens_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr_reb.p'
gen_file_ssp_noreb = cart_in + 'out_eceens_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr.p'

area = 'EAT'
ssp = 'ssp585'

reshist, resref = ctl.load_wrtool(file_hist.format(area))
#reshist_re, _ = ctl.load_wrtool(file_hist.format(area))
resssp, _ = ctl.load_wrtool(gen_file_ssp.format(ssp, area))
resssp_noreb, _ = ctl.load_wrtool(gen_file_ssp_noreb.format(ssp, area))

histbases = []
for mod in reshist.keys():
    histbases.append(np.mean(reshist[mod]['climate_mean'][:, 50:70, -8], axis = 0))

lat = reshist[mod]['lat'][50:70]
lon = reshist[mod]['lon'][50:70]

sspbases = []
for mod in resssp_noreb.keys():
    sspbases.append(np.mean(resssp_noreb[mod]['climate_mean'][:, 50:70, -8], axis = 0))

fig = plt.figure()
for ci in histbases:
    plt.plot(lat, ci, color = 'black', linestyle = ':', linewidth = 0.2)
plt.plot(lat, np.mean(histbases, axis = 0), color = 'black', linewidth = 2, label = 'hist')
for ci in sspbases:
    plt.plot(lat, ci, color = 'red', linestyle = ':', linewidth = 0.2)
plt.plot(lat, np.mean(sspbases, axis = 0), color = 'red', linewidth = 2, label = 'ssp585')
plt.legend()
fig.savefig(cart_out_orig + 'check_climatemean_vlat.pdf')

histbases = []
for mod in reshist.keys():
    histbases.append(np.mean(reshist[mod]['climate_mean'][:, 50:70, -8], axis = 1))

dates = reshist[mod]['climate_mean_dates']

sspbases = []
for mod in resssp_noreb.keys():
    sspbases.append(np.mean(resssp_noreb[mod]['climate_mean'][:, 50:70, -8], axis = 1))


fig = plt.figure()
for ci in histbases:
    plt.plot_date(dates, ci, color = 'black', linestyle = ':', linewidth = 0.2)
plt.plot_date(dates, np.mean(histbases, axis = 0), color = 'black', linewidth = 2, label = 'hist')
for ci in sspbases:
    plt.plot_date(dates, ci, color = 'red', linestyle = ':', linewidth = 0.2)
plt.plot_date(dates, np.mean(sspbases, axis = 0), color = 'red', linewidth = 2, label = 'ssp585')
plt.legend()
fig.savefig(cart_out_orig + 'check_climatemean_vtime.pdf')


histbases = []
for mod in reshist.keys():
    histbases.append(np.mean(reshist[mod]['climate_mean'], axis = 0))
histcoso = np.mean(histbases, axis = 0)

lat = reshist[mod]['lat']
lon = reshist[mod]['lon']

sspbases = []
for mod in resssp_noreb.keys():
    sspbases.append(np.mean(resssp_noreb[mod]['climate_mean'], axis = 0))
sspcoso = np.mean(sspbases, axis = 0)

ctl.plot_double_sidebyside(sspcoso, histcoso, lat, lon, filename = cart_out_orig + 'mapdiff_rebase.pdf', visualization = 'standard', central_lat_lon = None, cmap = 'RdBu_r', title = None, xlabel = None, ylabel = None, cb_label = None, stitle_1 = 'ssp585', stitle_2 = 'hist', cbar_range = None, plot_anomalies = True, n_color_levels = 21, draw_contour_lines = False, n_lines = 5, color_percentiles = (0,100), use_different_grids = False, bounding_lat = 30, plot_margins = 'EAT', add_rectangles = None, draw_grid = True, plot_type = 'filled_contour', verbose = False, lw_contour = 0.5)

ctl.plot_map_contour(sspcoso-histcoso, lat, lon, filename = cart_out_orig + 'mapdiff_rebase_diff.pdf', visualization = 'standard', central_lat_lon = None, cmap = 'RdBu_r', title = None, xlabel = None, ylabel = None, cb_label = None, cbar_range = None, plot_anomalies = True, n_color_levels = 21, draw_contour_lines = False, n_lines = 5, color_percentiles = (0,100), bounding_lat = 30, plot_margins = 'EAT', add_rectangles = None, draw_grid = True, plot_type = 'filled_contour', verbose = False, lw_contour = 0.5)


fig = plt.figure()
noreb = ctl.running_mean(resssp_noreb['EC-Earth3_r1i1p1f1']['pcs'][:, 0], 10)
reb = ctl.running_mean(resssp['EC-Earth3_r1i1p1f1']['pcs'][:, 0], 10)
plt.plot(noreb, label = 'rebase_ssp')
plt.plot(reb, label = 'rebase_hist')
plt.plot(reb-noreb, label = 'diff')
fig.savefig(cart_out_orig + 'check_first_eof_vs_rebase.pdf')

fig = plt.figure()
noreb = ctl.running_mean(resssp_noreb['EC-Earth3_r1i1p1f1']['pcs'][:300, 0], 10)
reb = ctl.running_mean(resssp['EC-Earth3_r1i1p1f1']['pcs'][:300, 0], 10)
plt.plot(noreb, label = 'rebase_ssp')
plt.plot(reb, label = 'rebase_hist')
plt.plot(reb-noreb, label = 'diff')
fig.savefig(cart_out_orig + 'check_first_eof_vs_rebase_zoom.pdf')
