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
import xarray as xr
import glob
import xclim

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 28
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.axisbelow'] = True

#############################################################################
cart_out = '/home/fabiano/Research/lavori/dani_eofs_tos/'

#date = [datetime.strptime('{}15'.format(int(da)), '%Y%m%d') for da in tos.time.values]

#allcos = []
pino = xr.load_dataset('/data-archimede/ORAS4/tos_Omon_ORAS4_opa0_195709-201412_r360x180.nc', use_cftime = True)
pino = pino.drop_vars('latitude')
pino = pino.drop_vars('longitude')

# only november
pino = pino.sel(time = slice('1960-01-01', '2014-12-31'))
pino11 = pino.sel(time = pino['time.month'] == 11)

filexp = '/data-archimede/historical/ecearth/a1tn/tos/r360x180/tos_Omon_EC-Earth3_hist*nc'
cose = glob.glob(filexp)

gigi = xr.open_mfdataset(cose, use_cftime = True)
gigi = gigi.rename({'latitude' : 'lat'})
gigi = gigi.rename({'longitude' : 'lon'})
# gigi = gigi.drop_vars('latitude')
# gigi = gigi.drop_vars('longitude')

# only november
gigi = gigi.sel(time = slice('1960-01-01', '2014-12-31'))
gigi11 = gigi.sel(time = gigi['time.month'] == 11)

# common mask

mask_obs = np.isnan(pino11.tos.values[0])
mask_exp = np.isnan(gigi11.tos.values[0])
mask = (mask_obs) | (mask_exp)

pi11tos = pino11.tos.values[~mask[np.newaxis, ...]]
gi11tos = gigi11.tos.values[~mask[np.newaxis, ...]]
lat = pino.lat.values
lon = pino.lon.values

############################################################

solver = ctl.eof_computation(pi11tos, latitude=pino.lat.values)

filout = cart_out + 'tos_eofs_opa0_withtrend.pdf'
ctl.plot_multimap_contour(solver.eofs(eofscaling=2)[:12], lat, lon, filout, plot_anomalies=True, cbar_range=(-1,1), subtitles= ['eof {}'.format(i) for i in range(12)], cb_label='T (K)')

### detrended
pinko, coeffs, var_reg, dats = ctl.remove_global_polytrend(lat, lon, pi11tos, pino11.time.values, None, deg = 1)

plt.ion()
plt.figure()
plt.plot(var_reg)
solver_dtr = ctl.eof_computation(pinko, latitude=pino.lat.values)

filout2 = cart_out + 'tos_eofs_opa0_detrended.pdf'
ctl.plot_multimap_contour(solver_dtr.eofs(eofscaling=2)[:12], lat, lon, filout2, plot_anomalies=True, cbar_range=(-1,1), subtitles= ['eof {}'.format(i) for i in range(12)], cb_label='T (K)')

# match
#eof2 = ctl.match_pc_sets(solver.eofs(eofscaling=2)[:12], solver_dtr.eofs(eofscaling=2)[:12], latitude = lat)
#filout3 = cart + 'tos_eofs_opa0_diff_dtr-wtr.pdf'
#ctl.plot_multimap_contour(eof2-solver.eofs(eofscaling=2)[:12], pino.lat.values, pino.lon.values, filout3, plot_anomalies=True, cbar_range=(-0.2,0.2), subtitles= ['eof {}'.format(i) for i in range(12)], cb_label='T (K)')

#######################################

solver_exp = ctl.eof_computation(gi11tos, latitude=gigi.lat.values)

filout = cart_out + 'tos_eofs_atn1_withtrend.pdf'
ctl.plot_multimap_contour(solver_exp.eofs(eofscaling=2)[:12], lat, lon, filout, plot_anomalies=True, cbar_range=(-1,1), subtitles= ['eof {}'.format(i) for i in range(12)], cb_label='T (K)')

### detrended
ginko, coeffs, var_reg, dats = ctl.remove_global_polytrend(lat, lon, gi11tos, gigi11.time.values, None, deg = 1)

plt.ion()
plt.figure()
plt.plot(var_reg)
solver_exp_dtr = ctl.eof_computation(ginko, latitude=gigi.lat.values)

filout2 = cart_out + 'tos_eofs_atn1_detrended.pdf'
ctl.plot_multimap_contour(solver_exp_dtr.eofs(eofscaling=2)[:12], lat, lon, filout2, plot_anomalies=True, cbar_range=(-1,1), subtitles= ['eof {}'.format(i) for i in range(12)], cb_label='T (K)')

pcs_ref = []
pcs_ref_dtr = []
for i in range(12):
    pcs_ref.append(solver.projectField(solver_exp.eofs(eofscaling=2)[i], neofs=12, eofscaling=0, weighted=True))
    pcs_ref_dtr.append(solver_dtr.projectField(solver_exp_dtr.eofs(eofscaling=2)[i], neofs=12, eofscaling=0, weighted=True))
