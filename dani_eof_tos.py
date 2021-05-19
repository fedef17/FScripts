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
pino11 = pino.sel(time = pino['time.month'] == 11)

solver = ctl.eof_computation(pino11.tos.values, latitude=pino.lat.values)

filout = cart_out + 'tos_eofs_opa0_withtrend.pdf'
ctl.plot_multimap_contour(solver.eofs(eofscaling=2)[:12], pino.lat.values, pino.lon.values, filout, plot_anomalies=True, cbar_range=(-1,1), subtitles= ['eof {}'.format(i) for i in range(12)], cb_label='T (K)')

### detrended
pinko, coeffs, var_reg, dats = ctl.remove_global_polytrend(pino11.lat.values, pino11.lon.values, pino11.tos.values, pino11.time.values, None, deg = 1)

plt.ion()
plt.figure()
plt.plot(var_reg)
solver_dtr = ctl.eof_computation(pinko, latitude=pino.lat.values)

filout2 = cart_out + 'tos_eofs_opa0_detrended.pdf'
ctl.plot_multimap_contour(solver_dtr.eofs(eofscaling=2)[:12], pino.lat.values, pino.lon.values, filout2, plot_anomalies=True, cbar_range=(-1,1), subtitles= ['eof {}'.format(i) for i in range(12)], cb_label='T (K)')

# match
eof2 = ctl.match_pc_sets(solver.eofs(eofscaling=2)[:12], solver_dtr.eofs(eofscaling=2)[:12], latitude = lat)

filout3 = cart + 'tos_eofs_opa0_diff_dtr-wtr.pdf'
ctl.plot_multimap_contour(eof2-solver.eofs(eofscaling=2)[:12], pino.lat.values, pino.lon.values, filout3, plot_anomalies=True, cbar_range=(-0.2,0.2), subtitles= ['eof {}'.format(i) for i in range(12)], cb_label='T (K)')
