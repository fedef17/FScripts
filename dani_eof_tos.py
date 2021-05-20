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

n_ref = 8
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


# ok, qui rileggo tutte le obs e tutti gli exp
pi11tos = []
pi11tos_dtr = []
gi11tos = []
gi11tos_dtr = []

for opas in range(5):
    pino = xr.load_dataset('/data-archimede/ORAS4/tos_Omon_ORAS4_opa{}_195709-201412_r360x180.nc'.format(opas), use_cftime = True)
    pino = pino.drop_vars('latitude')
    pino = pino.drop_vars('longitude')

    # only november
    pino = pino.sel(time = slice('1960-01-01', '2014-12-31'))
    pino11 = pino.sel(time = pino['time.month'] == 11)

    lat = pino.lat.values
    lon = pino.lon.values

    picoso = pino11.tos.values
    picoso[:, mask] = np.nan

    pi11tos.append(picoso)

    pinko, coeffs, var_reg, dats = ctl.remove_global_polytrend(lat, lon, picoso, pino11.time.values, None, deg = 1)

    pi11tos_dtr.append(pinko)

pi11tos = np.concatenate(pi11tos, axis = 0)
pi11tos_dtr = np.concatenate(pi11tos_dtr, axis = 0)

expnams = os.listdir('/data-archimede/historical/ecearth/')
filexp = '/data-archimede/historical/ecearth/{}/tos/r360x180/tos_Omon_EC-Earth3_hist*nc'

for nam in expnams:
    cose = glob.glob(filexp.format(nam))

    gigi = xr.open_mfdataset(cose, use_cftime = True)
    gigi = gigi.rename({'latitude' : 'lat'})
    gigi = gigi.rename({'longitude' : 'lon'})
    # gigi = gigi.drop_vars('latitude')
    # gigi = gigi.drop_vars('longitude')

    # only november
    gigi = gigi.sel(time = slice('1960-01-01', '2014-12-31'))
    gigi11 = gigi.sel(time = gigi['time.month'] == 11)

    gicoso = gigi11.tos.values
    gicoso[:, mask] = np.nan
    gi11tos.append(gicoso)

    ginko, coeffs, var_reg, dats = ctl.remove_global_polytrend(lat, lon, gicoso, pino11.time.values, None, deg = 1)

    gi11tos_dtr.append(ginko)

gi11tos = np.concatenate(gi11tos, axis = 0)
gi11tos_dtr = np.concatenate(gi11tos_dtr, axis = 0)

pickle.dump([pi11tos, pi11tos_dtr, gi11tos, gi11tos_dtr], open(cart_out + 'tos_data.nc', 'wb'))
############################################################

solver = ctl.eof_computation(pi11tos, latitude=pino.lat.values)

filout = cart_out + 'tos_eofs_obs_withtrend.pdf'
ctl.plot_multimap_contour(solver.eofs(eofscaling=2)[:n_ref], lat, lon, filout, plot_anomalies=True, cbar_range=(-1,1), subtitles= ['eof {}'.format(i) for i in range(n_ref)], cb_label='T (K)')

### detrended
solver_dtr = ctl.eof_computation(pi11tos_dtr, latitude=pino.lat.values)

filout2 = cart_out + 'tos_eofs_obs_detrended.pdf'
ctl.plot_multimap_contour(solver_dtr.eofs(eofscaling=2)[:n_ref], lat, lon, filout2, plot_anomalies=True, cbar_range=(-1,1), subtitles= ['eof {}'.format(i) for i in range(n_ref)], cb_label='T (K)')

# match
#eof2 = ctl.match_pc_sets(solver.eofs(eofscaling=2)[:n_ref], solver_dtr.eofs(eofscaling=2)[:n_ref], latitude = lat)
#filout3 = cart + 'tos_eofs_opa0_diff_dtr-wtr.pdf'
#ctl.plot_multimap_contour(eof2-solver.eofs(eofscaling=2)[:n_ref], pino.lat.values, pino.lon.values, filout3, plot_anomalies=True, cbar_range=(-0.2,0.2), subtitles= ['eof {}'.format(i) for i in range(n_ref)], cb_label='T (K)')

#######################################

obseofs = solver.eofs(eofscaling=2)[:n_ref]
obseofs_dtr = solver_dtr.eofs(eofscaling=2)[:n_ref]

solver_exp = ctl.eof_computation(gi11tos, latitude=gigi.lat.values)

okmatch = ctl.match_patterns(obseofs, solver_exp.eofs(eofscaling=2)[:n_ref+10], latitude = lat)
#expeofs = solver_exp.eofs(eofscaling=2)[:n_ref+10][okmatch]
expeofs = solver_exp.eofs(eofscaling=2)[:n_ref]

filout = cart_out + 'tos_eofs_exp_withtrend.pdf'
ctl.plot_multimap_contour(expeofs, lat, lon, filout, plot_anomalies=True, cbar_range=(-1,1), subtitles= ['eof {}'.format(i) for i in range(n_ref)], cb_label='T (K)')

### detrended
solver_exp_dtr = ctl.eof_computation(gi11tos_dtr, latitude=gigi.lat.values)

okmatch_dtr = ctl.match_patterns(obseofs_dtr, solver_exp_dtr.eofs(eofscaling=2)[:n_ref+10], latitude = lat)
#expeofs_dtr = solver_exp_dtr.eofs(eofscaling=2)[:n_ref+10][okmatch_dtr]
expeofs_dtr = solver_exp_dtr.eofs(eofscaling=2)[:n_ref]

filout2 = cart_out + 'tos_eofs_exp_detrended.pdf'
ctl.plot_multimap_contour(expeofs_dtr, lat, lon, filout2, plot_anomalies=True, cbar_range=(-1,1), subtitles= ['eof {}'.format(i) for i in range(n_ref)], cb_label='T (K)')

#### matched diffs

expeofs = solver_exp.eofs(eofscaling=2)[:n_ref+10][okmatch]
expeofs_dtr = solver_exp_dtr.eofs(eofscaling=2)[:n_ref+10][okmatch_dtr]

signs = np.array([np.sign(ctl.Rcorr(ob, ex, latitude = lat)) for ob,ex in zip(obseofs, expeofs)])
filout3 = cart + 'tos_eofs_diff_obs-exp_withtrend.pdf'
ctl.plot_multimap_contour(expeofs-signs[:, np.newaxis, np.newaxis]*obseofs, pino.lat.values, pino.lon.values, filout3, plot_anomalies=True, cbar_range=(-0.2,0.2), subtitles= ['eof {}'.format(i) for i in range(n_ref)], cb_label='T (K)')

signs = np.array([np.sign(ctl.Rcorr(ob, ex, latitude = lat)) for ob,ex in zip(obseofs_dtr, expeofs_dtr)])
filout3 = cart + 'tos_eofs_diff_obs-exp_detrended.pdf'
ctl.plot_multimap_contour(expeofs_dtr-signs[:, np.newaxis, np.newaxis]*obseofs_dtr, pino.lat.values, pino.lon.values, filout3, plot_anomalies=True, cbar_range=(-0.2,0.2), subtitles= ['eof {}'.format(i) for i in range(n_ref)], cb_label='T (K)')

####################################################

fig = plt.figure()
plt.plot(np.cumsum(solver.varianceFraction()[:20]), label = 'obs eofs')
plt.plot(np.cumsum(solver_exp.varianceFraction()[:20]), label = 'exp eofs')
plt.plot(np.cumsum(solver_exp.varianceFraction()[okmatch]), color = 'orange', linestyle = '--')
plt.axhline(0.5, color = 'grey', linestyle = '--')
plt.title('Explained variance')
fig.savefig(cart_out + 'explained_variance.pdf')

fig = plt.figure()
plt.plot(np.cumsum(solver_dtr.varianceFraction()[:20]), label = 'obs eofs')
plt.plot(np.cumsum(solver_exp_dtr.varianceFraction()[:20]), label = 'exp eofs')
plt.plot(np.cumsum(solver_exp_dtr.varianceFraction()[okmatch_dtr]), color = 'orange', linestyle = '--')
plt.axhline(0.5, color = 'grey', linestyle = '--')
plt.title('Explained variance (dtr)')
fig.savefig(cart_out + 'explained_variance_dtr.pdf')
# np.argwhere(np.cumsum(solver_exp_dtr.varianceFraction()) > 0.5)
# np.argwhere(np.cumsum(solver_exp_dtr.varianceFraction()) > 0.5)[0][0]
# np.argwhere(np.cumsum(solver_dtr.varianceFraction()) > 0.5)[0][0]



# pcs_ref = []
# pcs_ref_dtr = []
# for i in range(n_ref):
#     pcs_ref.append(solver.projectField(solver_exp.eofs(eofscaling=2)[i], neofs=n_ref, eofscaling=0, weighted=True))
#     pcs_ref_dtr.append(solver_dtr.projectField(solver_exp_dtr.eofs(eofscaling=2)[i], neofs=n_ref, eofscaling=0, weighted=True))
