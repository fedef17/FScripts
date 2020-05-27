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

#############################

cart_in = '/nas/reference/ERAInterim/daily/'
filnam = cart_in + '{}/{}_Aday_ERAInterim_201801_201812.nc'

cart_out = '/home/fabiano/Research/lavori/dynamics_exe/'

# 29 Dec 2018
u, ucoords, _ = ctl.readxDncfield(filnam.format('u', 'u'))
v, vcoords, _ = ctl.readxDncfield(filnam.format('v', 'v'))
ta, tacoords, _ = ctl.readxDncfield(filnam.format('ta', 'ta'))
zg, zgcoords, _ = ctl.readxDncfield(filnam.format('zg', 'zg'))
zg = 9.80665*zg
pv, pvcoords, _ = ctl.readxDncfield(filnam.format('pv', 'pv'))

#u.shape (365, 5, 241, 480)
lat_orig = ucoords['lat']
lon_orig = ucoords['lon']

levs = np.array([850, 500, 200, 50, 30])
area = (-180, 180, 10, 85)

#for da in np.arange(365):
da = -3
u_, lat, lon = ctl.sel_area(lat_orig, lon_orig, u[da, ::-1, ...], area)
v_, lat, lon = ctl.sel_area(lat_orig, lon_orig, v[da, ::-1, ...], area)
ta_, lat, lon = ctl.sel_area(lat_orig, lon_orig, ta[da, ::-1, ...], area)
zg_, lat, lon = ctl.sel_area(lat_orig, lon_orig, zg[da, ::-1, ...], area)
pv_, lat, lon = ctl.sel_area(lat_orig, lon_orig, pv[da, ::-1, ...], area)

Om = 2*np.pi/86400.
f = 2*Om*np.sin(np.deg2rad(lat))
f[f == 0.0] = 1.e-6
f = f[:, np.newaxis]

#f = np.reshape(np.repeat(f, len(lon)), zg[0,0].shape)
grad_zg = dict()
ug = dict()
vg = dict()
vortg = dict()
grad_ta = dict()


for i, lev in enumerate(levs):
    grad_zg[lev] = ctl.calc_gradient_2d(zg_[i], lat, lon)
    ug[lev] = -1/f * grad_zg[lev][1]
    vg[lev] = 1/f * grad_zg[lev][0]

    laplzg = ctl.calc_gradient_2d(grad_zg[lev][0], lat, lon)[0] + ctl.calc_gradient_2d(grad_zg[lev][1], lat, lon)[1]
    vortg[lev] = 1/f * laplzg
    vortg[lev][np.isinf(vortg[lev])] = np.nan

    grad_ta[lev] = ctl.calc_gradient_2d(ta_[i], lat, lon)

R = 8.31
ut = np.sum([R/f * np.log(lev1/lev2) * grad_ta[lev1][1] for lev1, lev2 in zip(levs[:-1], levs[1:])], axis = 0)
vt = np.sum([-R/f * np.log(lev1/lev2) * grad_ta[lev1][0] for lev1, lev2 in zip(levs[:-1], levs[1:])], axis = 0)

tam = np.mean(ta_, axis = 0)

quiver_scale = 500
vec_every = 10

figs = []
for i, lev in enumerate(levs):
    fig = ctl.plot_map_contour(zg_[i], lat, lon, add_vector_field = [u_[i], v_[i]], title = 'real winds - lev {} hPa'.format(lev), plot_anomalies = False, plot_margins = area, quiver_scale = quiver_scale, vec_every = vec_every)
    figs.append(fig)
    fig = ctl.plot_map_contour(zg_[i], lat, lon, add_vector_field = [ug[lev], vg[lev]], title = 'geostrophic winds - lev {} hPa'.format(lev), plot_anomalies = False, plot_margins = area, quiver_scale = quiver_scale, vec_every = vec_every)
    figs.append(fig)
    fig = ctl.plot_map_contour(zg_[i], lat, lon, add_vector_field = [u_[i] - ug[lev], v_[i] - vg[lev]], title = 'ageostrophic winds - lev {} hPa'.format(lev), plot_anomalies = False, plot_margins = area, quiver_scale = quiver_scale, vec_every = vec_every)
    figs.append(fig)
    fig = ctl.plot_map_contour(vortg[lev], lat, lon, add_contour_field = pv_[i], title = 'geostrophic vort and PV - lev {} hPa'.format(lev), plot_anomalies = False, plot_margins = area, vec_every = vec_every, plot_type = 'pcolormesh')
    figs.append(fig)

fig = ctl.plot_map_contour(tam, lat, lon, add_vector_field = [ut, vt], title = 'temperature and thermal wind', plot_anomalies = False, plot_margins = area, quiver_scale = quiver_scale, vec_every = vec_every)
figs.append(fig)

fig = ctl.plot_map_contour(tam, lat, lon, add_vector_field = [u_[-1]-ug[levs[0]]-ut, v_[-1]-vg[levs[0]]-vt], title = 'wind at height minus thermal wind', plot_anomalies = False, plot_margins = area, quiver_scale = quiver_scale, vec_every = vec_every)
figs.append(fig)

ctl.plot_pdfpages(cart_out + 'geostrophy_exe.pdf', figs)
