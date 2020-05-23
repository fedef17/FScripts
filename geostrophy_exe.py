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
u = u[-3, ...]
v, vcoords, _ = ctl.readxDncfield(filnam.format('v', 'v'))
v = v[-3, ...]
ta, tacoords, _ = ctl.readxDncfield(filnam.format('ta', 'ta'))
ta = ta[-3, ...]
zg, zgcoords, _ = ctl.readxDncfield(filnam.format('zg', 'zg'))
zg = 9.80665*zg[-3, ...]
pv, pvcoords, _ = ctl.readxDncfield(filnam.format('pv', 'pv'))
pv = pv[-3, ...]

#u.shape (365, 5, 241, 480)
lat = ucoords['lat']
lon = ucoords['lon']

Om = 2*np.pi/86400.
f = 2*Om*np.sin(np.deg2rad(lat))
f = np.reshape(np.repeat(f, len(lon)), zg[0].shape)

levs = np.array([850, 500, 200, 50, 30])

grad_zg = dict()
ug = dict()
vg = dict()
vortg = dict()
grad_ta = dict()

for i, lev in enumerate(levs):
    grad_zg[lev] = ctl.calc_gradient_2d(zg[i], lat, lon)
    ug[lev] = -1/f * grad_zg[lev][1]
    vg[lev] = 1/f * grad_zg[lev][0]

    laplzg = ctl.calc_gradient_2d(grad_zg[lev][0], lat, lon)[0] + ctl.calc_gradient_2d(grad_zg[lev][1], lat, lon)[1]
    vortg[lev] = 1/f * laplzg
    vortg[lev][np.isinf(vortg[lev])] = np.nan

    grad_ta[lev] = ctl.calc_gradient_2d(ta[i], lat, lon)

R = 8.31
ut = np.sum([ R/f * np.log(lev1/lev2) * grad_ta[lev][1] for lev1, lev2 in zip(levs[:-1], levs[1:])], axis = 0)
vt = np.sum([-R/f * np.log(lev1/lev2) * grad_ta[lev][0] for lev1, lev2 in zip(levs[:-1], levs[1:])], axis = 0)

tam = np.mean(ta, axis = 0)

figs = []
for i, lev in enumerate(levs):
    fig = ctl.plot_map_contour(zg[i], lat, lon, add_vector_field = [u[i], v[i]], title = 'real winds - lev {} hPa'.format(lev), plot_anomalies = False, plot_margins = (-100, 60, 20, 90), quiver_scale = 1000)
    figs.append(fig)
    fig = ctl.plot_map_contour(zg[i], lat, lon, add_vector_field = [ug[lev], vg[lev]], title = 'geostrophic winds - lev {} hPa'.format(lev), plot_anomalies = False, plot_margins = (-100, 60, 20, 90), quiver_scale = 1000)
    figs.append(fig)
    fig = ctl.plot_map_contour(zg[i], lat, lon, add_vector_field = [u[i] - ug[lev], v[i] - vg[lev]], title = 'ageostrophic winds - lev {} hPa'.format(lev), plot_anomalies = False, plot_margins = (-100, 60, 20, 90), quiver_scale = 1000)
    figs.append(fig)
    fig = ctl.plot_map_contour(vortg[lev], lat, lon, add_contour_field = pv[i], title = 'geostrophic vort and PV - lev {} hPa'.format(lev), plot_anomalies = False, plot_margins = (-100, 60, 20, 90))
    figs.append(fig)

fig = ctl.plot_map_contour(tam, lat, lon, add_vector_field = [ut, vt], title = 'temperature and thermal wind', plot_anomalies = False, plot_margins = (-100, 60, 20, 90), quiver_scale = 1000)
figs.append(fig)

fig = ctl.plot_map_contour(tam, lat, lon, add_vector_field = [u[-1]-ug[levs[0]]-ut, v[-1]-vg[levs[0]]-vt], title = 'wind at height minus thermal wind', plot_anomalies = False, plot_margins = (-100, 60, 20, 90), quiver_scale = 1000)
figs.append(fig)

ctl.plot_pdfpages(cart_out + 'geostrophy_exe.pdf', figs)
