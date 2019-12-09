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

###############################################################################
cose = np.arange(-2, 2.01, 0.01)
lanczos = np.sinc(cose)*np.sinc(cose/2.)
lanc20 = lanczos[::20]
lanc20 = lanc20/np.sum(lanc20)

def dothing(wind, coords, lanc = lanc20):
    area = [-60., 0., 20., 70.]
    lat = coords['lat']
    lon = coords['lon']
    wind_area, latsel, lonsel = ctl.sel_area(lat, lon, wind, area)
    wind_low = np.zeros(wind_area.shape)
    for ila, la in enumerate(latsel):
        for ilo, lo in enumerate(lonsel):
            wind_low[:, ila, ilo] = np.convolve(lanc20, wind_area[:, ila, ilo], mode = 'same')
    #wind_low = ctl.running_mean(wind_area, 10)

    wind_low_djf, dates = ctl.sel_season(wind_low, coords['dates'], 'DJF')

    return wind_low_djf, dates


area = [-60., 0., 20., 70.]
cart_out = '/home/fabiano/Research/lavori/Jet_latitude/'

#orogfi = '/data-hobbes/reference/ERAInterim/geopot_vegcover_25.nc'
orogfi = '/data-hobbes/reference/ERAInterim/geopot_vegcover.nc'
orog, coords, aux_info = ctl.readxDncfield(orogfi, select_var = 'z')
#orog = orog/9.80665
orogmask = orog > 1300.0
orogarea, _, _ = ctl.sel_area(coords['lat'], coords['lon'], orogmask, area)
orogarea = orogarea[0]
#ctl.plot_map_contour(orogmask, plot_type = 'pcolormesh')

cart_in = '/nas/reference/ERA40/daily/u/'
file_in = 'u_Aday_ERA40_{}01_{}12.nc'
yea = 1957
file_in_57 = 'u_Aday_ERA40_195709_195712.nc'

winds = []
alldates = []

wind, coords, aux_info = ctl.readxDncfield(cart_in + file_in_57, extract_level = 850)
wind_lo, dates = dothing(wind, coords)
winds.append(wind_lo)
alldates.append(dates)

wind_area, latsel, lonsel = ctl.sel_area(coords['lat'], coords['lon'], wind, area)

for yea in range(1958, 1979):
    wind, coords, aux_info = ctl.readxDncfield(cart_in + file_in.format(yea, yea), extract_level = 850)
    wind_lo, dates = dothing(wind, coords)
    winds.append(wind_lo)
    alldates.append(dates)

cart_in = '/nas/reference/ERAInterim/daily/u/'
file_in = 'u_Aday_ERAInterim_{}01_{}12.nc'
for yea in range(1979, 2015):
    wind, coords, aux_info = ctl.readxDncfield(cart_in + file_in.format(yea, yea), extract_level = 850)
    wind_lo, dates = dothing(wind, coords)
    winds.append(wind_lo)
    alldates.append(dates)

winds = np.concatenate(winds, axis = 0)
dates = np.concatenate(alldates, axis = 0)

for co in range(winds.shape[0]):
    winds[co][orogarea] = np.nan

wind_zon = np.nanmean(winds, axis = 2)

jli_speed = np.max(wind_zon, axis = 1)
jli = latsel[np.argmax(wind_zon, axis = 1)]

pickle.dump([jli, jli_speed], open(cart_out + 'jli_ERA.p', 'wb'))

plt.ion()
pdf = ctl.calc_pdf(jli)
plt.plot(latsel, pdf(latsel))
