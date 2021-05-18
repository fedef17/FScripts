#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
import glob

import matplotlib
#matplotlib.use('Agg') # This is to avoid the code crash if no Xwindow is available
from matplotlib import pyplot as plt

import pickle
from scipy import io
import iris
import xarray as xr

import climtools_lib as ctl
import climdiags as cd
from copy import deepcopy as copy
import itertools as itt

########################################################

#zg = '/data-hobbes/fabiano/OBS/ERA/ERA40+Int_daily_1957-2018_zg500_remap25_meters.nc'
cart_out = '/home/fabiano/Research/lavori/WeatherRegimes/ERA5_tests/'

# zgfil = '/nas/reference/ERA5/daily/zg500/ERA5_zg500_1979-2019_meters.nc'
# ufils = '/nas/reference/ERA5/daily/u/ERA5_u_5lev_*nc'
#
# ufils_li = glob.glob(ufils)
# ufils_li.sort()
#
# ua, coords, auxi = ctl.read_xr(ufils_li, extract_level_hPa = 850., regrid_to_deg = 2.5)
# zg, coords_zg, auxi_zg = ctl.read_xr(zgfil, extract_level_hPa = 850., regrid_to_deg = 2.5)
# zg, dates = ctl.sel_time_range(zg, coords_zg['dates'], year_range = (1979,2018))
#
# pickle.dump([ua, zg, coords], open(cart_out + 'zg_e_ua_ERA5.p', 'wb'))
ua, zg, coords = pickle.load(open(cart_out + 'zg_e_ua_ERA5.p', 'rb'))


ua_climate_mean, ua_dates_climate_mean, _ = ctl.daily_climatology(ua, coords['dates'], window = 20)
zg_climate_mean, zg_dates_climate_mean, _ = ctl.daily_climatology(zg, coords['dates'], window = 20)

ua_low = ctl.butter_filter(ua, 10)
ua_low_djfm, dates_djfm = ctl.sel_season(ua_low, coords['dates'], 'DJFM')

### Test jli.
figs = []
jli, jspeed, jdates = cd.jetlatindex(ua, coords['lat'], coords['lon'], coords['dates'], filter = 'butter')
figs.append(cd.plot_jli_w_speed(jli, jspeed, jdates, title = 'butterworth w 0.22'))
jli, jspeed, jdates = cd.jetlatindex(ua, coords['lat'], coords['lon'], coords['dates'], filter = 'butter')
figs.append(cd.plot_jli_w_speed(jli, jspeed, jdates, title = 'butterworth w auto bnd'), bnd_width = None)
jli, jspeed, jdates = cd.jetlatindex(ua, coords['lat'], coords['lon'], coords['dates'], filter = 'butter_worog')
figs.append(cd.plot_jli_w_speed(jli, jspeed, jdates, title = 'butterworth no greenland'))
jli, jspeed, jdates = cd.jetlatindex(ua, coords['lat'], coords['lon'], coords['dates'], filter = 'lanczos')
figs.append(cd.plot_jli_w_speed(jli, jspeed, jdates, title = 'lanczos'))
ctl.plot_pdfpages(cart_out + 'test_jli_filters.pdf', figs)




#results = cd.WRtool_core(var_season, lat, lon, dates_season, area, climate_mean = climate_mean, dates_climate_mean = dates_climate_mean, numpcs = inputs['numpcs'], perc = inputs['perc'], numclus = inputs['numclus'])
