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

zg_zon = np.mean(zg, axis = -1)
zg_eddy = zg - zg_zon[..., np.newaxis]

ua_climate_mean, dates_climate_mean, _ = ctl.daily_climatology(ua, coords['dates'], window = 20)
zg_climate_mean, dates_climate_mean, _ = ctl.daily_climatology(zg, coords['dates'], window = 20)
zg_climate_mean_eddy, dates_climate_mean, _ = ctl.daily_climatology(zg_eddy, coords['dates'], window = 20)
zg_climate_mean_zon, dates_climate_mean, _ = ctl.daily_climatology(zg_zon, coords['dates'], window = 20)

ua_low = ctl.butter_filter(ua, 10)
ua_low_djfm, dates_djfm = ctl.sel_season(ua_low, coords['dates'], 'DJFM')

zg_djfm, dates_djfm = ctl.sel_season(zg, coords['dates'], 'DJFM')
zg_eddy_djfm, dates_djfm = ctl.sel_season(zg_eddy, coords['dates'], 'DJFM')

# zg_low = ctl.butter_filter(zg, 10)
# zg_low_djfm, dates_djfm = ctl.sel_season(zg_low, coords['dates'], 'DJFM')
#
# zg_low_eddy = ctl.butter_filter(zg_eddy, 10)
# zg_low_eddy_djfm, dates_djfm = ctl.sel_season(zg_low_eddy, coords['dates'], 'DJFM')

### Test jli.
figs = []
# jli, jspeed, jdates = cd.jetlatindex(ua, coords['lat'], coords['lon'], coords['dates'], filter = 'butter')
# figs.append(cd.plot_jli_w_speed(jli, jspeed, jdates, title = 'butterworth w 0.22'))
# jli, jspeed, jdates = cd.jetlatindex(ua, coords['lat'], coords['lon'], coords['dates'], filter = 'butter')
# figs.append(cd.plot_jli_w_speed(jli, jspeed, jdates, title = 'butterworth w auto bnd', bnd_width = None))
# jli, jspeed, jdates = cd.jetlatindex(ua, coords['lat'], coords['lon'], coords['dates'], filter = 'butter', remove_orog = True)
# figs.append(cd.plot_jli_w_speed(jli, jspeed, jdates, title = 'butterworth no greenland'))
# jli, jspeed, jdates = cd.jetlatindex(ua, coords['lat'], coords['lon'], coords['dates'], filter = 'lanczos')
# figs.append(cd.plot_jli_w_speed(jli, jspeed, jdates, title = 'lanczos'))
# ctl.plot_pdfpages(cart_out + 'test_jli_filters.pdf', figs)

lat = coords['lat']
lon = coords['lon']

# res_all = dict()
#
# area = 'NATL'
# numpcs = 5
# numclus = 5
# perc = None
#
# res_all['ua_low_full_5'] = cd.WRtool_core(ua_low_djfm, lat, lon, dates_djfm, area, numpcs = numpcs, perc = perc, numclus = numclus, remove_climate_mean = False)
#
# res_all['ua_low_anom_5'] = cd.WRtool_core(ua_low_djfm, lat, lon, dates_djfm, area, climate_mean = ua_climate_mean, dates_climate_mean = dates_climate_mean, numpcs = numpcs, perc = perc, numclus = numclus, remove_climate_mean = True)
#
# numclus = 4
# res_all['ua_low_full_4'] = cd.WRtool_core(ua_low_djfm, lat, lon, dates_djfm, area, numpcs = numpcs, perc = perc, numclus = numclus, remove_climate_mean = False)
#
# res_all['ua_low_anom_4'] = cd.WRtool_core(ua_low_djfm, lat, lon, dates_djfm, area, climate_mean = ua_climate_mean, dates_climate_mean = dates_climate_mean, numpcs = numpcs, perc = perc, numclus = numclus, remove_climate_mean = True)
#
#
# area = 'EAT'
# numpcs = 4
# numclus = 4
#
# res_all['zg_eddy_full'] = cd.WRtool_core(zg_eddy_djfm, lat, lon, dates_djfm, area, numpcs = numpcs, perc = perc, numclus = numclus, remove_climate_mean = False)
#
# res_all['zg_eddy_anom'] =  cd.WRtool_core(zg_eddy_djfm, lat, lon, dates_djfm, area, climate_mean = zg_climate_mean_eddy, dates_climate_mean = dates_climate_mean, numpcs = numpcs, perc = perc, numclus = numclus, remove_climate_mean = True)
#
# res_all['zg_full'] = cd.WRtool_core(zg_djfm, lat, lon, dates_djfm, area, numpcs = numpcs, perc = perc, numclus = numclus, remove_climate_mean = False)
#
# res_all['zg_anom'] = cd.WRtool_core(zg_djfm, lat, lon, dates_djfm, area, climate_mean = zg_climate_mean, dates_climate_mean = dates_climate_mean, numpcs = numpcs, perc = perc, numclus = numclus, remove_climate_mean = True)
#
# res_all['zg_zon_anom'] = cd.WRtool_core(zg_djfm, lat, lon, dates_djfm, area, climate_mean = zg_climate_mean_zon[..., np.newaxis], dates_climate_mean = dates_climate_mean, numpcs = numpcs, perc = perc, numclus = numclus, remove_climate_mean = True)
#
# pickle.dump(res_all, open(cart_out + 'res_all_regs.p', 'wb'))
res_all = pickle.load(open(cart_out + 'res_all_regs.p', 'rb'))

figs = []
for ke in res_all:
    koze = res_all[ke]
    if 'ua' in ke:
        cbar_range = (-30, 30.)
        cb_label = 'u (m/s)'
    else:
        cbar_range = (-160, 160.)
        cb_label = 'zg (m)'

    fix_subplots_shape = (2, 2)
    if len(koze['cluspattern']) == 5: fix_subplots_shape = (2, 3)

    patts = koze['cluspattern']
    if ke == 'zg_full':
        patts = patts - np.mean(koze['cluspattern_area'])
    fig = cd.plot_regimes(koze['lat'], koze['lon'], patts, None, cbar_range = cbar_range, cb_label = cb_label, reg_freq = koze['freq_clus'])
    fig[0].suptitle(ke)
    figs.append(fig[0])

#figs = np.concatenate(figs)
ctl.plot_pdfpages(cart_out + 'test_regimes_pattern.pdf', figs)
