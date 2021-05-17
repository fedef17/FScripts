#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
import glob

import matplotlib
matplotlib.use('Agg') # This is to avoid the code crash if no Xwindow is available
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

##### low pass filters for ua

from scipy import signal

def butter_lowpass(cutoff, fs = 1, order=5):
    """
    Creates the filter. fs is the sample frequency (defaults = 1).
    """
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs = 1, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = signal.lfilter(b, a, data)
    return y

# Filter requirements.
# Get the filter coefficients so we can check its frequency response.
b, a = butter_lowpass(0.1) # 10 days

# Plot the frequency response.
w, h = signal.freqz(b, a, worN = 1024)#8000)
plt.subplot(2, 1, 1)
plt.plot(0.5*fs*w/np.pi, np.abs(h), 'b')
plt.plot(cutoff, 0.5*np.sqrt(2), 'ko')
plt.axvline(cutoff, color='k')
plt.xlim(0, 0.5*fs)
plt.title("Lowpass Filter Frequency Response")
plt.xlabel('Frequency [Hz]')
plt.grid()


data = ua[:, 0, 50]

# Filter the data, and plot both the original and filtered signals.
y = butter_lowpass_filter(data, 10)

plt.subplot(2, 1, 2)
plt.plot(np.arange(len(data)), data, 'b-', label='data')
plt.plot(np.arange(len(data)), y, 'g-', linewidth=2, label='filtered data')
plt.xlabel('Time [sec]')
plt.grid()
plt.legend()

plt.subplots_adjust(hspace=0.35)
plt.show()


#results = cd.WRtool_core(var_season, lat, lon, dates_season, area, climate_mean = climate_mean, dates_climate_mean = dates_climate_mean, numpcs = inputs['numpcs'], perc = inputs['perc'], numclus = inputs['numclus'])
