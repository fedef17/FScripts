#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import netCDF4 as nc
import climtools_lib as ctl
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

plt.ion()

sys.exit()

filename = '/data-hobbes/fabiano/SPHINX/zg_daily/ridottone.nc'
var_day, lat_day, lon_day, dates_day, time_units, var_units = ctl.read3Dncfield(filename)

dates_pdh_day = pd.to_datetime(dates_day)

filename='/data-hobbes/fabiano/SPHINX/tas_mon/ridotto.nc'
var_mon, lat_mon, lon_mon, dates_mon, time_units, var_units = ctl.read3Dncfield(filename)

dates_pdh_mon = pd.to_datetime(dates_mon)
