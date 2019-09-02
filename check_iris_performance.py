#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
from matplotlib import pyplot as plt
from matplotlib import cm
import pickle
from datetime import datetime
import iris
import netCDF4 as nc

import climtools_lib as ctl
import climdiags as cd

from matplotlib.colors import LogNorm

#######################################
nr = 1
ifile = '/data-hobbes/fabiano/OBS/ERA/ERA40+Int_daily_1957-2018_zg500_remap25_meters.nc'

print('IRIS read\n')
for i in range(nr):
    t1 = datetime.now()
    fh = iris.load(ifile)
    data = fh[0].data
    print(data.shape, np.max(data))
    t2 = datetime.now()

    diff = (t2-t1).total_seconds()
    print('Data loaded in {:8.3f} s\n'.format(diff))

del data
print('NetCDF read\n')
for i in range(nr):
    t1 = datetime.now()
    fh = nc.Dataset(ifile, mode='r')
    data = list(fh.variables.values())[-1][:]
    print(data.shape, np.max(data))
    t2 = datetime.now()

    diff = (t2-t1).total_seconds()
    print('Data loaded in {:8.3f} s\n'.format(diff))

del data
#ifile = '/data-hobbes/fabiano/PRIMAVERA/highres_SST_19792014/zg500_Aday_EC-Earth3-HR_T511_regrid25_1979-2014.nc'
#ifile = '/data-hobbes/fabiano/PRIMAVERA/hist_1950/Stream1/EC-Earth-3-HR/EC-Earth-3-HR_s1_allyears_remap25.nc'
ifile = '/data-hobbes/fabiano/PRIMAVERA/hist_1950/Stream1/EC-Earth-3-HR/EC-Earth-3-HR_s1_allyears_remap25_rechunk.nc'
# ifile = '/data-hobbes/fabiano/PRIMAVERA/hist_1950/Stream1/EC-Earth-3-HR/prova.nc'
#ifile = '/home/fabiano/Research/temp/EC-Earth-3-HR_s1_allyears_remap25.nc'

print('NetCDF read\n')
for i in range(nr):
    t1 = datetime.now()
    fh = nc.Dataset(ifile, mode='r')
    data = list(fh.variables.values())[-1][:]
    print(data.shape, np.max(data))
    t2 = datetime.now()

    diff = (t2-t1).total_seconds()
    print('Data loaded in {:8.3f} s\n'.format(diff))
del data

print('IRIS read\n')
for i in range(nr):
    t1 = datetime.now()
    fh = iris.load(ifile)
    data = fh[0].data
    print(data.shape, np.max(data))
    t2 = datetime.now()

    diff = (t2-t1).total_seconds()
    print('Data loaded in {:8.3f} s\n'.format(diff))
del data

# OUTPUT:
# IRIS read
#
# ((22402, 1, 73, 144), '6043.098685788551')
# Data loaded in    4.003 s
#
# ((22402, 1, 73, 144), '6043.098685788551')
# Data loaded in    4.159 s
#
# ((22402, 1, 73, 144), '6043.098685788551')
# Data loaded in    4.157 s
#
# NetCDF read
#
# ((22402, 1, 73, 144), '6043.098685788551')
# Data loaded in    1.949 s
#
# ((22402, 1, 73, 144), '6043.098685788551')
# Data loaded in    1.911 s
#
# ((22402, 1, 73, 144), '6043.098685788551')
# Data loaded in    1.903 s
#
# NetCDF read
#
# ((13149, 1, 73, 144), 6018.4404)
# Data loaded in    3.557 s
#
# ((13149, 1, 73, 144), 6018.4404)
# Data loaded in    0.521 s
#
# ((13149, 1, 73, 144), 6018.4404)
# Data loaded in    0.519 s
#
# IRIS read
#
# ((13149, 1, 73, 144), 6018.4404)
# Data loaded in   89.475 s
#
# ((13149, 1, 73, 144), 6018.4404)
# Data loaded in   75.511 s
#
# ((13149, 1, 73, 144), 6018.4404)
# Data loaded in   75.911 s
