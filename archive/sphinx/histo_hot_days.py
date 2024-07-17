#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import netCDF4 as nc
import climtools_lib as ctl
import pandas as pd
import numpy as np
from numpy import linalg as LA
from matplotlib import pyplot as plt
from matplotlib import colors as mpl_colors
from matplotlib import cm
import pickle

import matplotlib.animation as animation
from matplotlib.animation import ImageMagickFileWriter

########################################################
# read masks
filo = '/home/fedefab/Scrivania/Research/Post-doc/SPHINX/masks.nc'
fh = nc.Dataset(filo)
land_mask = fh.variables['RnfO.msk'][:]
#ocean_mask = fh.variables['RnfA.msk'][:]
land_mask = land_mask[128:, :]

#lat_range = [37, 45]
#lon_range = [7, 18]
lat_range = [37, 45]
lon_range = [-79, -70]
plt.ion()
binzz = np.linspace(-20,25,50)

cart = '/home/fedefab/Scrivania/Research/Post-doc/SPHINX/Daily_temp/'
var_all = dict()

for ran in [1950,2020,2090]:
    fil = 'lcb0_tas_{}-{}.nc'.format(ran, str(ran+5)[-2:])
    var, lat, lon, dates, time_units, var_units = ctl.read3Dncfield(cart+fil)
    var = var - 273.15
    var_all[ran] = var

    var_ok, lat_area, lon_area = ctl.sel_area(lat, lon, var, lat_range+lon_range)
    mask_ok, lat_area, lon_area = ctl.sel_area(lat, lon, land_mask, lat_range+lon_range)

    listavals = []
    for cos in var_ok:
        cosok = np.min(cos[mask_ok])
        #cosok = np.percentile(cos[mask_ok], 20)
        listavals.append(cosok)
    plt.hist(listavals, bins = binzz, label = str(ran), alpha = 0.5)

plt.legend()
plt.xlim(-20,None)
