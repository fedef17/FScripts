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

#cart_out = '/home/fabiano/Research/lavori/BOTTINO/analisi/'
cart = '/home/federico/work/Tipes/tipes_hosing/'

fils = ['uaday_cntrl_1850-1999_850hPa_remap25.nc', 'uaday_ho03_1850-1899_850hPa_remap25.nc']#, 'uaday_c3r5_1900-1999_850hPa_remap25.nc']
#filna = '/nas/BOTTINO/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/{}/r1i1p1f1/day_r25/ua/ua*nc'

allru = ['pi', 'ho03']#, 'c3r5']
colors = ['steelblue', 'indianred']#, 'forestgreen']

#############################################################################
resdict = dict()

for filone, ru, col in zip(fils, allru, colors):
    jli, jspeed, dates = cd.jli_from_files(cart + filone, orogfile = '/home/federico/work/Tipes/geopot_vegcover_25.nc')

    resdict[(ru, 'jli')] = jli
    resdict[(ru, 'jspeed')] = jspeed
    resdict[(ru, 'dates')] = dates

with open(cart_out + 'res_jli_hosing.p', 'wb') as filox:
    pickle.dump(resdict, filox)
