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

plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
titlefont = 22
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.axisbelow'] = True
plt.rcParams['legend.fontsize'] = 18

#############################################################################

if os.uname()[1] == 'hobbes':
    cart_out = '/home/fabiano/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'xaru':
    cart_out = '/home/fedef/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'tintin':
    cart_out = '/home/fabiano/work/lavori/BOTTINO/'

cart_out = cart_out + 'albedo/'
ctl.mkdir(cart_out)


allru = ['b050', 'b100']

for ru in allru:
    filz = glob.glob(cart_out + 'alb_{}_*nc'.format(ru))
    filz.sort()

    fields_sep = []
    fields_mar = []
    yeas = []
    for fi in filz:
        pino = xr.load_dataset(fi, use_cftime = True)
        #ctl.plot_map_contour(pino['rsuscs'])
        fields_sep.append(pino['rsuscs'][8])
        fields_mar.append(pino['rsuscs'][2])
        yeas.append(str(pino.time.values[0].year))

    ctl.plot_multimap_contour(fields_sep, filename = cart_out + 'greenland_albedo_sep_{}.pdf'.format(ru), cbar_range=[0,1], plot_anomalies=False, plot_margins = (-80, -10, 60, 85), subtitles = yeas, cb_label = 'Albedo (sept)', figsize = (16,9))

    ctl.plot_multimap_contour(fields_mar, filename = cart_out + 'antarctic_albedo_mar_{}.pdf'.format(ru), cbar_range=[0,1], plot_anomalies=False, visualization = 'nearside', bounding_lat = -45, central_lat_lon = (-90, 0), subtitles = yeas, cb_label = 'Albedo (march)')
