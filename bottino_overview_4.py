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
#from tunlib import gregplot_on_ax

from matplotlib.colors import LogNorm
from datetime import datetime

from scipy import stats
import xarray as xr
import glob
import xclim

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 22
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['axes.axisbelow'] = True

#############################################################################

if os.uname()[1] == 'hobbes':
    cart = '/home/fabiano/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'xaru':
    cart = '/home/fedef/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'tintin':
    cart = '/home/fabiano/work/lavori/BOTTINO/'

cart_in = cart + 'seasmean/'
cart_out = cart + 'trends/'
ctl.mkdir(cart_out)

allru = ['pi', 'b025', 'b050', 'b100']
allnams = ['piControl', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']

colors = ['black', 'forestgreen', 'orange', 'violet']

####################################################################################################

glomeans, pimean, yeamean, mapmean = pickle.load(open(cart_in + 'bottino_seasmean_2D.p', 'rb'))

# century trends. or bicentury?

for var in ['tas', 'pr', 'clt', 'rlut']:
    trendz = []
    tits = []
    for ru in allru[1:]:
        for cent, start, end in zip(['1st', '2nd', '5th'], [0, 100, 200], [100, 200, 500]):
            coso = yeamean[(ru, var)][start:end]
            gtas = glomeans[(ru, var)][start:end]
            print(coso.shape, gtas.shape)

            var_trend, var_intercept, var_trend_err, var_intercept_err = ctl.calc_trend_climatevar(gtas, coso)
            trendz.append(var_trend)
            tits.append(ru + ' - ' + cent)

    ctl.plot_multimap_contour(trendz, coso.lat, coso.lon, filename = cart_out+var+'_trendz.pdf', fix_subplots_shape = (3,3), figsize = (16, 9), cbar_range = (0, 2), subtitles = tits)
