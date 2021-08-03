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
    cart_out = '/home/fabiano/Research/lavori/BOTTINO/yearmean/'
elif os.uname()[1] == 'xaru':
    cart_out = '/home/fedef/Research/lavori/BOTTINO/yearmean/'
elif os.uname()[1] == 'tintin':
    cart_out = '/home/fabiano/work/lavori/BOTTINO/yearmean/'

ctl.mkdir(cart_out)

filna = '/nas/BOTTINO/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/{}/{}i1p1f1/{}/{}/*nc'

allru = ['pi', 'b025', 'b050', 'b100']
allnams = ['piControl', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']

colors = ['black', 'forestgreen', 'orange', 'violet']

####################################################################################################

miptab = 'Amon'

figs_glob = []
axs_glob = []

seamean = dict()

allseasons = ['NDJFM', 'MJJAS', 'year']

# for na, ru, col in zip(allnams, allru, colors):
#     print(ru)
#     mem = 'r1'
#
#     fils = glob.glob(filna.format(na, mem, miptab, 'pr'))
#
#     kose = xr.open_mfdataset(fils, use_cftime = True)
#     kose = kose.drop_vars('time_bnds')
#
#     for sea in allseasons:
#         seamean[(ru, sea)] = ctl.seasonal_set(kose['pr'], season = sea, seasonal_stat = 'sum')
#
# pickle.dump(seamean, open(cart_out + 'bottino_seamean_pr.p', 'wb'))

seamean = pickle.load(open(cart_out + 'bottino_seamean_pr.p', 'rb'))

secuse = dict()
for ru in allru:
    for sea in allseasons:
        for iiy in range(5):
            pino = seamean[(ru, sea)].isel(year = slice(100*iiy, 100*(iiy+1) )).mean('year')
            secuse[(ru, sea, iiy)] = pino

        pino = seamean[(ru, sea)].isel(year = slice(-200, None)).mean('year')
        secuse[(ru, sea, 'stab')] = pino

        if 'pi' in ke:
            secuse[(ru, sea, 'stab')] = seamean[(ru, sea)].mean('year')

pickle.dump(secuse, open(cart_out + 'bottino_seamean_pr_secular.p', 'wb'))
