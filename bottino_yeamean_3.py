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
#allvars_2D = 'clt pr psl rlut rsut tas uas'.split()
#allvars_2D = 'pr tas'.split()
allvars_2D = ['tas']
#allvars_3D = 'ta ua'.split()

#var_map_200 = 'clt pr psl tas rlut uas'.split()  # plot last 200 mean map, stddev, low/high var wrt pi
allnams2 = allnams + ['ssp585']
allru2 = allru + ['ssp585']
colors2 = colors + ['indianred']

figs_glob = []
axs_glob = []

yeamean = dict()

allseasons = ['year']

for na, ru, col in zip(allnams, allru, colors):
#for na, ru, col in zip(allnams2, allru2, colors2):
    print(ru)
    mem = 'r1'
    if na == 'ssp585': mem = 'r4'

    fils = np.concatenate([glob.glob(filna.format(na, mem, miptab, var)) for var in allvars_2D])

    kose = xr.open_mfdataset(fils, use_cftime = True)
    kose = kose.drop_vars('time_bnds')

    for var in allvars_2D:
        print(var)
        if var not in kose:
            continue

        yeamean[(ru, var, 'max')] = ctl.seasonal_set(kose[var], season = 'year', seasonal_stat = 'max')

        if var != 'pr':
            cosoye = kose[var].groupby("time.year").mean().compute()
            yeamean[(ru, var, 'mean')] = cosoye
            yeamean[(ru, var, 'min')] = ctl.seasonal_set(kose[var], season = 'year', seasonal_stat = 'min')
        else:
            yeamean[(ru, var, 'sum')] = ctl.seasonal_set(kose[var], season = 'year', seasonal_stat = 'sum')

pickle.dump(yeamean, open(cart_out + 'bottino_yeastat_tas.p', 'wb'))
