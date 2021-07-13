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

miptab = 'Amon_r25'
allvars_2D = 'clt pr psl rlut rsut tas uas'.split()
allvars_3D = 'ta ua'.split()

var_map_200 = 'clt pr psl tas rlut uas'.split()  # plot last 200 mean map, stddev, low/high var wrt pi
allnams2 = allnams + ['ssp585']
allru2 = allru + ['ssp585']
colors2 = colors + ['indianred']

figs_glob = []
axs_glob = []
pimean = dict()
glomeans = dict()
yeamean = dict()
mapmean = dict()

for na, ru, col in zip(allnams, allru, colors):
#for na, ru, col in zip(allnams2, allru2, colors2):
    print(ru)
    mem = 'r1'
    if na == 'ssp585': mem = 'r4'

    fils = np.concatenate([glob.glob(filna.format(na, mem, miptab, var)) for var in allvars_2D[:-1]])

    kose = xr.open_mfdataset(fils, use_cftime = True)
    kose = kose.drop_vars('time_bnds')

    # Separate for uas
    fils = glob.glob(filna.format(na, mem, miptab, allvars_2D[-1]))
    if len(fils) > 0:
        kosettt = xr.open_mfdataset(fils, use_cftime = True)
        kosettt = kosettt.drop_vars('time_bnds')
        kosettt = kosettt.drop_vars('height')
        kose = kose.assign(uas = kosettt.uas)

    for var in allvars_2D:
        print(var)
        if var not in kose:
            continue

        cosoye = kose[var].groupby("time.year").mean().compute()
        yeamean[(ru, var)] = cosoye


# 3D vars
for na, ru, col in zip(allnams, allru, colors):
    mem = 'r1'
    if na == 'ssp585': mem = 'r4'

    fils = np.concatenate([glob.glob(filna.format(na, mem, miptab, var)) for var in allvars_3D])

    kose = xr.open_mfdataset(fils, use_cftime = True)
    kose = kose.drop_vars('time_bnds')

    for var in allvars_3D:
        print(var)
        cosoye = kose[var].groupby("time.year").mean().compute()
        yeamean[(ru, var)] = cosoye


pickle.dump(yeamean, open(cart_out + 'bottino_yeamean.p', 'wb'))
