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

filcart = '/nas/BOTTINO/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/{}/'
filna = filcart + '{}/{}/{}/*nc'

####################################################################################################

miptab = 'Amon'
allvars_2D = 'pr tas'.split()

allseasons = ['DJFM', 'JJAS']

#for na, ru, col in zip(allnams, allru, colors):
for var in allvars_2D:
    yeamean = dict()
    seamean = dict()
    print(var)
    #for na, ru, col in zip(allnams2, allru2, colors2):
    for exp in ['historical', 'ssp585']:
        print(exp)
        members = os.listdir(filcart.format(exp))
        memok = []
        for mem in members:
            print(mem)
            fils = glob.glob(filna.format(exp, mem, miptab, var))

            if len(fils) == 0:
                print('NO data for {} {}'.format(var, exp, mem))
                continue
            memok.append(mem)

            kose = xr.open_mfdataset(fils, use_cftime = True)
            kose = kose.drop_vars('time_bnds')

            cosoye = kose[var].groupby("time.year").mean().compute()
            yeamean[(exp, mem, var)] = cosoye

            for sea in allseasons:
                seamean[(exp, mem, var, sea)] = ctl.seasonal_set(kose[var], season = sea, seasonal_stat = 'mean')

        cosoye = np.mean([yeamean[(exp, mem, var)] for mem in memok], axis = 0)
        yeamean[(exp, 'ensmean', var)] = cosoye
        cosoye = np.std([yeamean[(exp, mem, var)] for mem in memok], axis = 0)
        yeamean[(exp, 'ensstd', var)] = cosoye
        yeamean[(exp, 'members')] = memok

        pickle.dump([yeamean, seamean], open(cart_out + 'bottino_yeamean_3_{}_{}.p'.format(exp, var), 'wb'))
