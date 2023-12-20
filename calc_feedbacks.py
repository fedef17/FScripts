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
import pandas as pd
import glob

#import tunlib as tl

from scipy.optimize import Bounds, minimize, least_squares
import itertools as itt
from multiprocessing import Process, Queue, Pool
import random
#from multiprocessing import set_start_method
#from multiprocessing import get_context
#set_start_method("spawn")

import climtools_lib as ctl
import xarray as xr


plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 18

################################################################################

cart_out = '/home/fabiano/Research/lavori/TunECS/results/feedbacks/'
ctl.mkdir(cart_out)



exps = ['fml1', 'amc5', 'amc9']
allvars = ['tcc', 'ttr', 'ttrc', 'tsr', 'tsrc']
longnames = ['Cloud cover', 'OLR (W/m2)', 'clear-sky OLR (W/m2)', 'OSR (W/m2)', 'clear-sky OSR (W/m2)']


/data-hobbes/fabiano/TunECS/coupled/pic5/cmorized/cmor_1898/CMIP6/CMIP/EC-Earth-Consortium/EC-Earth3/piControl/r1i1p5f1/Amon/rlutcs/gr/v20210204/

cart_in = '/data-hobbes/fabiano/TunECS/coupled/'
finam = '{}_clim00-14_{}_r1.nc'

cose = dict()

for exp in exps:
    for varnam in allvars:
        var, coords, auxv = ctl.read_xr(cart_in.format(exp) + finam.format(exp, varnam))

        mean_field = np.mean(var, axis = 0)
        # zon_mean, zon_std = ctl.zonal_seas_climatology(var, coords['dates'], 'year')
        # glob_mean, glob_std = ctl.global_seas_climatology(var, coords['lat'], coords['dates'], 'year')

        cose[(exp, varnam)] = mean_field