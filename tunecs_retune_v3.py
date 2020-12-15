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

import tunlib as tl

from scipy.optimize import Bounds, minimize, least_squares
import itertools as itt

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 18

################################################################################

cart_in = '/data-hobbes/fabiano/TunECS/AMIP_exps/'
cart_out = '/home/fabiano/Research/lavori/TunECS/tuning/experiments/analysis/'

mname = '{:2s}{:1s}{:1s}'
fil = cart_in + '{:4s}/post/mon/Post_{:4d}/{:4s}_{:4d}_{}.nc'

testparams = ['ENTRORG', 'RPRCON', 'DETRPEN', 'RMFDEPS', 'RVICE', 'RSNOWLIN2', 'RCLDIFF', 'RLCRIT_UPHYS']
nums = np.arange(8)
letts = 'a b c d e f g h'.split()

uff_params = dict()
uff_params['RPRCON'] = 1.34E-3
uff_params['RVICE'] = 0.137
uff_params['RLCRITSNOW'] = 4.0E-5
uff_params['RSNOWLIN2'] = 0.035
uff_params['ENTRORG'] = 1.70E-4
uff_params['DETRPEN'] = 0.75E-4
uff_params['ENTRDD'] = 3.0E-4
uff_params['RMFDEPS'] = 0.3
uff_params['RCLDIFF'] = 3.E-6
uff_params['RCLDIFFC'] = 5.0
uff_params['RLCRIT_UPHYS'] = 0.875e-5

bounds = (np.array([1e-4, 9e-4, 0.4e-4, 0.1, 0.08, 0.015, 1.e-6, 0.7e-5]), np.array([2.4e-4, 2.e-3, 1.3e-4, 0.5, 0.2, 0.055, 5.e-6, 1.1e-5]))

diseqs = dict()
diseqs['m'] = -2
diseqs['n'] = -1
diseqs['p'] = 1
diseqs['q'] = 2
diseqs['l'] = -0.5
diseqs['r'] = 0.5

valchange = dict()
valchange['ENTRORG'] = np.array([1.05, 1.35, 2.05, 2.35, 1.53, 1.87])*1e-4 # Ollinaho: 1.8 +/- 0.3. bounds: 0.9-2.6
valchange['RPRCON'] = np.array([2.45, 1.9, 0.8, 0.25, 1.62, 1.07])*1e-3 # Ollinaho: 1.5 +/- 0.2. bounds: 0.7-2.1
valchange['DETRPEN'] = np.array([0.1, 0.25, 1.25, 1.75, 0.5, 1.0])*1e-4 # Ollinaho: 0.78 +/- 0.07. bounds: 0.4-1.1
valchange['RMFDEPS'] = np.array([0.02, 0.16, 0.44, 0.58, 0.23, 0.37])
valchange['RVICE'] = np.array([0.06, 0.098, 0.176, 0.214, 0.118, 0.157])
valchange['RSNOWLIN2'] = np.array([0.079, 0.057, 0.013, 0.001, 0.046, 0.024])
valchange['RCLDIFF'] = np.array([5, 4, 2, 1, 3.5, 2.5])*1e-6
valchange['RLCRIT_UPHYS'] = np.array([1.02, 0.95, 0.8, 0.73, 0.91, 0.84])*1e-5

allforc = ['pi', 'c4']
allchan = 'm n p q l r'.split()

def val_ok(param, change):
    iok = allchan.index(change)
    return valchange[param][iok]
##########################################

##### OK. This scripts estimates the best values for all parameters, once one or more of them are set.
xtol = 1.e-12
gtol = 1.e-12

allvars = ['toa_net', 'srf_net', 'cp', 'lsp', 'tcc']

climvars = dict()
for var in allvars:
    glo, zon = tl.clim_var('pi', var)
    climvars[var] = glo
    climvars[(var, 'zonal')] = zon
    print(var, glo)

print('\n\n ------------------------------ \n\n')

allpi = []
allcha = []

facs = np.arange(0.7, 1.4, 0.1)
perms = list(itt.combinations_with_replacement(list(facs), len(testparams)))

print(len(perms))
i = 0
uffpars = np.array([uff_params[par] for par in testparams])

for perm in perms:
    i+=1
    if i%1000 == 0:
        print(perm)
    newvals = np.array(perm)*uffpars
    #parset = dict(zip(testparams, newvals))

    pichan = tl.delta_pi_glob(newvals, testparams, var = 'toa_net', method = 'spline')
    c4pichan = tl.delta_c4pi_glob(newvals, testparams, var = 'toa_net', method = 'spline')

    allpi.append(pichan)
    allcha.append(c4pichan)

allpi = np.array(allpi)
allcha = np.array(allcha)
pickle.dump([perms, allpi, allcha], open(cart_out + 'paramspace_v1.p', 'wb'))

oks = np.abs(pichan) < 0.05
okchan = allcha[oks]

plt.ion()
plt.hist(okchan)
