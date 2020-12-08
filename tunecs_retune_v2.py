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

allsets = []
for rprval, c4pi_change in zip([0.0008, 0.001, 0.0012, 0.0015, 0.0017, 0.0019], [1., 0.7, 0.5, -0.5, -1, -1.7]):
    print('\n\n ------------------------------ \n\n')
    parset = {'RPRCON' : rprval}
    print('\n PARSET: \n')
    print(parset)
    if rprval < 0.00134:
        print('HIGH ECS! \n')
    else:
        print('LOW ECS! \n')

    okparams = ['ENTRORG', 'DETRPEN', 'RMFDEPS', 'RVICE', 'RSNOWLIN2', 'RCLDIFF', 'RLCRIT_UPHYS']
    start = [uff_params[par] for par in okparams]
    okbounds_lo = np.array([bo for bo, par in zip(bounds[0], testparams) if par in okparams])
    okbounds_hi = np.array([bo for bo, par in zip(bounds[1], testparams) if par in okparams])
    okbounds = (okbounds_lo, okbounds_hi)

    result = least_squares(tl.delta_maxmin_glob, start, jac = tl.jac_delta_maxmin_glob, args = (okparams, parset, 'toa_net', c4pi_change, 'deriv_edge', ), verbose=1, method = 'trf', bounds = okbounds, xtol = xtol, gtol = gtol)
    nuvals = result.x
    nudic = dict(zip(okparams, nuvals))
    parset.update(nudic)

    print('\n PARSET: \n')
    print(parset)

    print('\nChange in pi\n')
    for var in allvars:
        cglob, czon = tl.calc_change_var_allparams('pi', var, parset)
        print('{:8s}: {:6.3e}  {:6.3f}'.format(var, cglob, cglob/climvars[var]))
        print(('{:8s}: '+5*' {:6.3e}').format('zonal', *czon))

    print('\nChange in sensitivity\n')
    for var in allvars:
        cglob, czon = tl.calc_change_c4pi_allparams(var, parset)
        print('{:8s}: {:6.3e}  {:6.3f}'.format(var, cglob, cglob/climvars[var]))
        print(('{:8s}: '+5*' {:6.3e}').format('zonal', *czon))

    allsets.append(parset)


fig, ax = plt.subplots(figsize=(16,12))

colors = ctl.color_set(len(allsets))
for parset, col in zip(allsets, colors):
    parvals = [parset[param]/uff_params[param] for param in testparams]
    ax.scatter(np.arange(len(parvals)), parvals, color = col, s = 100, label = 'RPRCON = {:6.4f}'.format(parset['RPRCON']))

ax.set_xticks(np.arange(len(parvals)))
ax.set_xticklabels(testparams, size = 'large', rotation = 60)
ax.grid()
ax.set_ylabel('change wrt uff. value')
ax.legend()

fig.savefig(cart_out + 'tunecs_retune_v2.pdf')


# Warming: RPRCON- : 0.0010, RSNOWLIN2+ : 0.05, ENTRORG- : 0.00014 (or less),
# Cooling: RPRCON+ : 0.0018, RSNOWLIN2+ : 0.02,  ENTRORG+ : 0.00020, DETRPEN+/-, RMFDEPS +/-, RVICE +,  RCLDIFF-, RLCRIT_UPHYS +/-
# More prec:

# DETRPEN and RMFDEPS may be used to balance precipitation, in response to RPRCON/ENTRORG changes. Also, this might produce extreme wetting/drying in the future

#cglob, czon = tl.calc_change_var_allparams(forc, var, new_set, method = 'deriv_edge')

#cglob, czon = tl.calc_change_var(forc, param, var, valchange[param][iic], method = 'deriv_edge')
