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

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 15
###############################################################################

cart = '/home/fabiano/Research/lavori/TunECS/tuning/ifs-tuning-params/'

# IFS tuning
# comment

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

allparams = ['RPRCON', 'RVICE', 'RLCRITSNOW', 'RSNOWLIN2', 'ENTRORG', 'DETRPEN', 'ENTRDD', 'RMFDEPS', 'RCLDIFF', 'RCLDIFFC', 'RLCRIT_UPHYS']

diseqs = dict()
diseqs['m'] = -2
diseqs['n'] = -1
diseqs['p'] = 1
diseqs['q'] = 2
diseqs['l'] = -0.5
diseqs['r'] = 0.5


valchange = dict()
valchange['ENTRORG'] = np.array([1.05, 1.35, 2.05, 2.35, 1.53, 1.87])*1e-4
valchange['RPRCON'] = np.array([2.45, 1.9, 0.8, 0.25, 1.62, 1.07])*1e-3
valchange['DETRPEN'] = np.array([0.1, 0.25, 1.25, 1.75, 0.5, 1.0])*1e-4
valchange['RMFDEPS'] = np.array([0.02, 0.16, 0.44, 0.58, 0.23, 0.37])
valchange['RVICE'] = np.array([0.06, 0.098, 0.176, 0.214, 0.118, 0.157])
valchange['RSNOWLIN2'] = np.array([0.079, 0.057, 0.013, 0.001, 0.046, 0.024])
valchange['RCLDIFF'] = np.array([5, 4, 2, 1, 3.5, 2.5])*1e-6
valchange['RLCRIT_UPHYS'] = np.array([1.02, 0.95, 0.8, 0.73, 0.91, 0.84])*1e-5

testparams = ['ENTRORG', 'RPRCON', 'DETRPEN', 'RMFDEPS', 'RVICE', 'RSNOWLIN2', 'RCLDIFF', 'RLCRIT_UPHYS']
nums = np.arange(2,10)
letts = 'a b c d e f g h'.split()

filnom = cart + 'ifstun_{:1s}{:1s}.sh'

for nu, let, param in zip(nums, letts, testparams):
    for iic, change in enumerate(['m', 'n', 'p', 'q', 'l', 'r']):
        filout = filnom.format(change, let)
        filok = open(filout, 'w')
        filok.write('# IFS tuning\n')
        filok.write('# Change of param {} to get a disequilibrium of {} W/m2 at TOA\n'.format(param, diseqs[change]))

        for par in allparams:
            if par != param:
                filok.write('{}={:.3e}\n'.format(par, uff_params[par]))
            else:
                nuval = valchange[param][iic]
                filok.write('{}={:.3e}\n'.format(param, nuval))

        filok.close()
