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

cart = '/home/fabiano/Research/lavori/TunECS/tuning/sobol_params/'
ctl.mkdir(cart)

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


testparams = ['ENTRORG', 'RPRCON', 'DETRPEN', 'RMFDEPS', 'RVICE', 'RSNOWLIN2', 'RCLDIFF', 'RLCRIT_UPHYS']

#####
sob = stats.qmc.Sobol(d=8)
gino = sob.random_base2(m=6)
pickle.dump([sob, gino], open(cart + 'sobolset.p', 'wb'))

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(gino[:,0], gino[:,1], gino[:,2])
ax.scatter(gino[:,3], gino[:,4], gino[:,5])
ax.scatter(gino[:,6], gino[:,7], gino[:,0])
fig.savefig(cart + 'sobol_3d.pdf')

maxpert = 0.3 # maximum perturbation of 30%

filnom = cart + 'ifstun_sob{:02d}.sh'
for co in range(len(gino)):
    filout = filnom.format(co)
    filok = open(filout, 'w')
    filok.write('# IFS tuning\n')
    filok.write('# Perturbed params (sobol d=7, m=6; max pert 30%), combination {} \n'.format(co))

    indpa = 0
    for par in allparams:
        if par not in testparams:
            filok.write('{}={:.3e}\n'.format(par, uff_params[par]))
        else:
            gik = gino[co, indpa]
            pert = (gik - 0.5)*maxpert
            nuval = uff_params[par] * (1 + pert)
            # if testparams[indpa] != par:
            #     raise ValueError('sbajato')

            filok.write('{}={:.3e}\n'.format(par, nuval))
            indpa += 1

    filok.close()
