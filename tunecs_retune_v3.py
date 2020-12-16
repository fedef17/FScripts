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
from multiprocessing import Process, Queue, Pool
import random
#from multiprocessing import set_start_method
#from multiprocessing import get_context
#set_start_method("spawn")


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

meto = 'spline'

allpi = []
allcha = []
okperms = []
zonchan = []

#facs = np.arange(0.7, 1.4, 0.1)
#perms = list(itt.combinations_with_replacement(list(facs), len(testparams)))
facs = np.arange(10)
#perms = list(itt.product(list(facs), repeat = len(testparams)))
#random.shuffle(perms)

i = 0
uffpars = np.array([uff_params[par] for par in testparams])

arrcos = np.array([False, True, True, False, True, True])
range_ok = dict()
for par in testparams:
    print(par, np.min(valchange[par][arrcos]/uff_params[par]), np.max(valchange[par][arrcos]/uff_params[par]))
    range_ok[par] = np.linspace(np.min(valchange[par][arrcos]), np.max(valchange[par][arrcos]), 10)
    if np.min(range_ok[par]/uff_params[par]) < 0.7 or np.max(range_ok[par]/uff_params[par]) > 1.3:
        range_ok[par] = np.linspace(0.7*uff_params[par], 1.3*uff_params[par], 10)


########## parallel
def doforproc(q, i1, i2, meto = 'spline', facs = facs, testparams = testparams, range_ok = range_ok):
    perms = itt.product(list(facs), repeat = len(testparams))

    allpi = []
    allcha = []
    okperms = []
    zonchan = []
    i = 0
    for perm in perms:
        if i < i1:
            i+=1
            continue
        elif i > i2:
            break

        if (i-i1)%10000 == 0: print(1.0*(i-i1)/(i2-i1))
        newvals = [range_ok[par][p] for par, p in zip(testparams, perm)]

        pichan = tl.delta_pi_glob(newvals, testparams, var = 'toa_net', method = meto)
        if np.abs(pichan) < 0.1:
            c4pichan = tl.delta_c4pi_glob(newvals, testparams, var = 'toa_net', method = meto)

            allpi.append(pichan)
            allcha.append(c4pichan)
            okperms.append(perm)

            parset = dict(zip(testparams, newvals))
            var_change_glob, var_change_zonal = tl.calc_change_var_allparams('pi', 'toa_net', parset, method = meto)
            zonchan.append(var_change_zonal)

        i+=1

    q.put([allpi, allcha, okperms, zonchan])

    return

n_threads = 8

processi = []
coda = []
outputs = []

#ctx = get_context('spawn')

n_ok = int(len(perms)/n_threads)
for i in range(n_threads):
    q = Queue()
    coda.append(q)
    #perms_sp = perms[(i*n_ok):(i*n_ok)+n_ok]
    processi.append(Process(target=doforproc,args=(q, i*n_ok, i*n_ok+n_ok, )))
    #processi.append(ctx.Process(target=doforproc,args=(perms_sp, )))
    processi[i].start()

for i in range(n_threads):
    outputs.append(coda[i].get())

for i in range(n_threads):
    processi[i].join()

# for perm in perms:
#     i+=1
#     if i%1000 == 0:
#         print(1.0*i/len(perms))
#     #newvals = np.array(perm)*uffpars
#     newvals = [range_ok[par][p] for par, p in zip(testparams, perm)]
#     #parset = dict(zip(testparams, newvals))
#
#     pichan = tl.delta_pi_glob(newvals, testparams, var = 'toa_net', method = meto)
#     if np.abs(pichan) < 0.1:
#         c4pichan = tl.delta_c4pi_glob(newvals, testparams, var = 'toa_net', method = meto)
#
#         allpi.append(pichan)
#         allcha.append(c4pichan)
#         okperms.append(perm)
#
#         parset = dict(zip(testparams, newvals))
#         var_change_glob, var_change_zonal = tl.calc_change_var_allparams('pi', 'toa_net', parset, method = meto)
#         zonchan.append(var_change_zonal)

allpi = np.append([out[0] for out in outputs])
allcha = np.append([out[1] for out in outputs])
okperms = np.append([out[2] for out in outputs])
zonchan = np.stack(np.append([out[3] for out in outputs]))

allpi = np.array(allpi)
allcha = np.array(allcha)
zonchan = np.stack(zonchan)

pickle.dump([okperms, allpi, allcha, zonchan], open(cart_out + 'paramspace_v3_{}.p'.format(meto), 'wb'))
[okperms, allpi, allcha, zonchan] = pickle.load(open(cart_out + 'paramspace_v3_{}.p'.format(meto), 'rb'))

zonchan_mean = np.mean(np.abs(zonchan), axis = 1)

plt.ion()
fig = plt.figure()

for co in [1.0, 0.5, 0.4, 0.25]:
    oks = zonchan_mean < co
    okchan = allcha[oks]
    plt.hist(okchan, bins = 30)
    print(okchan.min(), okchan.max())

plt.xlabel('delta sensitivity 4xCO2-pi (W/m2)')
plt.title('Alternative parametrizations with < 0.1 W/m2 effect on toa_net in pi (< 0.5, 0.4, 0.25 in lat band)', fontsize = 'small')

fig.savefig(cart_out + 'altparam_sens_v3a.pdf')

newvals = np.stack([[range_ok[par][p] for par, p in zip(testparams, perm)] for perm in okperms])
okvals = np.all(0.8 <= newvals/uffpars[np.newaxis, :] <= 1.2, axis = 1)


fig = plt.figure()

oks = (zonchan_mean < 0.4) & okvals
okchan = allcha[oks]
plt.hist(okchan, bins = 30)
print(okchan.min(), okchan.max())

plt.xlabel('delta sensitivity 4xCO2-pi (W/m2)')
plt.title('Alt. par.: < 0.4 toa_net change in lat band, < 20% param change', fontsize = 'small')

fig.savefig(cart_out + 'altparam_sens_v3b.pdf')


fig = plt.figure()

oks = (zonchan_mean < 0.3) & okvals
okchan = allcha[oks]
plt.hist(okchan, bins = 30)
print(okchan.min(), okchan.max())

plt.xlabel('delta sensitivity 4xCO2-pi (W/m2)')
plt.title('Alt. par.: < 0.3 toa_net change in lat band, < 20% param change', fontsize = 'small')

fig.savefig(cart_out + 'altparam_sens_v3c.pdf')
