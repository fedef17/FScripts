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
cart_out_2 = '/home/fabiano/Research/lavori/TunECS/tuning/coupled/'

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
facs = np.arange(11)
#perms = list(itt.product(list(facs), repeat = len(testparams)))
#random.shuffle(perms)
ncosi = len(facs)**len(testparams)

i = 0
uffpars = np.array([uff_params[par] for par in testparams])

arrcos = np.array([False, True, True, False, True, True])
range_ok = dict()
for par in testparams:
    print(par, np.min(valchange[par][arrcos]/uff_params[par]), np.max(valchange[par][arrcos]/uff_params[par]))
    range_ok[par] = np.linspace(np.min(valchange[par][arrcos]), np.max(valchange[par][arrcos]), 11)
    if np.min(range_ok[par]/uff_params[par]) < 0.7 or np.max(range_ok[par]/uff_params[par]) > 1.3:
        range_ok[par] = np.linspace(0.7*uff_params[par], 1.3*uff_params[par], 11)


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

# processi = []
# coda = []
# outputs = []
#
# n_ok = int(ncosi/n_threads)
# print(n_ok)
# for i in range(n_threads):
#     q = Queue()
#     coda.append(q)
#     processi.append(Process(target=doforproc,args=(q, i*n_ok, i*n_ok+n_ok, )))
#     processi[i].start()
#
# for i in range(n_threads):
#     outputs.append(coda[i].get())
#
# for i in range(n_threads):
#     processi[i].join()
#
# allpi = np.concatenate([out[0] for out in outputs])
# allcha = np.concatenate([out[1] for out in outputs])
# okperms = np.concatenate([out[2] for out in outputs])
# zonchan = np.stack(np.concatenate([out[3] for out in outputs]))
#
# pickle.dump([okperms, allpi, allcha, zonchan], open(cart_out + 'paramspace_v4_{}.p'.format(meto), 'wb'))
[okperms, allpi, allcha, zonchan] = pickle.load(open(cart_out + 'paramspace_v4_{}.p'.format(meto), 'rb'))

zonchan_mean = np.mean(np.abs(zonchan), axis = 1)

#plt.ion()
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
pino = newvals/uffpars[np.newaxis, :]
okvals = np.all([np.all(pino <= 1.2, axis = 1), np.all(pino >= 0.8, axis = 1)], axis = 0)

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

### Now. Selecting the two params.
for thres in [1.0, 0.5, 0.4, 0.3]:
    print('\n\n ------------------------------------------------ \n')
    print('\n Threshold on zonal change in toa_net: {} W/m2 \n'.format(thres))
    oks = (zonchan_mean < thres)
    okpinok = pino[oks]
    oknewval = newvals[oks]
    okchan = allcha[oks]
    okzon = zonchan[oks, :]
    okpi = allpi[oks]
    okokper = okperms[oks]

    #### COLD PARAM
    chak = np.percentile(okchan, 0.1)
    zup = okchan <= chak
    print('\nselecting change larger than {} W/m2'.format(chak))

    fig, ax = plt.subplots(figsize=(16,12))
    #for co in okokper[zup]:
    for pin in okpinok[zup]:
        ax.plot(np.arange(8), pin, color = 'blue', linewidth = 0.1)

    zullo = np.argmin(okchan)
    ax.plot(np.arange(8), okpinok[zullo], color = 'orange', label = 'coldest', linewidth = 3)
    print('\ncoldest', okchan[zullo])
    print(oknewval[zullo])

    gnizi = np.array([np.sum(np.abs(pin-1)) for pin in okpinok[zup]])
    ziko = np.argmin(gnizi)
    ax.plot(np.arange(8), okpinok[zup][ziko], color = 'green', label = 'lowest param variation', linewidth = 3)
    print('\nlowest param variation', okchan[zup][ziko])
    print(oknewval[zup][ziko])

    ziko = np.argmin(zonchan_mean[oks][zup])
    ax.plot(np.arange(8), okpinok[zup][ziko], color = 'violet', label = 'lowest zonal variation', linewidth = 3)
    print('\nlowest zonal variation', okchan[zup][ziko])
    print(oknewval[zup][ziko])

    ax.legend()
    ax.set_xticks(np.arange(8))
    ax.set_xticklabels(testparams, size = 'large', rotation = 60)
    ax.set_ylabel('Relative change of param')

    fig.savefig(cart_out + 'cold_params_zon{:02d}.pdf'.format(int(thres*10)))


    #### WARM PARAM
    chak = np.percentile(okchan, 99.9)
    print('\n\nselecting change larger than {} W/m2'.format(chak))
    zup = okchan >= chak

    fig, ax = plt.subplots(figsize=(16,12))
    #for co in okokper[zup]:
    for pin in okpinok[zup]:
        ax.plot(np.arange(8), pin, color = 'blue', linewidth = 0.1)

    zullo = np.argmax(okchan)
    ax.plot(np.arange(8), okpinok[zullo], color = 'orange', label = 'warmest', linewidth = 3)
    print('\nwarmest', okchan[zullo])
    print(oknewval[zullo])

    gnizi = np.array([np.sum(np.abs(pin-1)) for pin in okpinok[zup]])
    ziko = np.argmin(gnizi)
    ax.plot(np.arange(8), okpinok[zup][ziko], color = 'green', label = 'lowest param variation', linewidth = 3)
    print('\nlowest param variation', okchan[zup][ziko])
    print(oknewval[zup][ziko])

    ziko = np.argmin(zonchan_mean[oks][zup])
    ax.plot(np.arange(8), okpinok[zup][ziko], color = 'violet', label = 'lowest zonal variation', linewidth = 3)
    print('\nlowest zonal variation', okchan[zup][ziko])
    print(oknewval[zup][ziko])

    ax.legend()
    ax.set_xticks(np.arange(8))
    ax.set_xticklabels(testparams, size = 'large', rotation = 60)
    ax.set_ylabel('Relative change of param')

    fig.savefig(cart_out + 'warm_params_zon{:02d}.pdf'.format(int(thres*10)))



##########################################################################
## campiono la distribuzione con 9 param alternative
okpinok = pino
oknewval = newvals
okchan = allcha
okzon = zonchan
okpi = allpi
okokper = okperms

print('PIOOOOOOOOO', np.min(okpi), np.max(okpi))

newsets = []
figs = []
deltas = []
for i, delta in enumerate(np.linspace(-3.2, 0.8, 10)):
    if np.abs(delta) > 0.1:
        print('\n\nselecting change close to {} W/m2'.format(delta))
        zup = np.abs(okchan - delta) < 0.05
    else:
        print('\n\nselecting change close to {} W/m2'.format(0.0))
        zup = np.abs(okchan) < 0.05

    if delta < -0.1:
        num = i+1
    elif delta > 0.1:
        num = i
    else:
        num = 0

    fig, ax = plt.subplots(figsize=(16,12))

    if np.sum(zup) > 1000:
        randcho = np.random.choice(np.sum(zup), 1000)
        for pin in okpinok[zup][randcho]:
            ax.plot(np.arange(8), pin, color = 'blue', linewidth = 0.1)
    else:
        for pin in okpinok[zup]:
            ax.plot(np.arange(8), pin, color = 'blue', linewidth = 0.1)

    gnizi1 = np.array([np.max(np.abs(pin-1)) for pin in okpinok[zup]])
    ziko1 = gnizi1 == np.min(gnizi1)
    gnizi = np.array([np.sum(np.abs(pin-1)) for pin in okpinok[zup][ziko1]])
    ziko = np.argmin(gnizi)
    ax.plot(np.arange(8), okpinok[zup][ziko1][ziko], color = 'green', label = 'lowest param variation', linewidth = 5)
    print('\nlowest param variation', okchan[zup][ziko1][ziko])
    print(oknewval[zup][ziko1][ziko])
    chak = okchan[zup][ziko1][ziko]
    if num == 0:
        newsets.insert(0, oknewval[zup][ziko1][ziko])
        deltas.insert(0, chak)
    else:
        deltas.append(chak)
        newsets.append(oknewval[zup][ziko1][ziko])

    ziko = np.argmin(zonchan_mean[zup])
    ax.plot(np.arange(8), okpinok[zup][ziko], color = 'violet', label = 'lowest zonal variation', linewidth = 5)
    print('\nlowest zonal variation', okchan[zup][ziko])
    print(oknewval[zup][ziko])

    ax.legend()
    ax.set_xticks(np.arange(8))
    ax.set_xticklabels(testparams, size = 'large', rotation = 60)
    ax.set_ylabel('Relative change of param')

    ax.set_ylim(0.67, 1.33)
    fig.suptitle('Alternative param: c{} -> {:6.2f} W/m2'.format(num, chak))
    #fig.savefig(cart_out + 'altparams_c{}.pdf'.format(i))
    if num == 0:
        figs.insert(0, fig)
    else:
        figs.append(fig)

ctl.plot_pdfpages(cart_out + 'altparams_9set.pdf', figs)

ctl.mkdir(cart_out_2 + 'tun_params')
for i, (nuvals, delta) in enumerate(zip(newsets, deltas)):
    parset = dict(zip(testparams, nuvals))
    filename = cart_out_2 + 'tun_params/ifstun_c{}.sh'.format(i)
    title = 'Alternative param: c{}. Should give a delta toa_net of {:6.2f} W/m2 in 4xCO2'.format(i, delta)
    tl.write_tunparams(filename, parset, title)


pich = []
c4ch = []
parsets = []
for nuvals in newsets:
    parset = dict(zip(testparams, nuvals))
    var_change_pi, _ = tl.calc_change_var_allparams('pi', 'toa_net', parset, method = meto)
    var_change_c4, _ = tl.calc_change_var_allparams('c4', 'toa_net', parset, method = meto)
    pich.append(var_change_pi)
    c4ch.append(var_change_c4)
    parsets.append(parset)

pickle.dump([parsets, pich, c4ch], open(cart_out_2 + 'altparams.p', 'wb'))
