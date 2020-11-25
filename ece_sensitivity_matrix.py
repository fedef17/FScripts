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

import tunlib as tl

plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 28
plt.rcParams['axes.labelsize'] = 22

#############################################################################
if os.uname()[1] == 'hobbes':
    cart = '/home/fabiano/Research/lavori/'
    cart_out = cart + 'TunECS/tuning/experiments/analysis/'
    cart_in = cart + 'TunECS/tuning/experiments/table/'
elif os.uname()[1] == 'ff-clevo':
    cart = '/home/fedefab/Scrivania/Research/Post-doc/lavori/'
    cart_out = cart + 'TunECS/tuning/analysis/'
    cart_in = cart + 'TunECS/tuning/table/'

ctl.mkdir(cart_out)
fil = cart_in + 'gregory_{:2s}{:1s}{:1s}.txt'

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

resdic = dict()
for forc in ['pi', 'c4']:
    for nu, let, param in zip(nums, letts, testparams):
        for iic, change in enumerate(['m', 'n', 'p', 'q', 'l', 'r']):
            filnam = fil.format(forc, change, let)
            print(filnam)
            if tl.check_file(filnam):
                print('ok!', filnam)
                anni, toa_net, srf_net, tas = tl.read_gregory(filnam)
                resdic[(forc, change, let, 'toa_net')] = toa_net
                resdic[(forc, change, let, 'srf_net')] = srf_net
                resdic[(forc, change, let, 'years')] = anni
                resdic[(forc, change, let, 'tas')] = tas


filnam = cart_in + 'gregory_tpa1.txt'
print(filnam)
if tl.check_file(filnam):
    print('ok!', filnam)
    anni, toa_net, srf_net, tas = tl.read_gregory(filnam)
    forc = 'pi'
    change = '0'
    let = '0'
    resdic[(forc, change, let, 'toa_net')] = toa_net
    resdic[(forc, change, let, 'srf_net')] = srf_net
    resdic[(forc, change, let, 'years')] = anni
    resdic[(forc, change, let, 'tas')] = tas

filnam = cart_in + 'gregory_t4a1.txt'
print(filnam)
if tl.check_file(filnam):
    print('ok!', filnam)
    anni, toa_net, srf_net, tas = tl.read_gregory(filnam)
    forc = 'c4'
    change = '0'
    let = '0'
    resdic[(forc, change, let, 'toa_net')] = toa_net
    resdic[(forc, change, let, 'srf_net')] = srf_net
    resdic[(forc, change, let, 'years')] = anni
    resdic[(forc, change, let, 'tas')] = tas

# Faccio uno stupidissimo scatter plot
allforc = ['pi', 'c4']
allchan = 'm n p q l r'.split()

forcsym = dict()
forcsym['pi'] = 'o'
forcsym['c4'] = 'x'

forccol = dict()
forccol['pi'] = 'lightseagreen'
forccol['c4'] = 'indianred'

changecol = dict()
changecol['m'] = 'darkblue'
changecol['n'] = 'steelblue'
changecol['p'] = 'darkorange'
changecol['q'] = 'indianred'
changecol['l'] = 'pink'
changecol['r'] = 'violet'

for tip in ['diff', 'abs']:
    for var in ['toa_net', 'srf_net']:
        fig, ax = plt.subplots(figsize=(16,12))
        for forc, shift in zip(allforc, [-0.05, 0.05]):

            change = '0'
            let = '0'
            scatt = []
            errs = []
            vals = resdic[(forc, change, let, var)]
            ctrl = np.mean(vals[1:])
            scatt.append(ctrl)
            errs.append(np.std(vals[1:]))

            ax.errorbar([shift], scatt, yerr = errs, fmt = 'none', color = 'black', capsize = 2, elinewidth = 1)
            ax.scatter([shift], scatt, color = 'black', marker = forcsym[forc], s = 100)

            for change in allchan:
                scatt = []
                errs = []
                for nu, let, param in zip(nums, letts, testparams):
                    if (forc, change, let, var) in resdic.keys():
                        print((forc, change, let, var))
                        vals = resdic[(forc, change, let, var)]
                        if tip == 'diff':
                            scatt.append(np.mean(vals[1:]) - ctrl)
                        else:
                            scatt.append(np.mean(vals[1:]))
                        errs.append(np.std(vals[1:]))
                    else:
                        scatt.append(np.nan)
                        errs.append(np.nan)

                ax.errorbar(nums+shift, scatt, yerr = errs, fmt = 'none', color = changecol[change], capsize = 2, elinewidth = 1)
                ax.scatter(nums+shift, scatt, color = changecol[change], marker = forcsym[forc], s = 100)


        ax.set_xticks(nums)
        ax.set_xticklabels(testparams, size = 'large', rotation = 60)
        ax.grid()
        ax.set_ylabel(var + ' (W/m2)')
        if tip == 'diff':
            fig.suptitle('Radiation difference relative to control')
        else:
            fig.suptitle('Net radiation')
        fig.savefig(cart_out + var+'_scattplot_{}.pdf'.format(tip))

## Delta
for var in ['toa_net', 'srf_net']:
    fig, ax = plt.subplots(figsize=(16,12))

    change = '0'
    let = '0'
    forc = 'pi'
    vals = resdic[(forc, change, let, var)]
    ctrl_pi = np.mean(vals[1:])
    forc = 'c4'
    vals = resdic[(forc, change, let, var)]
    ctrl_c4 = np.mean(vals[1:])

    for change in ['n', 'p']:
        scatt = []
        errs = []
        for nu, let, param in zip(nums, letts, testparams):
            vals_pi = resdic[('pi', change, let, var)]
            vals_c4 = resdic[('c4', change, let, var)]
            scatt.append((np.mean(vals_c4[1:]) - ctrl_c4)-(np.mean(vals_pi[1:]) - ctrl_pi))
            errs.append(np.mean([np.std(vals_c4[1:]),np.std(vals_pi[1:])]))

        ax.errorbar(nums, scatt, yerr = errs, fmt = 'none', color = changecol[change], capsize = 2, elinewidth = 1)
        ax.scatter(nums, scatt, color = changecol[change], marker = 'D', s = 100)


    ax.set_xticks(nums)
    ax.set_xticklabels(testparams, size = 'large', rotation = 60)
    ax.grid()
    ax.set_ylabel('delta '+ var + ' (W/m2)')

    fig.suptitle('Delta 4xco2 change - pi change')

    fig.savefig(cart_out + var+'_scattplot_{}.pdf'.format('delta'))

## Derivata con parametro normalizzato
derdic = dict()
derdic_err = dict()
for var in ['toa_net', 'srf_net']:
    fig, ax = plt.subplots(figsize=(16,12))

    ctrl = dict()
    change = '0'
    let = '0'
    forc = 'pi'
    vals = resdic[(forc, change, let, var)]
    ctrl['pi'] = np.mean(vals[1:])
    forc = 'c4'
    vals = resdic[(forc, change, let, var)]
    ctrl['c4'] = np.mean(vals[1:])

    for forc, shift in zip(allforc, [-0.05, 0.05]):
        ders = []
        err_ders = []
        for nu, let, param in zip(nums, letts, testparams):
            cose = []
            errs = []
            xval = []
            for ii, change in zip([1, -1, 2], ['n', '0', 'p']):
                if change == '0':
                    cose.append(ctrl[forc])
                else:
                    vals = resdic[(forc, change, let, var)]
                    cose.append(np.mean(vals[1:]))
                    errs.append(np.std(vals[1:]))

                if ii >= 0:
                    xval.append(valchange[param][ii])
                else:
                    xval.append(uff_params[param])

            deriv = np.gradient(np.array(cose), np.array(xval))

            print(param, deriv)
            ders.append(uff_params[param]*deriv[1])
            derdic[(forc, param, var)] = deriv[1]
            if tl.check_increasing(xval):
                derdic[(forc, param, var, 'left')] = deriv[0]
                derdic[(forc, param, var, 'right')] = deriv[-1]
            elif tl.check_decreasing(xval):
                derdic[(forc, param, var, 'left')] = deriv[-1]
                derdic[(forc, param, var, 'right')] = deriv[0]
            else:
                print(xval)
                raise ValueError('problema tenemos')

            errs[1] = -errs[1]
            errs.insert(1, 0.)
            deriv_err = np.gradient(np.array(errs), np.array(xval))
            err_ders.append(uff_params[param]*np.abs(deriv_err[1]))
            derdic_err[(forc, param, var)] = np.abs(deriv_err[1])

        ax.errorbar(nums+shift, ders, yerr = err_ders, fmt = 'none', color = forccol[forc], capsize = 2, elinewidth = 1)
        ax.scatter(nums+shift, ders, color = forccol[forc], marker = forcsym[forc], s = 100, label = forc)


    ax.legend()
    ax.set_xticks(nums)
    ax.set_xticklabels(testparams, size = 'large', rotation = 60)
    ax.grid()
    ax.axhline(0., color = 'black')
    ax.set_ylabel('uff_param * derivative of '+ var + ' (W/m2)')

    fig.suptitle('Derivative of {} with respect to params'.format(var))

    fig.savefig(cart_out + var+'_scattplot_{}.pdf'.format('deriv'))

with open(cart_out + 'der_sensmat_global.p', 'wb') as filox:
    pickle.dump([resdic, derdic, derdic_err], filox)
