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

# cosa voglio fare: per ogni membro, faccio la media del toa_long e del toa_short annuale. e metto a punto una funzione che dato il cambio dei parametri, mi da l'effetto atteso sulla toa_net. Quindi forse voglio anche qui le derivate. yes, derivate normalizzate per ogni parametro.
# prima leggo, poi faccio la media zonale?, meglio in bande? SI. [-90, -65, -40, -20, 20, 40, 65, 90]
# Salvo la media nelle bande per ogni exp
# poi faccio le derivate, banda per banda, e le plotto, normalizzate. un plot per pi e uno per c4, per ogni var.

lats = [-90, -65, -40, -20, 20, 40, 65, 90]
bands = [(la1, la2) for la1, la2 in zip(lats[:-1], lats[1:])]

allvars = ['ttr', 'tsr', 'tnr', 'tcc', 'cp', 'lsp']

resdic = dict()
resdic_err = dict()
for varnam in allvars:
    print(varnam)
    for forc in ['pi', 'c4']:
        for nu, let, param in zip(nums, letts, testparams):
            for iic, change in enumerate(['m', 'n', 'p', 'q', 'l', 'r']):
                mok = mname.format(forc, change, let)
                print(mok)

                listafil = [fil.format(mok, ye, mok, ye, varnam) for ye in range(1851, 1855)] # skipping the first year

                if tl.check_file(listafil[0]):
                    var, coords, aux_info = ctl.read_ensemble_iris(listafil)
                    varm, varstd = ctl.seasonal_climatology(var, coords['dates'], 'year')
                    varm = ctl.zonal_mean(varm)
                    varstd = ctl.zonal_mean(varstd)

                    for band in bands:
                        laok = (coords['lat'] > band[0]) & (coords['lat'] <= band[1])
                        resdic[(forc, change, let, varnam, band)] = np.mean(varm[laok])
                        resdic_err[(forc, change, let, varnam, band)] = np.mean(varstd[laok])

        # Loading the control
        mok = mname.format(forc, 0, 0)
        print(mok)

        listafil = [fil.format(mok, ye, mok, ye, varnam) for ye in range(1851, 1855)] # skipping the first year
        if tl.check_file(listafil[0]):
            var, coords, aux_info = ctl.read_ensemble_iris(listafil)
            varm, varstd = ctl.seasonal_climatology(var, coords['dates'], 'year')
            varm = ctl.zonal_mean(varm)
            varstd = ctl.zonal_mean(varstd)

            for band in bands:
                laok = (coords['lat'] > band[0]) & (coords['lat'] <= band[1])
                resdic[(forc, 0, 0, varnam, band)] = np.mean(varm[laok])
                resdic_err[(forc, 0, 0, varnam, band)] = np.mean(varstd[laok])


derdic = dict()
derdic_err = dict()
## Derivata con parametro normalizzato
for var in allvars:
    figs = []
    for band in bands:
        fig, ax = plt.subplots(figsize=(16,12))

        ctrl = dict()
        ctrl['pi'] = resdic[('pi', 0, 0, var, band)]
        ctrl['c4'] = resdic[('c4', 0, 0, var, band)]

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
                        cose.append(resdic[(forc, change, let, var, band)])
                        errs.append(resdic_err[(forc, change, let, var, band)])

                    if ii >= 0:
                        xval.append(valchange[param][ii])
                    else:
                        xval.append(uff_params[param])

                deriv = np.gradient(np.array(cose), np.array(xval))
                derdic[(forc, change, let, varnam, band)] = deriv[1]
                print(param, deriv)
                ders.append(uff_params[param]*deriv[1])
                errs[1] = -errs[1]
                errs.insert(1, 0.)
                deriv_err = np.gradient(np.array(errs), np.array(xval))
                err_ders.append(uff_params[param]*np.abs(deriv_err[1]))
                derdic_err[(forc, change, let, varnam, band)] = deriv_err[1]

            ax.errorbar(nums+shift, ders, yerr = err_ders, fmt = 'none', color = forccol[forc], capsize = 2, elinewidth = 1)
            ax.scatter(nums+shift, ders, color = forccol[forc], marker = forcsym[forc], s = 100, label = forc)

        ax.legend()
        ax.set_xticks(nums)
        ax.set_xticklabels(testparams, size = 'large', rotation = 60)
        ax.grid()
        ax.axhline(0., color = 'black')
        ax.set_ylabel('uff_param * derivative of '+ var + ' (W/m2)')

        fig.suptitle('Derivative of {} in lat band ({}, {})'.format(var, band[0], band[1]))
        figs.append(fig)

        #fig.savefig(cart_out + var+'_scattplot_{}.pdf'.format('deriv'))

    ctl.plot_pdfpages(cart_out + '{}_sensmat_zonal.pdf'.format(var), figs)

with open(cart_out + 'der_sensmat_zonal.p', 'wb') as filox:
    pickle.dump([resdic, resdic_err, derdic, derdic_err], filox)
