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

plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 28
plt.rcParams['axes.labelsize'] = 22

#############################################################################
cart_out = '/home/fabiano/Research/lavori/TunECS/tuning/experiments/analysis/'

# Loading the derivatives of the vars wrt params
with open(cart_out + 'der_sensmat_zonal.p', 'rb') as filox:
    _, _, derdic, _ = pickle.load(filox)

with open(cart_out + 'der_sensmat_global.p', 'rb') as filox:
    _, derdic_glo, _ = pickle.load(filox)

derdic.update(derdic_glo)

spldic = None

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

lats = [-90, -65, -40, -20, 20, 40, 65, 90]
bands = [(la1, la2) for la1, la2 in zip(lats[:-1], lats[1:])]
lacen = np.array([np.mean(laol) for laol in bands])

##################################################################################

def calc_change_var(forc, param, var, newpar_val, method = 'deriv', derdic = derdic, spldic = spldic, uff_params = uff_params):
    """
    Calculates the change in the global and zonal vars for a single parameter change, using the derivatives or the splines.

    Available methods: deriv, deriv_edge, spline
    """

    if method == 'deriv':
        der_g = derdic[(forc, param, var)]
        var_change_glob = (newpar_val-uff_params[param])*der_g

        der_z = np.array([derdic[(forc, param, var, band)] for band in bands])
        var_change_zonal = (newpar_val-uff_params[param])*der_z
    elif method == 'deriv_edge':
        if newpar_val < uff_params[param]:
            edge = 'left'
        else:
            edge = 'right'
        der_g = derdic[(forc, param, var, edge)]
        var_change_glob = (newpar_val-uff_params[param])*der_g

        der_z = np.array([derdic[(forc, param, var, band, edge)] for band in bands])
        var_change_zonal = (newpar_val-uff_params[param])*der_z
    elif method == 'spline':
        pass
    else:
        raise ValueError('method {} not available'.format(method))

    return var_change_glob, var_change_zonal


def calc_change_var_allparams(forc, var, newpar_set, method = 'deriv'):
    """
    Simply sums the change for each param.
    """
    var_change_glob = 0.
    var_change_zonal = np.zeros(len(bands))

    for param in newpar_set:
        cglob, czon = calc_change_var(forc, param, var, newpar_set[param], method = method)
        var_change_glob += cglob
        var_change_zonal += czon

    return var_change_glob, var_change_zonal


def read_gregory(filnam):
    with open(filnam, 'r') as filoz:
        filoz.readline()
        linee = filoz.readlines()
        anni = np.array([int(lin.split(':')[0].split('(')[1][:4]) for lin in linee])

        cose = np.stack([lin.rstrip().split(':')[1].split() for lin in linee])

        toa_net = np.array([float(co) for co in cose[:, 0]])
        srf_net = np.array([float(co) for co in cose[:, 1]])
        tas = np.array([float(co) for co in cose[:, 2]])

        gigi = np.argsort(anni)
        anni = anni[gigi]
        toa_net = toa_net[gigi]
        srf_net = srf_net[gigi]
        tas = tas[gigi]

    return anni, toa_net, srf_net, tas

def check_file(filnam):
    if os.path.exists(filnam):
        return True
    else:
        return False

def check_increasing(arr):
    return np.all(np.diff(arr) > 0)

def check_decreasing(arr):
    return np.all(np.diff(arr) < 0)
