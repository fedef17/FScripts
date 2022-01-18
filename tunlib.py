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
from scipy.interpolate import UnivariateSpline as spline

plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 28
plt.rcParams['axes.labelsize'] = 22

#############################################################################


cart_out = '/home/fabiano/Research/lavori/TunECS/tuning/experiments/analysis/'

def precalc_splines(cart_out):
    """
    Loads the splines for interpolation.
    """
    # Loading the derivatives of the vars wrt params
    with open(cart_out + 'der_sensmat_zonal.p', 'rb') as filox:
        resdic, _, derdic, _, linder, _, chandic, _ = pickle.load(filox)

    spldic = None
    spldic_der = None

    precalc_splines = False
    if precalc_splines:
        print('precalculating splines....')
        spldic = dict()
        spldic_der = dict()
        for ke in chandic:
            xs, ys = chandic[ke]
            spldic[ke] = spline(xs, ys, k = 2, s = 0, ext = 2) # interpolating spline of order 2, no extrapolation
            spldic_der[ke] = spldic[ke].derivative()

    return spldic, spldic_der

derdic = None
linder = None
spldic = None
chandic = None
spldic_der = None
# with open(cart_out + 'der_sensmat_global.p', 'rb') as filox:
#     gigi = pickle.load(filox)
#     derdic_glo = gigi[-2]
# derdic.update(derdic_glo)

allparams = ['RPRCON', 'RVICE', 'RLCRITSNOW', 'RSNOWLIN2', 'ENTRORG', 'DETRPEN', 'ENTRDD', 'RMFDEPS', 'RCLDIFF', 'RCLDIFFC', 'RLCRIT_UPHYS']
testparams = ['ENTRORG', 'RPRCON', 'DETRPEN', 'RMFDEPS', 'RVICE', 'RSNOWLIN2', 'RCLDIFF', 'RLCRIT_UPHYS']

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

#lats = [-90, -65, -40, -20, 20, 40, 65, 90]
lats = [-90, -60, -30, 30, 60, 90]
bands = [(la1, la2) for la1, la2 in zip(lats[:-1], lats[1:])]
lacen = np.array([np.mean(laol) for laol in bands])

# climvars = dict()
# for var in allvars:
#     glo, zon = clim_var('pi', var)
#     climvars[var] = glo
#     climvars[(var, 'zonal')] = zon
#     print(var, glo)

# var = 'toa_net'
# for var in ['toa_net', 'srf_net', 'ttr', 'tsr']:
#     filok = cart_out + 'sensmat_{}.dat'.format(var)
#     with open(filok, 'w') as flk:
#         flk.write('---- Sensitivity of {} to IFS parameters (ece 3.3.1) ---- \n\n'.format(var))
#         strfo = '{:8s}' + 5*'{:12.3e}' + '\n'
#         for forc in ['pi', 'c4']:
#             flk.write('\n\n\n ----------------- Forcing: ' + forc + '  -----------------\n')
#             for param in testparams:
#                 flk.write('\n' + param + '\n')
#                 xs = chandic[(forc, param, var, 'glob')][0]
#                 flk.write(strfo.format('val:', *xs))
#                 for band in bands+['glob']:
#                     chans = chandic[(forc, param, var, band)][1]
#                     if type(band) is str:
#                         flk.write(strfo.format(band, *chans))
#                     else:
#                         bandco = '{:2d}/{:2d}'.format(int(band[0]), int(band[1]))
#                         flk.write(strfo.format(bandco, *chans))



def gregplot_on_ax(ax, tas, toa, color = None, label = None, marker = 'D', nfirst = None, nlast = 0, calc_ERF = True, calc_ECS = True):
    """
    Plots on a gregory plot and calculates ERF (using first nfirst points) and ECS (using last nlast points).
    """

    toa5 = []
    tas5 = []
    for i in range(0, len(tas), 5):
        if len(tas) - i < 3: continue
        toa5.append(np.mean(toa[i:i+5]))
        tas5.append(np.mean(tas[i:i+5]))

    toa5 = np.array(toa5)
    tas5 = np.array(tas5)

    # toa5 = ctl.running_mean(toa, 5, remove_nans = True)
    # tas5 = ctl.running_mean(tas, 5, remove_nans = True)

    #ax.scatter(tas5, toa5, color = col, marker = mar, label = exp)
    ax.scatter(tas5[1:-1], toa5[1:-1], color = color, marker = marker, label = label)
    ax.scatter(tas5[0], toa5[0], color = color, marker = '>')
    ax.scatter(tas5[-1], toa5[-1], color = color, marker = '<')
    ax.plot(tas5, toa5, color = color, linewidth = 0.5)

    if calc_ERF:
        ax.scatter(tas[:nfirst], toa[:nfirst], s = 2, color = color)
        m, c, err_m, err_c = ctl.linear_regre_witherr(tas[:nfirst], toa[:nfirst])
        xino = np.array([0]+list(tas[:nfirst]))
        ax.plot(xino, c+m*xino, color = color, linestyle = '--', linewidth = 0.5)
        print('ERF: {} -> {:6.3f} +/- {:6.3f} W/m2'.format(label, c/2., err_c/2.))

    if calc_ECS:
        ax.scatter(tas[-nlast:], toa[-nlast:], s = 2, color = color)
        m, c, err_m, err_c = ctl.linear_regre_witherr(tas[-nlast:], toa[-nlast:])
        xino = np.array(list(tas[-nlast:])+[-c/m])
        ax.plot(xino, c+m*xino, color = color, linestyle = '--', linewidth = 0.5)
        print('ECS: {} -> {:6.3f} +/- {:6.3f} K'.format(label, -0.5*c/m, 0.5*(np.abs(err_c/c)+np.abs(err_m/m))*(-c/m)))

    return


##################################################################################
def write_tunparams(filename, new_parset, title):
    filok = open(filename, 'w')
    filok.write('# IFS tuning\n')
    filok.write('# ' + title + '\n')

    for par in allparams:
        if par in new_parset:
            filok.write('{}={:.3e}\n'.format(par, new_parset[par]))
        else:
            filok.write('{}={:.3e}\n'.format(par, uff_params[par]))

    filok.close()

    return

def clim_var(forc, var):
    """
    Gives the climatological values of var with the official param.
    """
    global_var = resdic[(forc, 0, 0, var, 'glob')]
    zonal_var = np.array([resdic[(forc, 0, 0, var, band)] for band in bands])

    return global_var, zonal_var


def calc_change_var(forc, param, var, newpar_val, method = 'deriv_edge', derdic = derdic, linder = linder, spldic = spldic, spldic_der = spldic_der, uff_params = uff_params, calc_zonal = True):
    """
    Calculates the change in the global and zonal vars for a single parameter change, using the derivatives or the splines.

    Available methods: deriv, deriv_edge, spline
    """

    if method == 'deriv':
        der_g = derdic[(forc, param, var, 'glob')]
        var_change_glob = (newpar_val-uff_params[param])*der_g

        if calc_zonal:
            der_z = np.array([derdic[(forc, param, var, band)] for band in bands])
            var_change_zonal = (newpar_val-uff_params[param])*der_z
    elif method == 'deriv_edge':
        if newpar_val < uff_params[param]:
            edge = 'left'
        else:
            edge = 'right'
        der_g = linder[(forc, param, var, 'glob', edge)]
        var_change_glob = (newpar_val-uff_params[param])*der_g

        if calc_zonal:
            der_z = np.array([linder[(forc, param, var, band, edge)] for band in bands])
            var_change_zonal = (newpar_val-uff_params[param])*der_z
    elif method == 'spline':
        spl = spldic[(forc, param, var, 'glob')]
        var_change_glob = spl(newpar_val)

        if calc_zonal:
            var_change_zonal = np.array([spldic[(forc, param, var, band)](newpar_val) for band in bands])
    else:
        raise ValueError('method {} not available'.format(method))

    if calc_zonal:
        return var_change_glob, var_change_zonal
    else:
        return var_change_glob, None


def jac_calc_change_var(forc, param, var, newpar_val, method = 'deriv_edge', derdic = derdic, linder = linder, spldic = spldic, uff_params = uff_params, calc_zonal = True):
    """
    Calculates the change in the global and zonal vars for a single parameter change, using the derivatives or the splines.

    Available methods: deriv, deriv_edge, spline
    """

    if method == 'deriv':
        der_g = derdic[(forc, param, var, 'glob')]
        var_change_glob = der_g

        if calc_zonal:
            der_z = np.array([derdic[(forc, param, var, band)] for band in bands])
            var_change_zonal = der_z
    elif method == 'deriv_edge':
        if newpar_val < uff_params[param]:
            edge = 'left'
        else:
            edge = 'right'
        der_g = linder[(forc, param, var, 'glob', edge)]
        var_change_glob = der_g

        if calc_zonal:
            der_z = np.array([linder[(forc, param, var, band, edge)] for band in bands])
            var_change_zonal = der_z
    elif method == 'spline':
        spl = spldic_der[(forc, param, var, 'glob')]
        var_change_glob = spl(newpar_val)

        if calc_zonal:
            var_change_zonal = np.array([spldic_der[(forc, param, var, band)](newpar_val) for band in bands])
    else:
        raise ValueError('method {} not available'.format(method))

    if calc_zonal:
        return var_change_glob, var_change_zonal
    else:
        return var_change_glob, None


def calc_change_var_allparams(forc, var, newpar_set, method = 'deriv_edge', calc_zonal = True):
    """
    Simply sums the change for each param.
    """
    var_change_glob = 0.
    if calc_zonal:
        var_change_zonal = np.zeros(len(bands))
    else:
        var_change_zonal = None

    for param in newpar_set:
        cglob, czon = calc_change_var(forc, param, var, newpar_set[param], method = method, calc_zonal = calc_zonal)
        var_change_glob += cglob
        if calc_zonal:
            var_change_zonal += czon

    return var_change_glob, var_change_zonal


def delta_pi_glob(newpars, okparams, fix_parset = None, var = 'toa_net', method = 'deriv_edge'):
    newpar_set = dict(zip(okparams, newpars))
    if fix_parset is not None:
        newpar_set.update(fix_parset)
    var_change_glob, var_change_zonal = calc_change_var_allparams('pi', var, newpar_set, method = method, calc_zonal = False)
    return var_change_glob


def jac_delta_pi_glob(newpars, okparams, fix_parset = None, var = 'toa_net', method = 'deriv_edge'):
    newpar_set = dict(zip(okparams, newpars))
    jac = []
    for param in newpar_set:
        der, _ = jac_calc_change_var('pi', param, var, newpar_set[param], method = method, calc_zonal = False)
        jac.append(der)

    return np.array(jac)


def delta_c4pi_glob(newpars, okparams, fix_parset = None, var = 'toa_net', method = 'deriv_edge'):
    newpar_set = dict(zip(okparams, newpars))
    if fix_parset is not None:
        newpar_set.update(fix_parset)
    var_change_glob, var_change_zonal = calc_change_c4pi_allparams(var, newpar_set, method = method, calc_zonal = False)
    return var_change_glob


def jac_delta_c4pi_glob(newpars, okparams, fix_parset = None, var = 'toa_net', method = 'deriv_edge'):
    newpar_set = dict(zip(okparams, newpars))
    jac = []
    for param in newpar_set:
        der, _ = jac_calc_change_c4pi(param, var, newpar_set[param], method = method, calc_zonal = False)
        jac.append(der)

    return np.array(jac)


def delta_maxmin_glob(newpars, okparams, fix_parset, var = 'toa_net', c4pi_change = 1.0, method = 'deriv_edge', pi_weight = 10.):
    """
    Minimizes the change for pi and gets the change in c4 closer to c4pi_change.
    """
    newpar_set = dict(zip(okparams, newpars))
    newpar_set.update(fix_parset)
    var_c4pi_glob, var_c4pi_zonal = calc_change_c4pi_allparams(var, newpar_set, method = method)
    var_pi_glob, var_pi_zonal = calc_change_var_allparams('pi', var, newpar_set, method = method)

    print(var_pi_glob, var_c4pi_glob, c4pi_change)
    delta = np.array([pi_weight*var_pi_glob, var_c4pi_glob - c4pi_change])

    return delta


def jac_delta_maxmin_glob(newpars, okparams, fix_parset, var = 'toa_net', c4pi_change = 1.0, method = 'deriv_edge'):
    """
    Minimizes the change for pi and gets the change in c4 closer to c4pi_change.
    """

    newpar_set = dict(zip(okparams, newpars))

    jac1 = []
    for param in newpar_set:
        der1, _ = jac_calc_change_c4pi(param, var, newpar_set[param], method = method)
        jac1.append(der1)

    jac2 = []
    for param in newpar_set:
        der2, _ = jac_calc_change_c4pi(param, var, newpar_set[param], method = method)
        jac2.append(der2)

    jac = np.stack([jac1, jac2])

    return jac


def calc_change_c4pi_allparams(var, newpar_set, method = 'deriv_edge', calc_zonal = True):
    """
    Gives the change in var in the c4 run wrt pi.
    """

    var_change_glob = 0.
    if calc_zonal:
        var_change_zonal = np.zeros(len(bands))
    else:
        var_change_zonal = None

    for param in newpar_set:
        cglob, czon = calc_change_c4pi(param, var, newpar_set[param], method = method, calc_zonal = calc_zonal)
        var_change_glob += cglob
        if calc_zonal:
            var_change_zonal += czon

    return var_change_glob, var_change_zonal


def calc_change_c4pi(param, var, newpar_val, method = 'deriv_edge', calc_zonal = True):
    """
    Gives the change in var in the c4 run wrt pi.
    """

    cglob_c4, czon_c4 = calc_change_var('c4', param, var, newpar_val, method = method, calc_zonal = calc_zonal)
    cglob_pi, czon_pi = calc_change_var('pi', param, var, newpar_val, method = method, calc_zonal = calc_zonal)

    var_change_glob = cglob_c4 - cglob_pi
    if calc_zonal:
        var_change_zonal = czon_c4 - czon_pi
    else:
        var_change_zonal = None

    return var_change_glob, var_change_zonal


def jac_calc_change_c4pi(param, var, newpar_val, method = 'deriv_edge', derdic = derdic, linder = linder, spldic = spldic, uff_params = uff_params):
    """
    Calculates the change in the global and zonal vars for a single parameter change, using the derivatives or the splines.

    Available methods: deriv, deriv_edge, spline
    """

    der_glob_c4, der_zonal_c4 = jac_calc_change_var('c4', param, var, newpar_val, method = method)
    der_glob_pi, der_zonal_pi = jac_calc_change_var('pi', param, var, newpar_val, method = method)

    der_glob = der_glob_c4-der_glob_pi
    der_zonal = der_zonal_c4-der_zonal_pi

    return der_glob, der_zonal


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


def read_toa_net(tablecart, expnam):
    fil = tablecart + '{}/gregory_{}.txt'.format(expnam, expnam)
    anni, toa_net, srf_net, tas = read_gregory(fil)
    toa_ok = np.mean(toa_net[1:])
    err_toa_ok = np.std(toa_net[1:])

    return toa_ok, err_toa_ok


def read_PI(tablecart, expnam, amip = True):
    """
    Reads the PI index for exps.
    """

    cart = tablecart + expnam + '/'
    allfils = [fi for fi in os.listdir(cart) if fi[:3] == 'PI2']

    pis = []
    anni = []
    for fi in allfils:
        with open(cart + fi, 'r') as filoz:
            linee = filoz.readlines()
            anno = int(linee[2].split()[-1])
            if amip:
                pi = float(linee[-1].split()[-1])
            else:
                pi = float(linee[-2].split()[-1])
            anni.append(anno)
            pis.append(pi)

    gigi = np.argsort(anni)
    pis = np.array(pis)[gigi]

    return np.mean(pis[1:])


def check_file(filnam):
    if os.path.exists(filnam):
        return True
    else:
        return False

def order_increasing(arr, *args):
    """
    Orders all args (and arr) according to the increasing order of arr.
    """
    nuord = np.argsort(arr)
    nuarr = np.array(arr)[nuord]

    nulist = []
    nulist.append(nuarr)

    for co in args:
        nulist.append(np.array(co)[nuord])

    return nulist

def check_increasing(arr):
    return np.all(np.diff(arr) > 0)

def check_decreasing(arr):
    return np.all(np.diff(arr) < 0)
