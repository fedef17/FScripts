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
import xarray as xr

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 18

################################################################################

cart_in = '/data-hobbes/fabiano/TunECS/AMIP_exps/'
cart_out = '/home/fabiano/Research/lavori/TunECS/tuning/experiments/param_fingerprint/'
ctl.mkdir(cart_out)

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

allforc = ['pi', 'c4']
allchan = 'm n p q l r'.split()

def val_ok(param, change):
    iok = allchan.index(change)
    return valchange[param][iok]

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

#lats = [-90, -65, -40, -20, 20, 40, 65, 90]
lats = [-90, -60, -30, 30, 60, 90]
bands = [(la1, la2) for la1, la2 in zip(lats[:-1], lats[1:])]
lacen = np.array([np.mean(laol) for laol in bands])

allvars = ['ttr', 'tsr', 'tcc']#, 'str', 'ssr', 'sshf', 'slhf', 'tcc', 'cp', 'lsp']

# tsr: This parameter is the incoming solar radiation (also known as shortwave radiation) minus the outgoing solar radiation at the top of the atmosphere. It is the amount of radiation passing through a horizontal plane. The incoming solar radiation is the amount received from the Sun. The outgoing solar radiation is the amount reflected and scattered by the Earth's atmosphere and surface.
# ttr: The thermal (also known as terrestrial or longwave) radiation emitted to space at the top of the atmosphere is commonly known as the Outgoing Longwave Radiation (OLR). The top net thermal radiation (this parameter) is equal to the negative of OLR.
# tsr + ttr = tnr. Positive -> incoming! (downward)
# ssr: This parameter is the amount of solar radiation (also known as shortwave radiation) that reaches a horizontal plane at the surface of the Earth (both direct and diffuse) minus the amount reflected by the Earth's surface (which is governed by the albedo).
# str: This parameter is the difference between downward and upward thermal radiation at the surface of the Earth. It the amount passing through a horizontal plane.
# ssr + str = snr. Positive -> downward
# MA!!!! ATTENZIONE!!!! snr del hiresclim non è solo la radiazione, ci sono anche i flussi di calore: sshf (Sensible heat flux) e slhf (Latent heat flux), che sono già netti e definiti POSITIVI verso il basso. The ECMWF convention for vertical fluxes is positive downwards.
# srf_net = ssr + str + sshf + slhf

# energia assorbita da atm: tnr-snr

# mapz = dict()
# mapz_err = dict()
# for forc in ['pi', 'c4']:
#     print(forc)
#     for varnam in allvars:
#         print(varnam)
#         for nu, let, param in zip(nums, letts, testparams):
#             cosy = []
#             parval = []
#             for iic, change in enumerate(['m', 'n', 'p', 'q', 'l', 'r']):
#                 mok = mname.format(forc, change, let)
#                 print(mok)
#
#                 if forc == 'pi' and let == 'h' and change in ['p', 'n']:
#                     listafil = [fil.format(mok, ye, mok, ye, varnam) for ye in range(1851, 1860)] # ten years for these, instead of 5
#                 else:
#                     listafil = [fil.format(mok, ye, mok, ye, varnam) for ye in range(1851, 1855)] # skipping the first year
#
#                 if tl.check_file(listafil[0]):
#                     if not tl.check_file(listafil[-1]): # Many jobs did not terminate!!
#                         listafil = listafil[:-1]
#
#                     gigi = xr.open_mfdataset(listafil, use_cftime = True, engine = 'netcdf4')
#                     gigok = gigi[varnam].groupby("time.year").mean().compute().values
#
#                     cosy.append(gigok)
#                     parval.append(len(gigok)*[val_ok(param, change)])
#
#             fields = np.concatenate(cosy, axis = 0)
#             relcha = np.concatenate(parval)/uff_params[param]
#
#             resu = ctl.calc_trend_climatevar(relcha, fields)
#             fgprint = resu[0]
#             fgprint_err = resu[2]
#
#             mapz[(forc, varnam, param)] = fgprint
#             mapz_err[(forc, varnam, param)] = fgprint_err
#
#
# for forc in ['pi', 'c4']:
#     for param in testparams:
#         mapz[(forc, 'toa_net', param)] = mapz[(forc, 'ttr', param)] + mapz[(forc, 'tsr', param)]
#         mapz_err[(forc, 'toa_net', param)] = mapz_err[(forc, 'ttr', param)] + mapz_err[(forc, 'tsr', param)]
#
# pickle.dump([mapz, mapz_err], open(cart_out + 'mapz_fgprint.p', 'wb'))

mapz, mapz_err = pickle.load(open(cart_out + 'mapz_fgprint.p', 'rb'))

allvars = ['ttr', 'tsr', 'tcc']
allvars = allvars + ['toa_net']
longnams = ['OLR (W/m2)', 'OSR (W/m2)', 'Cloud cover', 'Incoming Net TOA (W/m2)']
for forc in ['pi', 'c4']:
    for varnam, lonam in zip(allvars, longnams):
        mapes = [mapz[(forc, varnam, param)]/10 for param in testparams]
        if varnam == 'ttr':
            mapes = [-mapz[(forc, varnam, param)]/10 for param in testparams]
        elif varnam == 'tcc':
            mapes = [100.*mapz[(forc, varnam, param)]/10 for param in testparams]

        mapes_err = [np.abs(mapz[(forc, varnam, param)]) > mapz_err[(forc, varnam, param)] for param in testparams]

        fil = 'fgprint_{}_{}.pdf'.format(forc, varnam)
        ctl.plot_multimap_contour(mapes, gigi.lat, gigi.lon, filename = cart_out + fil, plot_anomalies = True, cbar_range = (-3, 3), fix_subplots_shape = (3,3), add_hatching = mapes_err, subtitles = testparams, cb_label = lonam)
