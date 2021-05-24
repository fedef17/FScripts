#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
import glob
from matplotlib import pyplot as plt
from matplotlib import cm

import pickle
import netCDF4 as nc

import climtools_lib as ctl
import climdiags as cd

from matplotlib.colors import LogNorm
from datetime import datetime

from scipy import stats

from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.metrics import r2_score

import itertools as itt
#import pymannkendall as mk

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 18
#############################################################################

if os.uname()[1] == 'hobbes':
    cart_in = '/home/fabiano/Research/lavori/KJ_regimes/'
else:
    cart_in = '/home/fedef/Research/lavori/KJ_regimes/'

# leggo metriche
gigi = pickle.load(open(cart_in + 'metrics/cmip6hist_data_dict.pkl', 'rb'))

# gigi['ACCESS-CM2_historical_r1i1p1f1'].keys()
# Out[5]: dict_keys(['U2', 'U4', 'R3', 'jet'])
#
# In [6]: gigi['ACCESS-CM2_historical_r1i1p1f1']['U4'].keys()
# Out[6]: dict_keys(['mean', 'BLK', 'NAO-', 'NAO+', 'AR'])
#
# In [7]: gigi['ACCESS-CM2_historical_r1i1p1f1']['U4']['BLK'].keys()
# Out[7]: dict_keys(['era20c_fidelity', 'era5_fidelity', 'stability', 'regime_occurrence', 'regime_persistence', 'variance_ratio'])
#
# In [8]: gigi['ACCESS-CM2_historical_r1i1p1f1']['U4']['mean'].keys()
# Out[8]: dict_keys(['era20c_fidelity', 'era5_fidelity', 'stability', 'regime_occurrence', 'regime_persistence', 'variance_ratio'])
models = list(gigi.keys())
models_ok = [mod.split('_')[0] + '_' + mod.split('_')[-1] for mod in models]

regtip = ['R3', 'U4']
metrics = dict()
#weights = dict()
metrics['R3'] = [('mean', 'era20c_fidelity'), ('mean', 'stability'), ('mean', 'regime_persistence'), ('mean', 'variance_ratio')] + [(rnam, 'regime_occurrence') for rnam in ['AR', 'NAO-', 'BLK']] + [(rnam, 'regime_persistence') for rnam in ['AR', 'NAO-', 'BLK']]
#weights['R3'] = np.append(np.ones(4), np.ones(6)/3.)

metrics['U4'] = [('mean', 'era20c_fidelity'), ('mean', 'stability'), ('mean', 'regime_persistence'), ('mean', 'variance_ratio')] + [(rnam, 'regime_occurrence') for rnam in ['AR', 'NAO-', 'BLK', 'NAO+']] + [(rnam, 'regime_persistence') for rnam in ['AR', 'NAO-', 'BLK', 'NAO+']]
#weights['U4'] = np.append(np.ones(4), np.ones(8)/4.)

# leggo drivers
drivs = os.listdir(cart_in + 'drivers/')
drivs = ['NorthAtlantic', 'Arctic', 'Stratosphere', 'JetSpeed', 'Tropics',
 'ENSO', 'EddyForcing'] # excluding: 'JetLatPDF'

dridic = dict()
for dri in drivs:
    cart_dri = cart_in + 'drivers/' + dri +'/'
    allfis = glob.glob(cart_dri + '*_ACCESS-CM2_r1*')
    drivtips = [fi.split('/')[-1].split('_ACCESS-CM2')[0] for fi in allfis]
    dridic[dri] = drivtips

# best 3 drivers for EAT
# varall = ['GW', 'PVS', 'UTWrGW', 'AArGW', 'NAWrGW', 'PST']
# pio2 = list(itt.combinations(varall, 3))

allscores = dict()

# faccio cose
regr_models = dict()
for reg in regtip:
    print(reg)
    allscores[reg] = []

    # Construct the dependent variables
    Y = []
    for mod in models:
        Y.append([gigi[mod][reg][metr[0]][metr[1]] for metr in metrics[reg]])
    Y = np.stack(Y)
    # standardizzo le Y?
    # scaler = StandardScaler().fit(Y)
    # Y = scaler.transform(Y)

    #sklearn.metrics.r2_score(y_true, y_pred, sample_weight)

    # pesca i drivers, uno per tipo, e fai il model
    allinds = [np.arange(len(dridic[dri])) for dri in drivs]
    allcombs = list(itt.product(*allinds))

    for ii, comb in enumerate(allcombs):
        print(ii)
        X = []
        for mod in models_ok:
            Xmod = []
            for dri, co in zip(drivs, comb):
                try:
                    filo = glob.glob(cart_in + 'drivers/' + dri + '/' + dridic[dri][co] + '*' + mod + '*')[0]
                    kos = np.loadtxt(filo)
                    if kos.size > 1:
                        kos = kos[0]
                    else:
                        kos = float(kos)
                except:
                    print('Not found: ', dri, dridic[dri][co], mod)
                    kos = np.nan

                Xmod.append(kos)

            X.append(Xmod)

        X = np.stack(X)
        if np.any(np.all(np.isnan(X), axis = 0)):
            baubo = np.array([dri +' '+ dridic[dri][co] for dri, co in zip(drivs, comb)])[np.where(np.all(np.isnan(X), axis = 0))[0]]
            print('All nans: ', baubo)
            continue

        if np.any(np.sum(np.isnan(X), axis = 0) > 3):
            baubo = np.array([dri +' '+ dridic[dri][co] for dri, co in zip(drivs, comb)])[np.where(np.sum(np.isnan(X), axis = 0) > 3)[0]]
            print('Too many nans: ', baubo)
            continue

        if np.any(np.isnan(X)):
            print('Replacing {} Nans in features with MMmean'.format(np.sum(np.isnan(X))))
            imputer = SimpleImputer()
            X = imputer.fit_transform(X)

        # STANDARDIZZO LE FEATURES
        scaler = StandardScaler().fit(X)
        X = scaler.transform(X)

        model = LinearRegression().fit(X, Y)
        #allscores[reg].append(model.score(X, Y))
        allscores[reg].append(r2_score(Y, model.predict(X), multioutput='raw_values'))
        #print(pi, model.score(X, Y))

pickle.dump(allscores, open(cart_in + 'allscores_7drivs.p', 'wb'))
# fine.
