#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
import glob
from matplotlib import pyplot as plt
from matplotlib import cm, colors

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
from sklearn.feature_selection import f_regression
import statsmodels.api as sm

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
metrnam = dict()
#weights = dict()
metrics['R3'] = [('mean', 'era20c_fidelity'), ('mean', 'stability'), ('mean', 'variance_ratio')] + [(rnam, 'regime_occurrence') for rnam in ['AR', 'NAO-', 'BLK']] + [(rnam, 'regime_persistence') for rnam in ['AR', 'NAO-', 'BLK']]
metrnam['R3'] = ['fidel20', 'stabil', 'var_ratio'] + ['occ: '+rnam for rnam in ['AR', 'NAO-', 'BLK']] + ['pers: '+rnam for rnam in ['AR', 'NAO-', 'BLK']]
#weights['R3'] = np.append(np.ones(4), np.ones(6)/3.)

metrics['U4'] = [('mean', 'era20c_fidelity'), ('mean', 'stability'), ('mean', 'variance_ratio')] + [(rnam, 'regime_occurrence') for rnam in ['AR', 'NAO-', 'BLK', 'NAO+']] + [(rnam, 'regime_persistence') for rnam in ['AR', 'NAO-', 'BLK', 'NAO+']]
metrnam['U4'] = ['fidel20', 'stabil', 'var_ratio'] + ['occ: '+rnam for rnam in ['AR', 'NAO-', 'BLK', 'NAO+']] + ['pers: '+rnam for rnam in ['AR', 'NAO-', 'BLK', 'NAO+']]
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

alldrivs = dict()
for mod in models_ok:
    for dri in drivs:
        for dritip in dridic[dri]:
            try:
                filo = glob.glob(cart_in + 'drivers/' + dri + '/' + dritip + '*' + mod + '*')[0]
                kos = np.loadtxt(filo)
                if kos.size > 1:
                    kos = kos[0]
                else:
                    kos = float(kos)
            except:
                print('Not found: ', dri, dritip, mod)
                kos = np.nan

            alldrivs[(dri, dritip, mod)] = kos

# resolutions
# nam, resi = ctl.read_from_txt(cart_in + 'resolutions.txt', n_skip=1, sep = ',')
# namok = [na.strip("'") for na in nam]
# namok = [nam.split('_')[0] + '_' + nam.split('_')[2] for nam in namok]

allinds = [np.arange(len(dridic[dri])) for dri in drivs]
allcombs = list(itt.product(*allinds))

okcombs = dict()
allscores = dict()

#faccio cose
for reg in regtip:
    print(reg)
    allscores[reg] = []
    okcombs[reg] = []

    # Construct the dependent variables
    Y = []
    for mod in models:
        Y.append([gigi[mod][reg][metr[0]][metr[1]] for metr in metrics[reg]])
    Y = np.stack(Y)
    # standardizzo le Y?
    scaler = StandardScaler().fit(Y)
    Y = scaler.transform(Y)

    # pesca i drivers, uno per tipo, e fai il model
    for ii, comb in enumerate(allcombs):
        print(ii)
        X = []
        for mod in models_ok:
            Xmod = [alldrivs[(dri, dridic[dri][co], mod)] for dri, co in zip(drivs, comb)]
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
        okcombs[reg].append(comb)
        #print(pi, model.score(X, Y))

pickle.dump(allscores, open(cart_in + 'allscores_7drivs.p', 'wb'))
# fine.


def make_XY(reg, comb, standard_Y = True):
    Y = []
    for mod in models:
        Y.append([gigi[mod][reg][metr[0]][metr[1]] for metr in metrics[reg]])
    Y = np.stack(Y)

    # standardizzo le Y?
    if standard_Y:
        scaler = StandardScaler().fit(Y)
        Y = scaler.transform(Y)

    # pesca i drivers, uno per tipo, e fai il model
    X = []
    for mod in models_ok:
        Xmod = [alldrivs[(dri, dridic[dri][co], mod)] for dri, co in zip(drivs, comb)]
        X.append(Xmod)

    X = np.stack(X)
    if np.any(np.isnan(X)):
        print('Replacing {} Nans in features with MMmean'.format(np.sum(np.isnan(X))))
        imputer = SimpleImputer()
        X = imputer.fit_transform(X)

    # STANDARDIZZO LE FEATURES
    scaler = StandardScaler().fit(X)
    X = scaler.transform(X)

    return X, Y

allscores = pickle.load(open(cart_in + 'allscores_7drivs.p', 'rb'))

colo = '#a50026 #d73027 #f46d43 #fdae61 #fee090 #ffffff #e0f3f8 #abd9e9 #74add1 #4575b4 #313695'
colo = colo.split()
colo = colo[::-1]
cmappa = colors.ListedColormap(colo)
cmappa.set_over('#800026') #662506
cmappa.set_under('#023858') #542788

for reg in regtip:
    ctl.printsep()
    ctl.printsep()
    print(reg)
    ire = int(reg[-1])
    allscores[reg] = np.stack(allscores[reg])
    okysing = np.argmax(allscores[reg], axis = 0)
    okyall = np.argmax(np.mean(allscores[reg], axis = 1))
    okymean = np.argmax(np.mean(allscores[reg][:, :ire], axis = 1))
    okyregs = np.argmax(np.mean(allscores[reg][:, ire:], axis = 1))

    print(okysing)
    print(okyall, okymean, okyregs)

    # print('single pred.')
    # for ii, (met, oky) in enumerate(zip(metrics[reg], okysing)):
    #     ctl.printsep()
    #     comb = okcombs[reg][oky]
    #     top_comb = [dri +' '+ dridic[dri][co] for dri, co in zip(drivs, comb)]
    #     print(met, allscores[reg][oky][ii])
    #     print(top_comb)

    ctl.printsep()
    fig_score, axs = plt.subplots(figsize=(12,8))

    for indcomb, namti, col in zip([okyall, okymean, okyregs], ['all', 'pattmean', 'occpers'], ['indianred', 'forestgreen', 'orange']):
        scoall = allscores[reg][indcomb]
        print('all', np.mean(scoall), np.min(scoall), np.max(scoall))
        comb = okcombs[reg][indcomb]
        top_comb = [dridic[dri][co] for dri, co in zip(drivs, comb)] #[dri +' '+
        print(top_comb)
        X, Y = make_XY(reg, comb)

        Xgi = sm.add_constant(X)
        pvals = []
        params = []
        rsq = []
        for ko in Y.T:
            est = sm.OLS(ko, Xgi)
            est2 = est.fit()
            pvals.append(est2.pvalues)
            params.append(est2.params)
            rsq.append(est2.rsquared)

        pvals = np.stack(pvals)
        params = np.stack(params)
        rsq = np.stack(rsq)
        pvals = pvals[:, 1:]
        params = params[:, 1:]

        skmodel = LinearRegression().fit(X, Y)
        print('CHECK!!')
        #print(params)
        print(np.max(np.abs((params-skmodel.coef_)/params)))
        print(rsq)
        print(r2_score(Y, skmodel.predict(X), multioutput='raw_values'))

        # pvals_skl = []
        # for ko in Y.T:
        #     F, pvals = f_regression(X, ko)
        #     pvals_skl.append(ko)
        #
        # print(pvals_skl)

        fig, ax = plt.subplots(figsize=(16,12))

        ndriv = len(comb)
        nval = len(scoall)
        vmi = np.percentile(np.abs(params), 95)
        ext = [0, nval, 0, ndriv]

        gigifig = ax.imshow(params.T, vmin = -vmi, vmax = vmi, cmap = cmappa, origin = 'lower',  extent = ext, aspect = 1)

        pvaok = pvals.T
        for ix in range(nval):
            for iy in range(ndriv):
                if pvaok[iy, ix] < 0.01:
                    ax.scatter(ix+0.5, iy+0.5, s=60, edgecolors = 'black', facecolors='white')
                elif pvaok[iy, ix] < 0.05:
                    ax.scatter(ix+0.5, iy+0.5, s=20, edgecolors = 'black', facecolors='white')

        ax.xaxis.tick_top()
        ax.set_xticks(0.5+np.arange(len(metrnam[reg])), minor = False)
        ax.set_xticklabels(metrnam[reg], ha='center', rotation = 30)
        ax.set_yticks(0.5+np.arange(len(top_comb)), minor = False)
        ax.set_yticklabels(top_comb, va='center', rotation = 30)

        for co in np.arange(len(metrnam[reg])):
            ax.axvline(co, color = 'white', linewidth = 0.1)
        for co in np.arange(len(top_comb)):
            ax.axhline(co, color = 'white', linewidth = 0.1)


        #cax = fig.add_subplot(gs[6, :])
        cax = plt.axes([0.1, 0.1, 0.8, 0.05])
        cb = plt.colorbar(gigifig, cax=cax, orientation='horizontal', extend = 'both')
        cb.ax.tick_params(labelsize=18)
        cb.set_label('Regression coefficient', fontsize=20)
        plt.subplots_adjust(left=0.1, bottom=0.17, right=0.98, top=0.92, wspace=0.05, hspace=0.20)

        fig.savefig(cart_in + 'Sm_{}_{}_7driv.pdf'.format(reg, namti))

        axs.plot(np.arange(len(rsq)), rsq, label = namti, color = col)
        axs.plot(np.arange(len(rsq)), scoall, color = col, linestyle = '--')
        print(reg, namti, np.mean(rsq[:ire]), np.mean(rsq[ire:]))
        print(reg, namti, np.mean(scoall[:ire]), np.mean(scoall[ire:]))
        #axs.scatter(np.arange(len(rsq)), rsq, color = col)

        ctl.printsep()

    axs.set_xticks(np.arange(len(rsq)), minor = False)
    axs.set_xticklabels(metrnam[reg], ha='center', rotation = 30)
    axs.legend()
    axs.set_ylabel(r'$R^2$')
    fig_score.savefig(cart_in + 'Rsquared_{}_7driv.pdf'.format(reg))
