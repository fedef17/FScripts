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
from copy import deepcopy as dcp

from sklearn.linear_model import LinearRegression, Ridge
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
#gigi = pickle.load(open(cart_in + 'metrics/cmip6hist_data_dict.pkl', 'rb'))
gigi_old = pickle.load(open(cart_in + 'metrics/cmip6hist_data_dict.pkl', 'rb'))

gigi = pickle.load(open(cart_in + 'metrics/cmip6_r3n_metrics.pkl', 'rb'))

for mod in gigi.keys():
    gigi[mod].update(gigi_old[mod])
    gigi[mod]['R3n']['mean'] = gigi_old[mod]['R3']['mean']

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
models_cm6_full = list(gigi.keys())
models_cmip6 = [mod.split('_')[0] + '_' + mod.split('_')[-1] for mod in models_cm6_full]

gigiprim_old = pickle.load(open(cart_in + 'metrics/primavera_metrics.pkl', 'rb'))
gigiprim = pickle.load(open(cart_in + 'metrics/prmva_r3n_metrics.pkl', 'rb'))

#models_prim = list(gigiprim.keys())
models_prim = [mo for mo in gigiprim.keys() if 'r1' in mo] # only first member
for mod in models_prim:
    gigiprim[mod].update(gigiprim_old[mod])
    gigiprim[mod]['R3n']['mean'] = gigiprim_old[mod]['R3']['mean']

gigi.update(gigiprim)

models_prim = ['CMCC-CM2-HR4_r1i1p1f1', 'CMCC-CM2-VHR4_r1i1p1f1', 'CNRM-CM6-1-HR_r1i1p1f2', 'CNRM-CM6-1_r1i1p1f2', 'EC-Earth3P-HR_r1i1p2f1', 'EC-Earth3P_r1i1p2f1', 'ECMWF-IFS-HR_r1i1p1f1', 'ECMWF-IFS-LR_r1i1p1f1', 'ECMWF-IFS-MR_r1i1p1f1', 'HadGEM3-GC31-HH_r1i1p1f1', 'HadGEM3-GC31-HM_r1i1p1f1', 'HadGEM3-GC31-LL_r1i2p1f1', 'HadGEM3-GC31-MM_r1i1p1f1', 'MPI-ESM1-2-HR_r1i1p1f1', 'MPI-ESM1-2-XR_r1i1p1f1']

gigi5_old = pickle.load(open(cart_in + 'metrics/cmip5_metrics.pkl', 'rb'))
gigi5 = pickle.load(open(cart_in + 'metrics/cmip5_r3n_metrics.pkl', 'rb'))

models_5 = [mo for mo in gigi5.keys() if 'r1' in mo] # only first member
models_5.append('CCSM4-historical-r6i1p1')

models_5_ok = [mod.split('_')[0] + '_' + mod.split('_')[-1] for mod in models_5]

for mod in models_5:
    gigi5[mod].update(gigi5_old[mod])
    gigi5[mod]['R3n']['mean'] = gigi5_old[mod]['R3']['mean']

gigi.update(gigi5)

resolution_file = cart_in + 'drivers/resolutions.txt'
cose = ctl.read_from_txt(resolution_file, n_skip = 1, dtype = float, first_column_header = True, sep = ',')
cose[1][cose[1] < 0] = np.nan
cosdi = dict(zip(cose[0], cose[1]))

models_ok_cm6prim = models_cmip6 + models_prim
models_cm6prim = models_cm6_full + models_prim

models_ok_all = models_cmip6 + models_prim + models_5_ok
models_all = models_cm6_full + models_prim + models_5

cart_out_gen = cart_in + 'results_{}/'

tuttecose = dict()

regtip = ['U4', 'R3', 'R3n']
metrics = dict()
metrnam = dict()
#weights = dict()

metrics['R3'] = [('mean', 'era20c_fidelity'), ('mean', 'stability'), ('mean', 'variance_ratio')] + [(rnam, 'regime_occurrence') for rnam in ['AR', 'NAO-', 'BLK']] + [(rnam, 'regime_persistence') for rnam in ['AR', 'NAO-', 'BLK']]
metrnam['R3'] = ['fidel20', 'stabil', 'var_ratio'] + ['occ: '+rnam for rnam in ['AR', 'NAO-', 'BLK']] + ['pers: '+rnam for rnam in ['AR', 'NAO-', 'BLK']]
#weights['R3'] = np.append(np.ones(4), np.ones(6)/3.)

metrics['U4'] = [('mean', 'era20c_fidelity'), ('mean', 'stability'), ('mean', 'variance_ratio')] + [(rnam, 'regime_occurrence') for rnam in ['AR', 'NAO-', 'BLK', 'NAO+']] + [(rnam, 'regime_persistence') for rnam in ['AR', 'NAO-', 'BLK', 'NAO+']]
metrnam['U4'] = ['fidel20', 'stabil', 'var_ratio'] + ['occ: '+rnam for rnam in ['AR', 'NAO-', 'BLK', 'NAO+']] + ['pers: '+rnam for rnam in ['AR', 'NAO-', 'BLK', 'NAO+']]

metrics['R3n'] = [('mean', 'era20c_fidelity'), ('mean', 'stability'), ('mean', 'variance_ratio')] + [(rnam, 'regime_occurrence') for rnam in ['AR', 'NAO-', 'BLK', 'Neutral']] + [(rnam, 'regime_persistence') for rnam in ['AR', 'NAO-', 'BLK', 'Neutral']]
metrnam['R3n'] = ['fidel20', 'stabil', 'var_ratio'] + ['occ: '+rnam for rnam in ['AR', 'NAO-', 'BLK', 'Neutral']] + ['pers: '+rnam for rnam in ['AR', 'NAO-', 'BLK', 'Neutral']]

# leggo drivers
drivs = os.listdir(cart_in + 'simple_drivers_1/')
#drivs = ['NorthAtlantic', 'Arctic', 'Stratosphere', 'JetSpeed', 'Tropics', 'ENSO', 'EddyForcing'] # excluding: 'JetLatPDF'
#drivs = ['NorthAtlantic', 'Arctic', 'Stratosphere', 'JetSpeed', 'Tropics', 'EddyForcing'] # excluding: 'JetLatPDF'
drivs = ['NorthAtlantic', 'Arctic', 'Stratosphere', 'JetSpeed', 'EddyForcing'] # excluding: 'JetLatPDF'

dridic = dict()
for dri in drivs:
    cart_dri = cart_in + 'simple_drivers_1/' + dri +'/'
    allfis = glob.glob(cart_dri + '*_ACCESS-CM2_r1*')
    drivtips = [fi.split('/')[-1].split('_ACCESS-CM2')[0] for fi in allfis]
    dridic[dri] = drivtips

# best 3 drivers for EAT
# varall = ['GW', 'PVS', 'UTWrGW', 'AArGW', 'NAWrGW', 'PST']
# pio2 = list(itt.combinations(varall, 3))

n_exc = 0
exc_list = []

alldrivs = dict()
for mod in models_ok_all:
    for dri in drivs:
        for dritip in dridic[dri]:
            try:
                if mod in models_cmip6:
                    filo = glob.glob(cart_in + 'simple_drivers_1/' + dri + '/' + dritip + '*' + mod + '*')[0]
                elif mod in models_prim:
                    filo = glob.glob(cart_in + 'simple_drivers_1/' + dri + '/PRIMAVERA/' + dritip + '*' + mod + '*')[0]
                elif mod in models_5_ok:
                    filo = glob.glob(cart_in + 'simple_drivers_1/' + dri + '/CMIP5/' + dritip + '*' + mod + '*')[0]
                kos = np.loadtxt(filo)
                if kos.size > 1:
                    kos = kos[0]
                else:
                    kos = float(kos)
            except:
                # print('Not found: ', dri, dritip, mod)
                n_exc += 1
                exc_list.append(dritip)
                kos = np.nan

            alldrivs[(dri, dritip, mod)] = kos

    if mod in models_5_ok:
        alldrivs[('resolution', 'atm_res', mod)] = cosdi[mod[:mod.find('historical')-1]][0]
        alldrivs[('resolution', 'atm_lev', mod)] = cosdi[mod[:mod.find('historical')-1]][1]
        alldrivs[('resolution', 'oce_res', mod)] = cosdi[mod[:mod.find('historical')-1]][3]
    else:
        alldrivs[('resolution', 'atm_res', mod)] = cosdi[mod][0]
        alldrivs[('resolution', 'atm_lev', mod)] = cosdi[mod][1]
        alldrivs[('resolution', 'oce_res', mod)] = cosdi[mod][3]

ctl.printsep()

print('Found {} exceptions:\n'.format(n_exc))
for dri in np.unique(exc_list):
    print(dri, len([dr for dr in exc_list if dr == dri]))

ctl.printsep()

drilis_all = []
for ke in dridic:
    for pu in dridic[ke]:
        drilis_all.append((ke, pu))

drilis_all.append(('resolution', 'atm_res'))

tuttecose['drivers'] = alldrivs
tuttecose['drivnam'] = drilis_all

models_ok_all.remove('HadGEM3-GC31-LL_r1i1p1f3')
models_ok_all.remove('HadGEM3-GC31-MM_r1i1p1f3')
models_all.remove('HadGEM3-GC31-LL_historical_r1i1p1f3')
models_all.remove('HadGEM3-GC31-MM_historical_r1i1p1f3')

for reg in regtip:
    for models, ensmod in zip([models_cm6_full, models_prim, models_5, models_cm6prim, models_all], ['solocmip6', 'soloprim', 'solocmip5', 'cmip6prim', 'all']):
        Y = []
        for mod in models:
            Y.append([gigi[mod][reg][metr[0]][metr[1]] for metr in metrics[reg]])
        Y = np.stack(Y)
        tuttecose[(ensmod, reg, 'metrics')] = Y
        tuttecose[(ensmod, reg, 'metrnam')] = metrics[reg]


def make_XY(reg, comb, models = models_all, models_ok = models_ok_all, alldrivs = alldrivs, drilis = drilis_all, metrics = metrics, standard_Y = True):
    Y = []
    for mod in models:
        Y.append([gigi[mod][reg][metr[0]][metr[1]] for metr in metrics[reg]])
    Y = np.stack(Y)
    if np.any(np.isnan(Y)):
        #print('!!! Replacing {} Nans in metrics!! with MMmean'.format(np.sum(np.isnan(Y))))
        imputer = SimpleImputer()
        Y = imputer.fit_transform(Y)

    # standardizzo le Y?
    if standard_Y:
        scaler = StandardScaler().fit(Y)
        Y = scaler.transform(Y)

    # pesca i drivers, uno per tipo, e fai il model
    X = []
    for mod in models_ok:
        Xmod = [alldrivs[(drilis[co][0], drilis[co][1], mod)] for co in comb]
        X.append(Xmod)

    X = np.stack(X)
    if np.any(np.isnan(X)):
        #print('Replacing {} Nans in features with MMmean'.format(np.sum(np.isnan(X))))
        imputer = SimpleImputer()
        X = imputer.fit_transform(X)

    #print('Mean and std before standardizing:')
    #print(np.mean(X, axis = 0))
    #print(np.std(X, axis = 0))

    # STANDARDIZZO LE FEATURES
    scaler = StandardScaler().fit(X)
    X = scaler.transform(X)

    # print('Cross correlation after std:')
    # print(np.cov(X.T)/np.cov(X.T)[0,0], '\n')

    return X, Y


for models_ok, models, ensmod in zip([models_cmip6, models_prim, models_5_ok, models_ok_cm6prim, models_ok_all], [models_cm6_full, models_prim, models_5, models_cm6prim, models_all], ['solocmip6', 'soloprim', 'solocmip5', 'cmip6prim', 'all']):

    drilis = dcp(drilis_all)
    if ensmod not in ['solocmip5', 'all']:
        drilis.append(('resolution', 'atm_lev'))
        drilis.append(('resolution', 'oce_res'))

    ctl.printsep()
    print(ensmod)
    ctl.printsep()

    #weights['R3'] = np.append(np.ones(4), np.ones(6)/3.)

    # resolutions
    # nam, resi = ctl.read_from_txt(cart_in + 'resolutions.txt', n_skip=1, sep = ',')
    # namok = [na.strip("'") for na in nam]
    # namok = [nam.split('_')[0] + '_' + nam.split('_')[2] for nam in namok]

    # allinds = [np.arange(len(dridic[dri])) for dri in drivs]
    # allcombs = list(itt.product(*allinds))


    ctl.mkdir(cart_out_gen.format(ensmod))

    for nu in np.arange(4, 11):
        print('Finding best {} drivers'.format(nu))
        allcombs = list(itt.combinations(range(len(drilis)), nu))

        okcombs = dict()
        allscores = dict()

        #faccio cose
        for reg in regtip:
            print(reg)
            allscores[reg] = []
            okcombs[reg] = []

            # pesca i drivers, uno per tipo, e fai il model
            for ii, comb in enumerate(allcombs):
                # print(ii, comb)
                X, Y = make_XY(reg, comb, models, models_ok, drilis = drilis)

                #model = LinearRegression().fit(X, Y)
                model = Ridge().fit(X, Y)

                #allscores[reg].append(model.score(X, Y))
                allscores[reg].append(r2_score(Y, model.predict(X), multioutput='raw_values'))
                okcombs[reg].append(comb)
                #print(pi, model.score(X, Y))

        pickle.dump([okcombs, allscores], open(cart_out_gen.format(ensmod) + 'allscores_v2_{}drivs_{}.p'.format(nu, ensmod), 'wb'))
        # fine.

        okcombs, allscores = pickle.load(open(cart_out_gen.format(ensmod) + 'allscores_v2_{}drivs_{}.p'.format(nu, ensmod), 'rb'))

        colo = '#a50026 #d73027 #f46d43 #fdae61 #fee090 #ffffff #e0f3f8 #abd9e9 #74add1 #4575b4 #313695'
        colo = colo.split()
        colo = colo[::-1]
        cmappa = colors.ListedColormap(colo)
        cmappa.set_over('#800026') #662506
        cmappa.set_under('#023858') #542788

        for reg in regtip:
            cart_out = cart_out_gen.format(ensmod) + reg + '/'
            ctl.mkdir(cart_out)

            # ctl.printsep()
            # ctl.printsep()
            print(reg)
            ire = int(reg[1])
            if reg == 'R3':
                ire = 3
            else:
                ire = 4
            allscores[reg] = np.stack(allscores[reg])
            okysing = np.argmax(allscores[reg], axis = 0)
            okyall = np.argmax(np.mean(allscores[reg], axis = 1))
            okymean = np.argmax(np.mean(allscores[reg][:, :ire], axis = 1))
            okyregs = np.argmax(np.mean(allscores[reg][:, ire:], axis = 1))

            # print(okysing)
            # print(okyall, okymean, okyregs)

            # print('single pred.')
            # for ii, (met, oky) in enumerate(zip(metrics[reg], okysing)):
            #     ctl.printsep()
            #     comb = okcombs[reg][oky]
            #     top_comb = [dri +' '+ dridic[dri][co] for dri, co in zip(drivs, comb)]
            #     print(met, allscores[reg][oky][ii])
            #     print(top_comb)

            # ctl.printsep()
            fig_score, axs = plt.subplots(figsize=(12,8))

            for indcomb, namti, col, linst in zip([okyall, okymean, okyregs], ['all', 'pattmean', 'occpers'], ['indianred', 'forestgreen', 'orange'], ['-', '--', ':']):
                scoall = allscores[reg][indcomb]
                #print('all', np.mean(scoall), np.min(scoall), np.max(scoall))
                comb = okcombs[reg][indcomb]
                top_comb = [drilis[co][1] for co in comb] #[dri +' '+
                #print(namti, top_comb)
                X, Y = make_XY(reg, comb, models, models_ok, drilis = drilis)

                crosscor = np.cov(X.T)/np.cov(X.T)[0,0]

                ## cross correlation after std
                fig, ax = plt.subplots(figsize=(16,12))

                ndriv = len(comb)
                vmi = 0.6
                ext = [0, ndriv, 0, ndriv]

                gigifig = ax.imshow(crosscor, vmin = -vmi, vmax = vmi, cmap = cmappa, origin = 'upper',  extent = ext, aspect = 1)

                ax.xaxis.tick_top()
                ax.set_xticks(0.5+np.arange(len(top_comb)), minor = False)
                ax.set_xticklabels(top_comb, ha='center', rotation = 20)
                ax.set_yticks(0.5+np.arange(len(top_comb)), minor = False)
                ax.set_yticklabels(top_comb[::-1], va='center', rotation = 30)

                #cax = fig.add_subplot(gs[6, :])
                cax = plt.axes([0.1, 0.1, 0.8, 0.05])
                cb = plt.colorbar(gigifig, cax=cax, orientation='horizontal', extend = 'both')
                cb.ax.tick_params(labelsize=18)
                cb.set_label('Cross correlation', fontsize=20)
                plt.subplots_adjust(left=0.1, bottom=0.17, right=0.98, top=0.92, wspace=0.05, hspace=0.20)

                fig.savefig(cart_out + 'Crosscorr_{}_{}_v2_{}driv_{}.pdf'.format(reg, namti, nu, ensmod))


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
                # print('CHECK!!')
                # #print(params)
                # print(np.max(np.abs((params-skmodel.coef_)/params)))
                # print(rsq)
                # print(r2_score(Y, skmodel.predict(X), multioutput='raw_values'))

                # pvals_skl = []
                # for ko in Y.T:
                #     F, pvals = f_regression(X, ko)
                #     pvals_skl.append(ko)
                #
                # print(pvals_skl)

                fig, ax = plt.subplots(figsize=(16,12))

                ndriv = len(comb)
                nval = len(scoall)
                #vmi = np.percentile(np.abs(params), 95)
                vmi = 0.6
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
                ax.set_yticklabels(top_comb, va='center', rotation = 45)

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

                fig.savefig(cart_out + 'Sm_{}_{}_v2_{}driv_{}.pdf'.format(reg, namti, nu, ensmod))

                axs.plot(np.arange(len(rsq)), rsq, label = namti, color = col, linestyle = linst)
                #axs.plot(np.arange(len(rsq)), scoall, color = col, linestyle = '--')

                # print(reg, namti, np.mean(rsq[:ire]), np.mean(rsq[ire:]))
                # print(reg, namti, np.mean(scoall[:ire]), np.mean(scoall[ire:]))
                #axs.scatter(np.arange(len(rsq)), rsq, color = col)

                ctl.printsep()

            axs.set_xticks(np.arange(len(rsq)), minor = False)
            axs.set_xticklabels(metrnam[reg], ha='center', rotation = 30)
            axs.legend()
            axs.set_ylabel(r'$R^2$')
            fig_score.savefig(cart_out + 'Rsquared_{}_v2_{}driv_{}.pdf'.format(reg, nu, ensmod))

pickle.dump(tuttecose, open(cart_out + 'tuttecose_wcmip5.p', 'wb'))
