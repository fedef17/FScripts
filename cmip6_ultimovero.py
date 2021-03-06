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
import pymannkendall as mk

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 18
#############################################################################

cart_in = '/home/fedef/Research/lavori/CMIP6/'

cart_out_orig = cart_in + 'Results_ultimo_v2/'
ctl.mkdir(cart_out_orig)

yr10 = 10 # length of running mean

cart_cmip5 = '/home/fedef/Research/lavori/CMIP6/Results_cmip5/{}_NDJFM/'
cart_v5 = '/home/fedef/Research/lavori/CMIP6/Results_v5_rebase/{}_NDJFM/'
#cart_v5 = '/home/fedef/Research/lavori/CMIP6/Results_v5_ensrebase/{}_NDJFM/'

numclus = 4
reg_names_area = dict()
reg_names_area['EAT'] = ['NAO+', 'SBL', 'AR', 'NAO-']
reg_names_area['PNA'] = ['PT', 'PNA+', 'PNA-', 'BR']
#okmods = ['ACCESS-CM2_r1i1p1f1', 'BCC-CSM2-MR_r1i1p1f1', 'CanESM5_r1i1p1f1', 'CESM2-WACCM_r1i1p1f1','CNRM-CM6-1-HR_r1i1p1f2', 'CNRM-CM6-1_r1i1p1f2','CNRM-ESM2-1_r1i1p1f2', 'EC-Earth3_r1i1p1f1', 'FGOALS-g3_r1i1p1f1','INM-CM4-8_r1i1p1f1', 'INM-CM5-0_r1i1p1f1', 'IPSL-CM6A-LR_r1i1p1f1','MIROC6_r1i1p1f1', 'MPI-ESM1-2-HR_r1i1p1f1','MPI-ESM1-2-LR_r1i1p1f1', 'MRI-ESM2-0_r1i1p1f1','NorESM2-LM_r1i1p1f1', 'NorESM2-MM_r1i1p1f1', 'UKESM1-0-LL_r1i1p1f2']
# mancano = ['BCC-ESM1_r1i1p1f1', 'CESM2_r1i1p1f1', 'GFDL-CM4_r1i1p1f1', 'HadGEM3-GC31-LL_r1i1p1f3', 'KACE-1-0-G_r1i1p1f1', 'MPI-ESM-1-2-HAM_r1i1p1f1']

allssps = 'ssp126 ssp245 ssp370 ssp585'.split()
allsimcol = ['hist', 'ssp126', 'ssp245', 'bau', 'ssp370', 'ssp585', 'rcp85_cmip5']
coldic = dict(zip(allsimcol, ctl.color_set(7)))
colssp = [coldic[ssp] for ssp in allssps]

area = 'EAT'
ssp = 'ssp585'
cose = dict()
for area in ['EAT', 'PNA']:
    #results_hist_refEOF, results_ref = ctl.load_wrtool(file_hist_refEOF.format(area))
    trend_ssp, residtime_ssp = pickle.load(open(cart_v5.format(area) + 'trends_wrfreq_e_restime_{}.p'.format(area), 'rb'))
    seasfreq, runfreq = pickle.load(open(cart_v5.format(area) + 'seasfreqs_{}_v4.p'.format(area), 'rb'))
    freqs, residtimes, patterns = pickle.load(open(cart_v5.format(area) + 'allresults_dicts_{}_v3.p'.format(area), 'rb'))

    freqs_cmip5, trend_ssp_cmip5, residtimes_cmip5 = pickle.load(open(cart_cmip5.format(area) + 'freqs_cmip5_{}.p'.format(area), 'rb'))
    seasfreq_cmip5, runfreq_cmip5 = pickle.load(open(cart_cmip5.format(area) + 'seasfreqs_cmip5_{}.p'.format(area), 'rb'))
    freqs.update(freqs_cmip5)
    trend_ssp.update(trend_ssp_cmip5)
    residtimes.update(residtimes_cmip5)
    seasfreq.update(seasfreq_cmip5)
    runfreq.update(runfreq_cmip5)

    cose[(area, 'freq')] = freqs
    cose[(area, 'trend')] = trend_ssp
    cose[(area, 'residtimes')] = residtimes
    cose[(area, 'runfreq')] = runfreq
    cose[(area, 'seasfreq')] = seasfreq

    okmods = dict()
    okmods['ssp585'] = [ke[1] for ke in freqs if 'ssp585' in ke and 'tot50' in ke and 'all' not in ke and 'rel' not in ke]
    okmods['rcp85'] = [ke[1] for ke in freqs if 'rcp85_cmip5' in ke and 'tot50' in ke and 'all' not in ke and 'rel' not in ke]
    print(okmods)

    ### model performance
    # cose[(area, 'var_ratio')] = plocos[('var_ratio', ke, area)]
    # plocos[('freqbias', ke, area)]
    # var_ratio = [results_hist_refEOF[cos][mod]['var_ratio'] for mod in okmods[cos]]
    # cose[(area, 'var_ratio')] = dict(zip(okmods[cos], var_ratio))
    # patc = [np.mean(results_hist_refEOF[cos][mod]['patcor']) for mod in okmods[cos]]
    # cose[(area, 'patc')] = dict(zip(okmods[cos], patc))
    # freqbias = np.array([np.mean(np.abs(results_hist_refEOF[cos][mod]['freq_clus']-results_ref['freq_clus'])) for mod in okmods[cos]])
    # cose[(area, 'freqbias')] = dict(zip(okmods[cos], freqbias))

filbi = '/home/fedef/Research/lavori/CMIP6/Results_cmip6_vs_cmip5_refEOF/histbiases.p'
with open(filbi, 'rb') as pino:
    res_short = pickle.load(pino)
kesal = list(res_short.keys())
for ke in kesal:
    if 'cmip5_refEOF' in ke:
        res_short[tuple(['rcp85'] + list(ke[1:]))] = res_short[ke]
    if 'cmip6_refEOF' in ke:
        res_short[tuple(['ssp585'] + list(ke[1:]))] = res_short[ke]

# da Virna
vicar = '/home/fedef/Research/lavori/CMIP6/virnas_tas/new/'
fils = ['Indices_{}.txt', 'Slopes_{}_ndjfm.txt', 'Slopes_{}_temp.txt', 'Trends_NDJFM_{}.txt']
nams = ['temp', 'slope', 'slope', 'slope']
tip = [0,0,1,2]

resvi = dict()
for na, fi, ti in zip(nams, fils, tip):
    for ssp in ['rcp85', 'ssp585']:

        if ti == 0:
            vmods, datatemp = ctl.read_from_txt(vicar + fi.format(ssp), n_skip = 2, sep = '\t')
            vkeys = ['AA', 'UTW', 'PST', 'PVS']
            for mod in okmods[ssp]:
                if mod.split('_')[0] in vmods:
                    ind = np.where(vmods == mod.split('_')[0])[0][0]
                    lin = datatemp[ind]
                    for ke, te in zip(vkeys, lin):
                        resvi[(ssp, mod, na, ke)] = te
                    resvi[(ssp, mod, na, 'RUTAW')] = resvi[(ssp, mod, na, 'UTW')]/resvi[(ssp, mod, na, 'AA')]
        elif ti == 1:
            vmods, datatemp = ctl.read_from_txt(vicar + fi.format(ssp), n_skip = 2, sep = '\t')
            vkeys = ['GW', 'NAW_annual']
            for mod in okmods[ssp]:
                if mod.split('_')[0] in vmods:
                    ind = np.where(vmods == mod.split('_')[0])[0][0]
                    lin = datatemp[ind]
                    for ke, te in zip(vkeys, lin):
                        resvi[(ssp, mod, na, ke)] = te
        else:
            vmods, datatemp = ctl.read_from_txt(vicar + fi.format(ssp), n_skip = 3, sep = '\t')
            vkeys = ['GW_ndjfm', 'NAW', 'SPGW', 'TNAW1', 'TNAW']
            for mod in okmods[ssp]:
                if mod.split('_')[0] in vmods:
                    ind = np.where(vmods == mod.split('_')[0])[0][0]
                    lin = datatemp[ind]
                    for ke, te in zip(vkeys, lin):
                        resvi[(ssp, mod, na, ke)] = te

vkeys = ['AA', 'UTW', 'PST', 'PVS', 'RUTAW', 'GW', 'NAW', 'SPGW', 'TNAW1', 'TNAW', 'NAW_annual']

for ke in list(resvi.keys()):
    if 'temp' in ke: continue
    for kk in vkeys:
        if kk in ke:
            nuke = tuple(list(ke[:-1])+[kk+'rGW'])
            kegw = tuple(list(ke[:-1])+['GW'])
            resvi[nuke] = resvi[ke]/resvi[kegw]

kose = np.unique([ke[:-1] for ke in resvi if 'slope' in ke], axis = 0)
for ke in kose:
    resvi[tuple(ke)+tuple(['NAW2rGW'])] = resvi[tuple(ke)+tuple(['TNAWrGW'])]+resvi[tuple(ke)+tuple(['SPGWrGW'])]
    resvi[tuple(ke)+tuple(['NAW2'])] = resvi[tuple(ke)+tuple(['TNAW'])]+resvi[tuple(ke)+tuple(['SPGW'])]

vkeys = np.unique([ke[-1] for ke in resvi])

ctl.mkdir(cart_out_orig)
cart_corr_sig = cart_out_orig + 'corrplots_regs/'
ctl.mkdir(cart_corr_sig)

varli = dict()

for ssp in ['rcp85', 'ssp585']:
    cos = ''
    if ssp == 'rcp85':
        cos = '_cmip5'
    for mod in okmods[ssp]:
        for area in ['EAT', 'PNA']:
            trend_ssp = cose[(area, 'trend')]
            freqs = cose[(area, 'freq')]
            varli[(ssp, mod, area, 'trend')] = np.array([trend_ssp[(ssp+cos, mod, 'trend', 'freq10', reg)] for reg in range(4)])
            varli[(ssp, mod, area, 'fr_fin')] = np.array(freqs[(ssp+cos, mod, 'tot50')])
            varli[(ssp, mod, area, 'fr_diff')] = np.array(freqs[(ssp+cos, mod, 'tot50')])-np.array(freqs[('hist'+cos, mod, 'tot50')])

allpeas = dict()
livar = 'trend'
for reg1 in range(4):
    for reg2 in range(4):
        x = []
        y = []
        for ssp in ['ssp585', 'rcp85']:
            xss = []
            yss = []
            for mod in okmods[ssp]:
                xss.append(varli[(ssp, mod, 'EAT', livar)][reg1])
                yss.append(varli[(ssp, mod, 'PNA', livar)][reg2])
            x.append(xss)
            y.append(yss)

        pears, pval = stats.pearsonr(np.concatenate(x), np.concatenate(y))
        filnam = cart_corr_sig + 'corr_EAT{}_PNA{}.pdf'.format(reg1, reg2)
        pea = ctl.plotcorr_wgroups(x, y, filename = filnam, xlabel = reg_names_area['EAT'][reg1], ylabel = reg_names_area['PNA'][reg2], groups = ['rcp85', 'ssp585'], colors = ['blue', 'red'], single_group_corr = True)
        allpeas[(reg1, reg2)] = (pears, pval)
        if abs(pears) > 0.3:
            print('-------->', reg1, reg2, '{:5.2f}'.format(pears), '{:6.3e}'.format(pval))
        else:
            print(reg1, reg2, '{:5.2f}'.format(pears), '{:6.3e}'.format(pval))

#### THE REGRESSION

from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler


na = 'slope'
varok = 'trend'

import itertools as itt

# best 3 drivers for EAT
varall = ['GW', 'PVS', 'UTWrGW', 'AArGW', 'NAWrGW', 'PST']
#varall = ['GW', 'PVS', 'RUTAW', 'NAWrGW', 'PST']
pio2 = list(itt.combinations(varall, 2))
pio4 = list(itt.combinations(varall, 4))
pio3 = list(itt.combinations(varall, 3))

print('\n\n\n With GW \n')
# partiamo con model1
print('BEST 3 DRIVERS!')
allscores = dict()
bestxs = dict()
bestys = dict()
bestset = dict()
tip = 'conGW'
for num, pio in zip([2,3,4], [pio2, pio3, pio4]):
    for area in ['EAT', 'PNA']:
        allscores[(tip, num, area)] = []
        for pi in pio:
            # if 'PVSrGW' in pi and 'PSTrGW' in pi: continue
            X = []
            Y = []
            for ssp in ['rcp85', 'ssp585']:
                for mod in okmods[ssp]:
                    if (ssp, mod, na, 'AA') in resvi:
                        Xmod = np.array([resvi[(ssp, mod, na, ke)] for ke in pi])
                        X.append(Xmod)

                        Y.append(varli[(ssp, mod, area, varok)])

            Y = np.stack(Y)
            # STANDARDIZZO LE FEATURES
            X = np.stack(X)
            scaler = StandardScaler().fit(X)
            X = scaler.transform(X)

            model = LinearRegression().fit(X, Y)
            allscores[(tip, num, area)].append(model.score(X, Y))
            print(pi, model.score(X, Y))

        argma = np.argmax(allscores[(tip, num, area)])
        bestset[(tip, num, area)] = list(pio)[argma]

        X = []
        for ssp in ['rcp85', 'ssp585']:
            for mod in okmods[ssp]:
                if (ssp, mod, na, 'AA') in resvi:
                    Xmod = np.array([resvi[(ssp, mod, na, ke)] for ke in bestset[(tip, num, area)]])
                    X.append(Xmod)
        X = np.stack(X)
        scaler = StandardScaler().fit(X)
        X = scaler.transform(X)

        bestxs[(tip, num, area)] = X
        bestys[(tip, num, area)] = Y
        print(tip, num, area, list(pio)[argma], allscores[(tip, num, area)][argma], '\n')
        print(np.cov(X.T)/np.cov(X.T)[0,0], '\n')


print('\n\n\n With NO GW scaling of drivers \n')
# best 3 drivers for EAT
varall = ['PVS', 'UTW', 'AArGW', 'NAWrGW', 'PST']
pio2 = list(itt.combinations(varall, 2))
pio4 = list(itt.combinations(varall, 4))
pio3 = list(itt.combinations(varall, 3))

# partiamo con model1
print('BEST 3 DRIVERS!')
tip = 'unscal'
for num, pio in zip([2,3,4], [pio2, pio3, pio4]):
    for area in ['EAT', 'PNA']:
        allscores[(tip, num, area)] = []
        for pi in pio:
            # if 'PVSrGW' in pi and 'PSTrGW' in pi: continue
            X = []
            Y = []
            for ssp in ['rcp85', 'ssp585']:
                for mod in okmods[ssp]:
                    if (ssp, mod, na, 'AA') in resvi:
                        Xmod = np.array([resvi[(ssp, mod, na, ke)] for ke in pi])
                        X.append(Xmod)

                        Y.append(varli[(ssp, mod, area, varok)])

            Y = np.stack(Y)
            # STANDARDIZZO LE FEATURES
            X = np.stack(X)
            scaler = StandardScaler().fit(X)
            X = scaler.transform(X)

            model = LinearRegression().fit(X, Y)
            allscores[(tip, num, area)].append(model.score(X, Y))
            print(pi, model.score(X, Y))

        argma = np.argmax(allscores[(tip, num, area)])
        bestset[(tip, num, area)] = list(pio)[argma]

        X = []
        for ssp in ['rcp85', 'ssp585']:
            for mod in okmods[ssp]:
                if (ssp, mod, na, 'AA') in resvi:
                    Xmod = np.array([resvi[(ssp, mod, na, ke)] for ke in bestset[(tip, num, area)]])
                    X.append(Xmod)
        X = np.stack(X)
        scaler = StandardScaler().fit(X)
        X = scaler.transform(X)

        bestxs[(tip, num, area)] = X
        bestys[(tip, num, area)] = Y
        print(tip, num, area, list(pio)[argma], allscores[(tip, num, area)][argma], '\n')
        print(np.cov(X.T)/np.cov(X.T)[0,0], '\n')


# RESCALING THE FREQ TRENDS TO GW
varall = ['RUTAW', 'PVSrGW', 'UTWrGW', 'AArGW', 'NAWrGW', 'PSTrGW']
pio2 = list(itt.combinations(varall, 2))
varall = ['PVSrGW', 'UTWrGW', 'AArGW', 'NAWrGW', 'PSTrGW']
pio3 = list(itt.combinations(varall, 3))
pio4 = list(itt.combinations(varall, 4))

tip = 'senzaGW'
# partiamo con model1
print('\n\n\n Without GW!\n')
for num, pio in zip([2,3,4], [pio2, pio3, pio4]):
    for area in ['EAT', 'PNA']:
        allscores[(tip, num, area)] = []
        for pi in pio:
            # if 'PVSrGW' in pi and 'PSTrGW' in pi: continue
            X = []
            Y = []
            for ssp in ['rcp85', 'ssp585']:
                for mod in okmods[ssp]:
                    if (ssp, mod, na, 'AA') in resvi:
                        Xmod = np.array([resvi[(ssp, mod, na, ke)] for ke in pi])
                        X.append(Xmod)

                        Y.append(varli[(ssp, mod, area, varok)])

            Y = np.stack(Y)
            # STANDARDIZZO LE FEATURES
            X = np.stack(X)
            scaler = StandardScaler().fit(X)
            X = scaler.transform(X)

            model = LinearRegression().fit(X, Y)
            allscores[(tip, num, area)].append(model.score(X, Y))
            print(pi, model.score(X, Y))

        argma = np.argmax(allscores[(tip, num, area)])
        bestset[(tip, num, area)] = list(pio)[argma]
        X = []
        for ssp in ['rcp85', 'ssp585']:
            for mod in okmods[ssp]:
                if (ssp, mod, na, 'AA') in resvi:
                    Xmod = np.array([resvi[(ssp, mod, na, ke)] for ke in bestset[(tip, num, area)]])
                    X.append(Xmod)
        X = np.stack(X)
        scaler = StandardScaler().fit(X)
        X = scaler.transform(X)

        bestxs[(tip, num, area)] = X
        bestys[(tip, num, area)] = Y
        print(tip, num, area, list(pio)[argma], allscores[(tip, num, area)][argma], '\n')
        print(np.cov(X.T)/np.cov(X.T)[0,0], '\n')


# RESCALING THE FREQ TRENDS TO GW
#varall = ['UTWrGW', 'AArGW', 'NAWrGW', 'PSTrGW'] # FIRST version of the article
#varall = ['UTWrGW', 'AArGW', 'PSTrGW', 'TNAWrGW', 'SPGWrGW'] # SECOND version of the article
varall = ['UTWrGW', 'AArGW', 'PSTrGW', 'NAWrGW'] # SECOND version of the article
pio2 = list(itt.combinations(varall, 2))
pio3 = list(itt.combinations(varall, 3))
pio4 = list(itt.combinations(varall, 4))

tip = 'GWscaling'
# partiamo con model1
print('\n\n\n SCALING FREQ trends for GW!\n')
for num, pio in zip([2,3,4], [pio2, pio3, pio4]):
    for area in ['EAT', 'PNA']:
        allscores[(tip, num, area)] = []
        for pi in pio:
            # if 'PVSrGW' in pi and 'PSTrGW' in pi: continue
            X = []
            Y = []
            for ssp in ['rcp85', 'ssp585']:
                for mod in okmods[ssp]:
                    if (ssp, mod, na, 'AA') in resvi:
                        Xmod = np.array([resvi[(ssp, mod, na, ke)] for ke in pi])
                        X.append(Xmod)

                        Y.append(varli[(ssp, mod, area, varok)]/resvi[(ssp, mod, na, 'GW')])

            Y = np.stack(Y)
            # STANDARDIZZO LE FEATURES
            X = np.stack(X)
            scaler = StandardScaler().fit(X)
            X = scaler.transform(X)

            model = LinearRegression().fit(X, Y)
            allscores[(tip, num, area)].append(model.score(X, Y))
            print(pi, model.score(X, Y))

        argma = np.argmax(allscores[(tip, num, area)])
        bestset[(tip, num, area)] = list(pio)[argma]
        X = []
        for ssp in ['rcp85', 'ssp585']:
            for mod in okmods[ssp]:
                if (ssp, mod, na, 'AA') in resvi:
                    Xmod = np.array([resvi[(ssp, mod, na, ke)] for ke in bestset[(tip, num, area)]])
                    X.append(Xmod)
        X = np.stack(X)
        scaler = StandardScaler().fit(X)
        X = scaler.transform(X)

        bestxs[(tip, num, area)] = X
        bestys[(tip, num, area)] = Y
        print(tip, num, area, list(pio)[argma], allscores[(tip, num, area)][argma], '\n')
        print(np.cov(X.T)/np.cov(X.T)[0,0], '\n')


varall = ['UTWrGW', 'AArGW', 'PSTrGW', 'TNAWrGW', 'SPGWrGW'] # SECOND version of the article
pio2 = list(itt.combinations(varall, 2))
pio3 = list(itt.combinations(varall, 3))
pio4 = list(itt.combinations(varall, 4))
pio5 = [varall]
tip = 'GWscaling2'
# partiamo con model1
print('\n\n\n SCALING FREQ trends for GW!\n')
for num, pio in zip([2,3,4,5], [pio2, pio3, pio4, pio5]):
    for area in ['EAT', 'PNA']:
        allscores[(tip, num, area)] = []
        for pi in pio:
            # if area == 'EAT' and ('TNAWrGW' in pi or 'SPGWrGW' in pi):
            #     print('bauuu')
            #     allscores[(tip, num, area)].append(0.)
            #     continue
            # if 'PVSrGW' in pi and 'PSTrGW' in pi: continue
            X = []
            Y = []
            for ssp in ['rcp85', 'ssp585']:
                for mod in okmods[ssp]:
                    if (ssp, mod, na, 'AA') in resvi:
                        Xmod = np.array([resvi[(ssp, mod, na, ke)] for ke in pi])
                        X.append(Xmod)

                        Y.append(varli[(ssp, mod, area, varok)]/resvi[(ssp, mod, na, 'GW')])

            Y = np.stack(Y)
            # STANDARDIZZO LE FEATURES
            X = np.stack(X)
            scaler = StandardScaler().fit(X)
            X = scaler.transform(X)

            model = LinearRegression().fit(X, Y)
            allscores[(tip, num, area)].append(model.score(X, Y))
            print(pi, model.score(X, Y))

        argma = np.argmax(allscores[(tip, num, area)])
        bestset[(tip, num, area)] = list(pio)[argma]
        X = []
        for ssp in ['rcp85', 'ssp585']:
            for mod in okmods[ssp]:
                if (ssp, mod, na, 'AA') in resvi:
                    Xmod = np.array([resvi[(ssp, mod, na, ke)] for ke in bestset[(tip, num, area)]])
                    X.append(Xmod)
        X = np.stack(X)
        scaler = StandardScaler().fit(X)
        X = scaler.transform(X)

        bestxs[(tip, num, area)] = X
        bestys[(tip, num, area)] = Y
        print(tip, num, area, list(pio)[argma], allscores[(tip, num, area)][argma], '\n')
        print(np.cov(X.T)/np.cov(X.T)[0,0], '\n')


import statsmodels.api as sm
from scipy import stats
# 2 EAT ('AArGW', 'PSTrGW') 0.30286198026464933
# 2 PNA ('AArGW', 'NAWrGW') 0.23600295791630332
# 3 EAT ('UTWrGW', 'AArGW', 'PSTrGW') 0.3278009127949955
# 3 PNA ('PVSrGW', 'AArGW', 'NAWrGW') 0.26665645026370693
resall = dict()
for tip in ['conGW', 'unscal', 'senzaGW', 'GWscaling', 'GWscaling2']:
    for num in [2,3,4,5]:
        for area in ['EAT', 'PNA']:
            if (tip, num, area) not in bestxs:
                continue
            xii = bestxs[(tip, num, area)]
            varfii = bestset[(tip, num, area)]
            y = bestys[(tip, num, area)]

            Xgi = sm.add_constant(xii)
            pvals = []
            params = []
            rsq = []
            for ko, nam in zip(y.T, reg_names_area[area]):
                print('\n')
                print(nam)
                est = sm.OLS(ko, Xgi)
                est2 = est.fit()
                pvals.append(est2.pvalues)
                params.append(est2.params)
                rsq.append(est2.rsquared)

            pvals = np.stack(pvals)
            params = np.stack(params)
            rsq = np.stack(rsq)
            resall[(tip, num, area, 'pvals')] = pvals[:, 1:]
            resall[(tip, num, area, 'params')] = params[:, 1:]
            resall[(tip, num, area, 'rsq')] = rsq


#############################################################################
from matplotlib import colors
import matplotlib.gridspec as gridspec

# FIGURA
colo = '#a50026 #d73027 #f46d43 #fdae61 #fee090 #ffffff #e0f3f8 #abd9e9 #74add1 #4575b4 #313695'
colo = colo.split()
colo = colo[::-1]
cmappa = colors.ListedColormap(colo)
cmappa.set_over('#800026') #662506
cmappa.set_under('#023858') #542788

for tip, vlim in zip(['conGW', 'unscal', 'senzaGW', 'GWscaling', 'GWscaling2'], [(-0.05, 0.05), (-0.05, 0.05), (-0.05, 0.05), (-1.2, 1.2), (-1.2, 1.2)]):

    piolo = np.concatenate([np.concatenate([resall[(tip, num, area, 'params')].flatten() for area in ['EAT', 'PNA']]) for num in [2,3]])
    vmin = np.min(piolo)
    vmax = np.max(piolo)
    vmi = np.round(np.max(np.abs([vmin, vmax])), decimals=2)

    fig = plt.figure(figsize=(12,8))
    gs = gridspec.GridSpec(7, 9)
    axes = []
    axes.append(fig.add_subplot(gs[:3, 1:5]))
    axes.append(fig.add_subplot(gs[:3, 5:]))
    axes.append(fig.add_subplot(gs[4:6, 1:5]))
    axes.append(fig.add_subplot(gs[4:6, 5:]))

    bau = -1
    for num in [3,2]:
        for area in ['EAT', 'PNA']:
            varf = bestset[(tip, num, area)]
            varf = [va if 'rGW' not in va else va[:-3] for va in varf]
            nvars = len(varf)
            ext = [0, 4, 0, nvars]

            bau += 1
            ax = axes[bau]
            gigi = ax.imshow(resall[(tip, num, area, 'params')].T, vmin = -vmi, vmax = vmi, cmap = cmappa, origin = 'lower',  extent = ext, aspect = 1)

            pvaok = resall[(tip, num, area, 'pvals')].T
            for ix in range(4):
                for iy in range(nvars):
                    if pvaok[iy, ix] < 0.01:
                        ax.scatter(ix+0.5, iy+0.5, s=60, edgecolors = 'black', facecolors='white')
                    elif pvaok[iy, ix] < 0.05:
                        ax.scatter(ix+0.5, iy+0.5, s=20, edgecolors = 'black', facecolors='white')

            ax.xaxis.tick_top()
            ax.set_xticks(0.5+np.arange(4), minor = False)
            ax.set_xticklabels(reg_names_area[area], ha='center')
            ax.set_yticks(0.5+np.arange(len(varf)), minor = False)
            ax.set_yticklabels(varf, va='center')
            if num == 3:
                if area == 'EAT':
                    ax.set_title(area)
                else:
                    ax.set_title('PAC')

    #cax = fig.add_subplot(gs[6, :])
    cax = plt.axes([0.1, 0.1, 0.8, 0.05])
    cb = plt.colorbar(gigi, cax=cax, orientation='horizontal')
    cb.ax.tick_params(labelsize=18)
    cb.set_label(r'Regression coefficient ($K^{-1}$)', fontsize=20)
    plt.subplots_adjust(left=0.02, bottom=0.13, right=0.98, top=0.92, wspace=0.05, hspace=0.20)

    ax.text(0.05, 0.75, '3 drivers', horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 20)
    ax.text(0.05, 0.35, '2 drivers', horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 20)

    fig.savefig(cart_out_orig + 'Regr_model_optimal_{}.pdf'.format(tip))

# 4 drivers
    fig = plt.figure(figsize=(16,8))
    gs = gridspec.GridSpec(5, 2)
    axes = []
    axes.append(fig.add_subplot(gs[:4, 0]))
    axes.append(fig.add_subplot(gs[:4, 1]))

    bau = -1
    num = 4
    for area in ['EAT', 'PNA']:
        varf = bestset[(tip, num, area)]
        varf = [va if 'rGW' not in va else va[:-3] for va in varf]
        nvars = len(varf)
        ext = [0, 4, 0, nvars]

        bau += 1
        ax = axes[bau]
        gigi = ax.imshow(resall[(tip, num, area, 'params')].T, vmin = -vmi, vmax = vmi, cmap = cmappa, origin = 'lower',  extent = ext, aspect = 1)

        pvaok = resall[(tip, num, area, 'pvals')].T
        for ix in range(4):
            for iy in range(nvars):
                if pvaok[iy, ix] < 0.01:
                    ax.scatter(ix+0.5, iy+0.5, s=60, edgecolors = 'black', facecolors='white')
                elif pvaok[iy, ix] < 0.05:
                    ax.scatter(ix+0.5, iy+0.5, s=20, edgecolors = 'black', facecolors='white')

        ax.xaxis.tick_top()
        ax.set_xticks(0.5+np.arange(4), minor = False)
        ax.set_xticklabels(reg_names_area[area], ha='center')
        ax.set_yticks(0.5+np.arange(len(varf)), minor = False)
        ax.set_yticklabels(varf, va='center')
        ax.set_title(area)

    #cax = fig.add_subplot(gs[6, :])
    cax = plt.axes([0.1, 0.1, 0.8, 0.05])
    cb = plt.colorbar(gigi, cax=cax, orientation='horizontal')
    cb.ax.tick_params(labelsize=18)
    cb.set_label(r'Regression coefficient ($yr^{-1}$)', fontsize=20)
    plt.subplots_adjust(left=0.02, bottom=0.13, right=0.98, top=0.92, wspace=0.05, hspace=0.20)

    fig.savefig(cart_out_orig + 'Regr_model_optimal_{}_4drivers.pdf'.format(tip))


    if (tip, 5, area) in bestxs:
        bau = -1
        num = 5
        fig = plt.figure(figsize=(16,8))
        gs = gridspec.GridSpec(5, 2)
        axes = []
        axes.append(fig.add_subplot(gs[:4, 0]))
        axes.append(fig.add_subplot(gs[:4, 1]))

        for area in ['EAT', 'PNA']:
            varf = bestset[(tip, num, area)]
            varf = [va if 'rGW' not in va else va[:-3] for va in varf]
            nvars = len(varf)
            ext = [0, 4, 0, nvars]

            bau += 1
            ax = axes[bau]
            gigi = ax.imshow(resall[(tip, num, area, 'params')].T, vmin = -vmi, vmax = vmi, cmap = cmappa, origin = 'lower',  extent = ext, aspect = 1)

            pvaok = resall[(tip, num, area, 'pvals')].T
            for ix in range(4):
                for iy in range(nvars):
                    if pvaok[iy, ix] < 0.01:
                        ax.scatter(ix+0.5, iy+0.5, s=60, edgecolors = 'black', facecolors='white')
                    elif pvaok[iy, ix] < 0.05:
                        ax.scatter(ix+0.5, iy+0.5, s=20, edgecolors = 'black', facecolors='white')

            ax.xaxis.tick_top()
            ax.set_xticks(0.5+np.arange(4), minor = False)
            ax.set_xticklabels(reg_names_area[area], ha='center')
            ax.set_yticks(0.5+np.arange(len(varf)), minor = False)
            ax.set_yticklabels(varf, va='center')

            if area == 'EAT':
                ax.set_title(area)
            else:
                ax.set_title('PAC')

        #cax = fig.add_subplot(gs[6, :])
        cax = plt.axes([0.1, 0.1, 0.8, 0.05])
        cb = plt.colorbar(gigi, cax=cax, orientation='horizontal')
        cb.ax.tick_params(labelsize=18)
        cb.set_label(r'Regression coefficient ($yr^{-1}$)', fontsize=20)
        plt.subplots_adjust(left=0.02, bottom=0.13, right=0.98, top=0.92, wspace=0.05, hspace=0.20)

        fig.savefig(cart_out_orig + 'Regr_model_TOT_{}_5drivers.pdf'.format(tip))

    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(111)
    for num, col in zip([2,3,4], ['indianred', 'steelblue', 'forestgreen']):
        print(tip, num)
        print(np.mean(resall[(tip, num, 'EAT', 'rsq')]))
        print(np.mean(resall[(tip, num, 'PNA', 'rsq')]))
        vals = np.concatenate([resall[(tip, num, 'EAT', 'rsq')], resall[(tip, num, 'PNA', 'rsq')]])
        ax.plot(np.arange(8), vals, label = '{} drivers'.format(num), color = col)
        ax.scatter(np.arange(8), vals, color = col)

        ax.set_xticks(np.arange(8), minor = False)
        ax.set_xticklabels(reg_names_area['EAT']+reg_names_area['PNA'], ha='center')
        ax.legend()
        ax.set_ylabel('R squared')
    fig.savefig(cart_out_orig + 'Rsquared_optimal_{}.pdf'.format(tip))


#varall = ['GW', 'PSTrGW', 'UTWrGW', 'AArGW', 'NAWrGW'] # First version
varall = ['GW', 'UTWrGW', 'AArGW', 'PSTrGW', 'TNAWrGW', 'SPGWrGW']
X = []
for ssp in ['rcp85', 'ssp585']:
    for mod in okmods[ssp]:
        if (ssp, mod, na, 'AA') in resvi:
            Xmod = np.array([resvi[(ssp, mod, na, ke)] for ke in varall])
            X.append(Xmod)
X = np.stack(X)
np.set_printoptions(precision=3)
print(varall)
print(np.mean(X, axis = 0))
print(np.std(X, axis = 0))

scaler = StandardScaler().fit(X)
X = scaler.transform(X)
print(np.cov(X.T)/np.cov(X.T)[0,0], '\n')


#varall = ['GW', 'UTW', 'AA', 'NAW', 'PST']  # First version
varall = ['GW', 'UTW', 'AA', 'PST', 'TNAW', 'SPGW']
X = []
for ssp in ['rcp85', 'ssp585']:
    for mod in okmods[ssp]:
        if (ssp, mod, na, 'AA') in resvi:
            Xmod = np.array([resvi[(ssp, mod, na, ke)] for ke in varall])
            X.append(Xmod)
X = np.stack(X)
np.set_printoptions(precision=3)
print(varall)
print(np.mean(X, axis = 0))
print(np.std(X, axis = 0))

scaler = StandardScaler().fit(X)
X = scaler.transform(X)
print(np.cov(X.T)/np.cov(X.T)[0,0], '\n')
