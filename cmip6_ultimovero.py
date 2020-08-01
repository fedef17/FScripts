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

if os.uname()[1] == 'hobbes':
    cart_in = '/home/fabiano/Research/lavori/CMIP6/'
elif os.uname()[1] == 'ff-clevo':
    cart_in = '/home/fedefab/Scrivania/Research/Post-doc/lavori/CMIP6/'

cart_out_orig = cart_in + 'Results_ultimo/'
ctl.mkdir(cart_out_orig)

cart_data = '/data-hobbes/fabiano/WR_CMIP6/'

file_hist_refEOF = cart_data + 'out_NEW_cmip6_hist_NDJFM_{}_4clus_4pcs_1957-2005_refEOF_dtr.p'

file_hist = cart_data + 'out_NEW_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refCLUS_dtr_light.p'
gen_file_ssp = cart_data + 'out_NEW_cmip6_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr_histrebase.p'

yr10 = 10 # length of running mean

cart_cmip5 = '/home/fedefab/Scrivania/Research/Post-doc/lavori/CMIP6/Results_cmip5/{}_NDJFM/'
cart_v5 = '/home/fedefab/Scrivania/Research/Post-doc/lavori/CMIP6/Results_v5_rebase/{}_NDJFM/'

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

filbi = '/home/fedefab/Scrivania/Research/Post-doc/lavori/CMIP6/Results_cmip6_vs_cmip5_refEOF/histbiases.p'
with open(filbi, 'rb') as pino:
    res_short = pickle.load(pino)
kesal = list(res_short.keys())
for ke in kesal:
    if 'cmip5_refEOF' in ke:
        res_short[tuple(['rcp85'] + list(ke[1:]))] = res_short[ke]
    if 'cmip6_refEOF' in ke:
        res_short[tuple(['ssp585'] + list(ke[1:]))] = res_short[ke]

# da Virna
vicar = '/home/fedefab/Scrivania/Research/Post-doc/lavori/CMIP6/virnas_tas/new/'
fils = ['Indices_{}.txt', 'Slopes_{}.txt', 'Slopes_{}_temp.txt']
nams = ['temp', 'slope', 'slope']
tip = [0,0,1]

resvi = dict()
for na, fi, ti in zip(nams, fils, tip):
    for ssp in ['rcp85', 'ssp585']:
        vmods, datatemp = ctl.read_from_txt(vicar + fi.format(ssp), n_skip = 2, sep = '\t')

        if ti == 0:
            vkeys = ['AA', 'UTW', 'PST', 'PVS']
            for mod in okmods[ssp]:
                if mod.split('_')[0] in vmods:
                    ind = np.where(vmods == mod.split('_')[0])[0][0]
                    lin = datatemp[ind]
                    for ke, te in zip(vkeys, lin):
                        resvi[(ssp, mod, na, ke)] = te
                    resvi[(ssp, mod, na, 'RUTAW')] = resvi[(ssp, mod, na, 'UTW')]/resvi[(ssp, mod, na, 'AA')]
        else:
            vkeys = ['GW', 'NAW']
            for mod in okmods[ssp]:
                if mod.split('_')[0] in vmods:
                    ind = np.where(vmods == mod.split('_')[0])[0][0]
                    lin = datatemp[ind]
                    for ke, te in zip(vkeys, lin):
                        resvi[(ssp, mod, na, ke)] = te

vkeys = ['AA', 'UTW', 'PST', 'PVS', 'NAW', 'RUTAW', 'GW']

for ke in list(resvi.keys()):
    if 'temp' in ke: continue
    for kk in vkeys[:-2]:
        if kk in ke:
            nuke = tuple(list(ke[:-1])+[kk+'rGW'])
            kegw = tuple(list(ke[:-1])+['GW'])
            resvi[nuke] = resvi[ke]/resvi[kegw]

vkeys = np.unique([ke[-1] for ke in resvi])

ctl.mkdir(cart_out_orig)
cart_corr_sig = cart_out_orig + 'corrplots/'
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

#### THE REGRESSION

from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler


na = 'slope'
varok = 'trend'

varall = ['GW', 'RUTAW', 'PVSrGW', 'UTWrGW', 'AArGW', 'NAWrGW', 'PSTrGW']

# best 3 drivers for EAT
import itertools as itt
pio3 = list(itt.combinations(varall, 3))
pio2 = list(itt.combinations(varall, 2))

print('\n\n\n With GW \n')
# partiamo con model1
print('BEST 3 DRIVERS!')
allscores = dict()
for num, pio in zip([2,3], [pio2, pio3]):
    for area in ['EAT', 'PNA']:
        allscores[(num, area)] = []
        for pi in pio:
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
            allscores[(num, area)].append(model.score(X, Y))
            print(pi, model.score(X, Y))

        argma = np.argmax(allscores[(num, area)])
        print(num, area, list(pio)[argma], allscores[(num, area)][argma], '\n')


# RESCALING THE FREQ TRENDS TO GW
varall = ['RUTAW', 'PVSrGW', 'UTWrGW', 'AArGW', 'NAWrGW', 'PSTrGW']
pio2 = list(itt.combinations(varall, 2))
varall = ['PVSrGW', 'UTWrGW', 'AArGW', 'NAWrGW', 'PSTrGW']
pio3 = list(itt.combinations(varall, 3))

# partiamo con model1
print('\n\n\n Without GW!\n')
allscores = dict()
for num, pio in zip([2,3], [pio2, pio3]):
    for area in ['EAT', 'PNA']:
        allscores[(num, area)] = []
        for pi in pio:
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
            allscores[(num, area)].append(model.score(X, Y))
            print(pi, model.score(X, Y))

        argma = np.argmax(allscores[(num, area)])
        print(num, area, list(pio)[argma], allscores[(num, area)][argma], '\n')


# RESCALING THE FREQ TRENDS TO GW
varall = ['RUTAW', 'PVSrGW', 'UTWrGW', 'AArGW', 'NAWrGW', 'PSTrGW']
pio2 = list(itt.combinations(varall, 2))
varall = ['PVSrGW', 'UTWrGW', 'AArGW', 'NAWrGW', 'PSTrGW']
pio3 = list(itt.combinations(varall, 3))

# partiamo con model1
print('\n\n\n SCALING FREQ trends for GW!\n')
allscores = dict()
bestxs = dict()
bestys = dict()
bestset = dict()
for num, pio in zip([2,3], [pio2, pio3]):
    for area in ['EAT', 'PNA']:
        allscores[(num, area)] = []
        for pi in pio:
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
            allscores[(num, area)].append(model.score(X, Y))
            print(pi, model.score(X, Y))

        argma = np.argmax(allscores[(num, area)])
        bestxs[(num, area)] = X
        bestys[(num, area)] = Y
        bestset[(num, area)] = list(pio)[argma]
        print(num, area, list(pio)[argma], allscores[(num, area)][argma], '\n')


import statsmodels.api as sm
from scipy import stats
# 2 EAT ('AArGW', 'PSTrGW') 0.30286198026464933
# 2 PNA ('AArGW', 'NAWrGW') 0.23600295791630332
# 3 EAT ('UTWrGW', 'AArGW', 'PSTrGW') 0.3278009127949955
# 3 PNA ('PVSrGW', 'AArGW', 'NAWrGW') 0.26665645026370693
resall = dict()
for num, pio in zip([2,3], [pio2, pio3]):
    for area in ['EAT', 'PNA']:
        xii = bestxs[(num, area)]
        varfii = bestset[(num, area)]
        y = bestys[(num, area)]

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
        resall[(num, area, 'pvals')] = pvals[:, 1:]
        resall[(num, area, 'params')] = params[:, 1:]
        resall[(num, area, 'rsq')] = rsq


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


fig = plt.figure(figsize=(12,8))
gs = gridspec.GridSpec(7, 9)
axes = []
axes.append(fig.add_subplot(gs[:3, 1:5]))
axes.append(fig.add_subplot(gs[:3, 5:]))
axes.append(fig.add_subplot(gs[4:6, 1:5]))
axes.append(fig.add_subplot(gs[4:6, 5:]))

bau = -1
for num, pio in zip([3,2], [pio2, pio3]):
    for area in ['EAT', 'PNA']:
        varf = bestset[(num, area)]
        varf = [va if 'rGW' not in va else va[:-3] for va in varf]
        nvars = len(varf)
        ext = [0, 4, 0, nvars]

        bau += 1
        ax = axes[bau]
        gigi = ax.imshow(resall[(num, area, 'params')].T, vmin = -1.2, vmax = 1.2, cmap = cmappa, origin = 'lower',  extent = ext, aspect = 1)

        pvaok = resall[(num, area, 'pvals')].T
        for ix in range(4):
            for iy in range(nvars):
                if pvaok[iy, ix] < 0.05:
                    ax.scatter(ix+0.5, iy+0.5, s=60, edgecolors = 'black', facecolors='white')
                elif pvaok[iy, ix] < 0.1:
                    ax.scatter(ix+0.5, iy+0.5, s=20, edgecolors = 'black', facecolors='white')

        ax.xaxis.tick_top()
        ax.set_xticks(0.5+np.arange(4), minor = False)
        ax.set_xticklabels(reg_names_area[area], ha='center')
        ax.set_yticks(0.5+np.arange(len(varf)), minor = False)
        ax.set_yticklabels(varf, va='center')
        if num == 3: ax.set_title(area)

#cax = fig.add_subplot(gs[6, :])
cax = plt.axes([0.1, 0.08, 0.8, 0.03])
cb = plt.colorbar(gigi, cax=cax, orientation='horizontal')
cb.ax.tick_params(labelsize=18)
cb.set_label('Regression coefficient', fontsize=20)
plt.subplots_adjust(left=0.02, bottom=0.13, right=0.98, top=0.92, wspace=0.05, hspace=0.20)

ax.text(0.05, 0.8, '3 drivers', horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 20)
ax.text(0.05, 0.4, '2 drivers', horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 20)

fig.savefig(cart_out_orig + 'Regr_model_optimal.pdf')

fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
for num, col in zip([2,3], ['indianred', 'steelblue']):
    vals = np.concatenate([resall[(num, 'EAT', 'rsq')], resall[(num, 'PNA', 'rsq')]])
    ax.plot(np.arange(8), vals, label = '{} drivers'.format(num), color = col)
    ax.scatter(np.arange(8), vals, color = col)

    ax.set_xticks(np.arange(8), minor = False)
    ax.set_xticklabels(reg_names_area['EAT']+reg_names_area['PNA'], ha='center')
    ax.legend()
    ax.set_ylabel('R squared')
fig.savefig(cart_out_orig + 'Rsquared_optimal.pdf'.format(pio))
