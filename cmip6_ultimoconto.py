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
        vmods, datatemp = ctl.read_from_txt(vicar + fi.format(ssp), n_skip = 2, sep = '\t', debug = True)

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


# allpeas = dict()
#
# for area in ['EAT', 'PNA']:
#     for ke in vkeys:
#         for livar in ['trend', 'fr_fin', 'fr_diff']:
#             co = (livar, ke, area)
#             x = []
#             y = []
#             for ssp in ['ssp585', 'rcp85']:
#                 xss = []
#                 yss = []
#                 for mod in okmods[ssp]:
#                     if (ssp, mod, 'slope', ke) in resvi.keys() and (ssp, mod, area, livar) in varli.keys():
#                         xss.append(resvi[(ssp, mod, 'slope', ke)])
#                         yss.append(varli[(ssp, mod, area, livar)][0])
#                 x.append(xss)
#                 y.append(yss)
#
#             pears, pval = stats.pearsonr(np.concatenate(x), np.concatenate(y))
#             filnam = cart_corr_sig + 'corr_{}_{}_{}.pdf'.format(*co)
#             pea = ctl.plotcorr_wgroups(x, y, filename = filnam, xlabel = ke, ylabel = livar, groups = ['rcp85', 'ssp585'], colors = ['blue', 'red'], single_group_corr = True)
#             allpeas[co] = (pears, pval)
#             if abs(pears) > 0.3:
#                 print('-------->', co, '{:5.2f}'.format(pears), '{:6.3e}'.format(pval))
#             else:
#                 print(co, '{:5.2f}'.format(pears), '{:6.3e}'.format(pval))



from sklearn.linear_model import LinearRegression

#model1: AA UTW PVS
varfit1 = ['PVS', 'UTW', 'AA']
#model2: RUTAW PST
varfit2 = ['RUTAW', 'PSTrGW']
#varfit2 = ['RUTAW', 'PVSrGW']
#model3: best of quellisopra + weights da var_ratio
#varfit3 = ['PVSrGW', 'UTWrGW', 'AArGW']
varfit3 = ['PVSrGW', 'UTWrGW', 'AArGW']
#model4
varfit4 = ['NAWrGW', 'PVSrGW', 'UTWrGW', 'AArGW']

varfit5 = ['NAWrGW', 'UTWrGW', 'AArGW']

varfit6 = ['RUTAW', 'NAWrGW']


# partiamo con model1
print('MODEL 1!')

na = 'slope'
varok = 'trend'

X = []
X2 = []
X3 = []
X4 = []
X5 = []
X6 = []

Y = []
Yall = []
Yeat = []
Ypna = []
Yfirst = []
weights = []
weights2 = []
#for ssp in ['ssp585']:
#for ssp in ['rcp85']:
for ssp in ['rcp85', 'ssp585']:
    for mod in okmods[ssp]:
        if (ssp, mod, na, 'AA') in resvi:
            Xmod = np.array([resvi[(ssp, mod, na, ke)] for ke in varfit1])
            X.append(Xmod)
            Xmod = np.array([resvi[(ssp, mod, na, ke)] for ke in varfit2])
            X2.append(Xmod)
            Xmod = np.array([resvi[(ssp, mod, na, ke)] for ke in varfit3])
            X3.append(Xmod)
            Xmod = np.array([resvi[(ssp, mod, na, ke)] for ke in varfit4])
            X4.append(Xmod)
            Xmod = np.array([resvi[(ssp, mod, na, ke)] for ke in varfit5])
            X5.append(Xmod)
            Xmod = np.array([resvi[(ssp, mod, na, ke)] for ke in varfit6])
            X6.append(Xmod)

            Ymod = np.concatenate([varli[(ssp, mod, 'EAT', varok)][:-1], varli[(ssp, mod, 'PNA', varok)][np.array([0,2,3])]])
            Y.append(Ymod)
            Ymod = np.concatenate([varli[(ssp, mod, 'EAT', varok)], varli[(ssp, mod, 'PNA', varok)]])
            Yall.append(Ymod)
            Yeat.append(varli[(ssp, mod, 'EAT', varok)][:2])
            Ypna.append(varli[(ssp, mod, 'PNA', varok)][np.array([0, -1])])
            Yfirst.append(np.array([varli[(ssp, mod, 'EAT', varok)][0], varli[(ssp, mod, 'PNA', varok)][0]]))

            wei = np.mean(res_short[(ssp, mod, 'EAT', 'patcor')])
            weights.append(wei)

            wei = res_short[(ssp, mod, 'EAT', 'var_ratio')]
            weights2.append(wei)


from sklearn.preprocessing import StandardScaler

# STANDARDIZZO LE FEATURES
X = np.stack(X)
X2 = np.stack(X2)
X3 = np.stack(X3)
X4 = np.stack(X4)
X5 = np.stack(X5)
X6 = np.stack(X6)

standardize_features = True
if standardize_features:
    scaler = StandardScaler().fit(X)
    X = scaler.transform(X)
    print('pio', X)

    scaler2 = StandardScaler().fit(X2)
    X2 = scaler2.transform(X2)
    print('pio', X2)

    scaler3 = StandardScaler().fit(X3)
    X3 = scaler3.transform(X3)
    print('pio', X3)

    scaler4 = StandardScaler().fit(X4)
    X4 = scaler4.transform(X4)
    print('pio', X4)

    scaler5 = StandardScaler().fit(X5)
    X5 = scaler5.transform(X5)
    print('pio', X5)

    scaler6 = StandardScaler().fit(X6)
    X6 = scaler6.transform(X6)
    print('pio', X6)

Y = np.stack(Y)
Yall = np.stack(Yall)
Yeat = np.stack(Yeat)
Ypna = np.stack(Ypna)
Yfirst = np.stack(Yfirst)
weights = np.array(weights)
weights2 = np.array(weights2)

print('6 regimes')
model1 = LinearRegression().fit(X, Y)
print(model1.score(X, Y))
print(model1.coef_)
#>>> reg.predict(np.array([[3, 5]]))
print('8 regimes')
model1 = LinearRegression().fit(X, Yall)
print(model1.score(X, Yall))
print(model1.coef_)




print('only 2 eat')
model1_eat = LinearRegression().fit(X, Yeat)
print(model1_eat.score(X, Yeat))
print(model1_eat.coef_)

print('only 2 pna')
model1_pna = LinearRegression().fit(X, Ypna)
print(model1_pna.score(X, Ypna))
print(model1_pna.coef_)

print('only NAO+ and PT')
model1_first = LinearRegression().fit(X, Yfirst)
print(model1_first.score(X, Yfirst))
print(model1_first.coef_)

print('only NAO+ and PT, weighted for patcor')
model1_first_wpatcor = LinearRegression().fit(X, Yfirst, sample_weight = weights)
print(model1_first_wpatcor.score(X, Yfirst, sample_weight = weights))
print(model1_first_wpatcor.coef_)

print('only NAO+ and PT, weighted for var_ratio')
model1_first_wvarrat = LinearRegression().fit(X, Yfirst, sample_weight = weights2)
print(model1_first_wvarrat.score(X, Yfirst, sample_weight = weights2))
print(model1_first_wvarrat.coef_)

#model2
print('MODEL 2!')
print('6 regimes')
model2 = LinearRegression().fit(X2, Y)
print(model2.score(X2, Y))
print(model2.coef_)

print('8 regimes')
model2 = LinearRegression().fit(X2, Yall)
print(model2.score(X2, Yall))
print(model2.coef_)

print('only 2 eat')
model2_eat = LinearRegression().fit(X2, Yeat)
print(model2_eat.score(X2, Yeat))
print(model2_eat.coef_)

print('only 2 pna')
model2_pna = LinearRegression().fit(X2, Ypna)
print(model2_pna.score(X2, Ypna))
print(model2_pna.coef_)

print('only NAO+ and PT')
model2_first = LinearRegression().fit(X2, Yfirst)
print(model2_first.score(X2, Yfirst))
print(model2_first.coef_)

print('only NAO+ and PT, weighted for patcor')
model2_first_wpatcor = LinearRegression().fit(X2, Yfirst, sample_weight = weights)
print(model2_first_wpatcor.score(X2, Yfirst, sample_weight = weights))
print(model2_first_wpatcor.coef_)

print('only NAO+ and PT, weighted for var_ratio')
model2_first_wvarrat = LinearRegression().fit(X2, Yfirst, sample_weight = weights2)
print(model2_first_wvarrat.score(X2, Yfirst, sample_weight = weights2))
print(model2_first_wvarrat.coef_)


# MODEL 3
print('MODEL 3!')
print('8 regimes')
model3 = LinearRegression().fit(X3, Yall)
print(model3.score(X3, Yall))
print(model3.coef_)

print('only 2 eat')
model3_eat = LinearRegression().fit(X3, Yeat)
print(model3_eat.score(X3, Yeat))
print(model3_eat.coef_)

print('only NAO+ and PT')
model3_first = LinearRegression().fit(X3, Yfirst)
print(model3_first.score(X3, Yfirst))
print(model3_first.coef_)


# MODEL 4
print('MODEL 4!')
print('8 regimes')
model4 = LinearRegression().fit(X4, Yall)
print(model4.score(X4, Yall))
print(model4.coef_)

print('only 2 eat')
model4_eat = LinearRegression().fit(X4, Yeat)
print(model4_eat.score(X4, Yeat))
print(model4_eat.coef_)

print('only NAO+ and PT')
model4_first = LinearRegression().fit(X4, Yfirst)
print(model4_first.score(X4, Yfirst))
print(model4_first.coef_)

# MODEL 4
print('MODEL 5!')
print('8 regimes')
model5 = LinearRegression().fit(X5, Yall)
print(model5.score(X5, Yall))
print(model5.coef_)

# MODEL 4
print('MODEL 6!')
print('8 regimes')
model6 = LinearRegression().fit(X6, Yall)
print(model6.score(X6, Yall))
print(model6.coef_)


import statsmodels.api as sm
from scipy import stats

resall = dict()
for ii, (xii, varfii) in enumerate(zip([X, X2, X3, X4, X5, X6], [varfit1, varfit2, varfit3, varfit4, varfit5, varfit6])):
    Xgi = sm.add_constant(xii)
    pvals = []
    params = []
    rsq = []
    for ko, nam in zip(Yall.T, reg_names_area['EAT']+reg_names_area['PNA']):
        print('\n')
        print(nam)
        est = sm.OLS(ko, Xgi)
        est2 = est.fit()
        pvals.append(est2.pvalues)
        params.append(est2.params)
        rsq.append(est2.rsquared)
        #print(est2.summary())

    pvals = np.stack(pvals)
    params = np.stack(params)
    rsq = np.stack(rsq)
    resall[(ii, 'pvals')] = pvals[:, 1:]
    resall[(ii, 'params')] = params[:, 1:]
    resall[(ii, 'rsq')] = rsq


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
for j, mod, varf in zip([0,1], [model1, model2], [varfit1, varfit2]):
    if j == 0:
        ext = [0, 4, 0, 3]
    else:
        ext = [0, 4, 0, 2]
    i=1
    bau += 1
    ax = axes[bau]
    area = 'EAT'
    #ax = fig.add_subplot(2, 2, i + 2*j)
    gigi = ax.imshow(mod.coef_[:4, :].T, vmin = -0.05, vmax = 0.05, cmap = cmappa, origin = 'lower',  extent = ext, aspect = 1)
    ax.xaxis.tick_top()
    ax.set_xticks(0.5+np.arange(4), minor = False)
    #if j == 0:
    ax.set_xticklabels(reg_names_area[area], ha='center')
    ax.set_yticks(0.5+np.arange(len(varf)), minor = False)
    ax.set_yticklabels(varf, va='center')
    if j == 0: ax.set_title(area)

    bau += 1
    ax = axes[bau]
    i=2
    area = 'PNA'
    #ax = fig.add_subplot(2, 2, i + 2*j)
    gigi = ax.imshow(mod.coef_[4:, :].T, vmin = -0.05, vmax = 0.05, cmap = cmappa, origin = 'lower', extent = ext, aspect = 1)
    ax.xaxis.tick_top()
    ax.set_xticks(0.5+np.arange(4), minor = False)
    #if j == 0:
    ax.set_xticklabels(reg_names_area[area], ha='center')
    ax.set_yticks(0.5+np.arange(len(varf)), minor = False)
    ax.set_yticklabels(varf, va='center')
    if j == 0: ax.set_title(area)

#cax = fig.add_subplot(gs[6, :])
cax = plt.axes([0.1, 0.08, 0.8, 0.03])
cb = plt.colorbar(gigi, cax=cax, orientation='horizontal')
cb.ax.tick_params(labelsize=18)
cb.set_label('Regression coefficient', fontsize=20)
plt.subplots_adjust(left=0.02, bottom=0.13, right=0.98, top=0.92, wspace=0.05, hspace=0.20)

ax.text(0.05, 0.8, 'Regr. 1', horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 20)
ax.text(0.05, 0.4, 'Regr. 2', horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 20)

fig.savefig(cart_out_orig + 'Regr_model.pdf')


fig = plt.figure(figsize=(12,8))
gs = gridspec.GridSpec(7, 9)
axes = []
axes.append(fig.add_subplot(gs[:3, 1:5]))
axes.append(fig.add_subplot(gs[:3, 5:]))
axes.append(fig.add_subplot(gs[4:6, 1:5]))
axes.append(fig.add_subplot(gs[4:6, 5:]))

bau = -1
for j, mod, varf in zip([0,1], [model3, model2], [varfit3, varfit2]):
    nvars = len(varf)
    ext = [0, 4, 0, nvars]

    varf = [va if 'rGW' not in va else va[:-3] for va in varf]

    if j == 0:
        ii = 2
    else:
        ii = 1
    i=1
    bau += 1
    ax = axes[bau]
    area = 'EAT'
    #ax = fig.add_subplot(2, 2, i + 2*j)
    gigi = ax.imshow(mod.coef_[:4, :].T, vmin = -0.05, vmax = 0.05, cmap = cmappa, origin = 'lower',  extent = ext, aspect = 1)

    pvaok = resall[(ii, 'pvals')][:4, :].T
    for ix in range(4):
        for iy in range(nvars):
            if pvaok[iy, ix] < 0.05:
                ax.scatter(ix+0.5, iy+0.5, s=60, edgecolors = 'black', facecolors='white')
            elif pvaok[iy, ix] < 0.1:
                ax.scatter(ix+0.5, iy+0.5, s=20, edgecolors = 'black', facecolors='white')

    ax.xaxis.tick_top()
    ax.set_xticks(0.5+np.arange(4), minor = False)
    #if j == 0:
    ax.set_xticklabels(reg_names_area[area], ha='center')
    ax.set_yticks(0.5+np.arange(len(varf)), minor = False)
    ax.set_yticklabels(varf, va='center')
    if j == 0: ax.set_title(area)

    bau += 1
    ax = axes[bau]
    i=2
    area = 'PNA'
    #ax = fig.add_subplot(2, 2, i + 2*j)
    gigi = ax.imshow(mod.coef_[4:, :].T, vmin = -0.05, vmax = 0.05, cmap = cmappa, origin = 'lower', extent = ext, aspect = 1)
    pvaok = resall[(ii, 'pvals')][4:, :].T
    for ix in range(4):
        for iy in range(nvars):
            if pvaok[iy, ix] < 0.05:
                ax.scatter(ix+0.5, iy+0.5, s=60, edgecolors = 'black', facecolors='white')
            elif pvaok[iy, ix] < 0.1:
                ax.scatter(ix+0.5, iy+0.5, s=20, edgecolors = 'black', facecolors='white')

    ax.xaxis.tick_top()
    ax.set_xticks(0.5+np.arange(4), minor = False)
    #if j == 0:
    ax.set_xticklabels(reg_names_area[area], ha='center')
    ax.set_yticks(0.5+np.arange(len(varf)), minor = False)
    ax.set_yticklabels(varf, va='center')
    if j == 0: ax.set_title(area)

#cax = fig.add_subplot(gs[6, :])
cax = plt.axes([0.1, 0.08, 0.8, 0.03])
cb = plt.colorbar(gigi, cax=cax, orientation='horizontal')
cb.ax.tick_params(labelsize=18)
cb.set_label('Regression coefficient', fontsize=20)
plt.subplots_adjust(left=0.02, bottom=0.13, right=0.98, top=0.92, wspace=0.05, hspace=0.20)

ax.text(0.05, 0.8, 'Regr. 1', horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 20)
ax.text(0.05, 0.4, 'Regr. 2', horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 20)

fig.savefig(cart_out_orig + 'Regr_model_rGW.pdf')


fig = plt.figure(figsize=(12,8))
gs = gridspec.GridSpec(7, 9)
axes = []
axes.append(fig.add_subplot(gs[:3, 1:5]))
axes.append(fig.add_subplot(gs[:3, 5:]))
axes.append(fig.add_subplot(gs[4:6, 1:5]))
axes.append(fig.add_subplot(gs[4:6, 5:]))

bau = -1
for j, mod, varf in zip([0,1], [model5, model6], [varfit5, varfit6]):
    nvars = len(varf)
    ext = [0, 4, 0, nvars]

    varf = [va if 'rGW' not in va else va[:-3] for va in varf]

    if j == 0:
        ii = 4
    else:
        ii = 5
    i=1
    bau += 1
    ax = axes[bau]
    area = 'EAT'
    #ax = fig.add_subplot(2, 2, i + 2*j)
    gigi = ax.imshow(mod.coef_[:4, :].T, vmin = -0.05, vmax = 0.05, cmap = cmappa, origin = 'lower',  extent = ext, aspect = 1)

    pvaok = resall[(ii, 'pvals')][:4, :].T
    for ix in range(4):
        for iy in range(nvars):
            if pvaok[iy, ix] < 0.05:
                ax.scatter(ix+0.5, iy+0.5, s=60, edgecolors = 'black', facecolors='white')
            elif pvaok[iy, ix] < 0.1:
                ax.scatter(ix+0.5, iy+0.5, s=20, edgecolors = 'black', facecolors='white')

    ax.xaxis.tick_top()
    ax.set_xticks(0.5+np.arange(4), minor = False)
    #if j == 0:
    ax.set_xticklabels(reg_names_area[area], ha='center')
    ax.set_yticks(0.5+np.arange(len(varf)), minor = False)
    ax.set_yticklabels(varf, va='center')
    if j == 0: ax.set_title(area)

    bau += 1
    ax = axes[bau]
    i=2
    area = 'PNA'
    #ax = fig.add_subplot(2, 2, i + 2*j)
    gigi = ax.imshow(mod.coef_[4:, :].T, vmin = -0.05, vmax = 0.05, cmap = cmappa, origin = 'lower', extent = ext, aspect = 1)
    pvaok = resall[(ii, 'pvals')][4:, :].T
    for ix in range(4):
        for iy in range(nvars):
            if pvaok[iy, ix] < 0.05:
                ax.scatter(ix+0.5, iy+0.5, s=60, edgecolors = 'black', facecolors='white')
            elif pvaok[iy, ix] < 0.1:
                ax.scatter(ix+0.5, iy+0.5, s=20, edgecolors = 'black', facecolors='white')

    ax.xaxis.tick_top()
    ax.set_xticks(0.5+np.arange(4), minor = False)
    #if j == 0:
    ax.set_xticklabels(reg_names_area[area], ha='center')
    ax.set_yticks(0.5+np.arange(len(varf)), minor = False)
    ax.set_yticklabels(varf, va='center')
    if j == 0: ax.set_title(area)

#cax = fig.add_subplot(gs[6, :])
cax = plt.axes([0.1, 0.08, 0.8, 0.03])
cb = plt.colorbar(gigi, cax=cax, orientation='horizontal')
cb.ax.tick_params(labelsize=18)
cb.set_label('Regression coefficient', fontsize=20)
plt.subplots_adjust(left=0.02, bottom=0.13, right=0.98, top=0.92, wspace=0.05, hspace=0.20)

ax.text(0.05, 0.8, 'Regr. 1', horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 20)
ax.text(0.05, 0.4, 'Regr. 2', horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 20)

fig.savefig(cart_out_orig + 'Regr_model_rGW_v56.pdf')


pio = ''
for ii in [2,3]:
    if ii == 3: pio = '_wNAW'
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(111)
    ax.plot(resall[(ii, 'rsq')], label = 'Regr. 1', color = 'indianred')
    ax.plot(resall[(1, 'rsq')], label = 'Regr. 2', color = 'steelblue')
    ax.scatter(np.arange(8), resall[(ii, 'rsq')], color = 'indianred')
    ax.scatter(np.arange(8), resall[(1, 'rsq')], color = 'steelblue')
    ax.set_xticks(np.arange(8), minor = False)
    ax.set_xticklabels(reg_names_area['EAT']+reg_names_area['PNA'], ha='center')
    ax.legend()
    ax.set_ylabel('R squared')
    fig.savefig(cart_out_orig + 'Rsquared_rGW{}.pdf'.format(pio))

fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
ax.plot(resall[(4, 'rsq')], label = 'Regr. 1', color = 'indianred')
ax.plot(resall[(5, 'rsq')], label = 'Regr. 2', color = 'steelblue')
ax.scatter(np.arange(8), resall[(4, 'rsq')], color = 'indianred')
ax.scatter(np.arange(8), resall[(5, 'rsq')], color = 'steelblue')
ax.set_xticks(np.arange(8), minor = False)
ax.set_xticklabels(reg_names_area['EAT']+reg_names_area['PNA'], ha='center')
ax.legend()
ax.set_ylabel('R squared')
fig.savefig(cart_out_orig + 'Rsquared_rGW_v56.pdf'.format(pio))


fig = plt.figure(figsize=(12,8))
gs = gridspec.GridSpec(8, 9)
axes = []
axes.append(fig.add_subplot(gs[:4, 1:5]))
axes.append(fig.add_subplot(gs[:4, 5:]))
axes.append(fig.add_subplot(gs[5:7, 1:5]))
axes.append(fig.add_subplot(gs[5:7, 5:]))

bau = -1
for j, mod, varf in zip([0,1], [model4, model2], [varfit4, varfit2]):
    nvars = len(varf)
    ext = [0, 4, 0, nvars]

    varf = [va if 'rGW' not in va else va[:-3] for va in varf]

    if j == 0:
        ii = 3
    else:
        ii = 1
    i=1
    bau += 1
    ax = axes[bau]
    area = 'EAT'
    #ax = fig.add_subplot(2, 2, i + 2*j)
    gigi = ax.imshow(mod.coef_[:4, :].T, vmin = -0.05, vmax = 0.05, cmap = cmappa, origin = 'lower',  extent = ext, aspect = 1)

    pvaok = resall[(ii, 'pvals')][:4, :].T
    for ix in range(4):
        for iy in range(nvars):
            if pvaok[iy, ix] < 0.05:
                ax.scatter(ix+0.5, iy+0.5, s=60, edgecolors = 'black', facecolors='white')
            elif pvaok[iy, ix] < 0.1:
                ax.scatter(ix+0.5, iy+0.5, s=20, edgecolors = 'black', facecolors='white')

    ax.xaxis.tick_top()
    ax.set_xticks(0.5+np.arange(4), minor = False)
    #if j == 0:
    ax.set_xticklabels(reg_names_area[area], ha='center')
    ax.set_yticks(0.5+np.arange(len(varf)), minor = False)
    ax.set_yticklabels(varf, va='center')
    if j == 0: ax.set_title(area)

    bau += 1
    ax = axes[bau]
    i=2
    area = 'PNA'
    #ax = fig.add_subplot(2, 2, i + 2*j)
    gigi = ax.imshow(mod.coef_[4:, :].T, vmin = -0.05, vmax = 0.05, cmap = cmappa, origin = 'lower', extent = ext, aspect = 1)
    pvaok = resall[(ii, 'pvals')][4:, :].T
    for ix in range(4):
        for iy in range(nvars):
            if pvaok[iy, ix] < 0.05:
                ax.scatter(ix+0.5, iy+0.5, s=60, edgecolors = 'black', facecolors='white')
            elif pvaok[iy, ix] < 0.1:
                ax.scatter(ix+0.5, iy+0.5, s=20, edgecolors = 'black', facecolors='white')

    ax.xaxis.tick_top()
    ax.set_xticks(0.5+np.arange(4), minor = False)
    #if j == 0:
    ax.set_xticklabels(reg_names_area[area], ha='center')
    ax.set_yticks(0.5+np.arange(len(varf)), minor = False)
    ax.set_yticklabels(varf, va='center')
    if j == 0: ax.set_title(area)

#cax = fig.add_subplot(gs[6, :])
cax = plt.axes([0.1, 0.08, 0.8, 0.03])
cb = plt.colorbar(gigi, cax=cax, orientation='horizontal')
cb.ax.tick_params(labelsize=18)
cb.set_label('Regression coefficient', fontsize=20)
plt.subplots_adjust(left=0.02, bottom=0.13, right=0.98, top=0.92, wspace=0.05, hspace=0.20)

ax.text(0.05, 0.7, 'Regr. 1', horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 20)
ax.text(0.05, 0.3, 'Regr. 2', horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 20)

fig.savefig(cart_out_orig + 'Regr_model_rGW_wNAW.pdf')


for ii, (mod, varf) in enumerate(zip([model1, model2, model3, model4], [varfit1, varfit2, varfit3, varfit4])):
    fig = plt.figure(figsize=(12,8))
    gs = gridspec.GridSpec(7, 9)
    axes = []
    axes.append(fig.add_subplot(gs[:5, 1:5]))
    axes.append(fig.add_subplot(gs[:5, 5:]))
    axes.append(fig.add_subplot(gs[5:7, 1:5]))
    axes.append(fig.add_subplot(gs[5:7, 5:]))

    nvars = len(varf)

    j = 0
    bau = -1
    varf = [va if 'rGW' not in va else va[:-3] for va in varf]
    ext = [0, 4, 0, nvars]

    i=1
    bau += 1
    ax = axes[bau]
    area = 'EAT'
    gigi = ax.imshow(mod.coef_[:4, :].T, vmin = -0.05, vmax = 0.05, cmap = cmappa, origin = 'lower',  extent = ext, aspect = 1)
    ax.xaxis.tick_top()
    ax.set_xticks(0.5+np.arange(4), minor = False)

    ax.set_xticklabels(reg_names_area[area], ha='center')
    ax.set_yticks(0.5+np.arange(len(varf)), minor = False)
    ax.set_yticklabels(varf, va='center')
    if j == 0: ax.set_title(area)

    bau += 1
    ax = axes[bau]
    i=2
    area = 'PNA'

    gigi = ax.imshow(mod.coef_[4:, :].T, vmin = -0.05, vmax = 0.05, cmap = cmappa, origin = 'lower', extent = ext, aspect = 1)
    ax.xaxis.tick_top()
    ax.set_xticks(0.5+np.arange(4), minor = False)

    ax.set_xticklabels(reg_names_area[area], ha='center')
    ax.set_yticks(0.5+np.arange(len(varf)), minor = False)
    ax.set_yticklabels([])
    if j == 0: ax.set_title(area)

    cax = plt.axes([0.1, 0.08, 0.8, 0.03])
    cb = plt.colorbar(gigi, cax=cax, orientation='horizontal')
    cb.ax.tick_params(labelsize=18)
    cb.set_label('Regression coefficient', fontsize=20)
    plt.subplots_adjust(left=0.02, bottom=0.13, right=0.98, top=0.92, wspace=0.05, hspace=0.20)

    # ax.text(0.05, 0.8, 'Regr. 1', horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 20)
    # ax.text(0.05, 0.4, 'Regr. 2', horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 20)

    bau += 1
    ax = axes[bau]
    gigi = ax.imshow(resall[(ii, 'rsq')][np.newaxis, :4], vmin = 0., vmax = 0.5, cmap = 'Greys', aspect = 1)

    bau += 1
    ax = axes[bau]
    gigi = ax.imshow(resall[(ii, 'rsq')][np.newaxis, :4], vmin = 0., vmax = 0.5, cmap = 'Greys', aspect = 1)

    fig.savefig(cart_out_orig + 'Regr_model_{}.pdf'.format(ii))
