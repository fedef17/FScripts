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
reg_names_area['PNA'] = ['PT', 'PNA+', 'PNA-', 'AR']
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

# da Virna
vicar = '/home/fedefab/Scrivania/Research/Post-doc/lavori/CMIP6/virnas_tas/new/'
fils = ['Indices_{}.txt', 'Slopes_{}.txt']
nams = ['temp', 'slope']

resvi = dict()
for na, fi in zip(nams, fils):
    for ssp in ['rcp85', 'ssp585']:
        vmods, datatemp = ctl.read_from_txt(vicar + fi.format(ssp), n_skip = 2, sep = '\t', debug = True)

        vkeys = ['AA', 'UTW', 'PVt', 'PVu']
        for mod in okmods[ssp]:
            if mod.split('_')[0] in vmods:
                ind = np.where(vmods == mod.split('_')[0])[0][0]
                lin = datatemp[ind]
                for ke, te in zip(vkeys, lin):
                    resvi[(ssp, mod, na, ke)] = te
                resvi[(ssp, mod, na, 'RUTAW')] = resvi[(ssp, mod, na, 'UTW')]/resvi[(ssp, mod, na, 'AA')]

vkeys += ['RUTAW']

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


allpeas = dict()

for area in ['EAT', 'PNA']:
    for ke in vkeys:
        for livar in ['trend', 'fr_fin', 'fr_diff']:
            co = (livar, ke, area)
            x = []
            y = []
            for ssp in ['ssp585', 'rcp85']:
                xss = []
                yss = []
                for mod in okmods[ssp]:
                    if (ssp, mod, 'slope', ke) in resvi.keys() and (ssp, mod, area, livar) in varli.keys():
                        xss.append(resvi[(ssp, mod, 'slope', ke)])
                        yss.append(varli[(ssp, mod, area, livar)][0])
                x.append(xss)
                y.append(yss)

            pears, pval = stats.pearsonr(np.concatenate(x), np.concatenate(y))
            filnam = cart_corr_sig + 'corr_{}_{}_{}.pdf'.format(*co)
            pea = ctl.plotcorr_wgroups(x, y, filename = filnam, xlabel = ke, ylabel = livar, groups = ['rcp85', 'ssp585'], colors = ['blue', 'red'], single_group_corr = True)
            allpeas[co] = (pears, pval)
            if abs(pears) > 0.3:
                print('-------->', co, '{:5.2f}'.format(pears), '{:6.3e}'.format(pval))
            else:
                print(co, '{:5.2f}'.format(pears), '{:6.3e}'.format(pval))



from sklearn.linear_model import LinearRegression

#model1: AA UTW PVu
varfit1 = ['AA', 'UTW', 'PVu']
#model2: RUTAW PVt
varfit2 = ['RUTAW', 'PVt']
#model3: best of quellisopra + weights da var_ratio

# partiamo con model1
print('MODEL 1!')

na = 'slope'
varok = 'trend'

X = []
X2 = []
Y = []
Yeat = []
Ypna = []
Yfirst = []
weights = []
for ssp in ['rcp85', 'ssp585']:
    for mod in okmods[ssp]:
        if (ssp, mod, na, 'AA') in resvi:
            Xmod = np.array([resvi[(ssp, mod, na, ke)] for ke in varfit1])
            X.append(Xmod)
            Xmod = np.array([resvi[(ssp, mod, na, ke)] for ke in varfit2])
            X2.append(Xmod)
            Ymod = np.concatenate([varli[(ssp, mod, 'EAT', varok)][:-1], varli[(ssp, mod, 'PNA', varok)][:-1]])
            Y.append(Ymod)
            Yeat.append(varli[(ssp, mod, 'EAT', varok)][:-1])
            Ypna.append(varli[(ssp, mod, 'PNA', varok)][:-1])
            Yfirst.append(np.array([varli[(ssp, mod, 'EAT', varok)][0], varli[(ssp, mod, 'PNA', varok)][0]]))

            weights.append(wei)

X = np.stack(X)
X2 = np.stack(X2)
Y = np.stack(Y)
Yeat = np.stack(Yeat)
Ypna = np.stack(Ypna)
Yfirst = np.stack(Yfirst)

model1 = LinearRegression().fit(X, Y)
print(model1.score(X, Y))
print(model1.coef_)
#>>> reg.predict(np.array([[3, 5]]))

model1_eat = LinearRegression().fit(X, Yeat)
print(model1_eat.score(X, Yeat))
print(model1_eat.coef_)

model1_pna = LinearRegression().fit(X, Ypna)
print(model1_pna.score(X, Ypna))
print(model1_pna.coef_)

model1_first = LinearRegression().fit(X, Yfirst)
print(model1_first.score(X, Yfirst))
print(model1_first.coef_)

#model2
print('MODEL 2!')
model2 = LinearRegression().fit(X2, Y)
print(model2.score(X2, Y))
print(model2.coef_)

model2_eat = LinearRegression().fit(X2, Yeat)
print(model2_eat.score(X2, Yeat))
print(model2_eat.coef_)

model2_pna = LinearRegression().fit(X2, Ypna)
print(model2_pna.score(X2, Ypna))
print(model2_pna.coef_)

model2_first = LinearRegression().fit(X2, Yfirst)
print(model2_first.score(X2, Yfirst))
print(model2_first.coef_)
