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
for area in ['EAT']:#, 'PNA']:
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
    okmods['rcp85'] = [ke[1] for ke in freqs if 'rcp85' in ke and 'tot50' in ke and 'all' not in ke and 'rel' not in ke]
    print(okmods)

    ### model performance
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
        for area in ['EAT']:#, 'PNA']:
            trend_ssp = cose[(area, 'trend')]
            freqs = cose[(area, 'freq')]
            varli[(ssp, mod, area, 'trend')] = np.array([trend_ssp[(ssp+cos, mod, 'trend', 'freq10', reg)] for reg in range(4)])
            varli[(ssp, mod, area, 'fr_fin')] = np.array(freqs[(ssp+cos, mod, 'last20')])
            varli[(ssp, mod, area, 'fr_diff')] = np.array(freqs[(ssp+cos, mod, 'tot50')])-np.array(freqs[('hist'+cos, mod, 'tot50')])


allpeas = dict()

for area in ['EAT']:#, 'PNA']:
    for ke in vkeys:
        for livar in ['trend', 'fr_fin', 'fr_diff']:
            co = (livar, ke, area)
            x = []
            y = []
            for ssp in ['ssp585']:#'rcp85',
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
            pea = ctl.plotcorr_wgroups(x, y, filename = filnam, xlabel = ke, ylabel = livar, groups = ['ssp585'], colors = ['red'], single_group_corr = True)#['rcp85', 'ssp585'], colors = ['blue', 'red'], single_group_corr = True)
            allpeas[co] = (pears, pval)
            if abs(pears) > 0.3:
                print(co, '{:5.2f}'.format(pears))
