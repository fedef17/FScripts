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

cart_out_orig = cart_in + 'Results_paper_last_rebase/'
ctl.mkdir(cart_out_orig)

cart_data = '/data-hobbes/fabiano/WR_CMIP6/'

file_hist_refEOF = cart_data + 'out_NEW_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refEOF_dtr.p'
file_hist = cart_data + 'out_NEW_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refCLUS_dtr_light.p'
gen_file_ssp = cart_data + 'out_NEW_cmip6_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr_histrebase.p'

yr10 = 10 # length of running mean

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
for area in ['EAT', 'PNA']:
    cart_out = cart_out_orig + area + '/'
    ctl.mkdir(cart_out)

    results_hist, results_ref = ctl.load_wrtool(file_hist.format(area))
    results_hist_refEOF, _ = ctl.load_wrtool(file_hist_refEOF.format(area))

    cart_lui = cart_in + 'Results_v5_rebase/{}_NDJFM/'.format(area)
    freqs, residtimes, patterns = pickle.load(open(cart_lui + 'allresults_dicts_{}_v3.p'.format(area), 'rb'))
    trend_ssp, residtime_ssp = pickle.load(open(cart_lui + 'trends_wrfreq_e_restime_{}.p'.format(area), 'rb'))
    seasfreq, runfreq = pickle.load(open(cart_lui + 'seasfreqs_{}_v4.p'.format(area), 'rb'))

    resssp = dict()
    for ssp in allssps:
        print('SSP '+ssp)
        results_ssp, _ = ctl.load_wrtool(gen_file_ssp.format(ssp, area))

        okmods = [ke[1] for ke in freqs if ssp in ke and 'tot50' in ke and 'all' not in ke and 'rel' not in ke]
        if ssp == 'ssp126':
            okmods = [mod for mod in okmods if mod != 'FGOALS-g3_r1i1p1f1']
        print(okmods)

        #['BCC-CSM2-MR_r1i1p1f1', 'CanESM5_r1i1p1f1', 'CESM2-WACCM_r1i1p1f1\', 'CNRM-CM6-1_r1i1p1f2', 'CNRM-ESM2-1_r1i1p1f2', 'EC-Earth3_r1i1p1f1', 'FGOALS-g3_r1i1p1f1', 'INM-CM4-8_r1i1p1f1', 'INM-CM5-0_r1i1p1f1', 'IPSL-CM6A-LR_r1i1p1f1', 'MIROC6_r1i1p1f1', 'MPI-ESM1-2-HR_r1i1p1f1', 'MRI-ESM2-0_r1i1p1f1', 'UKESM1-0-LL_r1i1p1f2']

        # For each model write down: freq diff ssp585, trend, trend significance & year of emergence, 20-yr variability (and 50, 5, 1 yr)
        # variability defined as std dev of N-yr means, trend removed. With bootstrap? mmm questionable.
        resu = open(cart_out + 'results_short_{}_{}.dat'.format(area, ssp), 'w')

        allsims = ['hist', 'ssp126', 'ssp245', 'ssp370', 'ssp585']
        tit = 5*'{:10s}'+'\n'
        tit = tit.format(*allsims)
        stringa = 5*'{:10.1f}'+'\n'

        # da Virna
        vmods, datatemp = ctl.read_from_txt(cart_in + 'virnas_tas/DT_{}.txt'.format(ssp), n_skip = 2, sep = '\t', n_stop = 27)
        vkeys = ['DT_Global', 'DT_High', 'DT_Mid', 'DT_Low', 'DT_NAT']
        tempmods = dict()
        for mod, lin in zip(vmods, datatemp):
            tempmods[mod] = lin

        okmods_mod = [ke.split('_')[0] for ke in okmods]

        deltaT = np.array([tempmods[mod][0] for mod in okmods_mod])
        #AA = np.array([tempmods[mod][1]/tempmods[mod][0] for mod in okmods_mod])
        AA = np.array([tempmods[mod][1] for mod in okmods_mod])
        print('AA', np.mean(AA), np.max(AA), np.min(AA))
        Anat = np.array([tempmods[mod][4] for mod in okmods_mod])
        print('Anat', np.mean(Anat), np.max(Anat), np.min(Anat))
        PE_grad = np.array([tempmods[mod][1]-tempmods[mod][3] for mod in okmods_mod])
        Pole_NA_grad = np.array([tempmods[mod][1]-tempmods[mod][4] for mod in okmods_mod])
        NA_EQ_grad = np.array([tempmods[mod][4]-tempmods[mod][3] for mod in okmods_mod])

        # model performance
        var_ratio = [results_hist_refEOF[mod]['var_ratio'] for mod in okmods]
        effcen_refeof = [results_hist_refEOF[mod]['eff_centroids'] for mod in okmods]
        effcen = [results_hist[mod]['eff_centroids'] for mod in okmods]

        ssp585res = dict()
        for mod, delt, aaa, ana, vrat, cen_re, cen, peg, png, neg in zip(okmods, deltaT, AA, Anat, var_ratio, effcen_refeof, effcen, PE_grad, Pole_NA_grad, NA_EQ_grad):
            ctl.printsep(resu)
            ctl.printsep(resu)
            resu.write('\n\n' + mod + '\n')
            cose = dict()
            #cose['namcorr'] = nam
            cose['deltaT'] = delt
            cose['AA'] = aaa
            cose['ANAT'] = ana
            cose['PE_grad'] = peg
            cose['Pole_NA_grad'] = png
            cose['NA_EQ_grad'] = neg

            cose['var_ratio'] = vrat
            cose['cen_rcorr'] = np.mean([ctl.Rcorr(ce1, ce2) for ce1, ce2 in zip(cen_re, cen)])

            # Frequency
            for reg in range(4):
                resu.write('Regime {} frequency (50 and 20-yr period minus 50-yr reference)\n'.format(reg))
                allres50 = [freqs[(sim, mod, 'tot50')][reg] for sim in allsims]
                allres50 = [allres50[0]] + list(np.array(allres50[1:])-allres50[0])
                resu.write(stringa.format(*allres50))
                allres20 = [freqs[(sim, mod, 'last20')][reg] for sim in allsims]
                allres20 = list(np.array(allres20)-allres50[0])
                resu.write(stringa.format(*allres20))

            cose['frNAO'] = freqs[('hist', mod, 'tot50')][0]
            cose['frSBL'] = freqs[('hist', mod, 'tot50')][1]
            cose['fdNAO50'] = freqs[(ssp, mod, 'tot50')][0]-freqs[('hist', mod, 'tot50')][0]
            cose['fdNAO20'] = freqs[(ssp, mod, 'last20')][0]-freqs[('hist', mod, 'tot50')][0]
            cose['fdSBL50'] = freqs[(ssp, mod, 'tot50')][1]-freqs[('hist', mod, 'tot50')][1]
            cose['fdSBL20'] = freqs[(ssp, mod, 'last20')][1]-freqs[('hist', mod, 'tot50')][1]

            for reg in range(4):
                resu.write('Regime {} persistence (50 and 20-yr period minus 50-yr reference)\n'.format(reg))
                allres50 = [residtimes[(sim, mod, 'mean', reg)] for sim in allsims]
                allres50 = [allres50[0]] + list(np.array(allres50[1:])-allres50[0])
                resu.write(stringa.format(*allres50))
                allres20 = [residtimes[(sim, mod, 'mean_last20', reg)] for sim in allsims]
                allres20 = list(np.array(allres20)-allres50[0])
                resu.write(stringa.format(*allres20))

            cose['perNAO'] = residtimes[('hist', mod, 'mean', 0)]
            cose['perSBL'] = residtimes[('hist', mod, 'mean', 1)]
            cose['pdNAO50'] = residtimes[(ssp, mod, 'mean', 0)] - residtimes[('hist', mod, 'mean', 0)]
            cose['pdNAO20'] = residtimes[(ssp, mod, 'mean_last20', 0)] - residtimes[('hist', mod, 'mean', 0)]
            cose['pdSBL50'] = residtimes[(ssp, mod, 'mean', 1)] - residtimes[('hist', mod, 'mean', 1)]
            cose['pdSBL20'] = residtimes[(ssp, mod, 'mean_last20', 1)] - residtimes[('hist', mod, 'mean', 1)]

            for reg in range(4):
                resu.write('Trend Regime {} seasonal frequency (and error, second line)\n'.format(reg))
                allres = [np.nan] + [trend_ssp[(sim, mod, 'trend', 'seafreq', reg)] for sim in allsims[1:]]
                resu.write(stringa.format(*allres))
                allres = [np.nan] + [trend_ssp[(sim, mod, 'errtrend', 'seafreq', reg)] for sim in allsims[1:]]
                resu.write(stringa.format(*allres))
                #
                # resu.write('Trend Regime {} 10-yr running mean of seas. frequency (and error, second line)\n'.format(reg))
                # allres = [np.nan] + [trend_ssp[(sim, mod, 'trend', 'freq10', reg)] for sim in allsims[1:]]
                # resu.write(stringa.format(*allres))
                # allres = [np.nan] + [trend_ssp[(sim, mod, 'errtrend', 'freq10', reg)] for sim in allsims[1:]]
                # resu.write(stringa.format(*allres))

            cose['trendNAO'] = 10*trend_ssp[(ssp, mod, 'trend', 'freq10', 0)]
            cose['trendSBL'] = 10*trend_ssp[(ssp, mod, 'trend', 'freq10', 1)]
            cose['trendAR'] = 10*trend_ssp[(ssp, mod, 'trend', 'freq10', 2)]
            cose['trendNAOneg'] = 10*trend_ssp[(ssp, mod, 'trend', 'freq10', 3)]
            cose['trendNAO_divDT'] = 10*trend_ssp[(ssp, mod, 'trend', 'freq10', 0)]/cose['deltaT']
            cose['trendSBL_divDT'] = 10*trend_ssp[(ssp, mod, 'trend', 'freq10', 1)]/cose['deltaT']
            cose['trendAR_divDT'] = 10*trend_ssp[(ssp, mod, 'trend', 'freq10', 2)]/cose['deltaT']
            cose['trendNAOneg_divDT'] = 10*trend_ssp[(ssp, mod, 'trend', 'freq10', 3)]/cose['deltaT']

            cose['fdNAO50'] = freqs[(ssp, mod, 'tot50')][0]-freqs[('hist', mod, 'tot50')][0]
            cose['fdSBL50'] = freqs[(ssp, mod, 'tot50')][1]-freqs[('hist', mod, 'tot50')][1]
            cose['fdAR50'] = freqs[(ssp, mod, 'tot50')][2]-freqs[('hist', mod, 'tot50')][2]
            cose['fdNAOneg50'] = freqs[(ssp, mod, 'tot50')][3]-freqs[('hist', mod, 'tot50')][3]
            cose['fdNAO50_divDT'] = cose['fdNAO50']/cose['deltaT']
            cose['fdSBL50_divDT'] = cose['fdSBL50']/cose['deltaT']
            cose['fdAR50_divDT'] = cose['fdAR50']/cose['deltaT']
            cose['fdNAOneg50_divDT'] = cose['fdNAOneg50']/cose['deltaT']

            gigi = ctl.running_mean(seasfreq[(ssp, mod, 0)], yr10, remove_nans=True)
            pio = mk.original_test(gigi)
            cose['trendNAO_pval'] = pio.p

            if pio.h:
                trendo = pio.trend
                nutrend = pio.trend

                ima = -1
                while nutrend == trendo:
                    gigi = ctl.running_mean(seasfreq[(ssp, mod, 0)][:ima], yr10, remove_nans=True)
                    pio = mk.original_test(gigi)
                    nutrend = pio.trend
                    ima = ima - 1

                cose['trendNAO_yeme'] = 2100+ima
            else:
                cose['trendNAO_yeme'] = np.nan


            gigi = ctl.running_mean(seasfreq[(ssp, mod, 1)], yr10, remove_nans=True)
            pio = mk.original_test(gigi)
            cose['trendSBL_pval'] = pio.p

            if pio.h:
                trendo = pio.trend
                nutrend = pio.trend

                ima = -1
                while nutrend == trendo:
                    gigi = ctl.running_mean(seasfreq[(ssp, mod, 1)][:ima], yr10, remove_nans=True)
                    pio = mk.original_test(gigi)
                    nutrend = pio.trend
                    ima = ima - 1

                cose['trendSBL_yeme'] = 2100+ima
            else:
                cose['trendSBL_yeme'] = np.nan

            # variability
            years = np.arange(2015, 2100)
            seafr_dtr = seasfreq[(ssp, mod, 0)]-trend_ssp[(ssp, mod, 'trend', 'seafreq', 0)]*years
            cose['frvar_1yr'] = np.std(seafr_dtr)

            # for ye in [5, 10, 20, 50]:
            #     okse = ctl.running_mean(seafr_dtr, ye, remove_nans = True)
            #     cose['frvar_{}yr'.format(ye)] = np.std(okse)

            for ye in [5, 10, 20, 50]:
                cose['frvar_{}yr'.format(ye)] = ctl.simple_bootstrap_err(seafr_dtr, ye)

            ssp585res[mod] = cose
            ctl.printsep(resu)


        resu.close()

        resssp[ssp] = ssp585res

        coseall = dict()
        for ke in cose:
            coseall[ke] = np.mean([ssp585res[mod][ke] for mod in okmods])
        ssp585res['MMM'] = coseall

        cose = dict()
        cose['frNAO'] = results_ref['freq_clus'][0]
        cose['frSBL'] = results_ref['freq_clus'][1]
        cose['perNAO'] = np.mean(results_ref['resid_times'][0])
        cose['perSBL'] = np.mean(results_ref['resid_times'][1])

        seafr_dtr = seasfreq[('hist', 'ref', 0)]
        cose['frvar_1yr'] = np.std(seafr_dtr)
        for ye in [5, 10, 20, 50]:
            cose['frvar_{}yr'.format(ye)] = ctl.simple_bootstrap_err(seafr_dtr, ye)

        for ke in ssp585res[okmods[0]]:
            if ke not in cose:
                cose[ke] = np.nan

        ssp585res['ref'] = cose

    pickle.dump(resssp, open(cart_out + 'resu_allssps.p', 'wb'))
    #print(ssp585res)

    #freqhistNAO  freqdiff50NAO  freqdiff20NAO  trendNAO  yearofemergence?  20-yr variab  corrwithNAM DeltaTAS  finalAA
    allke = ['frNAO', 'fdNAO50', 'fdNAO20', 'trendNAO', 'trendNAO_pval', 'trendNAO_yeme', 'frvar_5yr', 'frvar_20yr', 'perNAO', 'pdNAO50', 'pdNAO20', 'deltaT', 'AA', 'ANAT']
    allnams = ['model', 'fr', 'fd50', 'fd20', 'trend', 'pval', 'yeme', 'var5', 'var20', 'pers', 'pd50', 'pd20', 'DT', 'AA', 'ANAT']
    stringa = '{:15s}'+ 4*'{:8.1f}' + '{:10.1e}' + '{:10.0f}' + 5*'{:8.1f}'+ 3*'{:6.2f}' + '\n'
    stringavirg = '{:15s},'+ 4*'{:8.1f},' + '{:10.1e},' + '{:10.0f},' + 5*'{:8.1f},' + 2*'{:6.2f},' +'{:6.2f}'
    stringatit = '{:15s}' + 4*'{:>8s}'+2*'{:>10s}'+5*'{:>8s}'+ 3*'{:>6s}'+'\n'
    stringatitvirg = '{:15s},' + 4*'{:>8s},'+2*'{:>10s},'+5*'{:>8s},'+ 2*'{:>6s},'+'{:>6s}'

    resu = open(cart_out + 'table_ssp585_{}.dat'.format(area), 'w')

    resu.write('First regime\n')
    resu.write(stringatit.format(*allnams))
    print(stringatitvirg.format(*allnams))

    for mod in okmods + ['MMM', 'ref']:
        bau = []
        bau.append(mod.split('_')[0])
        for ke in allke:
            bau.append(ssp585res[mod][ke])

        resu.write(stringa.format(*bau))
        print(stringavirg.format(*bau))

    ctl.printsep(resu)

    resu.write('Second regime\n')
    resu.write(stringatit.format(*allnams))

    allke = ['frSBL', 'fdSBL50', 'fdSBL20', 'trendSBL', 'trendSBL_pval', 'trendSBL_yeme', 'frvar_5yr', 'frvar_20yr', 'perSBL', 'pdSBL50', 'pdSBL20', 'deltaT', 'AA', 'ANAT']

    for mod in okmods + ['MMM', 'ref']:
        bau = []
        bau.append(mod.split('_')[0])
        for ke in allke:
            bau.append(ssp585res[mod][ke])

        resu.write(stringa.format(*bau))
        print(stringavirg.format(*bau))

    resu.close()

    # Guardiamo anche le correlazioni va là.
    cart_corr = cart_out + 'corrplots_allssps_v2/'
    ctl.mkdir(cart_corr)

    coppie = [('fdNAO50', 'deltaT'), ('fdNAO50', 'AA'), ('fdNAO50', 'ANAT'), ('fdSBL50', 'deltaT'), ('fdSBL50', 'AA'), ('fdSBL50', 'ANAT'), ('fdAR50', 'deltaT'), ('fdAR50', 'AA'), ('fdAR50', 'ANAT'), ('fdNAOneg50', 'deltaT'), ('fdNAOneg50', 'AA'), ('fdNAOneg50', 'ANAT'), ('trendNAO', 'deltaT'), ('trendSBL', 'deltaT'), ('trendAR', 'deltaT'), ('trendNAOneg', 'deltaT'), ('trendNAO', 'AA'), ('trendSBL', 'AA'), ('trendAR', 'AA'), ('trendNAOneg', 'AA'), ('trendNAO', 'ANAT'), ('trendSBL', 'ANAT'), ('trendAR', 'ANAT'), ('trendNAOneg', 'ANAT'), ('trendNAO', 'var_ratio'), ('deltaT', 'var_ratio'), ('trendNAO', 'cen_rcorr'), ('var_ratio', 'cen_rcorr'), ('frNAO', 'trendNAO'), ('AA', 'ANAT'), ('PE_grad', 'trendNAO'), ('PE_grad', 'trendNAOneg'), ('Pole_NA_grad', 'trendNAO'), ('Pole_NA_grad', 'trendNAOneg'), ('NA_EQ_grad', 'trendNAO'), ('NA_EQ_grad', 'trendNAOneg'), ('PE_grad', 'trendSBL'), ('PE_grad', 'trendAR'), ('Pole_NA_grad', 'trendSBL'), ('Pole_NA_grad', 'trendAR'), ('NA_EQ_grad', 'trendSBL'), ('NA_EQ_grad', 'trendAR'), ('fdNAO50', 'trendNAO'), ('fdNAOneg50', 'trendNAOneg'), ('fdSBL50', 'trendSBL'), ('fdAR50', 'trendAR')]

    allpeas = dict()
    for co in coppie:
        x = []
        y = []
        for ssp in allssps:
            okmods = [ke for ke in resssp[ssp].keys() if ke not in ['ref', 'MMM']]
            xss = [resssp[ssp][mod][co[0]] for mod in okmods]
            yss = [resssp[ssp][mod][co[1]] for mod in okmods]
            #print(ssp, np.mean(xss), np.mean(yss))
            x.append(xss)
            y.append(yss)
        filnam = cart_corr + 'corr_{}_{}.pdf'.format(co[0], co[1])
        pea = ctl.plotcorr_wgroups(x, y, filename = filnam, xlabel = co[0], ylabel = co[1], groups = allssps, colors = colssp)
        pears, pval = stats.pearsonr(np.concatenate(x), np.concatenate(y))
        allpeas[co] = (pears, pval)

    #print(allpeas)

    # cart_corr = cart_corr + 'senzaINMCM4/'
    # ctl.mkdir(cart_corr)
    # okmods_m1 = [mod for mod in okmods if 'INM-CM4-8' not in mod]
    # for co in coppie:
    #     x = [ssp585res[mod][co[0]] for mod in okmods_m1]
    #     y = [ssp585res[mod][co[1]] for mod in okmods_m1]
    #     filnam = cart_corr + 'corr_{}_{}_senzaINMCM4.pdf'.format(co[0], co[1])
    #     pea = ctl.plotcorr(x, y, filename = filnam, xlabel = co[0], ylabel = co[1])
    #     pears, pval = stats.pearsonr(x, y)
    #     allpeas[tuple(list(co)+['senzaINMCM4'])] = (pears, pval)

    #print(allpeas)
    pickle.dump(allpeas, open(cart_corr+'allcorrs.p', 'wb'))

    print('SIGNIFICANT CORRS?')
    for co in coppie:
        if allpeas[co][1] < 0.1:
            print(co, allpeas[co])

    # for co in coppie:
    #     if allpeas[tuple(list(co)+['senzaINMCM4'])][1] < 0.1:
    #         print(co, 'senzaINMCM4', allpeas[tuple(list(co)+['senzaINMCM4'])])

    print('HIGH CORRS?')
    for co in coppie:
        if np.abs(allpeas[co][0]) > 0.4:
            print(co, allpeas[co])

    # for co in coppie:
    #     if np.abs(allpeas[tuple(list(co)+['senzaINMCM4'])][0]) > 0.4:
    #         print(co, 'senzaINMCM4', allpeas[tuple(list(co)+['senzaINMCM4'])])

    #########################################################

    # distribution of seasonal frequency
    allsims = ['hist', 'ssp126', 'ssp245', 'ssp370', 'ssp585']
    colsim = ctl.color_set(5)

    frgri = np.arange(0, 80)

    for nye in [20, 50]:
        fig = plt.figure(figsize = (16,12))
        axes = []
        calcpdfs = dict()
        for reg in range(4):
            ax = fig.add_subplot(2, 2, reg+1)

            distrib = dict()
            for sim, col in zip(allsims, colsim):
                distrib[sim] = np.concatenate([seasfreq[(sim, ke, reg)][-nye:] for ke in okmods])

                pdfhi = ctl.calc_pdf(distrib[sim])
                calcpd = pdfhi(frgri)
                calcpd = calcpd/np.sum(calcpd)
                calcpdfs[(sim, reg)] = calcpd
                ax.plot(frgri, calcpd, color = col, linewidth = 2)
                axes.append(ax)

        ctl.adjust_ax_scale(axes)
        ctl.custom_legend(fig, colsim, allsims, ncol = 3)
        fig.savefig(cart_out + 'Freq_pdf_{}_{}y.pdf'.format(area, nye))

        fig = plt.figure(figsize = (16,12))
        axes = []
        for reg in range(4):
            ax = fig.add_subplot(2, 2, reg+1)
            for sim, col in zip(allsims[1:], colsim[1:]):
                ax.plot(frgri, calcpdfs[(sim, reg)]-calcpdfs[('hist', reg)], color = col, linewidth = 2)
                axes.append(ax)

        ctl.adjust_ax_scale(axes)

        ctl.custom_legend(fig, colsim[1:], allsims[1:], ncol = 3)
        fig.savefig(cart_out + 'Freq_pdfdiff_{}_{}y.pdf'.format(area, nye))
