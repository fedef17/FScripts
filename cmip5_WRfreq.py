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

#############################################################################
if os.uname()[1] == 'hobbes':
    cart_in = '/home/fabiano/Research/lavori/CMIP6/'
elif os.uname()[1] == 'ff-clevo':
    cart_in = '/home/fedefab/Scrivania/Research/Post-doc/lavori/CMIP6/'

cart_out_orig = cart_in + 'Results_cmip5/'
ctl.mkdir(cart_out_orig)

#file_hist_refEOF = cart_in + 'cmip6_hist/out_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refEOF_dtr.p'
file_hist = cart_in + 'cmip5_hist_reb/out_cmip5_hist_reb_NDJFM_{}_4clus_4pcs_allyrs_refCLUS_dtr.p'
#file_hist_refEOF = cart_in + 'cmip6_hist/out_cmip6_hist_NDJFM_EAT_4clus_4pcs_1964-2014_refEOF.p'
gen_file_ssp = cart_in + 'cmip5_{}/out_cmip5_{}_reb_NDJFM_{}_4clus_4pcs_2005-2100_refCLUS_dtr_reb.p'

numclus = 4
reg_names_area = dict()
reg_names_area['EAT'] = ['NAO+', 'SBL', 'AR', 'NAO-']
reg_names_area['PNA'] = ['PT', 'PNA+', 'PNA-', 'AR']

#allssps = 'ssp119 ssp126 ssp245 ssp370 ssp585'.split()
#allssps = 'ssp126 ssp245 ssp370 ssp585'.split()
allssps = ['rcp85']

for area in ['EAT']:#, 'PNA']:
    freqs = dict() # tot50 e last20
    trend_ssp = dict()
    residtimes = dict()

    cart_out = cart_out_orig + '{}_NDJFM/'.format(area)
    ctl.mkdir(cart_out)

    results = pickle.load(open(file_hist.format(area), 'rb'))
    results_hist = results['models']
    results_ref = results['reference']

    yr0 = 1950
    yr1 = 2005
    allyr = np.arange(1950, 2005)
    yr = allyr

    # Erasing incomplete runs
    avlen = np.median([len(results_hist[ke]['labels']) for ke in results_hist.keys()])
    for ke in tuple(results_hist.keys()):
        if len(results_hist[ke]['labels']) < avlen-1000:
            del results_hist[ke]
        elif len(results_hist[ke]['labels']) > avlen+100:
            # there is some duplicated year
            labs, dats = ctl.seasonal_set(results_hist[ke]['labels'], results_hist[ke]['dates'], None)
            pcs, dats = ctl.seasonal_set(results_hist[ke]['pcs'], results_hist[ke]['dates'], None)
            yeas = np.array([da[0].year for da in dats])
            labs_ok = []
            dats_ok = []
            pcs_ok = []
            for ye in np.arange(1950, 2005):
                okse = np.where(yeas == ye)[0][0]
                labs_ok.append(labs[okse])
                dats_ok.append(dats[okse])
                pcs_ok.append(pcs[okse])
            results_hist[ke]['labels'] = np.concatenate(labs_ok)
            results_hist[ke]['dates'] = np.concatenate(dats_ok)
            results_hist[ke]['pcs'] = np.concatenate(pcs_ok)

    okmods = [cos for cos in results_hist.keys()]
    print(okmods)
    print(len(okmods))

    runfreq = dict()
    seasfreq = dict()

    fig = plt.figure(figsize = (16,12))
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg+1)
        cosi = []
        for mem in okmods:
            labok, datok = ctl.sel_time_range(results_hist[mem]['labels'], results_hist[mem]['dates'], ctl.range_years(yr0, yr1))
            seasfr, yr = ctl.calc_seasonal_clus_freq(labok, datok, numclus)
            seasfreq[('hist', mem, reg)] = seasfr[reg, :]
            seas20 = np.array(ctl.running_mean(seasfr[reg, :], 20))
            print(mem, len(seas20))
            if len(seas20) == len(yr):
                ax.plot(yr, seas20)
                cosi.append(seas20)
            else:
                print(mem, len(seas20), 'too short')
        coso = np.mean(cosi, axis = 0)
        runfreq[('hist', reg)] = coso
        ax.plot(yr, coso, color = 'black', linewidth = 3)
        ax.set_title(reg_names_area[area][reg])

    fig.savefig(cart_out + 'models_run20_{}_hist.pdf'.format(area))

    reg_names = reg_names_area[area]

    seasfr, yr_ref = ctl.calc_seasonal_clus_freq(results_ref['labels'], results_ref['dates'], numclus)
    for reg in range(4):
        seasfreq[('hist', 'ref', reg)] = seasfr[reg, :]

    fig = plt.figure(figsize = (16,12))
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg+1)
        cosi = []
        for mem in okmods:
            seas20 = np.array(ctl.running_mean(seasfreq[('hist', mem, reg)], 20))
            if len(seas20) == len(allyr):
                cosi.append(seas20)
            else:
                print(mem, 'too short')
            cosi.append(seas20)
        coso = np.mean(cosi, axis = 0)
        runfreq[('run20', reg)] = coso
        coserr = np.std(cosi, axis = 0)
        runfreq[('run20err', reg)] = coserr/np.sqrt(len(okmods)-1)
        ax.fill_between(yr, coso-coserr, coso+coserr, color = 'steelblue', alpha = 0.3)
        ax.plot(yr, coso, color = 'black', linewidth = 3)

        seas20ref = np.array(ctl.running_mean(seasfreq[('hist', 'ref', reg)], 20))
        ax.plot(yr_ref, seas20ref, color = 'red', linewidth = 2, linestyle = '--')

        ax.set_title(reg_names_area[area][reg])

    fig.savefig(cart_out + 'long_run20_{}_hist.pdf'.format(area))

    # now for rcp85
    print('rcp85')
    results = pickle.load(open(gen_file_ssp.format('rcp85', 'rcp85', area), 'rb'))
    results_ssp = results['models']

    avlen = np.median([len(results_ssp[ke]['labels']) for ke in results_ssp.keys()])
    for ke in tuple(results_ssp.keys()):
        if len(results_ssp[ke]['labels']) < avlen-1000:
            del results_ssp[ke]
        elif len(results_ssp[ke]['labels']) > avlen+100:
            # there is some duplicated year
            labs, dats = ctl.seasonal_set(results_ssp[ke]['labels'], results_ssp[ke]['dates'], None)
            pcs, dats = ctl.seasonal_set(results_ssp[ke]['pcs'], results_ssp[ke]['dates'], None)
            yeas = np.array([da[0].year for da in dats])
            labs_ok = []
            dats_ok = []
            pcs_ok = []
            for ye in np.arange(2005, 2100):
                okse = np.where(yeas == ye)[0][0]
                labs_ok.append(labs[okse])
                dats_ok.append(dats[okse])
                pcs_ok.append(pcs[okse])
            results_ssp[ke]['labels'] = np.concatenate(labs_ok)
            results_ssp[ke]['dates'] = np.concatenate(dats_ok)
            results_ssp[ke]['pcs'] = np.concatenate(pcs_ok)

    okmods = [cos for cos in results_ssp.keys()]
    allyr = np.arange(2005, 2101)
    yr0 = 2005
    yr1 = 2101

    fig = plt.figure(figsize = (16,12))
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg+1)
        cosi = []
        for mem in okmods:
            labok, datok = ctl.sel_time_range(results_ssp[mem]['labels'], results_ssp[mem]['dates'], ctl.range_years(yr0, yr1))
            seasfr, yr = ctl.calc_seasonal_clus_freq(labok, datok, numclus)
            seasfreq[('rcp85', mem, reg)] = seasfr[reg, :]

            m, c, err_m, err_c = ctl.linear_regre_witherr(yr, seasfr[reg, :])
            trend_ssp[('rcp85_cmip5', mem, 'trend', 'seafreq', reg)] = m
            trend_ssp[('rcp85_cmip5', mem, 'errtrend', 'seafreq', reg)] = err_m

            seas10 = ctl.running_mean(seasfr[reg, :], 10)
            m, c, err_m, err_c = ctl.linear_regre_witherr(np.array(yr[~np.isnan(seas10)]), np.array(seas10[~np.isnan(seas10)]))
            trend_ssp[('rcp85_cmip5', mem, 'trend', 'freq10', reg)] = m
            trend_ssp[('rcp85_cmip5', mem, 'errtrend', 'freq10', reg)] = err_m

            seas20 = np.array(ctl.running_mean(seasfr[reg, :], 20))
            print(mem, len(seas20))
            if len(seas20) == len(allyr):
                ax.plot(allyr, seas20)
                cosi.append(seas20)
            else:
                print(mem, len(seas20), 'too short')
        coso = np.mean(cosi, axis = 0)
        runfreq[('rcp85', reg)] = coso
        ax.plot(yr, coso, color = 'black', linewidth = 3)
        ax.set_title(reg_names_area[area][reg])

    fig.savefig(cart_out + 'models_run20_{}_rcp85.pdf'.format(area))

    reg_names = reg_names_area[area]

    fig = plt.figure(figsize = (16,12))
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg+1)
        cosi = []
        for mem in okmods:
            seas20 = np.array(ctl.running_mean(seasfreq[('rcp85', mem, reg)], 20))
            if len(seas20) == len(allyr):
                cosi.append(seas20)
            else:
                print(mem, 'too short')
        coso = np.mean(cosi, axis = 0)
        runfreq[('run20', reg)] = coso
        coserr = np.std(cosi, axis = 0)
        runfreq[('run20err', reg)] = coserr/np.sqrt(len(okmods)-1)
        ax.fill_between(yr, coso-coserr, coso+coserr, color = 'steelblue', alpha = 0.3)
        ax.plot(yr, coso, color = 'black', linewidth = 3)

        ax.set_title(reg_names_area[area][reg])

    fig.savefig(cart_out + 'long_run20_{}_rcp85.pdf'.format(area))

    numclus = 4
    for ke in results_hist.keys():
        dat1 = pd.Timestamp('09-01-1964').to_pydatetime()
        dat2 = pd.Timestamp('04-01-2005').to_pydatetime()
        labs, dats = ctl.sel_time_range(results_hist[ke]['labels'], results_hist[ke]['dates'], (dat1, dat2))
        freqs[('hist_cmip5', ke, 'tot50')] = ctl.calc_clus_freq(labs, numclus)

        restim, _, _ = ctl.calc_regime_residtimes(labs, dats)
        for reg in range(numclus):
            residtimes[('hist_cmip5', mem, 'mean', reg)] = np.mean(restim[reg])
            residtimes[('hist_cmip5', mem, 'p90', reg)] = np.percentile(restim[reg], 90)

    freqs[('hist_cmip5', 'all', 'tot50')] = np.array([freqs[('hist_cmip5', ke, 'tot50')] for ke in results_hist.keys()])

    for ke in results_ssp.keys():
        dat1 = pd.Timestamp('09-01-2050').to_pydatetime()
        dat2 = pd.Timestamp('04-01-2100').to_pydatetime()
        labs, dats = ctl.sel_time_range(results_ssp[ke]['labels'], results_ssp[ke]['dates'], (dat1, dat2))
        freqs[('rcp85_cmip5', ke, 'tot50')] = ctl.calc_clus_freq(labs, numclus)

        restim, _, _ = ctl.calc_regime_residtimes(labs, dats)
        for reg in range(numclus):
            residtimes[('rcp85_cmip5', mem, 'mean', reg)] = np.mean(restim[reg])
            residtimes[('rcp85_cmip5', mem, 'p90', reg)] = np.percentile(restim[reg], 90)

    freqs[('rcp85_cmip5', 'all', 'tot50')] = np.array([freqs[('rcp85_cmip5', ke, 'tot50')] for ke in results_ssp.keys()])

    for cos in ['mean', 'p90']:
        residtimes[('hist_cmip5', 'all', cos, reg)] = np.array([residtimes[('hist_cmip5', mod, cos, reg)] for mod in okmods])
        residtimes[('rcp85_cmip5', 'all', cos, reg)] = np.array([residtimes[('rcp85_cmip5', mod, cos, reg)] for mod in okmods])

    ssp = 'rcp85_cmip5'
    for reg in range(4):
        trend_ssp[(ssp, 'all', 'trend', 'seafreq', reg)] = np.array([trend_ssp[(ssp, mem, 'trend', 'seafreq', reg)] for mem in okmods])
        trend_ssp[(ssp, 'all', 'errtrend', 'seafreq', reg)] = np.array([trend_ssp[(ssp, mem, 'errtrend', 'seafreq', reg)] for mem in okmods])
        trend_ssp[(ssp, 'all', 'trend', 'freq10', reg)] = np.array([trend_ssp[(ssp, mem, 'trend', 'freq10', reg)] for mem in okmods])
        trend_ssp[(ssp, 'all', 'errtrend', 'freq10', reg)] = np.array([trend_ssp[(ssp, mem, 'errtrend', 'freq10', reg)] for mem in okmods])

    pickle.dump([freqs, trend_ssp, residtimes], open(cart_out + 'freqs_cmip5_{}.p'.format(area), 'wb'))

    allsims = ['hist_cmip5', 'rcp85_cmip5']
    colsim = ctl.color_set(2)

    figall = plt.figure(figsize = (16,12))
    axes = []
    cos = 'tot50'
    for reg in range(4):
        ax = figall.add_subplot(2, 2, reg + 1)
        axes.append(ax)

        allpercs = dict()
        for nu in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(nu)] = [np.percentile(freqs[(ssp, 'all', cos)][:, reg], nu) for ssp in allsims]

        ax.axhline(allpercs['p50'][0], color = 'gray', linewidth = 0.5)

        ctl.boxplot_on_ax(ax, allpercs, allsims, colsim, plot_mean = False, plot_ensmeans = False, plot_minmax = False)
        ax.set_xticks([])
        ax.set_title(reg_names[reg])
        if reg == 0: ax.set_ylabel('Regime frequency')

        ax.scatter(0, results_ref['freq_clus'][reg], color = 'black', marker = '*', s = 5)

    ctl.custom_legend(figall, colsim, allsims, ncol = 4)
    figall.savefig(cart_out + 'WRfreq_cmip5_{}.pdf'.format(area))
