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

plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 28
plt.rcParams['axes.labelsize'] = 22

#############################################################################
cart_v5 = '/home/fabiano/Research/lavori/CMIP6/Results_v5_rebase/{}_NDJFM/'
cart_cmip5 = '/home/fabiano/Research/lavori/CMIP6/Results_cmip5/{}_NDJFM/'

yr10 = 10 # length of running mean
#dtrtyp = 'light'
dtrtyp = 'histrebase'
cart_out_orig = '/home/fabiano/Research/lavori/CMIP6/Results_v6_eceens/'
ctl.mkdir(cart_out_orig)

cart_in = '/data-hobbes/fabiano/WR_CMIP6/'
#file_hist = cart_in + 'out_NEW_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refCLUS_dtr_light.p'
file_hist = cart_in + 'out_eceens_hist_NDJFM_{}_4clus_4pcs_1964-2014_refCLUS_dtr.p'
file_hist_refEOF = cart_in + 'out_eceens_hist_NDJFM_{}_4clus_4pcs_1964-2014_refEOF_dtr.p'
gen_file_ssp = cart_in + 'out_eceens_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr_reb.p'

numclus = 4
reg_names_area = dict()
reg_names_area['EAT'] = ['NAO+', 'SBL', 'AR', 'NAO-']
reg_names_area['PNA'] = ['PT', 'PNA+', 'PNA-', 'BR']

clatlo = dict()
clatlo['EAT'] = (70., -20.)
clatlo['PNA'] = (70., -120.)

#allssps = 'ssp119 ssp126 ssp245 ssp370 ssp585'.split()
#allssps = 'ssp126 ssp245 ssp370 ssp585'.split()
allssps = ['ssp585']

ttests = dict()
for area in ['EAT', 'PNA']:
    cart_out = cart_out_orig + '{}_NDJFM/'.format(area)
    ctl.mkdir(cart_out)

    results_hist, results_ref = ctl.load_wrtool(file_hist.format(area))
    res_hist_refEOF, _ = ctl.load_wrtool(file_hist_refEOF.format(area))

    # Erasing incomplete runs
    for ke in tuple(results_hist.keys()):
        if len(results_hist[ke]['labels']) < 7000:
            del results_hist[ke]

    results_ssp = dict()
    for ssp in allssps:
        results_ssp[ssp], _ = ctl.load_wrtool(gen_file_ssp.format(ssp, area, dtrtyp))

        # Erasing incomplete runs
        for ke in tuple(results_ssp[ssp].keys()):
            if len(results_ssp[ssp][ke]['labels']) < 12000:
                del results_ssp[ssp][ke]
            elif len(results_ssp[ssp][ke]['labels']) > 13000:
                # there is some duplicated year
                labs, dats = ctl.seasonal_set(results_ssp[ssp][ke]['labels'], results_ssp[ssp][ke]['dates'], None)
                pcs, dats = ctl.seasonal_set(results_ssp[ssp][ke]['pcs'], results_ssp[ssp][ke]['dates'], None)
                yeas = np.array([da[0].year for da in dats])
                labs_ok = []
                dats_ok = []
                pcs_ok = []
                for ye in np.arange(2015, 2100):
                    okse = np.where(yeas == ye)[0][0]
                    labs_ok.append(labs[okse])
                    dats_ok.append(dats[okse])
                    pcs_ok.append(pcs[okse])
                results_ssp[ssp][ke]['labels'] = np.concatenate(labs_ok)
                results_ssp[ssp][ke]['dates'] = np.concatenate(dats_ok)
                results_ssp[ssp][ke]['pcs'] = np.concatenate(pcs_ok)

    mod_hist = np.unique([cos.split('_')[0] for cos in results_hist.keys()])
    mod_ssp = dict()
    for ssp in allssps:
        mod_ssp[ssp] = np.unique([cos.split('_')[0] for cos in results_ssp[ssp].keys()])
        print(ssp, mod_ssp[ssp])

    okmods_hist = list(results_hist.keys())
    okmods_ssp = list(results_ssp[ssp].keys())
    okmods = okmods_ssp

    runfreq = dict()
    seasfreq = dict()

    fig = plt.figure(figsize = (16,12))
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg+1)
        cosi = []
        for mem in okmods_hist:
            seasfr, yr = ctl.calc_seasonal_clus_freq(results_hist[mem]['labels'], results_hist[mem]['dates'], numclus)
            seasfreq[('hist', mem, reg)] = seasfr[reg, :]
            seas10 = np.array(ctl.running_mean(seasfr[reg, :], yr10))
            ax.plot(yr, seas10)
            cosi.append(seas10)
        coso = np.mean(cosi, axis = 0)
        runfreq[('hist', reg)] = coso
        ax.plot(yr, coso, color = 'black', linewidth = 3)
        ax.set_title(reg_names_area[area][reg])

    fig.savefig(cart_out + 'models_freq10_{}_hist.pdf'.format(area, ssp))

    trend_ssp = dict()
    residtime_ssp = dict()

    for ssp in allssps:
        fig = plt.figure(figsize = (16,12))
        for reg in range(4):
            ax = fig.add_subplot(2, 2, reg+1)
            cosi = []
            for mem in okmods_ssp:
                if mem not in results_ssp[ssp].keys(): continue
                seasfr, yr = ctl.calc_seasonal_clus_freq(results_ssp[ssp][mem]['labels'], results_ssp[ssp][mem]['dates'], numclus)
                seasfreq[(ssp, mem, reg)] = seasfr[reg, :]
                seas10 = np.array(ctl.running_mean(seasfr[reg, :], yr10))
                ax.plot(yr, seas10, label = mem)
                cosi.append(seas10)

            coso = np.mean(cosi, axis = 0)
            runfreq[(ssp, reg)] = coso
            ax.plot(yr, coso, color = 'black', linewidth = 3)
            print(yr[0], runfreq[('hist', reg)][-1])
            ax.scatter(yr[0], np.nanmean(runfreq[('hist', reg)]), s = 40, color = 'black')
            ax.scatter(yr[0], runfreq[('hist', reg)][~np.isnan(runfreq[('hist', reg)])][-1], s = 20, color = 'red')
            ax.set_xlim(2010, 2100)
            ax.set_title(reg_names_area[area][reg])

        fig.savefig(cart_out + 'models_freq10_{}_{}.pdf'.format(area, ssp))

        for mem in okmods_ssp:
            if mem not in results_ssp[ssp].keys(): continue
            seasfr, yr = ctl.calc_seasonal_clus_freq(results_ssp[ssp][mem]['labels'], results_ssp[ssp][mem]['dates'], numclus)

            for reg in range(4):
                m, c, err_m, err_c = ctl.linear_regre_witherr(yr, seasfr[reg, :])
                trend_ssp[(ssp, mem, 'trend', 'seafreq', reg)] = m
                trend_ssp[(ssp, mem, 'errtrend', 'seafreq', reg)] = err_m

                seas10 = ctl.running_mean(seasfr[reg, :], yr10)
                m, c, err_m, err_c = ctl.linear_regre_witherr(np.array(yr[~np.isnan(seas10)]), np.array(seas10[~np.isnan(seas10)]))
                trend_ssp[(ssp, mem, 'trend', 'freq10', reg)] = m
                trend_ssp[(ssp, mem, 'errtrend', 'freq10', reg)] = err_m

            # devo fare ogni dieci anni e selezionare
            restimem = dict()
            for reg in range(4):
                restimem[reg] = []

            for ye in np.arange(2015, 2091, 5):
                dat1 = pd.Timestamp('09-01-{}'.format(ye)).to_pydatetime()
                dat2 = pd.Timestamp('04-01-{}'.format(ye+10)).to_pydatetime()
                labs, dats = ctl.sel_time_range(results_ssp[ssp][mem]['labels'], results_ssp[ssp][mem]['dates'], (dat1, dat2))

                resti, _, _ = ctl.calc_regime_residtimes(labs, dats)
                for reg in range(4):
                    restimem[reg].append(np.mean(resti[reg]))

            for reg in range(4):
                m, c, err_m, err_c = ctl.linear_regre_witherr(np.arange(2015, 2091, 5), np.array(restimem[reg]))
                residtime_ssp[(ssp, mem, 'trend', reg)] = m
                residtime_ssp[(ssp, mem, 'errtrend', reg)] = err_m


    for ssp in allssps:
        for reg in range(4):
            trend_ssp[(ssp, 'all', 'trend', 'seafreq', reg)] = np.array([trend_ssp[(ssp, mem, 'trend', 'seafreq', reg)] for mem in okmods_ssp])
            trend_ssp[(ssp, 'all', 'errtrend', 'seafreq', reg)] = np.array([trend_ssp[(ssp, mem, 'errtrend', 'seafreq', reg)] for mem in okmods_ssp])
            trend_ssp[(ssp, 'all', 'trend', 'freq10', reg)] = np.array([trend_ssp[(ssp, mem, 'trend', 'freq10', reg)] for mem in okmods_ssp])
            trend_ssp[(ssp, 'all', 'errtrend', 'freq10', reg)] = np.array([trend_ssp[(ssp, mem, 'errtrend', 'freq10', reg)] for mem in okmods_ssp])
            residtime_ssp[(ssp, 'all', 'trend', reg)] = np.array([residtime_ssp[(ssp, mem, 'trend', reg)] for mem in okmods_ssp])
            residtime_ssp[(ssp, 'all', 'errtrend', reg)] = np.array([residtime_ssp[(ssp, mem, 'errtrend', reg)] for mem in okmods_ssp])


    pickle.dump([trend_ssp, residtime_ssp], open(cart_out + 'trends_wrfreq_e_restime_{}.p'.format(area), 'wb'))

    #allsims = ['hist', 'ssp126', 'ssp245', 'ssp370', 'ssp585']
    allsims = ['hist', 'ssp585']
    colsim = [ctl.color_set(5)[0], ctl.color_set(5)[-1]]

    reg_names = reg_names_area[area]

    kess = [ke for ke in trend_ssp.keys() if 'all' in ke]
    print(kess)

    for pio in ['trend', 'errtrend']:
        for cos in ['seafreq', 'freq10']:
            fig = plt.figure(figsize = (16,12))
            axes = []
            for reg in range(4):
                ax = fig.add_subplot(2, 2, reg + 1)
                axes.append(ax)

                allpercs = dict()
                allpercs['mean'] = [np.mean(trend_ssp[(ssp, 'all', pio, cos, reg)]) for ssp in allsims[1:]]

                for nu in [10, 25, 50, 75, 90]:
                    allpercs['p{}'.format(nu)] = [np.percentile(trend_ssp[(ssp, 'all', pio, cos, reg)], nu) for ssp in allsims[1:]]

                allpercs['ens_min'] = [np.min(trend_ssp[(ssp, 'all', pio, cos, reg)]) for ssp in allsims[1:]]
                allpercs['ens_max'] = [np.max(trend_ssp[(ssp, 'all', pio, cos, reg)]) for ssp in allsims[1:]]

                ctl.boxplot_on_ax(ax, allpercs, allsims[1:], colsim[1:], plot_mean = False, plot_ensmeans = False, plot_minmax = True)

                ax.axhline(0, color = 'gray', linewidth = 0.5)
                ax.set_xticks([])
                ax.set_title(reg_names[reg])

            ctl.adjust_ax_scale(axes)

            ctl.custom_legend(fig, colsim[1:], allsims[1:], ncol = 4)
            #fig.suptitle('30yr bsp of WR freq. in 2050-2100 wrt 1964-2014')
            fig.savefig(cart_out + '{}_WRfreq_allssp_{}_{}.pdf'.format(pio, area, cos))


        fig = plt.figure(figsize = (16,12))
        axes = []
        for reg in range(4):
            ax = fig.add_subplot(2, 2, reg + 1)
            axes.append(ax)

            allpercs = dict()
            allpercs['mean'] = [np.mean(residtime_ssp[(ssp, 'all', pio, reg)]) for ssp in allsims[1:]]

            for nu in [10, 25, 50, 75, 90]:
                allpercs['p{}'.format(nu)] = [np.percentile(residtime_ssp[(ssp, 'all', pio, reg)], nu) for ssp in allsims[1:]]

            allpercs['ens_min'] = [np.min(residtime_ssp[(ssp, 'all', pio, reg)]) for ssp in allsims[1:]]
            allpercs['ens_max'] = [np.max(residtime_ssp[(ssp, 'all', pio, reg)]) for ssp in allsims[1:]]

            ctl.boxplot_on_ax(ax, allpercs, allsims[1:], colsim[1:], plot_mean = False, plot_ensmeans = False, plot_minmax = True)

            ax.axhline(0, color = 'gray', linewidth = 0.5)
            ax.set_xticks([])
            ax.set_title(reg_names[reg])

        ctl.adjust_ax_scale(axes)

        ctl.custom_legend(fig, colsim[1:], allsims[1:], ncol = 4)
        #fig.suptitle('30yr bsp of WR freq. in 2050-2100 wrt 1964-2014')
        fig.savefig(cart_out + '{}_residtimes_allssp_{}_{}.pdf'.format(pio, area, cos))

    fig = plt.figure(figsize = (16,12))
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg+1)
        for col, ssp in zip(colsim[1:], allssps):
            ax.plot(yr, runfreq[(ssp, reg)], label = ssp, color = col)
        ax.set_title(reg_names_area[area][reg])

    ctl.custom_legend(fig, colsim[1:], allsims[1:], ncol = 4)
    fig.savefig(cart_out + 'allssps_freq10_{}.pdf'.format(area))


    seasfr, yr_ref = ctl.calc_seasonal_clus_freq(results_ref['labels'], results_ref['dates'], numclus)
    for reg in range(4):
        seasfreq[('hist', 'ref', reg)] = seasfr[reg, :]

    yr = np.arange(1965, 2100)
    for ssp in allssps:
        fig = plt.figure(figsize = (16,12))
        for reg in range(4):
            ax = fig.add_subplot(2, 2, reg+1)
            cosi = []
            for mem in okmods_hist:
                cosi.append(seasfreq[('hist', mem, reg)])
            coso1 = np.mean(cosi, axis = 0)
            coso1err = np.std(cosi, axis = 0)
            cosi = []
            for mem in okmods_ssp:
                cosi.append(seasfreq[(ssp, mem, reg)])
            coso2 = np.mean(cosi, axis = 0)
            coso2err = np.std(cosi, axis = 0)

            coso = np.concatenate([coso1, coso2])
            seas20 = np.array(ctl.running_mean(coso), 20)
            runfreq[(ssp, 'run20', reg)] = seas20
            coserr = np.std(cosi, axis = 0)
            cosoerr = np.concatenate([coso1err, coso2err])
            seas20err = np.array(ctl.running_mean(cosoerr), 20)
            runfreq[(ssp, 'run20err', reg)] = seas20err
            ax.fill_between(yr, seas20-seas20err, seas20+seas20err, color = 'steelblue', alpha = 0.3)
            ax.plot(yr, coso, color = 'black', linewidth = 3)

            seas20ref = np.array(ctl.running_mean(seasfreq[('hist', 'ref', reg)], 20))
            ax.plot(yr_ref, seas20ref, color = 'red', linewidth = 2, linestyle = '--')

            ax.set_title(reg_names_area[area][reg])
            ax.axvline(2015, color = 'lightslategray', linewidth = 0.2)

        fig.savefig(cart_out + 'long_freq20_{}_{}_e_hist.pdf'.format(area, ssp))

    plt.close('all')

    pickle.dump([seasfreq, runfreq], open(cart_out + 'seasfreqs_{}_v4.p'.format(area), 'wb'))

    #################### la parte di v2
    freqs = dict() # tot50 e last20
    residtimes = dict() # mean e p90
    patterns = dict()

    for mem in okmods:
        freqs[('hist', mem, 'tot50')] = ctl.calc_clus_freq(results_hist[mem]['labels'], numclus)

        for reg in range(numclus):
            alltimes = results_hist[mem]['resid_times'][reg]
            residtimes[('hist', mem, 'mean', reg)] = np.mean(alltimes)
            residtimes[('hist', mem, 'p90', reg)] = np.percentile(alltimes, 90)

        patterns[('hist', mem, 'tot50')] = results_hist[mem]['eff_centroids']

        dat1 = pd.Timestamp('09-01-1995').to_pydatetime()
        dat2 = pd.Timestamp('04-01-2014').to_pydatetime()
        alllabs_20, dats = ctl.sel_time_range(results_hist[mem]['labels'], results_hist[mem]['dates'], (dat1, dat2))
        pcs, dats = ctl.sel_time_range(results_hist[mem]['pcs'], results_hist[mem]['dates'], (dat1, dat2))

        alltimes_20, _, _ = ctl.calc_regime_residtimes(alllabs_20, dats)
        effcen = ctl.calc_effective_centroids(pcs, alllabs_20, numclus)

        for reg in range(numclus):
            residtimes[('hist', mem, 'mean_last20', reg)] = np.mean(alltimes_20[reg])
            residtimes[('hist', mem, 'p90_last20', reg)] = np.percentile(alltimes_20[reg], 90)

        patterns[('hist', mem, 'last20')] = effcen
        freqs[('hist', mem, 'last20')] = ctl.calc_clus_freq(alllabs_20, numclus)

    for ssp in allssps:
        for mem in okmods:
            if mem not in results_ssp[ssp].keys(): continue
            ## Attach all members labels
            dat1 = pd.Timestamp('09-01-2050').to_pydatetime()
            dat2 = pd.Timestamp('04-01-2100').to_pydatetime()
            labs, dats = ctl.sel_time_range(results_ssp[ssp][mem]['labels'], results_ssp[ssp][mem]['dates'], (dat1, dat2))

            restim, _, _ = ctl.calc_regime_residtimes(labs, dats)
            pcs, dats = ctl.sel_time_range(results_ssp[ssp][mem]['pcs'], results_ssp[ssp][mem]['dates'], (dat1, dat2))
            effcen = ctl.calc_effective_centroids(pcs, labs, numclus)

            patterns[(ssp, mem, 'tot50')] = effcen
            freqs[(ssp, mem, 'tot50')] = ctl.calc_clus_freq(labs, numclus)

            for reg in range(numclus):
                residtimes[(ssp, mem, 'mean', reg)] = np.mean(restim[reg])
                residtimes[(ssp, mem, 'p90', reg)] = np.percentile(restim[reg], 90)


            dat1 = pd.Timestamp('09-01-2081').to_pydatetime()
            dat2 = pd.Timestamp('04-01-2100').to_pydatetime()
            labs, dats = ctl.sel_time_range(results_ssp[ssp][mem]['labels'], results_ssp[ssp][mem]['dates'], (dat1, dat2))

            restim, _, _ = ctl.calc_regime_residtimes(labs, dats)

            pcs, dats = ctl.sel_time_range(results_ssp[ssp][mem]['pcs'], results_ssp[ssp][mem]['dates'], (dat1, dat2))
            effcen = ctl.calc_effective_centroids(pcs, labs, numclus)
            patterns[(ssp, mem, 'last20')] = effcen

            freqs[(ssp, mem, 'last20')] = ctl.calc_clus_freq(labs, numclus)

            for reg in range(numclus):
                residtimes[(ssp, mem, 'mean_last20', reg)] = np.mean(restim[reg])
                residtimes[(ssp, mem, 'p90_last20', reg)] = np.percentile(restim[reg], 90)

    ##### salvo le distrib per ogni ssp. Quelle assolute (con tutti i modelli) e quelle relative (solo con i modelli che ci sono anche in hist)
    ### questo è cambiato, uso gli stessi modelli per tutto
    freqs[('hist', 'all', 'tot50')] = np.stack([freqs[('hist', mod, 'tot50')] for mod in okmods])
    freqs[('hist', 'all', 'last20')] = np.stack([freqs[('hist', mod, 'last20')] for mod in okmods])

    for cos in ['last20', 'tot50']:
        patterns[('hist', 'mean', cos)] = np.mean([patterns[('hist', mod, cos)] for mod in okmods], axis = 0)
        patterns[('hist', 'std', cos)] = np.std([patterns[('hist', mod, cos)] for mod in okmods], axis = 0)


    for reg in range(numclus):
        for cos in ['mean', 'p90', 'mean_last20', 'p90_last20']:
            residtimes[('hist', 'all', cos, reg)] = np.array([residtimes[('hist', mod, cos, reg)] for mod in okmods])

    for ssp in allssps:
        freqs[(ssp, 'all', 'tot50')] = np.stack([freqs[(ssp, mod, 'tot50')] for mod in okmods])
        freqs[(ssp, 'all', 'last20')] = np.stack([freqs[(ssp, mod, 'last20')] for mod in okmods])

        for cos in ['last20', 'tot50']:
            patterns[(ssp, 'mean', cos)] = np.mean([patterns[(ssp, mod, cos)] for mod in okmods], axis = 0)
            patterns[(ssp, 'std', cos)] = np.std([patterns[(ssp, mod, cos)] for mod in okmods], axis = 0)

            patterns[(ssp, 'mean_diff', cos)] = np.mean([patterns[(ssp, mod, cos)]-patterns[('hist', mod, cos)] for mod in okmods], axis = 0)
            patterns[(ssp, 'std_diff', cos)] = np.std([patterns[(ssp, mod, cos)]-patterns[('hist', mod, cos)] for mod in okmods], axis = 0)

        # rel
        modoks = okmods
        freqs[(ssp, 'rel', 'tot50')] = np.stack([freqs[(ssp, mod, 'tot50')]-freqs[('hist', mod, 'tot50')] for mod in modoks])
        freqs[(ssp, 'rel', 'last20')] = np.stack([freqs[(ssp, mod, 'last20')]-freqs[('hist', mod, 'last20')] for mod in modoks])

        for reg in range(numclus):
            for cos in ['mean', 'p90', 'mean_last20', 'p90_last20']:
                residtimes[(ssp, 'all', cos, reg)] = np.array([residtimes[(ssp, mod, cos, reg)] for mod in okmods])

                residtimes[(ssp, 'rel', cos, reg)] = np.array([residtimes[(ssp, mod, cos, reg)]-residtimes[('hist', mod, cos, reg)] for mod in modoks])

    for ke in patterns.keys():
        gigi = patterns[ke][..., np.newaxis, np.newaxis] * results_ref['eofs_ref_pcs'][np.newaxis, ...]
        patterns[ke] = np.sum(gigi, axis = 1)

    pickle.dump([freqs, residtimes, patterns], open(cart_out + 'allresults_dicts_{}_v3.p'.format(area), 'wb'))

####

for area in ['EAT', 'PNA']:
    cart_out = cart_out_orig + '{}_NDJFM/'.format(area)
    ctl.mkdir(cart_out)
    reg_names = reg_names_area[area]

    results_hist, results_ref = ctl.load_wrtool(file_hist.format(area))

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

    # allssps = ['ssp126', 'ssp245', 'ssp370', 'ssp585', 'rcp85_cmip5']
    # allsims = ['hist', 'ssp126', 'ssp245', 'ssp370', 'ssp585', 'rcp85_cmip5']
    # allsims_wcmip5 = ['hist', 'ssp126', 'ssp245', 'ssp370', 'ssp585', 'hist_cmip5', 'rcp85_cmip5']
    # allsimcol = ['hist', 'ssp126', 'ssp245', 'hist_cmip5', 'ssp370', 'ssp585', 'rcp85_cmip5']
    allssps = ['ssp585']
    allsims = ['hist', 'ssp585']
    allsimcol = ['hist', 'ssp126', 'ssp245', 'hist_cmip5', 'ssp370', 'ssp585', 'rcp85_cmip5']
    coldic = dict(zip(allsimcol, ctl.color_set(7)))
    colsim = [coldic[ssp] for ssp in allsims]
    colssp = [coldic[ssp] for ssp in allssps]

    print('T-TEST for {}'.format(area))
    for reg in range(4):
        print('REGIME',reg)
        for ssp in allssps:
            if ssp == 'rcp85_cmip5':
                a = freqs[('hist_cmip5', 'all', cos)][:, reg]
            else:
                a = freqs[('hist', 'all', cos)][:, reg]
            print(ssp)
            b = freqs[(ssp, 'all', cos)][:, reg]
            ttests[('freq', area, reg, ssp)] = stats.ttest_ind(a, b, equal_var = False)
            print(ttests[('freq', area, reg, ssp)])

    figall = plt.figure(figsize = (16,12))
    axes = []
    cos = 'tot50'
    for reg in range(4):
        ax = figall.add_subplot(2, 2, reg + 1)
        axes.append(ax)

        histmean = np.mean(freqs[('hist', 'all', 'tot50')][:, reg])
        allpercs = dict()
        for nu in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(nu)] = [np.percentile(freqs[(ssp, 'all', cos)][:, reg], nu) - histmean for ssp in allsims]

        ax.axhline(allpercs['p50'][0], color = 'gray', linewidth = 0.5)

        positions = list(np.arange(len(allsims)-1)*0.7)
        positions.append(positions[-1]+0.3+0.7)
        ctl.boxplot_on_ax(ax, allpercs, allsims, colsim, plot_mean = False, plot_ensmeans = False, plot_minmax = False, positions = positions)
        ax.set_xticks([])
        ax.set_title(reg_names[reg])
        #ax.text(1.0, 1.0, na, horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 20)
        if reg == 0 or reg == 2: ax.set_ylabel('Regime frequency anomaly')

        #ax.scatter(0, results_ref['freq_clus'][reg], color = 'black', marker = '*', s = 5)
        for pos, ssp in zip(positions[1:], allsims[1:]):
            if ttests[('freq', area, reg, ssp)].pvalue < 0.05:
                ax.scatter(pos, -10, color = 'black', marker = '*', s = 30)

        ax.axvline(np.mean([positions[-1], positions[-2]]), color = 'lightslategray', linewidth = 0.2, linestyle = '--')
        xli = ax.get_xlim()

    ctl.adjust_ax_scale(axes)

    ctl.custom_legend(figall, colsim, allsims, ncol = 3)
    figall.savefig(cart_out + 'WRfreq_{}_{}_FINAL.pdf'.format(area, cos))


    figall = plt.figure(figsize = (16,12))
    axes = []
    for reg in range(4):
        ax = figall.add_subplot(2, 2, reg + 1)
        allpercs = dict()
        for nu in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(nu)] = [np.percentile(trend_ssp[(ssp, 'all', 'trend', 'freq10', reg)], nu) for ssp in allssps]

        ctl.boxplot_on_ax(ax, allpercs, allssps, colssp, plot_mean = False, plot_ensmeans = False, plot_minmax = False, positions = positions[1:])

        ax.axhline(0, color = 'gray', linewidth = 0.5)
        ax.set_xticks([])
        if reg == 0 or reg == 2: ax.set_ylabel('Trend in regime frequency (1/yr)')
        axes.append(ax)
        ax.set_xlim(xli)
        ax.axvline(np.mean([positions[-1], positions[-2]]), color = 'lightslategray', linewidth = 0.2, linestyle = '--')

    ctl.adjust_ax_scale(axes)

    ctl.custom_legend(figall, colsim, allsims, ncol = 3)
    figall.savefig(cart_out + 'Trends_{}_FINAL.pdf'.format(area))
