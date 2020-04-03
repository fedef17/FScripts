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

cart_out_orig = cart_in + 'Results_v4/'
ctl.mkdir(cart_out_orig)

file_hist_refEOF = cart_in + 'cmip6_hist/out_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refEOF_dtr.p'
file_hist = cart_in + 'cmip6_hist/out_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refCLUS_dtr.p'
#file_hist_refEOF = cart_in + 'cmip6_hist/out_cmip6_hist_NDJFM_EAT_4clus_4pcs_1964-2014_refEOF.p'
gen_file_ssp = cart_in + 'cmip6_{}/out_cmip6_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr.p'

numclus = 4
reg_names_area = dict()
reg_names_area['EAT'] = ['NAO+', 'SBL', 'AR', 'NAO-']
reg_names_area['PNA'] = ['PT', 'PNA+', 'PNA-', 'AR']


#allssps = 'ssp119 ssp126 ssp245 ssp370 ssp585'.split()
allssps = 'ssp126 ssp245 ssp370 ssp585'.split()

area = 'EAT'
for area in ['EAT', 'PNA']:
    pdfssp = dict()

    cart_out = cart_out_orig + '{}_NDJFM/'.format(area)
    ctl.mkdir(cart_out)

    results_hist, results_ref = pickle.load(open(file_hist.format(area), 'rb'))
    res_hist_refEOF, _ = pickle.load(open(file_hist_refEOF.format(area), 'rb'))

    # Erasing incomplete runs
    for ke in tuple(results_hist.keys()):
        if len(results_hist[ke]['labels']) < 7000:
            del results_hist[ke]

    results_ssp = dict()
    for ssp in allssps:
        results_ssp[ssp] = pickle.load(open(gen_file_ssp.format(ssp, ssp, area), 'rb'))['models']

        # Erasing incomplete runs
        for ke in tuple(results_ssp[ssp].keys()):
            if len(results_ssp[ssp][ke]['labels']) < 12000:
                del results_ssp[ssp][ke]

    mod_hist = np.unique([cos.split('_')[0] for cos in results_hist.keys()])
    mod_ssp = dict()
    for ssp in allssps:
        mod_ssp[ssp] = np.unique([cos.split('_')[0] for cos in results_ssp[ssp].keys()])

    print('keeping only models used in ssps')
    print(mod_hist)
    mod_hist = [mod for mod in mod_hist if np.any([mod in mod_ssp[ssp] for ssp in allssps])]
    print(mod_hist)

    okmods = [mod for mod in results_hist.keys() if np.all([mod in results_ssp[ssp].keys() for ssp in allssps[1:]])]
    okmods.remove('EC-Earth3_r4i1p1f1')
    print(okmods)

    runfreq = dict()
    seasfreq = dict()

    fig = plt.figure(figsize = (16,12))
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg+1)
        cosi = []
        for mem in okmods:
            seasfr, yr = ctl.calc_seasonal_clus_freq(results_hist[mem]['labels'], results_hist[mem]['dates'], numclus)
            seasfreq[('hist', mem, reg)] = seasfr[reg, :]
            seas10 = np.array(ctl.running_mean(seasfr[reg, :], 10))
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
            for mem in okmods:
                if mem not in results_ssp[ssp].keys(): continue
                seasfr, yr = ctl.calc_seasonal_clus_freq(results_ssp[ssp][mem]['labels'], results_ssp[ssp][mem]['dates'], numclus)
                seasfreq[(ssp, mem, reg)] = seasfr[reg, :]
                seas10 = np.array(ctl.running_mean(seasfr[reg, :], 10))
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

        for mem in okmods:
            if mem not in results_ssp[ssp].keys(): continue
            seasfr, yr = ctl.calc_seasonal_clus_freq(results_ssp[ssp][mem]['labels'], results_ssp[ssp][mem]['dates'], numclus)

            for reg in range(4):
                m, c, err_m, err_c = ctl.linear_regre_witherr(yr, seasfr[reg, :])
                trend_ssp[(ssp, mem, 'trend', 'seafreq', reg)] = m
                trend_ssp[(ssp, mem, 'errtrend', 'seafreq', reg)] = err_m

                seas10 = ctl.running_mean(seasfr[reg, :], 10)
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
            trend_ssp[(ssp, 'all', 'trend', 'seafreq', reg)] = np.array([trend_ssp[(ssp, mem, 'trend', 'seafreq', reg)] for mem in okmods])
            trend_ssp[(ssp, 'all', 'errtrend', 'seafreq', reg)] = np.array([trend_ssp[(ssp, mem, 'errtrend', 'seafreq', reg)] for mem in okmods])
            trend_ssp[(ssp, 'all', 'trend', 'freq10', reg)] = np.array([trend_ssp[(ssp, mem, 'trend', 'freq10', reg)] for mem in okmods])
            trend_ssp[(ssp, 'all', 'errtrend', 'freq10', reg)] = np.array([trend_ssp[(ssp, mem, 'errtrend', 'freq10', reg)] for mem in okmods])
            residtime_ssp[(ssp, 'all', 'trend', reg)] = np.array([residtime_ssp[(ssp, mem, 'trend', reg)] for mem in okmods])
            residtime_ssp[(ssp, 'all', 'errtrend', reg)] = np.array([residtime_ssp[(ssp, mem, 'errtrend', reg)] for mem in okmods])

    allsims = ['hist', 'ssp126', 'ssp245', 'ssp370', 'ssp585']
    colsim = ctl.color_set(5)

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

            ctl.custom_legend(fig, colsim[1:], allsims[1:], ncol = 3)
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

        ctl.custom_legend(fig, colsim[1:], allsims[1:], ncol = 3)
        #fig.suptitle('30yr bsp of WR freq. in 2050-2100 wrt 1964-2014')
        fig.savefig(cart_out + '{}_residtimes_allssp_{}_{}.pdf'.format(pio, area, cos))

    fig = plt.figure(figsize = (16,12))
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg+1)
        for col, ssp in zip(colsim[1:], allssps):
            ax.plot(yr, runfreq[(ssp, reg)], label = ssp, color = col)
        ax.set_title(reg_names_area[area][reg])

    ctl.custom_legend(fig, colsim[1:], allsims[1:], ncol = 3)
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
            for mem in okmods:
                #seas20 = np.array(ctl.running_mean(np.concatenate([seasfreq[('hist', mem, reg)], seasfreq[(ssp, mem, reg)]]), 20))
                seas20 = np.array(ctl.lowpass_lanczos(np.concatenate([seasfreq[('hist', mem, reg)], seasfreq[(ssp, mem, reg)]]), 15))
                #ax.plot(yr, seas10)
                cosi.append(seas20)
            coso = np.mean(cosi, axis = 0)
            runfreq[(ssp, 'lanc20', reg)] = coso
            coserr = np.std(cosi, axis = 0)
            runfreq[(ssp, 'lanc20err', reg)] = coserr/np.sqrt(len(okmods)-1)
            ax.fill_between(yr, coso-coserr, coso+coserr, color = 'steelblue', alpha = 0.3)
            ax.plot(yr, coso, color = 'black', linewidth = 3)

            #seas20ref = np.array(ctl.running_mean(seasfreq[('hist', 'ref', reg)], 20))
            seas20ref = np.array(ctl.lowpass_lanczos(seasfreq[('hist', 'ref', reg)], 15))
            ax.plot(yr_ref, seas20ref, color = 'red', linewidth = 2, linestyle = '--')

            ax.set_title(reg_names_area[area][reg])
            ax.axvline(2015, color = 'lightslategray', linewidth = 0.2)

        fig.savefig(cart_out + 'long_lanc20_{}_{}_e_hist.pdf'.format(area, ssp))

    yr = np.arange(1965, 2100)
    for ssp in allssps:
        fig = plt.figure(figsize = (16,12))
        for reg in range(4):
            ax = fig.add_subplot(2, 2, reg+1)
            cosi = []
            for mem in okmods:
                seas20 = np.array(ctl.running_mean(np.conatenate([seasfreq[('hist', mem, reg)], seasfreq[(ssp, mem, reg)]]), 20))
                #ax.plot(yr, seas10)
                cosi.append(seas20)
            coso = np.mean(cosi, axis = 0)
            runfreq[(ssp, 'run20', reg)] = coso
            coserr = np.std(cosi, axis = 0)
            runfreq[(ssp, 'run20err', reg)] = coserr/np.sqrt(len(okmods)-1)
            ax.fill_between(yr, coso-coserr, coso+coserr, color = 'steelblue', alpha = 0.3)
            ax.plot(yr, coso, color = 'black', linewidth = 3)

            seas20ref = np.array(ctl.running_mean(seasfreq[('hist', 'ref', reg)], 20))
            ax.plot(yr_ref, seas20ref, color = 'red', linewidth = 2, linestyle = '--')

            ax.set_title(reg_names_area[area][reg])
            ax.axvline(2015, color = 'lightslategray', linewidth = 0.2)

        fig.savefig(cart_out + 'long_freq20_{}_{}_e_hist.pdf'.format(area, ssp))


    fig = plt.figure(figsize = (16,12))
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg+1)
        for col, ssp in zip(colsim[1:], allssps):
            coso = runfreq[(ssp, 'lanc20', reg)]
            coserr = runfreq[(ssp, 'lanc20err', reg)]
            ax.fill_between(yr, coso-coserr, coso+coserr, color = col, alpha = 0.05)
            ax.plot(yr, coso, label = ssp, color = col, linewidth = 2)
        ax.set_title(reg_names_area[area][reg])
        ax.axvline(2015, color = 'lightslategray', linewidth = 0.2)

    ctl.custom_legend(fig, colsim[1:], allsims[1:], ncol = 3)
    fig.savefig(cart_out + 'allssps_lanc20_{}.pdf'.format(area))

    fig = plt.figure(figsize = (16,12))
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg+1)
        for col, ssp in zip(colsim[1:], allssps):
            coso = runfreq[(ssp, 'run20', reg)]
            coserr = runfreq[(ssp, 'run20err', reg)]
            ax.fill_between(yr, coso-coserr, coso+coserr, color = col, alpha = 0.05)
            ax.plot(yr, coso, label = ssp, color = col, linewidth = 2)
        ax.set_title(reg_names_area[area][reg])
        ax.axvline(2015, color = 'lightslategray', linewidth = 0.2)

    ctl.custom_legend(fig, colsim[1:], allsims[1:], ncol = 3)
    fig.savefig(cart_out + 'allssps_freq20_{}.pdf'.format(area))

    pickle.dump([seasfreq, runfreq], open(cart_out + 'seasfreqs_{}_v4.p'.format(area), 'wb'))

    cart_out_nc = cart_out + 'clus_freq_index/'
    ctl.mkdir(cart_out_nc)
    cart_out_nc_ok = cart_out_nc + '{}_NDJFM/'.format(area)
    ctl.mkdir(cart_out_nc_ok)
    for sim in allsims:
        for mod in okmods:
            outfil = cart_out_nc_ok + 'clus_freq_NDJFM_{}_{}_{}.nc'.format(area, sim, mod)

            var = [seasfreq[(sim, mem, reg)] for reg in np.arange(4)]
            if sim == 'hist':
                yrs = np.arange(1964, 2014)
            else:
                yrs = np.arange(2015, 2100)

            dates = np.array([pd.to_datetime('{}1101'.format(year), format='%Y%m%d').to_pydatetime() for year in yrs])

            long_names = []
            for i, fre in enumerate(var):
                long_names.append('clus {} frequency'.format(i))

            ctl.save_iris_N_timeseries(outfil, var, dates = dates, time_units = 'day since 1850-01-01 00:00:00 UTC', time_cal = 'proleptic_gregorian', long_names = long_names)

    #################### la parte di v2
    freqs = dict() # tot50 e last20
    residtimes = dict() # mean e p90

    for mem in okmods:
        freqs[('hist', mem, 'tot50')] = ctl.calc_clus_freq(results_hist[mem]['labels'], numclus)

        for reg in range(numclus):
            alltimes = results_hist[mem]['resid_times'][reg]
            residtimes[('hist', mem, 'mean', reg)] = np.mean(alltimes)
            residtimes[('hist', mem, 'p90', reg)] = np.percentile(alltimes, 90)

        dat1 = pd.Timestamp('09-01-1995').to_pydatetime()
        dat2 = pd.Timestamp('04-01-2014').to_pydatetime()
        alllabs_20, dats = ctl.sel_time_range(results_hist[mem]['labels'], results_hist[mem]['dates'], (dat1, dat2))

        alltimes_20, _, _ = ctl.calc_regime_residtimes(labs, dats)

        for reg in range(numclus):
            residtimes[('hist', mem, 'mean_last20', reg)] = np.mean(alltimes_20[reg])
            residtimes[('hist', mem, 'p90_last20', reg)] = np.percentile(alltimes_20[reg], 90)

        freqs[('hist', mem, 'last20')] = ctl.calc_clus_freq(alllabs_20, numclus)

    for ssp in allssps:
        for mem in okmods:
            if mem not in results_ssp[ssp].keys(): continue
            ## Attach all members labels
            dat1 = pd.Timestamp('09-01-2050').to_pydatetime()
            dat2 = pd.Timestamp('04-01-2100').to_pydatetime()
            labs, dats = ctl.sel_time_range(results_ssp[ssp][mem]['labels'], results_ssp[ssp][mem]['dates'], (dat1, dat2))

            restim, _, _ = ctl.calc_regime_residtimes(labs, dats)

            freqs[(ssp, mem, 'tot50')] = ctl.calc_clus_freq(labs, numclus)

            for reg in range(numclus):
                residtimes[(ssp, mem, 'mean', reg)] = np.mean(restim[reg])
                residtimes[(ssp, mem, 'p90', reg)] = np.percentile(restim[reg], 90)


            dat1 = pd.Timestamp('09-01-2081').to_pydatetime()
            dat2 = pd.Timestamp('04-01-2100').to_pydatetime()
            labs, dats = ctl.sel_time_range(results_ssp[ssp][mem]['labels'], results_ssp[ssp][mem]['dates'], (dat1, dat2))

            restim, _, _ = ctl.calc_regime_residtimes(labs, dats)

            freqs[(ssp, mem, 'last20')] = ctl.calc_clus_freq(labs, numclus)

            for reg in range(numclus):
                residtimes[(ssp, mem, 'mean_last20', reg)] = np.mean(restim[reg])
                residtimes[(ssp, mem, 'p90_last20', reg)] = np.percentile(restim[reg], 90)

    ##### salvo le distrib per ogni ssp. Quelle assolute (con tutti i modelli) e quelle relative (solo con i modelli che ci sono anche in hist)
    freqs[('hist', 'all', 'tot50')] = np.stack([freqs[('hist', mod, 'tot50')] for mod in okmods])
    freqs[('hist', 'all', 'last20')] = np.stack([freqs[('hist', mod, 'last20')] for mod in okmods])
    for reg in range(numclus):
        for cos in ['mean', 'p90', 'mean_last20', 'p90_last20']:
            residtimes[('hist', 'all', cos, reg)] = np.array([residtimes[('hist', mod, cos, reg)] for mod in okmods])

    for ssp in allssps:
        freqs[(ssp, 'all', 'tot50')] = np.stack([freqs[(ssp, mod, 'tot50')] for mod in okmods])
        freqs[(ssp, 'all', 'last20')] = np.stack([freqs[(ssp, mod, 'last20')] for mod in okmods])

        # rel
        modoks = okmods
        freqs[(ssp, 'rel', 'tot50')] = np.stack([freqs[(ssp, mod, 'tot50')]-freqs[('hist', mod, 'tot50')] for mod in modoks])
        freqs[(ssp, 'rel', 'last20')] = np.stack([freqs[(ssp, mod, 'last20')]-freqs[('hist', mod, 'last20')] for mod in modoks])

        for reg in range(numclus):
            for cos in ['mean', 'p90', 'mean_last20', 'p90_last20']:
                residtimes[(ssp, 'all', cos, reg)] = np.array([residtimes[(ssp, mod, cos, reg)] for mod in okmods])

                residtimes[(ssp, 'rel', cos, reg)] = np.array([residtimes[(ssp, mod, cos, reg)]-residtimes[('hist', mod, cos, reg)] for mod in modoks])


    pickle.dump([freqs, residtimes], open(cart_out + 'allresults_dicts_{}_v3.p'.format(area), 'wb'))
    freqs, residtimes = pickle.load(open(cart_out + 'allresults_dicts_{}_v3.p'.format(area), 'rb'))

    #### Grafico con tutti gli ssp
    allsims = ['hist', 'ssp126', 'ssp245', 'ssp370', 'ssp585']
    colsim = ctl.color_set(5)

    reg_names = reg_names_area[area]
    for cos in ['last20', 'tot50']:
        fig = plt.figure(figsize = (16,12))
        axes = []
        for reg in range(4):
            ax = fig.add_subplot(2, 2, reg + 1)
            axes.append(ax)

            allpercs = dict()
            allpercs['mean'] = [np.mean(freqs[(ssp, 'rel', cos)][:, reg]) for ssp in allsims[1:]]

            for nu in [10, 25, 50, 75, 90]:
                allpercs['p{}'.format(nu)] = [np.percentile(freqs[(ssp, 'rel', cos)][:, reg], nu) for ssp in allsims[1:]]

            allpercs['ens_min'] = [np.min(freqs[(ssp, 'rel', cos)][:, reg]) for ssp in allsims[1:]]
            allpercs['ens_max'] = [np.max(freqs[(ssp, 'rel', cos)][:, reg]) for ssp in allsims[1:]]

            ctl.boxplot_on_ax(ax, allpercs, allsims[1:], colsim[1:], plot_mean = True, plot_ensmeans = False, plot_minmax = True)

            ax.axhline(0, color = 'gray', linewidth = 0.5)
            ax.set_xticks([])
            ax.set_title(reg_names[reg])

        ctl.adjust_ax_scale(axes)

        ctl.custom_legend(fig, colsim[1:], allsims[1:], ncol = 3)
        #fig.suptitle('30yr bsp of WR freq. in 2050-2100 wrt 1964-2014')
        fig.savefig(cart_out + 'WRfreqchange_allssp_{}_{}.pdf'.format(area, cos))


    for cos in ['mean', 'p90', 'mean_last20', 'p90_last20']:
        fig = plt.figure(figsize = (16,12))
        axes = []
        for reg in range(4):
            ax = fig.add_subplot(2, 2, reg + 1)
            axes.append(ax)

            allpercs = dict()
            allpercs['mean'] = [np.mean(residtimes[(ssp, 'rel', cos, reg)]) for ssp in allsims[1:]]

            for nu in [10, 25, 50, 75, 90]:
                allpercs['p{}'.format(nu)] = [np.percentile(residtimes[(ssp, 'rel', cos, reg)], nu) for ssp in allsims[1:]]

            allpercs['ens_min'] = [np.min(residtimes[(ssp, 'rel', cos, reg)]) for ssp in allsims[1:]]
            allpercs['ens_max'] = [np.max(residtimes[(ssp, 'rel', cos, reg)]) for ssp in allsims[1:]]

            ctl.boxplot_on_ax(ax, allpercs, allsims[1:], colsim[1:], plot_mean = True, plot_ensmeans = False, plot_minmax = True)

            ax.axhline(0, color = 'gray', linewidth = 0.5)
            ax.set_xticks([])
            ax.set_title(reg_names[reg])

        ctl.adjust_ax_scale(axes)

        ctl.custom_legend(fig, colsim[1:], allsims[1:], ncol = 3)
        #fig.suptitle('30yr bsp of WR freq. in 2050-2100 wrt 1964-2014')
        fig.savefig(cart_out + 'Restime_change_allssp_{}_{}.pdf'.format(area, cos))

    #### Absolute values
    dat1 = pd.Timestamp('09-01-1995').to_pydatetime()
    dat2 = pd.Timestamp('04-01-2014').to_pydatetime()
    labs, dats = ctl.sel_time_range(results_ref['labels'], results_ref['dates'], (dat1, dat2))
    results_ref['freq_clus_last20'] = ctl.calc_clus_freq(labs, numclus)

    figall = plt.figure(figsize = (28,12))
    axescos = dict()
    for cos in ['tot50', 'last20']:
        fig = plt.figure(figsize = (16,12))
        axes = []
        axescos[cos] = []
        for reg in range(4):
            ax = fig.add_subplot(2, 2, reg + 1)
            if cos == 'tot50':
                ax2 = figall.add_subplot(2, 4, reg + 1)
            else:
                ax2 = figall.add_subplot(2, 4, 4 + reg + 1)

            axes.append(ax)
            axescos[cos].append(ax2)

            allpercs = dict()
            allpercs['mean'] = [np.mean(freqs[(ssp, 'all', cos)][:, reg]) for ssp in allsims]

            for nu in [10, 25, 50, 75, 90]:
                allpercs['p{}'.format(nu)] = [np.percentile(freqs[(ssp, 'all', cos)][:, reg], nu) for ssp in allsims]

            allpercs['ens_min'] = [np.min(freqs[(ssp, 'all', cos)][:, reg]) for ssp in allsims]
            allpercs['ens_max'] = [np.max(freqs[(ssp, 'all', cos)][:, reg]) for ssp in allsims]

            ax.axhline(allpercs['p50'][0], color = 'gray', linewidth = 0.5)
            ctl.boxplot_on_ax(ax, allpercs, allsims, colsim, plot_mean = True, plot_ensmeans = False, plot_minmax = True)
            ax.set_xticks([])
            ax.set_title(reg_names[reg])

            ax2.axhline(allpercs['p50'][0], color = 'gray', linewidth = 0.5)
            ctl.boxplot_on_ax(ax2, allpercs, allsims, colsim, plot_mean = True, plot_ensmeans = False, plot_minmax = True)
            ax2.set_xticks([])
            if cos == 'tot50':
                ax2.set_title(reg_names[reg])
                ax2.scatter(0, results_ref['freq_clus'][reg], color = 'black', marker = '*', s = 5)
            else:
                ax2.scatter(0, results_ref['freq_clus_last20'][reg], color = 'black', marker = '*', s = 5)

        #ctl.adjust_ax_scale(axes)

        ctl.custom_legend(fig, colsim, allsims, ncol = 3)
        #fig.suptitle('30yr bsp of WR freq. in 2050-2100 wrt 1964-2014')
        fig.savefig(cart_out + 'WRfreq_allssp_{}_{}.pdf'.format(area, cos))

    ctl.custom_legend(figall, colsim, allsims, ncol = 3)
    figall.savefig(cart_out + 'WRfreq_allssp_{}_8box.pdf'.format(area, cos))

    fig = plt.figure(figsize = (16,12))
    axes = []

    allsims_dub = np.concatenate([[si,si] for si in allsims])
    colsims_dub = np.concatenate([[co,co] for co in colsim])
    positions = [0]
    positions.append(positions[-1]+0.1+0.7)
    for co in colsim[1:]:
        positions.append(positions[-1]+0.4+0.7)
        positions.append(positions[-1]+0.1+0.7)

    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg + 1)
        axes.append(ax)

        allpercs = dict()
        allpercs['mean'] = np.concatenate([[np.mean(freqs[(ssp, 'all', cos)][:, reg]) for cos in ['tot50', 'last20']] for ssp in allsims])

        for nu in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(nu)] = np.concatenate([[np.percentile(freqs[(ssp, 'all', cos)][:, reg], nu) for cos in ['tot50', 'last20']] for ssp in allsims])

        allpercs['ens_min'] = np.concatenate([[np.min(freqs[(ssp, 'all', cos)][:, reg]) for cos in ['tot50', 'last20']] for ssp in allsims])
        allpercs['ens_max'] = np.concatenate([[np.max(freqs[(ssp, 'all', cos)][:, reg]) for cos in ['tot50', 'last20']] for ssp in allsims])

        ax.axhline(allpercs['p50'][0], color = 'gray', linewidth = 0.5)
        ctl.boxplot_on_ax(ax, allpercs, allsims_dub, colsims_dub, plot_mean = True, plot_ensmeans = False, plot_minmax = False, positions = positions)

        ax.set_xticks([])
        ax.set_title(reg_names[reg])
        ax.scatter(positions[0], results_ref['freq_clus'][reg], color = 'black', marker = '*', s = 5)
        ax.scatter(positions[1], results_ref['freq_clus_last20'][reg], color = 'black', marker = '*', s = 2)

    #ctl.adjust_ax_scale(axes)

    ctl.custom_legend(fig, colsim, allsims, ncol = 3)
    #fig.suptitle('30yr bsp of WR freq. in 2050-2100 wrt 1964-2014')
    fig.savefig(cart_out + 'WRfreq_allssp_{}_bothrefs.pdf'.format(area, cos))

    for pio in ['mean', 'p90']:
        figall = plt.figure(figsize = (28,12))
        axescos = dict()
        for gigi in ['', '_last20']:
            cos = pio+gigi

            fig = plt.figure(figsize = (16,12))
            axes = []
            axescos[gigi] = []
            for reg in range(4):
                ax = fig.add_subplot(2, 2, reg + 1)

                if gigi == '':
                    ax2 = figall.add_subplot(2, 4, reg + 1)
                else:
                    ax2 = figall.add_subplot(2, 4, 4 + reg + 1)

                axes.append(ax)
                axescos[gigi].append(ax2)

                allpercs = dict()
                allpercs['mean'] = [np.mean(residtimes[(ssp, 'all', cos, reg)]) for ssp in allsims]

                for nu in [10, 25, 50, 75, 90]:
                    allpercs['p{}'.format(nu)] = [np.percentile(residtimes[(ssp, 'all', cos, reg)], nu) for ssp in allsims]

                allpercs['ens_min'] = [np.min(residtimes[(ssp, 'all', cos, reg)]) for ssp in allsims]
                allpercs['ens_max'] = [np.max(residtimes[(ssp, 'all', cos, reg)]) for ssp in allsims]

                ax.axhline(allpercs['p50'][0], color = 'gray', linewidth = 0.5)
                ctl.boxplot_on_ax(ax, allpercs, allsims, colsim, plot_mean = True, plot_ensmeans = False, plot_minmax = True)
                # ax.axhline(0, color = 'gray', linewidth = 0.5)
                ax.set_xticks([])
                ax.set_title(reg_names[reg])

                ax2.axhline(allpercs['p50'][0], color = 'gray', linewidth = 0.5)
                ctl.boxplot_on_ax(ax2, allpercs, allsims, colsim, plot_mean = True, plot_ensmeans = False, plot_minmax = True)
                # ax.axhline(0, color = 'gray', linewidth = 0.5)
                ax2.set_xticks([])
                if gigi == '':
                    ax2.set_title(reg_names[reg])

            #ctl.adjust_ax_scale(axes)

            ctl.custom_legend(fig, colsim, allsims, ncol = 3)
            #fig.suptitle('30yr bsp of WR freq. in 2050-2100 wrt 1964-2014')
            fig.savefig(cart_out + 'Restime_allssp_{}_{}.pdf'.format(area, cos))

        ctl.custom_legend(figall, colsim, allsims, ncol = 3)
        #fig.suptitle('30yr bsp of WR freq. in 2050-2100 wrt 1964-2014')
        figall.savefig(cart_out + 'Restime_allssp_{}_{}_8box.pdf'.format(area, pio))
