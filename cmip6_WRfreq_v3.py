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

cart_in = '/home/fabiano/Research/lavori/CMIP6/'
cart_in = '/home/fedefab/Scrivania/Research/Post-doc/lavori/CMIP6/'
cart_out_orig = cart_in + 'Results_v3/'
ctl.mkdir(cart_out_orig)

file_hist_refEOF = cart_in + 'cmip6_hist/out_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refEOF_dtr.p'
file_hist = cart_in + 'cmip6_hist/out_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refCLUS_dtr.p'
#file_hist_refEOF = cart_in + 'cmip6_hist/out_cmip6_hist_NDJFM_EAT_4clus_4pcs_1964-2014_refEOF.p'
gen_file_ssp = cart_in + 'cmip6_{}/out_cmip6_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr.p'

numclus = 4
reg_names_area = dict()
reg_names_area['EAT'] = ['NAO+', 'SBL', 'AR', 'NAO-']
reg_names_area['PNA'] = ['PT', 'PNA+', 'PNA-', 'AR']


allssps = 'ssp119 ssp126 ssp245 ssp370 ssp585'.split()

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
    runfreq = dict()

    fig = plt.figure(figsize = (16,12))
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg+1)
        cosi = []
        for mod in mod_hist:
            allmems = [cos for cos in results_hist.keys() if cos.split('_')[0] == mod]
            cosimem = []
            print(allmems)
            for mem in allmems:
                seasfr, yr = ctl.calc_seasonal_clus_freq(results_hist[mem]['labels'], results_hist[mem]['dates'], numclus)
                seas10 = np.array(ctl.running_mean(seasfr[reg, :], 10))
                ax.plot(yr, seas10)
                cosimem.append(seas10)

            if len(allmems) > 1:
                cosi.append(np.mean(cosimem, axis = 0))
            else:
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
            for mod in mod_ssp[ssp]:
                mem = [cos for cos in results_ssp[ssp].keys() if cos.split('_')[0] == mod][0]
                seasfr, yr = ctl.calc_seasonal_clus_freq(results_ssp[ssp][mem]['labels'], results_ssp[ssp][mem]['dates'], numclus)
                seas10 = np.array(ctl.running_mean(seasfr[reg, :], 10))
                ax.plot(yr, seas10, label = mod)
                cosi.append(seas10)
                print(mem, len(seas10))
            coso = np.mean(cosi, axis = 0)
            runfreq[(ssp, reg)] = coso
            ax.plot(yr, coso, color = 'black', linewidth = 3)
            print(yr[0], runfreq[('hist', reg)][-1])
            ax.scatter(yr[0], np.nanmean(runfreq[('hist', reg)]), s = 40, color = 'black')
            ax.scatter(yr[0], runfreq[('hist', reg)][~np.isnan(runfreq[('hist', reg)])][-1], s = 20, color = 'red')
            ax.set_xlim(2010, 2100)
            ax.set_title(reg_names_area[area][reg])

        fig.savefig(cart_out + 'models_freq10_{}_{}.pdf'.format(area, ssp))

        continue

        for mod in mod_ssp[ssp]:
            allmems = [cos for cos in results_ssp[ssp].keys() if cos.split('_')[0] == mod]
            ## Attach all members labels

            for reg in range(4):
                trendall = []
                for mem in allmems:
                    seasfr, yr = ctl.calc_seasonal_clus_freq(results_ssp[ssp][mem]['labels'], results_ssp[ssp][mem]['dates'], numclus)
                    m, c, err_m, err_c = ctl.linear_regre_witherr(yr, seasfr[reg, :])
                    trendall.append(m)

                if len(allmems) == 1:
                    trend_ssp[(ssp, mod, 'trend', 'seafreq', reg)] = m
                    trend_ssp[(ssp, mod, 'errtrend', 'seafreq', reg)] = err_m
                else:
                    trend_ssp[(ssp, mod, 'trend', 'seafreq', reg)] = np.mean(trendall)
                    trend_ssp[(ssp, mod, 'errtrend', 'seafreq', reg)] = np.std(trendall)

            for reg in range(4):
                trendall = []
                for mem in allmems:
                    seasfr, yr = ctl.calc_seasonal_clus_freq(results_ssp[ssp][mem]['labels'], results_ssp[ssp][mem]['dates'], numclus)
                    seas10 = ctl.running_mean(seasfr[reg, :], 10)
                    m, c, err_m, err_c = ctl.linear_regre_witherr(np.array(yr[~np.isnan(seas10)]), np.array(seas10[~np.isnan(seas10)]))
                    trendall.append(m)

                if len(allmems) == 1:
                    trend_ssp[(ssp, mod, 'trend', 'freq10', reg)] = m
                    trend_ssp[(ssp, mod, 'errtrend', 'freq10', reg)] = err_m
                else:
                    trend_ssp[(ssp, mod, 'trend', 'freq10', reg)] = np.mean(trendall)
                    trend_ssp[(ssp, mod, 'errtrend', 'freq10', reg)] = np.std(trendall)

            # devo fare ogni dieci anni e selezionare
            restimem = dict()
            mem = [cos for cos in results_ssp[ssp].keys() if cos.split('_')[0] == mod][0]
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
                residtime_ssp[(ssp, mod, 'trend', reg)] = m
                residtime_ssp[(ssp, mod, 'errtrend', reg)] = err_m

    continue

    for ssp in allssps:
        for reg in range(4):
            trend_ssp[(ssp, 'all', 'trend', 'seafreq', reg)] = np.array([trend_ssp[(ssp, mod, 'trend', 'seafreq', reg)] for mod in mod_ssp[ssp]])
            trend_ssp[(ssp, 'all', 'errtrend', 'seafreq', reg)] = np.array([trend_ssp[(ssp, mod, 'errtrend', 'seafreq', reg)] for mod in mod_ssp[ssp]])
            trend_ssp[(ssp, 'all', 'trend', 'freq10', reg)] = np.array([trend_ssp[(ssp, mod, 'trend', 'freq10', reg)] for mod in mod_ssp[ssp]])
            trend_ssp[(ssp, 'all', 'errtrend', 'freq10', reg)] = np.array([trend_ssp[(ssp, mod, 'errtrend', 'freq10', reg)] for mod in mod_ssp[ssp]])
            residtime_ssp[(ssp, 'all', 'trend', reg)] = np.array([residtime_ssp[(ssp, mod, 'trend', reg)] for mod in mod_ssp[ssp]])
            residtime_ssp[(ssp, 'all', 'errtrend', reg)] = np.array([residtime_ssp[(ssp, mod, 'errtrend', reg)] for mod in mod_ssp[ssp]])

    allsims = ['hist', 'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']
    colsim = ctl.color_set(6)

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
