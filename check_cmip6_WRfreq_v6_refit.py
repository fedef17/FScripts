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
cart_out_orig = '/home/fabiano/Research/lavori/CMIP6/Check_results_v6_refit/'
ctl.mkdir(cart_out_orig)

cart_in = '/data-hobbes/fabiano/WR_CMIP6/'
file_hist = cart_in + 'out_NEW_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refCLUS_dtr_light.p'
file_ssp = cart_in + 'out_NEW_cmip6_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr_light.p'
file_ssp_rebase = cart_in + 'out_NEW_cmip6_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr_histrebase.p'
fil_ece_ssp = cart_in + 'out_eceens_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr.p'
fil_ece_ssp_rebase = cart_in + 'out_eceens_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr_reb.p'

file_refit = cart_in + 'out_NEW_cmip6_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr_refit.p'
file_refit2 = cart_in + 'out_NEW_cmip6_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr_refit_rebasetot.p'

ssp = 'ssp585'
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

    # Erasing incomplete runs
    for ke in tuple(results_hist.keys()):
        if len(results_hist[ke]['labels']) < 7000:
            del results_hist[ke]

    results_ssp, _ = ctl.load_wrtool(file_ssp.format(ssp, area))
    results_ssp_rebase, _ = ctl.load_wrtool(file_ssp_rebase.format(ssp, area))

    results_refit, _ = ctl.load_wrtool(file_refit.format(ssp, area))
    results_refit2, _ = ctl.load_wrtool(file_refit2.format(ssp, area))

    # Adding the ece ensemble
    ece_ssp, _ = ctl.load_wrtool(fil_ece_ssp.format(ssp, area))
    ece_ssp_rebase, _ = ctl.load_wrtool(fil_ece_ssp_rebase.format(ssp, area))
    for mod in ece_ssp.keys():
        results_hist[mod] = results_hist['EC-Earth3_r1i1p1f1']

    del ece_ssp['EC-Earth3_r1i1p1f1']
    del ece_ssp_rebase['EC-Earth3_r1i1p1f1']
    results_ssp.update(ece_ssp)
    results_ssp_rebase.update(ece_ssp_rebase)

    # Erasing incomplete runs
    for ke in tuple(results_ssp.keys()):
        if len(results_ssp[ke]['labels']) < 12000:
            del results_ssp[ke]
        elif len(results_ssp[ke]['labels']) > 13000:
            # there is some duplicated year
            for cosone in [results_ssp, results_ssp_rebase]:
                labs, dats = ctl.seasonal_set(cosone[ke]['labels'], cosone[ke]['dates'], None)
                pcs, dats = ctl.seasonal_set(cosone[ke]['pcs'], cosone[ke]['dates'], None)
                yeas = np.array([da[0].year for da in dats])
                labs_ok = []
                dats_ok = []
                pcs_ok = []
                for ye in np.arange(2015, 2100):
                    okse = np.where(yeas == ye)[0][0]
                    labs_ok.append(labs[okse])
                    dats_ok.append(dats[okse])
                    pcs_ok.append(pcs[okse])
                cosone[ke]['labels'] = np.concatenate(labs_ok)
                cosone[ke]['dates'] = np.concatenate(dats_ok)
                cosone[ke]['pcs'] = np.concatenate(pcs_ok)

            # same for results_refit
            for cosone in [results_refit, results_refit2]:
                labs, dats = ctl.seasonal_set(cosone[ke]['labels'], cosone[ke]['dates'], None)
                pcs, dats = ctl.seasonal_set(cosone[ke]['pcs'], cosone[ke]['dates'], None)
                yeas = np.array([da[0].year for da in dats])
                labs_ok = []
                dats_ok = []
                pcs_ok = []
                for ye in np.arange(1964, 2100):
                    okse = np.where(yeas == ye)[0][0]
                    labs_ok.append(labs[okse])
                    dats_ok.append(dats[okse])
                    pcs_ok.append(pcs[okse])
                cosone[ke]['labels'] = np.concatenate(labs_ok)
                cosone[ke]['dates'] = np.concatenate(dats_ok)
                cosone[ke]['pcs'] = np.concatenate(pcs_ok)

    mod_hist = np.unique([cos.split('_')[0] for cos in results_hist.keys()])
    mod_ssp = np.unique([cos.split('_')[0] for cos in results_ssp.keys()])
    print(ssp, mod_ssp)

    okmods_hist = list(results_hist.keys())
    okmods_ssp = list(results_ssp.keys())
    okmods = [mod for mod in okmods_ssp if mod in okmods_hist]

    # appiccico hist e ssp
    for cosone in [results_ssp, results_ssp_rebase]:
        for mem in okmods:
            labs = np.concatenate([results_hist[mem]['labels'], cosone[mem]['labels']])
            dats = np.concatenate([results_hist[mem]['dates'], cosone[mem]['dates']])
            pcs = np.concatenate([results_hist[mem]['pcs'], cosone[mem]['pcs']], axis = 0)
            cosone[mem]['labels'] = labs
            cosone[mem]['dates'] = dats
            cosone[mem]['pcs'] = pcs

    resdict = dict()
    resdict['noreb'] = results_ssp
    resdict['histrebase'] = results_ssp_rebase
    resdict['refit'] = results_refit
    resdict['refit_rbtot'] = results_refit2
    alltips = tuple(resdict.keys())
    colorz = ctl.color_set(len(alltips))

    runfreq = dict()
    seasfreq = dict()

    figs = []
    for mem in okmods:
        fig = plt.figure(figsize = (16,12))
        plt.title(mem)
        for reg in range(4):
            ax = fig.add_subplot(2, 2, reg+1)
            cosi = []
            for tip, col in zip(alltips, colorz):
                seasfr, yr = ctl.calc_seasonal_clus_freq(resdict[tip][mem]['labels'], resdict[tip][mem]['dates'], numclus)
                seasfreq[(tip, mem, reg)] = seasfr[reg, :]
                seas1 = np.array(ctl.running_mean(seasfr[reg, :], 20))
                runfreq[(tip, mem, reg)] = seas1
                ax.plot(yr, seas1, label = tip, color = col)
            ax.set_title(reg_names_area[area][reg])
            ax.legend()
        figs.append(fig)

    ctl.plot_pdfpages(cart_out + 'check_models_refit_{}.pdf'.format(area), figs)

    okmods = [mod for mod in okmods if 'EC-Earth3' not in mod]
    okmods.append('EC-Earth3_r1i1p1f1')

    fig = plt.figure(figsize = (16,12))
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg+1)
        for tip, col in zip(alltips, colorz):
            cosi = [runfreq[(tip, mem, reg)] for mem in okmods]
            coso = np.mean(cosi, axis = 0)
            runfreq[(tip, reg)] = coso
            coserr = np.std(cosi, axis = 0)
            ax.fill_between(yr, coso-coserr, coso+coserr, color = col, alpha = 0.3)
            ax.plot(yr, coso, label = tip, color = col)

        ax.set_title(reg_names_area[area][reg])
        ax.legend()

    fig.savefig(cart_out + 'check_all_refit_{}.pdf'.format(area))

    trend_ssp = dict()

    for tip in alltips:
        for reg in range(4):
            for mem in okmods:
                seasfr = seasfreq[(tip, mem, reg)]
                seas10 = np.array(ctl.running_mean(seasfr, 10))

                okcose = ~(np.isnan(seas10)) & (yr >= 2015)
                m, c, err_m, err_c = ctl.linear_regre_witherr(np.array(yr[okcose]), np.array(seas10[okcose]))
                trend_ssp[(mem, tip, reg)] = m
                trend_ssp[(mem, tip, reg, 'err')] = err_m

            trend_ssp[(tip, reg)] = np.array([trend_ssp[(mem, tip, reg)] for mem in okmods])
            trend_ssp[(tip, reg, 'err')] = np.array([trend_ssp[(mem, tip, reg, 'err')] for mem in okmods])


    freqs = dict()
    for tip in alltips:
        for mem in okmods:
            dat1 = pd.Timestamp('09-01-1964').to_pydatetime()
            dat2 = pd.Timestamp('04-01-2014').to_pydatetime()
            labs, dats = ctl.sel_time_range(resdict[tip][mem]['labels'], resdict[tip][mem]['dates'], (dat1, dat2))
            freqs[('hist', mem, tip)] = ctl.calc_clus_freq(labs, numclus)

            dat1 = pd.Timestamp('09-01-2050').to_pydatetime()
            dat2 = pd.Timestamp('04-01-2100').to_pydatetime()
            labs, dats = ctl.sel_time_range(resdict[tip][mem]['labels'], resdict[tip][mem]['dates'], (dat1, dat2))

            # restim, _, _ = ctl.calc_regime_residtimes(labs, dats)
            freqs[(ssp, mem, tip)] = ctl.calc_clus_freq(labs, numclus)

        freqs[('hist', tip)] = np.stack([freqs[('hist', mem, tip)] for mem in okmods])
        freqs[(ssp, tip)] = np.stack([freqs[(ssp, mem, tip)] for mem in okmods])

    pickle.dump([seasfreq, runfreq, freqs, trend_ssp], open(cart_out + 'check_refit_{}.p'.format(area), 'wb'))

    reg_names = reg_names_area[area]

    figall = plt.figure(figsize = (16,12))
    axes = []
    for reg in range(4):
        ax = figall.add_subplot(2, 2, reg + 1)
        axes.append(ax)

        histmean = dict()
        for tip in alltips:
            histmean[tip] = np.mean(freqs[('hist', tip)][:, reg])

        allpercs = dict()
        for nu in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(nu)] = [np.percentile(freqs[(ssp, tip)][:, reg], nu) - histmean[tip] for tip in alltips]

        ctl.boxplot_on_ax(ax, allpercs, alltips, colorz, plot_mean = False, plot_ensmeans = False, plot_minmax = False)
        ax.axhline(0, color = 'gray', linewidth = 0.5)
        ax.set_xticks([])
        ax.set_title(reg_names[reg])

        if reg == 0 or reg == 2: ax.set_ylabel('Regime frequency anomaly')

    ctl.adjust_ax_scale(axes)

    ctl.custom_legend(figall, colorz, alltips, ncol = 2)
    figall.savefig(cart_out + 'check_WRfreq_alltips_{}.pdf'.format(area))


    figall = plt.figure(figsize = (16,12))
    axes = []
    for reg in range(4):
        ax = figall.add_subplot(2, 2, reg + 1)

        allpercs = dict()
        for nu in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(nu)] = [np.percentile(trend_ssp[(tip, reg)], nu) for tip in alltips]

        ctl.boxplot_on_ax(ax, allpercs, alltips, colorz, plot_mean = False, plot_ensmeans = False, plot_minmax = False)

        ax.axhline(0, color = 'gray', linewidth = 0.5)
        ax.set_xticks([])
        if reg == 0 or reg == 2: ax.set_ylabel('Trend in regime frequency (1/yr)')
        axes.append(ax)

    ctl.adjust_ax_scale(axes)

    ctl.custom_legend(figall, colorz, alltips, ncol = 3)
    figall.savefig(cart_out + 'check_trends_alltips_{}.pdf'.format(area))


    figall = plt.figure(figsize = (16,12))
    axes = []
    for reg in range(4):
        ax = figall.add_subplot(2, 2, reg + 1)
        axes.append(ax)

        histmean = dict()
        for tip in alltips:
            histmean[tip] = np.mean(freqs[('hist', tip)][:, reg])

        data = [freqs[(ssp, tip)][:, reg]-histmean[tip] for tip in alltips]

        parts = ax.violinplot(data, positions = None, widths=0.4, showmeans=True, showextrema=True, showmedians=True, quantiles=[0.25, 0.75], bw_method=0.5)
        for pc, col in zip(parts['bodies'], colorz):
            pc.set_facecolor(col)
            pc.set_edgecolor(col)
            pc.set_alpha(0.5)

        ax.axhline(0, color = 'gray', linewidth = 0.5)
        ax.set_xticks([])
        ax.set_title(reg_names[reg])

        if reg == 0 or reg == 2: ax.set_ylabel('Regime frequency anomaly')

    ctl.adjust_ax_scale(axes)

    ctl.custom_legend(figall, colorz, alltips, ncol = 2)
    figall.savefig(cart_out + 'check_WRfreq_alltips_violin_{}.pdf'.format(area))
