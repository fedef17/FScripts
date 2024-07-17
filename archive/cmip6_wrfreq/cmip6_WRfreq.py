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

gen_file_hist = cart_in + 'cmip6_hist/out_cmip6_hist_{}_{}_4clus_4pcs_1964-2014_{}_dtr.p'
#file_hist_refEOF = cart_in + 'cmip6_hist/out_cmip6_hist_NDJFM_EAT_4clus_4pcs_1964-2014_refEOF.p'
gen_file_ssp585 = cart_in + 'cmip6_{}/out_cmip6_{}_{}_{}_4clus_4pcs_2050-2100_{}_dtr.p'

numclus = 4
reg_names_area = dict()
reg_names_area['EAT'] = ['NAO+', 'SBL', 'AR', 'NAO-']
reg_names_area['PNA'] = ['PT', 'PNA+', 'PNA-', 'AR']

#kes = [('ssp585', 'EAT', 'NDJFM', 'refEOF'), ('ssp585', 'EAT', 'NDJFM', 'refCLUS'), ('ssp585', 'PNA', 'NDJFM', 'refCLUS'), ('ssp126', 'EAT', 'NDJFM', 'refCLUS'), ('ssp245', 'EAT', 'NDJFM', 'refCLUS'), ('ssp370', 'EAT', 'NDJFM', 'refCLUS'), ('ssp119', 'EAT', 'NDJFM', 'refCLUS'), ('ssp126', 'PNA', 'NDJFM', 'refCLUS'), ('ssp245', 'PNA', 'NDJFM', 'refCLUS'), ('ssp370', 'PNA', 'NDJFM', 'refCLUS'), ('ssp119', 'PNA', 'NDJFM', 'refCLUS')]
allssps = 'ssp119 ssp126 ssp245 ssp370 ssp585'.split()
kes = []
for seas in ['NDJFM', 'DJFM']:
    for area in ['EAT', 'PNA']:
        for ssp in allssps:
            kes.append((ssp, area, seas, 'refCLUS'))


allpercs_ssp = dict()
allfreqs = dict()

for ke in kes:
    ssp, area, season, tip = ke
    file_hist = gen_file_hist.format(season, area, tip)
    file_ssp585 = gen_file_ssp585.format(ssp, ssp, season, area, tip)

    cart_out = cart_in + 'FutWR_{}_{}_{}_{}/'.format(*ke)
    ctl.mkdir(cart_out)

    results_hist, results_ref = pickle.load(open(file_hist, 'rb'))
    results_ssp585, results_ref = pickle.load(open(file_ssp585, 'rb'))
    if tip == 'refCLUS':
        res_hist_refEOF, _ = pickle.load(open(gen_file_hist.format(season, area, 'refEOF'), 'rb'))

    ### plot 0: WR freq difference tra ssp e hist, per ogni WR
    reg_names = reg_names_area[area]
    modoks = [mod for mod in results_ssp585.keys() if mod in results_hist.keys()]
    coloks = ctl.color_set(len(modoks))

    fig = plt.figure(figsize = (16,12))
    axes = []
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg + 1)
        axes.append(ax)

        i = 0
        allfreqs[(ssp, area, season, tip, reg)] = []
        for mod, col in zip(modoks, coloks):
            frdiff = results_ssp585[mod]['freq_clus'][reg] - results_hist[mod]['freq_clus'][reg]
            ax.scatter(i, frdiff, marker = 'D', color = col, s = 25)
            allfreqs[(ssp, area, season, tip, reg)].append(results_ssp585[mod]['freq_clus'][reg]- results_hist[mod]['freq_clus'][reg])
            # allfreqs[('hist', tip, reg)].append(results_hist[mod]['freq_clus'][reg])
            i += 1

        ax.axhline(0, color = 'gray', linewidth = 0.5)
        ax.set_xticks([])
        ax.set_title(reg_names[reg])

    ctl.adjust_ax_scale(axes)
    ctl.custom_legend(fig, coloks, modoks, ncol = 4)
    fig.suptitle('Change of WR freq. in 2050-2100 wrt 1964-2014')

    fig.savefig(cart_out + 'WR_freq_change_dtr_{}_{}_{}_{}.pdf'.format(*ke))

    # Same but only for long-lasting (> 5 days) regimes
    for mod in modoks:
        results_ssp585[mod]['av_res'] = np.array([np.mean(results_ssp585[mod]['resid_times'][reg]) for reg in range(4)])
        results_hist[mod]['av_res'] = np.array([np.mean(results_hist[mod]['resid_times'][reg]) for reg in range(4)])

        labs = ctl.regime_filter_long(results_ssp585[mod]['labels'], results_ssp585[mod]['dates'], days_thres = 5)
        results_ssp585[mod]['freq_clus_long5'] = ctl.calc_clus_freq(labs[labs >= 0], numclus)
        labs = ctl.regime_filter_long(results_hist[mod]['labels'], results_hist[mod]['dates'], days_thres = 5)
        results_hist[mod]['freq_clus_long5'] = ctl.calc_clus_freq(labs[labs >= 0], numclus)

    fig = plt.figure(figsize = (16,12))
    axes = []
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg + 1)
        axes.append(ax)
        allfreqs[(ssp, area, season, 'refCLUS_long5', reg)] = []

        i = 0
        for mod, col in zip(modoks, coloks):
            frdiff = results_ssp585[mod]['freq_clus_long5'][reg] - results_hist[mod]['freq_clus_long5'][reg]
            allfreqs[(ssp, area, season, 'refCLUS_long5', reg)].append(frdiff)
            ax.scatter(i, frdiff, marker = 'D', color = col, s = 25)
            i += 1

        ax.axhline(0, color = 'gray', linewidth = 0.5)
        ax.set_xticks([])
        ax.set_title(reg_names[reg])

    ctl.adjust_ax_scale(axes)
    ctl.custom_legend(fig, coloks, modoks, ncol = 4)
    fig.suptitle('Change of WR freq. in 2050-2100 wrt 1964-2014 (> 5 days persistence)')

    fig.savefig(cart_out + 'WR_freq_change_dtr_long5_{}_{}_{}_{}.pdf'.format(*ke))


    ### WR freq change for 2081-2100 vs 1995-2014
    for mod in modoks:
        dat1 = pd.Timestamp('09-01-2081').to_pydatetime()
        dat2 = pd.Timestamp('04-01-2100').to_pydatetime()
        labs, dats = ctl.sel_time_range(results_ssp585[mod]['labels'], results_ssp585[mod]['dates'], (dat1, dat2))
        results_ssp585[mod]['freq_clus_last20'] = ctl.calc_clus_freq(labs, numclus)

        dat1 = pd.Timestamp('09-01-1995').to_pydatetime()
        dat2 = pd.Timestamp('04-01-2014').to_pydatetime()
        labs, dats = ctl.sel_time_range(results_hist[mod]['labels'], results_hist[mod]['dates'], (dat1, dat2))
        results_hist[mod]['freq_clus_last20'] = ctl.calc_clus_freq(labs, numclus)

    fig = plt.figure(figsize = (16,12))
    axes = []
    for reg in range(4):
        allfreqs[(ssp, area, season, 'refCLUS_last20', reg)] = []
        ax = fig.add_subplot(2, 2, reg + 1)
        axes.append(ax)

        i = 0
        for mod, col in zip(modoks, coloks):
            frdiff = results_ssp585[mod]['freq_clus_last20'][reg] - results_hist[mod]['freq_clus_last20'][reg]
            allfreqs[(ssp, area, season, 'refCLUS_last20', reg)].append(frdiff)
            ax.scatter(i, frdiff, marker = 'D', color = col, s = 25)
            i += 1

        ax.axhline(0, color = 'gray', linewidth = 0.5)
        ax.set_xticks([])
        ax.set_title(reg_names[reg])

    ctl.adjust_ax_scale(axes)
    ctl.custom_legend(fig, coloks, modoks, ncol = 4)
    fig.suptitle('Change of WR freq. in 2081-2100 wrt 1995-2014')

    fig.savefig(cart_out + 'WR_freq_change_dtr_last20_{}_{}_{}_{}.pdf'.format(*ke))


    fig = plt.figure(figsize = (16,12))
    axes = []
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg + 1)
        axes.append(ax)
        allfreqs[(ssp, area, season, 'refCLUS_avres', reg)] = []

        i = 0
        for mod, col in zip(modoks, coloks):
            frdiff = results_ssp585[mod]['av_res'][reg] - results_hist[mod]['av_res'][reg]
            allfreqs[(ssp, area, season, 'refCLUS_avres', reg)].append(frdiff)
            ax.scatter(i, frdiff, marker = 'D', color = col, s = 25)
            i += 1

        ax.axhline(0, color = 'gray', linewidth = 0.5)
        ax.set_xticks([])
        ax.set_title(reg_names[reg])

    ctl.adjust_ax_scale(axes)
    ctl.custom_legend(fig, coloks, modoks, ncol = 4)
    fig.suptitle('Change of avg. regime residence time in 2050-2100 wrt 1964-2014')

    fig.savefig(cart_out + 'WR_av_restime_dtr_{}_{}_{}_{}.pdf'.format(*ke))


    # regime frequency change for subseasons: ND, DJF, FM
    allseas = ['ND', 'DJF', 'FM']

    for sea in allseas:
        for mod in modoks:
            lab_seas, dates_seas = ctl.sel_season(results_ssp585[mod]['labels'], results_ssp585[mod]['dates'], sea)
            results_ssp585[mod]['freq_clus_'+sea] = ctl.calc_clus_freq(lab_seas, numclus)
            lab_seas, dates_seas = ctl.sel_season(results_hist[mod]['labels'], results_hist[mod]['dates'], sea)
            results_hist[mod]['freq_clus_'+sea] = ctl.calc_clus_freq(lab_seas, numclus)

        fig = plt.figure(figsize = (16,12))
        axes = []
        for reg in range(4):
            ax = fig.add_subplot(2, 2, reg + 1)
            axes.append(ax)

            i = 0
            for mod, col in zip(modoks, coloks):
                frdiff = results_ssp585[mod]['freq_clus_'+sea][reg] - results_hist[mod]['freq_clus_'+sea][reg]
                ax.scatter(i, frdiff, marker = 'D', color = col, s = 25)
                i += 1

            ax.axhline(0, color = 'gray', linewidth = 0.5)
            ax.set_xticks([])
            ax.set_title(reg_names[reg])

        ctl.adjust_ax_scale(axes)
        ctl.custom_legend(fig, coloks, modoks, ncol = 4)
        fig.suptitle('Change of WR freq. in 2050-2100 wrt 1964-2014 - {}'.format(sea))

        fig.savefig(cart_out + 'WR_freq_change_dtr_{}_{}_{}_{}_subsea{}.pdf'.format(*ke, sea))


    #### Prima posso guardare la cosa con le seasonal frequencies. O monthly? No seasonal.
    ## Devo. calcolare le seasonal freq. calcolare i percentili per hist e per ssp. Fare differenze di ogni percentile. Plottare. Ok vado.
    for mod in modoks:
        seasfr, yr = ctl.calc_seasonal_clus_freq(results_hist[mod]['labels'], results_hist[mod]['dates'], numclus)
        results_hist[mod]['freq_clus_seasonal'] = seasfr
        seasfr, yr = ctl.calc_seasonal_clus_freq(results_ssp585[mod]['labels'], results_ssp585[mod]['dates'], numclus)
        results_ssp585[mod]['freq_clus_seasonal'] = seasfr

    allpercs_regs = []
    for reg in range(numclus):
        allpercs = dict()
        allpercs['mean'] = []
        for pp in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(pp)] = []

        for mod in modoks:
            hist_distr = results_hist[mod]['freq_clus_seasonal'][reg, ]
            ssp_distr = results_ssp585[mod]['freq_clus_seasonal'][reg, ]
            # allmedian = np.median(np.concatenate([hist_distr, ssp_distr]))
            histmedian = np.median(hist_distr)

            allpercs['mean'].append(np.mean(hist_distr)-histmedian)
            allpercs['mean'].append(np.mean(ssp_distr)-histmedian)

            for pp in [10, 25, 50, 75, 90]:
                allpercs['p{}'.format(pp)].append(np.percentile(hist_distr, pp)-histmedian)
                allpercs['p{}'.format(pp)].append(np.percentile(ssp_distr, pp)-histmedian)

        allpercs_regs.append(allpercs)

    modoks_dub = []
    coloks_dub = []
    positions = []
    versions = []
    pos = 0
    for mod, col in zip(modoks, coloks):
        modoks_dub.append(mod)
        modoks_dub.append(mod)
        coloks_dub.append(col)
        coloks_dub.append(col)
        positions.append(pos)
        positions.append(pos+0.7)
        pos+=1.7
        versions.append('hist')
        versions.append('ssp585')

    edge_colors = ['steelblue', 'indianred']*len(modoks)

    positions.append(pos+0.5)
    positions.append(pos+1.5)

    fig = plt.figure(figsize = (16,12))
    axes = []
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg + 1)
        axes.append(ax)

        ctl.boxplot_on_ax(ax, allpercs_regs[reg], modoks_dub, coloks_dub, edge_colors = edge_colors, positions = positions, versions = versions, ens_colors = ['steelblue', 'indianred'], ens_names = ['hist', 'ssp585'], plot_mean = False)

        ax.axhline(0, color = 'gray', linewidth = 0.5)
        ax.set_xticks([])
        ax.set_title(reg_names[reg])

    ctl.adjust_ax_scale(axes)
    ctl.custom_legend(fig, coloks + ['steelblue', 'indianred'], modoks + ['hist', 'ssp585'], ncol = 4)
    fig.suptitle('distribution of seasonal WR freq. in 2050-2100 wrt 1964-2014')

    fig.savefig(cart_out + 'Dist_seasonal_WRfreq_dtr_{}_{}_{}_{}.pdf'.format(*ke))


    ### 30-yr bootstraps of the WR freq.
    for mod in modoks:
        clusfreqs = ctl.bootstrap(results_hist[mod]['labels'], results_hist[mod]['dates'], None, apply_func = ctl.calc_clus_freq, func_args = [numclus], n_choice = 30, n_bootstrap = 500)
        results_hist[mod]['freq_clus_30yrboot'] = np.stack(clusfreqs)

        clusfreqs = ctl.bootstrap(results_ssp585[mod]['labels'], results_ssp585[mod]['dates'], None, apply_func = ctl.calc_clus_freq, func_args = [numclus], n_choice = 30, n_bootstrap = 500)
        results_ssp585[mod]['freq_clus_30yrboot'] = np.stack(clusfreqs)

    allpercs_regs = []
    for reg in range(numclus):
        allpercs = dict()
        allpercs['mean'] = []
        for pp in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(pp)] = []

        for mod in modoks:
            hist_distr = results_hist[mod]['freq_clus_30yrboot'][:,reg]
            ssp_distr = results_ssp585[mod]['freq_clus_30yrboot'][:,reg]
            #allmedian = np.median(np.concatenate([hist_distr, ssp_distr]))
            histmedian = np.median(hist_distr)

            allpercs['mean'].append(np.mean(hist_distr)-histmedian)
            allpercs['mean'].append(np.mean(ssp_distr)-histmedian)

            for pp in [10, 25, 50, 75, 90]:
                allpercs['p{}'.format(pp)].append(np.percentile(hist_distr, pp)-histmedian)
                allpercs['p{}'.format(pp)].append(np.percentile(ssp_distr, pp)-histmedian)

        allpercs_regs.append(allpercs)

    for reg in range(numclus):
        allpercs_ssp[(ssp, area, season, tip, reg)] = dict()
        allpercs_ssp[('hist', area, season, tip, reg)] = dict()
        for kiav in allpercs_regs[reg].keys():
            allpercs_ssp[(ssp, area, season, tip, reg)][kiav] = np.mean(allpercs_regs[reg][kiav][1::2])
            allpercs_ssp[('hist', area, season, tip, reg)][kiav] = np.mean(allpercs_regs[reg][kiav][0::2])

        allpercs_ssp[(ssp, area, season, tip, reg)]['ens_min'] = np.min(allpercs_regs[reg]['p50'][1::2])
        allpercs_ssp[(ssp, area, season, tip, reg)]['ens_max'] = np.max(allpercs_regs[reg]['p50'][1::2])
        allpercs_ssp[('hist', area, season, tip, reg)]['ens_min'] = np.nan
        allpercs_ssp[('hist', area, season, tip, reg)]['ens_max'] = np.nan

    fig = plt.figure(figsize = (16,12))
    axes = []
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg + 1)
        axes.append(ax)

        ctl.boxplot_on_ax(ax, allpercs_regs[reg], modoks_dub, coloks_dub, edge_colors = edge_colors, positions = positions, versions = versions, ens_colors = ['steelblue', 'indianred'], ens_names = ['hist', 'ssp585'], plot_mean = False)

        ax.axhline(0, color = 'gray', linewidth = 0.5)
        ax.set_xticks([])
        ax.set_title(reg_names[reg])

    ctl.adjust_ax_scale(axes)
    ctl.custom_legend(fig, coloks + ['steelblue', 'indianred'], modoks + ['hist', 'ssp585'], ncol = 4)
    fig.suptitle('30yr bootstrap of WR freq. in 2050-2100 wrt 1964-2014')

    fig.savefig(cart_out + 'Dist_30yr_WRfreq_dtr_{}_{}_{}_{}.pdf'.format(*ke))


    ####### RECALCULATING THE FUTURE LABELS BASED ON THE HIST CENTROIDS
    ### 30-yr bootstraps of the WR freq.
    if tip == 'refCLUS':
        for reg in range(numclus):
            allfreqs[(ssp, area, season, 'refMODCLUS', reg)] = []

        for mod in modoks:
            clusfreqs = ctl.bootstrap(res_hist_refEOF[mod]['labels'], res_hist_refEOF[mod]['dates'], None, apply_func = ctl.calc_clus_freq, func_args = [numclus], n_choice = 30, n_bootstrap = 500)
            res_hist_refEOF[mod]['freq_clus_30yrboot'] = np.stack(clusfreqs)

            centroids = res_hist_refEOF[mod]['centroids']

            pcs = results_ssp585[mod]['pcs']
            nulabels = ctl.assign_to_closest_cluster(pcs, centroids)

            nufreqs = ctl.calc_clus_freq(nulabels, numclus)

            for reg in range(numclus):
                allfreqs[(ssp, area, season, 'refMODCLUS', reg)].append(nufreqs[reg] - res_hist_refEOF[mod]['freq_clus'][reg])

            clusfreqs = ctl.bootstrap(nulabels, results_ssp585[mod]['dates'], None, apply_func = ctl.calc_clus_freq, func_args = [numclus], n_choice = 30, n_bootstrap = 500)
            results_ssp585[mod]['freq_clus_30yrboot_refMODCLUS'] = np.stack(clusfreqs)

        allpercs_regs = []
        for reg in range(numclus):
            allpercs = dict()
            allpercs['mean'] = []
            for pp in [10, 25, 50, 75, 90]:
                allpercs['p{}'.format(pp)] = []

            for mod in modoks:
                hist_distr = res_hist_refEOF[mod]['freq_clus_30yrboot'][:,reg]
                ssp_distr = results_ssp585[mod]['freq_clus_30yrboot_refMODCLUS'][:,reg]
                # allmedian = np.median(np.concatenate([hist_distr, ssp_distr]))
                histmedian = np.median(hist_distr)

                allpercs['mean'].append(np.mean(hist_distr)-histmedian)
                allpercs['mean'].append(np.mean(ssp_distr)-histmedian)

                for pp in [10, 25, 50, 75, 90]:
                    allpercs['p{}'.format(pp)].append(np.percentile(hist_distr, pp)-histmedian)
                    allpercs['p{}'.format(pp)].append(np.percentile(ssp_distr, pp)-histmedian)

            allpercs_regs.append(allpercs)

        for reg in range(numclus):
            allpercs_ssp[(ssp, area, season, 'refMODCLUS', reg)] = dict()
            allpercs_ssp[('hist', area, season, 'refMODCLUS', reg)] = dict()
            for kiav in allpercs_regs[reg].keys():
                allpercs_ssp[(ssp, area, season, 'refMODCLUS', reg)][kiav] = np.mean(allpercs_regs[reg][kiav][1::2])
                allpercs_ssp[('hist', area, season, 'refMODCLUS', reg)][kiav] = np.mean(allpercs_regs[reg][kiav][0::2])

            allpercs_ssp[(ssp, area, season, 'refMODCLUS', reg)]['ens_min'] = np.min(allpercs_regs[reg]['p50'][1::2])
            allpercs_ssp[(ssp, area, season, 'refMODCLUS', reg)]['ens_max'] = np.max(allpercs_regs[reg]['p50'][1::2])
            allpercs_ssp[('hist', area, season, 'refMODCLUS', reg)]['ens_max'] = np.nan
            allpercs_ssp[('hist', area, season, 'refMODCLUS', reg)]['ens_min'] = np.nan


        fig = plt.figure(figsize = (16,12))
        axes = []
        for reg in range(4):
            ax = fig.add_subplot(2, 2, reg + 1)
            axes.append(ax)

            ctl.boxplot_on_ax(ax, allpercs_regs[reg], modoks_dub, coloks_dub, edge_colors = edge_colors, positions = positions, versions = versions, ens_colors = ['steelblue', 'indianred'], ens_names = ['hist', 'ssp585'], plot_mean = False)

            ax.axhline(0, color = 'gray', linewidth = 0.5)
            ax.set_xticks([])
            ax.set_title(reg_names[reg])

        ctl.adjust_ax_scale(axes)
        ctl.custom_legend(fig, coloks + ['steelblue', 'indianred'], modoks + ['hist', 'ssp585'], ncol = 4)
        fig.suptitle('30yr bootstrap of WR freq. in 2050-2100 wrt 1964-2014')

        fig.savefig(cart_out + 'Dist_30yr_WRfreq_dtr_{}_{}_{}_refMODCLUS.pdf'.format(*ke[:-1]))


    #### il maxgrads
    if 'maxgrad' in results_hist[mod]:
        allpercs_regs = []
        for reg in range(numclus):
            allpercs = dict()
            allpercs['mean'] = []
            for pp in [10, 25, 50, 75, 90]:
                allpercs['p{}'.format(pp)] = []

            for mod in modoks:
                okcosi = results_hist[mod]['labels'] == reg
                hist_distr = results_hist[mod]['maxgrad'][okcosi]

                okcosi = results_ssp585[mod]['labels'] == reg
                ssp_distr = results_ssp585[mod]['maxgrad'][okcosi]
                # allmedian = np.median(np.concatenate([hist_distr, ssp_distr]))
                histmedian = np.median(hist_distr)

                allpercs['mean'].append(np.mean(hist_distr)-histmedian)
                allpercs['mean'].append(np.mean(ssp_distr)-histmedian)

                for pp in [10, 25, 50, 75, 90]:
                    allpercs['p{}'.format(pp)].append(np.percentile(hist_distr, pp)-histmedian)
                    allpercs['p{}'.format(pp)].append(np.percentile(ssp_distr, pp)-histmedian)

            allpercs_regs.append(allpercs)


        fig = plt.figure(figsize = (16,12))
        axes = []
        for reg in range(4):
            ax = fig.add_subplot(2, 2, reg + 1)
            axes.append(ax)

            ctl.boxplot_on_ax(ax, allpercs_regs[reg], modoks_dub, coloks_dub, edge_colors = edge_colors, positions = positions, versions = versions, ens_colors = ['steelblue', 'indianred'], ens_names = ['hist', 'ssp585'], plot_mean = False)

            ax.axhline(0, color = 'gray', linewidth = 0.5)
            ax.set_xticks([])
            ax.set_title(reg_names[reg])

        ctl.adjust_ax_scale(axes)
        ctl.custom_legend(fig, coloks + ['steelblue', 'indianred'], modoks + ['hist', 'ssp585'], ncol = 4)
        fig.suptitle('Dist of max gradient in 2050-2100 wrt 1964-2014')

        fig.savefig(cart_out + 'Dist_WRmaxgrad_dtr_{}_{}_{}_{}.pdf'.format(*ke))

    ### Guardo le clouds. Faccio differenze tra clouds? O plotto tutto? Boh. Posso intanto vedere se i centroids si spostano sistematicamente. Non credo.

pickle.dump(allpercs_ssp, open(cart_in + 'freqdist_30yr_allssps.p', 'wb'))
pickle.dump(allfreqs, open(cart_in + 'allfreq_allssps.p', 'wb'))


#### Grafico con tutti gli ssp
area = 'EAT'
season = 'NDJFM'
tip = 'refCLUS'

allpercke = allpercs_ssp[(ssp, area, season, tip, reg)].keys()

allsims = ['hist', 'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']
colsim = ctl.color_set(6)

figs = dict()
axes = []
for season in ['NDJFM', 'DJFM']:
    for area in ['EAT', 'PNA']:
        reg_names = reg_names_area[area]
        for tip in ['refCLUS', 'refMODCLUS']:
            figs[(area, tip)] = plt.figure(figsize = (16,12))
            axes = []
            for reg in range(4):
                ax = figs[(area, tip)].add_subplot(2, 2, reg + 1)
                axes.append(ax)

                allpercs = dict()
                for ke in allpercke:
                    allpercs[ke] = [allpercs_ssp[(ssp, area, season, tip, reg)][ke] for ssp in allsims]

                ctl.boxplot_on_ax(ax, allpercs, allsims, colsim, plot_mean = False, plot_ensmeans = False, plot_minmax = True)

                ax.axhline(0, color = 'gray', linewidth = 0.5)
                ax.set_xticks([])
                ax.set_title(reg_names[reg])

            ctl.adjust_ax_scale(axes)

            ctl.custom_legend(figs[(area, tip)], colsim, allsims, ncol = 3)
            #fig.suptitle('30yr bsp of WR freq. in 2050-2100 wrt 1964-2014')
            figs[(area, tip)].savefig(cart_in + 'Dist_30yr_WRfreq_allssp_{}_{}_{}.pdf'.format(area, season, tip))

            figs[(area, tip)] = plt.figure(figsize = (16,12))
            axes = []
            for reg in range(4):
                ax = figs[(area, tip)].add_subplot(2, 2, reg + 1)
                axes.append(ax)

                allpercs = dict()
                allpercs['mean'] = [np.mean(allfreqs[(ssp, area, season, tip, reg)]) for ssp in allsims[1:]]

                for nu in [10, 25, 50, 75, 90]:
                    allpercs['p{}'.format(nu)] = [np.percentile(allfreqs[(ssp,  area, season, tip, reg)], nu) for ssp in allsims[1:]]

                allpercs['ens_min'] = [np.min(allfreqs[(ssp,  area, season, tip, reg)]) for ssp in allsims[1:]]
                allpercs['ens_max'] = [np.max(allfreqs[(ssp,  area, season, tip, reg)]) for ssp in allsims[1:]]

                ctl.boxplot_on_ax(ax, allpercs, allsims[1:], colsim[1:], plot_mean = False, plot_ensmeans = False, plot_minmax = True)

                ax.axhline(0, color = 'gray', linewidth = 0.5)
                ax.set_xticks([])
                ax.set_title(reg_names[reg])

            ctl.adjust_ax_scale(axes)

            ctl.custom_legend(figs[(area, tip)], colsim[1:], allsims[1:], ncol = 3)
            #fig.suptitle('30yr bsp of WR freq. in 2050-2100 wrt 1964-2014')
            figs[(area, tip)].savefig(cart_in + 'WRfreqchange_allssp_{}_{}_{}.pdf'.format(area, season, tip))

            fig = plt.figure(figsize = (16,12))
            axes = []
            for reg in range(4):
                ax = fig.add_subplot(2, 2, reg + 1)
                axes.append(ax)

                allpercs = dict()
                allpercs['mean'] = [np.mean(allfreqs[(ssp, area, season, tip, reg)]) for ssp in allsims[1:]]

                for nu in [10, 25, 50, 75, 90]:
                    allpercs['p{}'.format(nu)] = [np.percentile(allfreqs[(ssp,  area, season, tip, reg)], nu) for ssp in allsims[1:]]

                allpercs['ens_min'] = [np.min(allfreqs[(ssp,  area, season, tip, reg)]) for ssp in allsims[1:]]
                allpercs['ens_max'] = [np.max(allfreqs[(ssp,  area, season, tip, reg)]) for ssp in allsims[1:]]

                ctl.boxplot_on_ax(ax, allpercs, allsims[1:], colsim[1:], plot_mean = False, plot_ensmeans = False, plot_minmax = True)

                ax.axhline(0, color = 'gray', linewidth = 0.5)
                ax.set_xticks([])
                ax.set_title(reg_names[reg])

            ctl.adjust_ax_scale(axes)

            ctl.custom_legend(fig, colsim[1:], allsims[1:], ncol = 3)
            #fig.suptitle('30yr bsp of WR freq. in 2050-2100 wrt 1964-2014')
            fig.savefig(cart_in + 'WRfreqchange_allssp_{}_{}_{}.pdf'.format(area, season, tip))

            if tip == 'refCLUS':
                # Last 20
                fig = plt.figure(figsize = (16,12))
                axes = []
                for reg in range(4):
                    ax = fig.add_subplot(2, 2, reg + 1)
                    axes.append(ax)

                    allpercs = dict()
                    allpercs['mean'] = [np.mean(allfreqs[(ssp, area, season, 'refCLUS_last20', reg)]) for ssp in allsims[1:]]

                    for nu in [10, 25, 50, 75, 90]:
                        allpercs['p{}'.format(nu)] = [np.percentile(allfreqs[(ssp,  area, season, 'refCLUS_last20', reg)], nu) for ssp in allsims[1:]]

                    allpercs['ens_min'] = [np.min(allfreqs[(ssp,  area, season, 'refCLUS_last20', reg)]) for ssp in allsims[1:]]
                    allpercs['ens_max'] = [np.max(allfreqs[(ssp,  area, season, 'refCLUS_last20', reg)]) for ssp in allsims[1:]]

                    ctl.boxplot_on_ax(ax, allpercs, allsims[1:], colsim[1:], plot_mean = False, plot_ensmeans = False, plot_minmax = True)

                    ax.axhline(0, color = 'gray', linewidth = 0.5)
                    ax.set_xticks([])
                    ax.set_title(reg_names[reg])

                ctl.adjust_ax_scale(axes)

                ctl.custom_legend(fig, colsim[1:], allsims[1:], ncol = 3)
                #fig.suptitle('30yr bsp of WR freq. in 2050-2100 wrt 1964-2014')
                fig.savefig(cart_in + 'WRfrch_LAST20_allssp_{}_{}_{}.pdf'.format(area, season, tip))


                # Long 5
                fig = plt.figure(figsize = (16,12))
                axes = []
                for reg in range(4):
                    ax = fig.add_subplot(2, 2, reg + 1)
                    axes.append(ax)

                    allpercs = dict()
                    allpercs['mean'] = [np.mean(allfreqs[(ssp, area, season, 'refCLUS_long5', reg)]) for ssp in allsims[1:]]

                    for nu in [10, 25, 50, 75, 90]:
                        allpercs['p{}'.format(nu)] = [np.percentile(allfreqs[(ssp,  area, season, 'refCLUS_long5', reg)], nu) for ssp in allsims[1:]]

                    allpercs['ens_min'] = [np.min(allfreqs[(ssp,  area, season, 'refCLUS_long5', reg)]) for ssp in allsims[1:]]
                    allpercs['ens_max'] = [np.max(allfreqs[(ssp,  area, season, 'refCLUS_long5', reg)]) for ssp in allsims[1:]]

                    ctl.boxplot_on_ax(ax, allpercs, allsims[1:], colsim[1:], plot_mean = False, plot_ensmeans = False, plot_minmax = True)

                    ax.axhline(0, color = 'gray', linewidth = 0.5)
                    ax.set_xticks([])
                    ax.set_title(reg_names[reg])

                ctl.adjust_ax_scale(axes)

                ctl.custom_legend(fig, colsim[1:], allsims[1:], ncol = 3)
                #fig.suptitle('30yr bsp of WR freq. in 2050-2100 wrt 1964-2014')
                fig.savefig(cart_in + 'WRfrch_long5_allssp_{}_{}_{}.pdf'.format(area, season, tip))


                # Last 20
                fig = plt.figure(figsize = (16,12))
                axes = []
                for reg in range(4):
                    ax = fig.add_subplot(2, 2, reg + 1)
                    axes.append(ax)

                    allpercs = dict()
                    allpercs['mean'] = [np.mean(allfreqs[(ssp, area, season, 'refCLUS_avres', reg)]) for ssp in allsims[1:]]

                    for nu in [10, 25, 50, 75, 90]:
                        allpercs['p{}'.format(nu)] = [np.percentile(allfreqs[(ssp,  area, season, 'refCLUS_avres', reg)], nu) for ssp in allsims[1:]]

                    allpercs['ens_min'] = [np.min(allfreqs[(ssp,  area, season, 'refCLUS_avres', reg)]) for ssp in allsims[1:]]
                    allpercs['ens_max'] = [np.max(allfreqs[(ssp,  area, season, 'refCLUS_avres', reg)]) for ssp in allsims[1:]]

                    ctl.boxplot_on_ax(ax, allpercs, allsims[1:], colsim[1:], plot_mean = False, plot_ensmeans = False, plot_minmax = True)

                    ax.axhline(0, color = 'gray', linewidth = 0.5)
                    ax.set_xticks([])
                    ax.set_title(reg_names[reg])

                ctl.adjust_ax_scale(axes)

                ctl.custom_legend(fig, colsim[1:], allsims[1:], ncol = 3)
                #fig.suptitle('30yr bsp of WR freq. in 2050-2100 wrt 1964-2014')
                fig.savefig(cart_in + 'WR_avres_allssp_{}_{}_{}.pdf'.format(area, season, tip))
