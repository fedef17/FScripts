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
cart_out_orig = cart_in + 'Results_v2/'

file_hist_refEOF = cart_in + 'cmip6_hist/out_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refEOF_dtr.p'
file_hist = cart_in + 'cmip6_hist/out_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refCLUS_dtr.p'
#file_hist_refEOF = cart_in + 'cmip6_hist/out_cmip6_hist_NDJFM_EAT_4clus_4pcs_1964-2014_refEOF.p'
gen_file_ssp = cart_in + 'cmip6_{}/out_cmip6_{}_NDJFM_{}_4clus_4pcs_2050-2100_refCLUS_dtr.p'

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
        results_ssp[ssp], _ = pickle.load(open(gen_file_ssp.format(ssp, ssp, area), 'rb'))

        # Erasing incomplete runs
        for ke in tuple(results_ssp[ssp].keys()):
            if len(results_ssp[ssp][ke]['labels']) < 7000:
                del results_ssp[ssp][ke]

    mod_hist = np.unique([cos.split('_')[0] for cos in results_hist.keys()])
    mod_ssp = dict()
    for ssp in allssps:
        mod_ssp[ssp] = np.unique([cos.split('_')[0] for cos in results_ssp[ssp].keys()])

    ### Voglio: freqs, resid_time, resid_time_90, eff_centroids_ssp, centroids_hist (e patcor, rms)
    freqs = dict() # tot50 e last20
    residtimes = dict() # mean e p90
    eff_centroids = dict()

    patterns_refEOF = dict()

    for mod in mod_hist:
        allmems = [cos for cos in results_hist.keys() if cos.split('_')[0] == mod]
        ## Attach all members labels
        alllabs = np.concatenate([results_hist[mem]['labels'] for mem in allmems])
        freqs[('hist', mod, 'tot50')] = ctl.calc_clus_freq(alllabs, numclus)

        for reg in range(numclus):
            alltimes = np.concatenate([results_hist[mem]['resid_times'][reg] for mem in allmems])
            residtimes[('hist', mod, 'mean', reg)] = np.mean(alltimes)
            residtimes[('hist', mod, 'p90', reg)] = np.percentile(alltimes, 90)

        alllabs_20 = []
        alltimes_20 = []
        for mem in allmems:
            dat1 = pd.Timestamp('09-01-1995').to_pydatetime()
            dat2 = pd.Timestamp('04-01-2014').to_pydatetime()
            labs, dats = ctl.sel_time_range(results_hist[mem]['labels'], results_hist[mem]['dates'], (dat1, dat2))
            alllabs_20.append(labs)
            print(mem, len(labs))

            restim, _, _ = ctl.calc_regime_residtimes(labs, dats)
            alltimes_20.append(restim)

        alllabs_20 = np.concatenate(alllabs_20)
        alltimes_20 = [np.concatenate([cos[reg] for cos in alltimes_20]) for reg in range(numclus)]
        for reg in range(numclus):
            residtimes[('hist', mod, 'mean_last20', reg)] = np.mean(alltimes_20[reg])
            residtimes[('hist', mod, 'p90_last20', reg)] = np.percentile(alltimes_20[reg], 90)

        freqs[('hist', mod, 'last20')] = ctl.calc_clus_freq(alllabs_20, numclus)

        for reg in range(numclus):
            patterns_refEOF[('centroids', mod, reg)] = np.mean([res_hist_refEOF[mem]['centroids'][reg] for mem in allmems], axis = 0)
            patterns_refEOF[('patcor', mod, reg)] = np.mean([res_hist_refEOF[mem]['patcor'][reg] for mem in allmems])
            patterns_refEOF[('centdist', mod, reg)] = np.mean([ctl.distance(res_hist_refEOF[mem]['centroids'][reg], results_ref['centroids'][reg]) for mem in allmems])

        patterns_refEOF[('var_ratio', mod)] = np.mean([res_hist_refEOF[mem]['var_ratio'] for mem in allmems])

        for reg in range(numclus):
            eff_centroids[('hist', mod, reg)] = np.mean([results_hist[mem]['eff_centroids'][reg] for mem in allmems], axis = 0)


    for ssp in allssps:
        for mod in mod_ssp[ssp]:
            allmems = [cos for cos in results_ssp[ssp].keys() if cos.split('_')[0] == mod]
            ## Attach all members labels
            alllabs = np.concatenate([results_ssp[ssp][mem]['labels'] for mem in allmems])
            freqs[(ssp, mod, 'tot50')] = ctl.calc_clus_freq(alllabs, numclus)

            alllabs_20 = []
            alltimes_20 = []
            for mem in allmems:
                dat1 = pd.Timestamp('09-01-2081').to_pydatetime()
                dat2 = pd.Timestamp('04-01-2100').to_pydatetime()
                labs, dats = ctl.sel_time_range(results_ssp[ssp][mem]['labels'], results_ssp[ssp][mem]['dates'], (dat1, dat2))
                alllabs_20.append(labs)
                restim, _, _ = ctl.calc_regime_residtimes(labs, dats)
                alltimes_20.append(restim)

            alllabs_20 = np.concatenate(alllabs_20)
            alltimes_20 = [np.concatenate([cos[reg] for cos in alltimes_20]) for reg in range(numclus)]
            for reg in range(numclus):
                residtimes[(ssp, mod, 'mean_last20', reg)] = np.mean(alltimes_20[reg])
                residtimes[(ssp, mod, 'p90_last20', reg)] = np.percentile(alltimes_20[reg], 90)

            freqs[(ssp, mod, 'last20')] = ctl.calc_clus_freq(alllabs_20, numclus)

            for reg in range(numclus):
                alltimes = np.concatenate([results_ssp[ssp][mem]['resid_times'][reg] for mem in allmems])
                residtimes[(ssp, mod, 'mean', reg)] = np.mean(alltimes)
                residtimes[(ssp, mod, 'p90', reg)] = np.percentile(alltimes, 90)

            for reg in range(numclus):
                eff_centroids[(ssp, mod, reg)] = np.mean([results_ssp[ssp][mem]['eff_centroids'][reg] for mem in allmems], axis = 0)

    ##### salvo le distrib per ogni ssp. Quelle assolute (con tutti i modelli) e quelle relative (solo con i modelli che ci sono anche in hist)
    freqs[('hist', 'all', 'tot50')] = np.stack([freqs[('hist', mod, 'tot50')] for mod in mod_hist])
    freqs[('hist', 'all', 'last20')] = np.stack([freqs[('hist', mod, 'last20')] for mod in mod_hist])
    for reg in range(numclus):
        for cos in ['mean', 'p90', 'mean_last20', 'p90_last20']:
            residtimes[('hist', 'all', cos, reg)] = np.array([residtimes[('hist', mod, cos, reg)] for mod in mod_hist])

    for ssp in allssps:
        freqs[(ssp, 'all', 'tot50')] = np.stack([freqs[(ssp, mod, 'tot50')] for mod in mod_ssp[ssp]])
        freqs[(ssp, 'all', 'last20')] = np.stack([freqs[(ssp, mod, 'last20')] for mod in mod_ssp[ssp]])

        # rel
        modoks = [mod for mod in mod_hist if mod in mod_ssp[ssp]]
        freqs[(ssp, 'rel', 'tot50')] = np.stack([freqs[(ssp, mod, 'tot50')]-freqs[('hist', mod, 'tot50')] for mod in modoks])
        freqs[(ssp, 'rel', 'last20')] = np.stack([freqs[(ssp, mod, 'last20')]-freqs[('hist', mod, 'last20')] for mod in modoks])

        for reg in range(numclus):
            for cos in ['mean', 'p90', 'mean_last20', 'p90_last20']:
                residtimes[(ssp, 'all', cos, reg)] = np.array([residtimes[(ssp, mod, cos, reg)] for mod in mod_ssp[ssp]])

                residtimes[(ssp, 'rel', cos, reg)] = np.array([residtimes[(ssp, mod, cos, reg)]-residtimes[('hist', mod, cos, reg)] for mod in modoks])


    pickle.dump([freqs, residtimes, eff_centroids, patterns_refEOF], open(cart_out + 'allresults_dicts_{}.p'.format(area), 'wb'))
    freqs, residtimes, eff_centroids, patterns_refEOF = pickle.load(open(cart_out + 'allresults_dicts_{}.p'.format(area), 'rb'))

    #### Grafico con tutti gli ssp
    allsims = ['hist', 'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']
    colsim = ctl.color_set(6)

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

            ctl.boxplot_on_ax(ax, allpercs, allsims[1:], colsim[1:], plot_mean = False, plot_ensmeans = False, plot_minmax = True)

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

            ctl.boxplot_on_ax(ax, allpercs, allsims[1:], colsim[1:], plot_mean = False, plot_ensmeans = False, plot_minmax = True)

            ax.axhline(0, color = 'gray', linewidth = 0.5)
            ax.set_xticks([])
            ax.set_title(reg_names[reg])

        ctl.adjust_ax_scale(axes)

        ctl.custom_legend(fig, colsim[1:], allsims[1:], ncol = 3)
        #fig.suptitle('30yr bsp of WR freq. in 2050-2100 wrt 1964-2014')
        fig.savefig(cart_out + 'Restime_change_allssp_{}_{}.pdf'.format(area, cos))


    #### Absolute values
    for cos in ['last20', 'tot50']:
        fig = plt.figure(figsize = (16,12))
        axes = []
        for reg in range(4):
            ax = fig.add_subplot(2, 2, reg + 1)
            axes.append(ax)

            allpercs = dict()
            allpercs['mean'] = [np.mean(freqs[(ssp, 'all', cos)][:, reg]) for ssp in allsims]

            for nu in [10, 25, 50, 75, 90]:
                allpercs['p{}'.format(nu)] = [np.percentile(freqs[(ssp, 'all', cos)][:, reg], nu) for ssp in allsims]

            allpercs['ens_min'] = [np.min(freqs[(ssp, 'all', cos)][:, reg]) for ssp in allsims]
            allpercs['ens_max'] = [np.max(freqs[(ssp, 'all', cos)][:, reg]) for ssp in allsims]

            ax.axhline(allpercs['p50'][0], color = 'gray', linewidth = 0.5)
            ctl.boxplot_on_ax(ax, allpercs, allsims, colsim, plot_mean = False, plot_ensmeans = False, plot_minmax = True)

            ax.set_xticks([])
            ax.set_title(reg_names[reg])

        #ctl.adjust_ax_scale(axes)

        ctl.custom_legend(fig, colsim, allsims, ncol = 3)
        #fig.suptitle('30yr bsp of WR freq. in 2050-2100 wrt 1964-2014')
        fig.savefig(cart_out + 'WRfreq_allssp_{}_{}.pdf'.format(area, cos))

    dat1 = pd.Timestamp('09-01-1995').to_pydatetime()
    dat2 = pd.Timestamp('04-01-2014').to_pydatetime()
    labs, dats = ctl.sel_time_range(results_ref['labels'], results_ref['dates'], (dat1, dat2))
    results_ref['freq_clus_last20'] = ctl.calc_clus_freq(labs, numclus)

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
        ctl.boxplot_on_ax(ax, allpercs, allsims_dub, colsims_dub, plot_mean = False, plot_ensmeans = False, plot_minmax = False, positions = positions)

        ax.set_xticks([])
        ax.set_title(reg_names[reg])
        ax.scatter(positions[0], results_ref['freq_clus'][reg], color = 'black', marker = '*', s = 5)
        ax.scatter(positions[1], results_ref['freq_clus_last20'][reg], color = 'black', marker = '*', s = 2)

    #ctl.adjust_ax_scale(axes)

    ctl.custom_legend(fig, colsim, allsims, ncol = 3)
    #fig.suptitle('30yr bsp of WR freq. in 2050-2100 wrt 1964-2014')
    fig.savefig(cart_out + 'WRfreq_allssp_{}_bothrefs.pdf'.format(area, cos))


    for cos in ['mean', 'p90', 'mean_last20', 'p90_last20']:
        fig = plt.figure(figsize = (16,12))
        axes = []
        for reg in range(4):
            ax = fig.add_subplot(2, 2, reg + 1)
            axes.append(ax)

            allpercs = dict()
            allpercs['mean'] = [np.mean(residtimes[(ssp, 'all', cos, reg)]) for ssp in allsims]

            for nu in [10, 25, 50, 75, 90]:
                allpercs['p{}'.format(nu)] = [np.percentile(residtimes[(ssp, 'all', cos, reg)], nu) for ssp in allsims]

            allpercs['ens_min'] = [np.min(residtimes[(ssp, 'all', cos, reg)]) for ssp in allsims]
            allpercs['ens_max'] = [np.max(residtimes[(ssp, 'all', cos, reg)]) for ssp in allsims]

            ax.axhline(allpercs['p50'][0], color = 'gray', linewidth = 0.5)
            ctl.boxplot_on_ax(ax, allpercs, allsims, colsim, plot_mean = False, plot_ensmeans = False, plot_minmax = True)

            # ax.axhline(0, color = 'gray', linewidth = 0.5)
            ax.set_xticks([])
            ax.set_title(reg_names[reg])

        #ctl.adjust_ax_scale(axes)

        ctl.custom_legend(fig, colsim, allsims, ncol = 3)
        #fig.suptitle('30yr bsp of WR freq. in 2050-2100 wrt 1964-2014')
        fig.savefig(cart_out + 'Restime_allssp_{}_{}.pdf'.format(area, cos))


    ######### PATTERNS
    colormods = ctl.color_set(len(mod_hist))
    fig = plt.figure(figsize = (16,12))
    axes = []
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg + 1, polar = True)
        axes.append(ax)
        # patterns_refEOF[('patcor', mod, reg)] = np.mean([res_hist_refEOF[mem]['patcor'][reg] for mem in allmems])
        # patterns_refEOF[('centdist', mod, reg)]

        models = [patterns_refEOF[('centroids', mod, reg)] for mod in mod_hist]
        observation = results_ref['centroids'][reg]
        ctl.Taylor_plot(models, observation, ax = ax, colors = colormods, only_first_quarter = True)

        ax.set_title(reg_names[reg])

    ctl.adjust_ax_scale(axes)

    ctl.custom_legend(fig, colormods, mod_hist, ncol = 5)
    #fig.suptitle('30yr bsp of WR freq. in 2050-2100 wrt 1964-2014')
    fig.savefig(cart_out + 'Taylor_hist_refEOF_{}.pdf'.format(area))


    ## scatter plots var_ratio contro fr change
    colini = [col for col,mod in zip(colormods, mod_hist) if mod in mod_ssp['ssp585']]
    fig = plt.figure(figsize = (16,12))
    axes = []
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg + 1)
        axes.append(ax)
        #varrat = [patterns_refEOF[('var_ratio', mod)] for mod in mod_hist if mod in mod_ssp['ssp585']]
        varrat = [patterns_refEOF[('centdist', mod, reg)] for mod in mod_hist if mod in mod_ssp['ssp585']]
        frchan = freqs[('ssp585', 'rel', 'tot50')][:, reg]

        ax.scatter(varrat, frchan, c = colini)
        ax.set_title(reg_names[reg])

    ctl.adjust_ax_scale(axes)

    ctl.custom_legend(fig, colini, mod_hist, ncol = 5)
    #fig.suptitle('30yr bsp of WR freq. in 2050-2100 wrt 1964-2014')
    fig.savefig(cart_out + 'Scatter_hist_centdist_{}.pdf'.format(area))

    ### All centroids
    fig = plt.figure(figsize = (16,12))
    ax = fig.add_subplot(111)
    for reg in range(numclus):
        models = [patterns_refEOF[('centroids', mod, reg)] for mod in mod_hist]
        xs = [cos[0] for cos in models]
        ys = [cos[1] for cos in models]

        ax.scatter(xs, ys, c = colormods, s = 10)

    observation = results_ref['centroids']
    ax.scatter(observation[:, 0], observation[:, 1], c = 'black', s = 50)

    ctl.custom_legend(fig, colormods, mod_hist, ncol = 5)
    fig.savefig(cart_out + 'Phasespace_hist_refEOF_{}.pdf'.format(area))

    ### All centroids refCLUS
    fig = plt.figure(figsize = (16,12))
    ax = fig.add_subplot(111)
    for reg in range(numclus):
        models = [eff_centroids[('hist', mod, reg)] for mod in mod_hist]
        xs = [cos[0] for cos in models]
        ys = [cos[1] for cos in models]

        ax.scatter(xs, ys, c = colsim[0], s = 10)

        for col, ssp in zip(colsim[1:], allsims[1:]):
            models = [eff_centroids[(ssp, mod, reg)] for mod in mod_ssp[ssp]]
            xs = [cos[0] for cos in models]
            ys = [cos[1] for cos in models]

            ax.scatter(xs, ys, c = col, s = 10)

    observation = results_ref['centroids']
    ax.scatter(observation[:, 0], observation[:, 1], c = 'black', s = 50)

    ctl.custom_legend(fig, colsim, allsims, ncol = 3)
    fig.savefig(cart_out + 'Phasespace_allssp_refCLUS_{}.pdf'.format(area))


    print('Building pdfs...')
    fig = plt.figure(figsize = (16,12))
    ax = fig.add_subplot(111)

    xss = np.linspace(-2000., 2000., 201)
    xi_grid, yi_grid = np.meshgrid(xss, xss)

    print('hist')
    okpcs_all = []
    for modmem in results_hist.keys():
        dat1 = pd.Timestamp('09-01-1995').to_pydatetime()
        dat2 = pd.Timestamp('04-01-2014').to_pydatetime()
        okpc = results_hist[modmem]['pcs'][:, :2]
        okpcok, dats = ctl.sel_time_range(okpc, results_hist[modmem]['dates'], (dat1, dat2))
        okpcs_all.append(okpcok)

    okpc = np.concatenate(okpcs_all, axis = 0)
    cent = np.mean(okpc, axis = 0)

    kufu = ctl.calc_pdf(okpc.T)

    zi = kufu(np.vstack([xi_grid.flatten(), yi_grid.flatten()]))
    zi = zi/np.max(zi)
    pdfssp[('hist', 'all')] = zi

    for reg in range(numclus):
        print(reg)
        okpcs_all = []
        for modmem in results_hist.keys():
            okclus = results_hist[modmem]['labels'] == reg
            okpc = results_hist[modmem]['pcs'][okclus, :2]
            okpcs_all.append(okpc)

        okpc = np.concatenate(okpcs_all, axis = 0)
        cent = np.mean(okpc, axis = 0)

        kufu = ctl.calc_pdf(okpc.T)
        #pdfssp[('hist', reg)] = kufu

        cmappa = ctl.custom_alphagradient_cmap(colsim[0])

        zi = kufu(np.vstack([xi_grid.flatten(), yi_grid.flatten()]))
        zi = zi/np.max(zi)
        pdfssp[('hist', reg)] = zi

        cont = ax.contour(xi_grid, yi_grid, zi.reshape(xi_grid.shape), [0.5, 0.8], cmap = cmappa)#, linewidths = lw)
        ax.scatter(cent[0], cent[1], color = colsim[0], s = 10, marker = 'x')

    for ssp, col in zip(allssps, colsim[1:]):
        okpcs_all = []
        for modmem in results_ssp[ssp].keys():
            dat1 = pd.Timestamp('09-01-2081').to_pydatetime()
            dat2 = pd.Timestamp('04-01-2100').to_pydatetime()
            okpc = results_ssp[ssp][modmem]['pcs'][:, :2]
            okpcok, dats = ctl.sel_time_range(okpc, results_ssp[ssp][modmem]['dates'], (dat1, dat2))
            okpcs_all.append(okpcok)

        okpc = np.concatenate(okpcs_all, axis = 0)
        cent = np.mean(okpc, axis = 0)

        kufu = ctl.calc_pdf(okpc.T)

        zi = kufu(np.vstack([xi_grid.flatten(), yi_grid.flatten()]))
        zi = zi/np.max(zi)
        pdfssp[(ssp, 'all')] = zi

        for reg in range(numclus):
            print(ssp, reg)
            okpcs_all = []
            for modmem in results_ssp[ssp].keys():
                okclus = results_ssp[ssp][modmem]['labels'] == reg
                okpc = results_ssp[ssp][modmem]['pcs'][okclus, :2]
                okpcs_all.append(okpc)

            okpc = np.concatenate(okpcs_all, axis = 0)
            cent = np.mean(okpc, axis = 0)

            kufu = ctl.calc_pdf(okpc.T)

            cmappa = ctl.custom_alphagradient_cmap(col)

            zi = kufu(np.vstack([xi_grid.flatten(), yi_grid.flatten()]))
            zi = zi/np.max(zi)
            pdfssp[(ssp, reg)] = zi

            cont = ax.contour(xi_grid, yi_grid, zi.reshape(xi_grid.shape), [0.5, 0.8], cmap = cmappa)
            ax.scatter(cent[0], cent[1], color = col, s = 10, marker = 'x')

    ctl.custom_legend(fig, colsim, allsims, ncol = 3)
    fig.savefig(cart_out + 'Clouds_allssp_refCLUS_{}.pdf'.format(area))

    pickle.dump(pdfssp, open(cart_out + 'pdfs_refCLUS_last20_{}.p'.format(area), 'wb'))
