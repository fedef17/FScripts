#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
from matplotlib import pyplot as plt
from matplotlib import cm

import pickle

import climtools_lib as ctl
import climdiags as cd

from matplotlib.colors import LogNorm
from datetime import datetime

from scipy import stats

#######################################
vtag = 'v7'

#cart_out = '/home/fabiano/Research/articoli/Papers/primavera_regimes/figures/figures_{}/'.format(vtag)
cart_out = '/home/fabiano/Research/articoli/Papers/primavera_regimes/figures/figures_{}_refCLUS/'.format(vtag)
if not os.path.exists(cart_out): os.mkdir(cart_out)

cart = '/home/fabiano/Research/lavori/WeatherRegimes/prima_coup_{}/'.format(vtag)
#filogen = cart + 'out_prima_coup_{}_DJF_EAT_4clus_4pcs_1957-2014_refEOF.p'.format(vtag)
filogen = cart + 'out_prima_coup_{}_DJF_EAT_4clus_4pcs_1957-2014_refCLUS.p'.format(vtag)

model_names = ['AWI-CM-1-1-LR', 'AWI-CM-1-1-HR', 'CMCC-CM2-HR4', 'CMCC-CM2-VHR4', 'CNRM-CM6-1', 'CNRM-CM6-1-HR', 'EC-Earth3P', 'EC-Earth3P-HR', 'ECMWF-IFS-LR', 'ECMWF-IFS-MR', 'ECMWF-IFS-HR', 'MPI-ESM1-2-HR', 'MPI-ESM1-2-XR', 'HadGEM3-GC31-LL', 'HadGEM3-GC31-MM', 'HadGEM3-GC31-HM', 'HadGEM3-GC31-HH']

vers = ['LR', 'HR', 'LR', 'HR', 'LR', 'HR', 'LR', 'HR', 'LR', 'MR', 'HR', 'LR', 'HR', 'LR', 'MR', 'MR', 'HR'] + ['OBS']
model_coups = ['AWI-CM-1-1', 'CMCC-CM2', 'CNRM-CM6-1', 'EC-Earth3P', 'ECMWF-IFS', 'MPI-ESM1-2', 'HadGEM3-GC31']

model_names_all = model_names + ['ERA']

colors = ctl.color_set(len(model_coups)*2, sns_palette = 'Paired')
colorscoup = [np.mean([col1, col2], axis = 0) for col1, col2 in zip(colors[:-1:2], colors[1::2])]

mr_cos = np.where(np.array(vers) == 'MR')[0]
for gi in mr_cos:
    colors.insert(gi, np.mean(colors[gi-1:gi+1], axis = 0))

colors_wERA = colors + ['black']

# # Ora tutte le statistiche
# atm_resolution = np.array([250, 100, 100, 25, 250, 50, 100, 50, 50, 50, 25, 100, 50, 250, 100, 50, 50])
# oce_resolution = np.array([50, 25, 25, 25, 100, 25, 100, 25, 100, 25, 25, 40, 40, 100, 25, 25, 8])
atm_resolution = np.array([245, 94, 97, 23, 250, 44, 100, 46, 48, 50, 27, 103, 52, 255, 106, 54, 56])
oce_resolution = np.array([50, 29, 21, 22, 96, 23, 98, 24, 102, 25, 26, 38, 42, 104, 27, 28, 8])

if len(atm_resolution) != len(model_names): sys.exit()
if len(oce_resolution) != len(model_names): sys.exit()

sym_all = []
for cos in vers:
    if cos == 'LR':
        sym_all.append('o')
    elif cos == 'HR':
        sym_all.append('P')
    else:
        sym_all.append('d')

lw = 2
wi = 0.6

#regnam = ['NAO +', 'Sc. BL', 'NAO -', 'AR']
regnam = ['NAO +', 'Sc. BL', 'AR', 'NAO -']

results, results_ref = pickle.load(open(filogen, 'rb'))
results['ERA_0'] = results_ref

results['ERA_0']['significance'] = 99.8

for mod in results:
    #results[mod]['significance'] = significance[mod]
    results[mod]['varopt'] = ctl.calc_varopt_molt(results[mod]['pcs'], results[mod]['centroids'], results[mod]['labels'])

    #results[mod]['autocorr_wlag'] = ctl.calc_autocorr_wlag(results[mod]['pcs'], results[mod]['dates'])

results.pop('HadGEM3-GC31-LL_r1i1p2f1')
results.pop('EC-Earth3P_r1i1p1f1')
results.pop('EC-Earth3P-HR_r1i1p1f1')

# Qui ci metto: RMS, patcor, significance, optimal_ratio, clus_frequency, clus_persistence
allnams = ['significance', 'varopt', 'RMS', 'patcor', 'freq_clus', 'resid_times_mean', 'resid_times_90']
longnams = ['Sharpness', 'Optimal ratio', 'RMS', '1 - Pattern correlation', 'Regime frequency bias', 'Average regime persistence', '90th percentile regime persistence']

figures = []
for cos, tit in zip(allnams[:2], longnams[:2]):
    if cos not in results['EC-Earth3P_r1i1p2f1'].keys():
        continue
    # WITH ATM RESOLUTION
    fig = plt.figure(figsize = (16, 12))
    ax = plt.subplot(111)

    allvals_mods = []
    allerr_mods = []
    for mod, res, col, mark in zip(model_names, atm_resolution, colors, sym_all):
        allvals = []
        for ke in results.keys():
            if mod == ke.split('_')[0]:
                allvals.append(results[ke][cos])

        if len(allvals) == 1:
            val = allvals[0]
            ax.scatter(res, val, color = col, marker = mark, s = 100, linewidth = lw)
            allvals_mods.append(val)
            allerr_mods.append(None)
        else:
            val = np.mean(allvals)
            errval = np.std(allvals)
            ax.scatter(res, val, color = col, marker = mark, s = 100, linewidth = lw)
            ax.errorbar(res, val, yerr = errval, color = col, linewidth = lw, capsize = 4)
            allvals_mods.append(val)
            allerr_mods.append(errval)

    ax.axhline(results['ERA_0'][cos], color = 'grey', alpha = 0.6)

    ax.set_ylabel(tit)
    ax.set_xlabel('Eff. atm. resolution (km)')

    fig = ctl.custom_legend(fig, colors, model_names, ncol = 6)
    fig.savefig(cart_out + cos+'_{}_atmresolution.pdf'.format(vtag))
    figures.append(fig)

    # WITH OCEAN RESOLUTION
    fig = plt.figure(figsize = (16, 12))
    ax = plt.subplot(111)

    for mod, res, col, mark in zip(model_names, oce_resolution, colors, sym_all):
        allvals = []
        for ke in results.keys():
            if mod == ke.split('_')[0]:
                allvals.append(results[ke][cos])

        if len(allvals) == 1:
            val = allvals[0]
            ax.scatter(res, val, color = col, marker = mark, s = 100, linewidth = lw)
        else:
            val = np.mean(allvals)
            errval = np.std(allvals)
            ax.scatter(res, val, color = col, marker = mark, s = 100, linewidth = lw)
            ax.errorbar(res, val, yerr = errval, color = col, linewidth = lw, capsize = 4)

    ax.axhline(results['ERA_0'][cos], color = 'grey', alpha = 0.6)

    ax.set_ylabel(tit)
    ax.set_xlabel('Eff. ocean resolution (km)')

    fig = ctl.custom_legend(fig, colors, model_names, ncol = 6)
    fig.savefig(cart_out + cos+'_{}_oceresolution.pdf'.format(vtag))
    figures.append(fig)


    # Normal bar plot
    fig = plt.figure(figsize = (16, 12))
    ax = plt.subplot(111)

    ax.set_ylabel(tit)
    ax.set_xticks([])

    i = 0
    wi = 0.6
    allvals = []
    allerr = []
    for mod, col, vv in zip(model_names, colors, vers):
        rmsall = []
        for ke in results.keys():
            if mod == ke.split('_')[0]:
                rmsall.append(results[ke][cos])

        if len(rmsall) == 0: continue

        val = np.mean(rmsall)
        allvals.append(val)

        if len(rmsall) == 1:
            ax.bar(i, val, width = wi, color = col)
            allerr.append(np.nan)
        else:
            print('{} has {} ens members'.format(mod, len(rmsall)))
            okrms = np.mean(rmsall)
            rmserr_lo = np.min(rmsall)
            rmserr_hi = np.max(rmsall)
            print(okrms, rmserr_lo, rmserr_hi)
            yerr = np.array([[okrms-rmserr_lo, rmserr_hi-okrms]]).T
            allerr.append(np.std(rmsall))

            ax.bar(i, okrms, width = wi, color = col, yerr = yerr, capsize = 5)

        i += 0.7
        if vv == 'HR': i += 0.2

    all_LR = np.mean([val for val, vv in zip(allvals, vers) if vv == 'LR'])
    all_HR = np.mean([val for val, vv in zip(allvals, vers) if vv == 'HR'])
    all_LR_std = np.std([val for val, vv in zip(allvals, vers) if vv == 'LR'])
    all_HR_std = np.std([val for val, vv in zip(allvals, vers) if vv == 'HR'])

    i+=0.2
    col_LR = 'grey'
    col_HR = 'black'
    ax.bar(i, all_LR, width = wi, color = col_LR, yerr = all_LR_std, capsize = 5)
    i += 0.7
    ax.bar(i, all_HR, width = wi, color = col_HR, yerr = all_HR_std, capsize = 5)

    ax.axhline(results['ERA_0'][cos], color = 'grey', alpha = 0.6)

    ctl.custom_legend(fig, colors, model_names, ncol = 6)
    fig.savefig(cart_out + '{}_{}_barplot.pdf'.format(cos, vtag))
    figures.append(fig)

    # HR bias - LR bias
    fig = plt.figure(figsize = (16, 12))

    ax = plt.subplot(111)
    eraval = results['ERA_0'][cos]
    allvalsdiff = []
    allerrdiff = []
    for modcou in model_coups:
        v_hr = [val for val, mod, vv in zip(allvals, model_names, vers) if modcou in mod and vv == 'HR'][0]
        v_lr = [val for val, mod, vv in zip(allvals, model_names, vers) if modcou in mod and vv == 'LR'][0]
        vdiff = v_hr-v_lr
        allvalsdiff.append(vdiff)

        v_err_hr = [val for val, mod, vv in zip(allerr, model_names, vers) if modcou in mod and vv == 'HR'][0]
        v_err_lr = [val for val, mod, vv in zip(allerr, model_names, vers) if modcou in mod and vv == 'LR'][0]
        errdiff = np.nanmax([v_err_hr, v_err_lr])
        allerrdiff.append(errdiff)

        # if v_err_hr is None and v_err_lr is None:
        #     errdiff = None
        # elif v_err_hr is not None and

    i = 0
    for val, err, col, modcou in zip(allvalsdiff, allerrdiff, colorscoup, model_coups):
        if np.isnan(err):
            ax.bar(i, val, color = col, width = wi)
        else:
            ax.bar(i, val, color = col, width = wi, yerr = err, capsize = 5)
        i+=0.7

    i+=0.3
    ax.bar(i, np.mean(allvalsdiff), color = 'grey', width = wi, yerr = np.std(allvalsdiff), capsize = 5)

    ax.set_xticks([])
    ax.axhline(0., color = 'grey', alpha = 0.6)

    ax.set_ylabel(tit + ' (difference: HR - LR)')
    #ax.set_xlabel('Eff. atm. resolution (km)')

    fig = ctl.custom_legend(fig, colorscoup+['grey'], model_coups+['avg'])
    fig.savefig(cart_out + cos +'_{}_1vs1.pdf'.format(vtag))
    figures.append(fig)


nuvars = dict()
for mod in model_names_all:
    nuvars[mod] = dict()
    for cos in allnams[2:5]:
        rmsall = []
        for ke in results.keys():
            if mod == ke.split('_')[0] and cos in results[ke].keys():
                rmsall.append(results[ke][cos])

        if len(rmsall) == 0:
            print('aFJAFJAJFAFJAJFAJFJFJAFJFAJ!!!!')
            print(mod)
            continue

        val = np.mean(rmsall, axis = 0)
        nuvars[mod][cos] = val

        nuvars[mod][cos+'_err'] = np.std(rmsall, axis = 0)



for cos, tit in zip(allnams[2:5], longnams[2:5]):
    for tipres, resolution in zip(['atm', 'oce'], [atm_resolution, oce_resolution]):
        fig = plt.figure(figsize = (16, 12))
        ind = 0
        axes = []
        for reg in range(4):
            ind += 1
            ax = plt.subplot(2, 2, ind)

            if cos == 'patcor':
                allvals = np.array([(1.-nuvars[mod][cos][reg]) for mod in model_names])
                allerr = np.array([nuvars[mod][cos+'_err'][reg] for mod in model_names])
            elif cos == 'RMS':
                allvals = np.array([nuvars[mod][cos][reg] for mod in model_names])
                allerr = np.array([nuvars[mod][cos+'_err'][reg] for mod in model_names])
            elif cos == 'freq_clus':
                allvals = np.array([nuvars[mod][cos][reg]-nuvars['ERA'][cos][reg] for mod in model_names])
                allerr = np.array([nuvars[mod][cos+'_err'][reg] for mod in model_names])

            for res, err, val, col, mark in zip(resolution, allerr, allvals, colors, sym_all):
                ax.scatter(res, val, color = col, marker = mark, s = 100, linewidth = lw)
                if not np.isnan(err):
                    ax.errorbar(res, val, yerr = err, color = col, linewidth = lw, capsize = 5)

            ax.set_ylabel(tit)
            #ax.set_xlabel('Eff. atm. resolution (km)')

            ax.axhline(0., color = 'grey', alpha = 0.6)

            # ax.set_title('Reg {}'.format(reg))
            ax.set_title(regnam[reg])
            axes.append(ax)
        ctl.adjust_ax_scale(axes)

        fig.suptitle(tit)
        plt.subplots_adjust(top = 0.9)
        fig = ctl.custom_legend(fig, colors, model_names, ncol = 6)
        fig.savefig(cart_out + cos+'_{}_{}resolution.pdf'.format(vtag, tipres))
        figures.append(fig)

    fig = plt.figure(figsize = (16, 12))
    ind = 0
    axes = []
    for reg in range(4):
        ind += 1
        ax = plt.subplot(2, 2, ind)

        if cos == 'patcor':
            allvals = np.array([(1.-nuvars[mod][cos][reg]) for mod in model_names])
            allerr = np.array([nuvars[mod][cos+'_err'][reg] for mod in model_names])
        elif cos == 'RMS':
            allvals = np.array([nuvars[mod][cos][reg] for mod in model_names])
            allerr = np.array([nuvars[mod][cos+'_err'][reg] for mod in model_names])
        elif cos == 'freq_clus':
            allvals = np.array([nuvars[mod][cos][reg]-nuvars['ERA'][cos][reg] for mod in model_names])
            allerr = np.array([nuvars[mod][cos+'_err'][reg] for mod in model_names])

        allvalsdiff = []
        allerrdiff = []
        for modcou in model_coups:
            v_hr = [val for val, mod, vv in zip(allvals, model_names, vers) if modcou in mod and vv == 'HR'][0]
            v_lr = [val for val, mod, vv in zip(allvals, model_names, vers) if modcou in mod and vv == 'LR'][0]
            #vdiff = abs(v_lr-eraval)-abs(v_hr-eraval)
            vdiff = v_hr-v_lr
            allvalsdiff.append(vdiff)

            v_err_hr = [val for val, mod, vv in zip(allerr, model_names, vers) if modcou in mod and vv == 'HR'][0]
            v_err_lr = [val for val, mod, vv in zip(allerr, model_names, vers) if modcou in mod and vv == 'LR'][0]
            errdiff = np.nanmax([v_err_hr, v_err_lr])
            allerrdiff.append(errdiff)

        i = 0
        for val, err, col, modcou in zip(allvalsdiff, allerrdiff, colorscoup, model_coups):
            if np.isnan(err):
                ax.bar(i, val, color = col, width = wi)
            else:
                ax.bar(i, val, color = col, width = wi, yerr = err, capsize = 5)
            i+=0.7

        i+=0.3
        ax.bar(i, np.mean(allvalsdiff), color = 'grey', width = wi, yerr = np.std(allvalsdiff), capsize = 5)

        ax.set_xticks([])
        ax.axhline(0., color = 'grey', alpha = 0.6)

        ax.set_ylabel(tit + ' (difference: HR - LR)')

        ax.set_title(regnam[reg])
        axes.append(ax)
    ctl.adjust_ax_scale(axes)

    fig.suptitle(tit)
    plt.subplots_adjust(top = 0.9)
    fig = ctl.custom_legend(fig, colorscoup+['grey'], model_coups+['avg'])

    fig.savefig(cart_out + cos+'_{}_1vs1.pdf'.format(vtag))
    figures.append(fig)


### mo ci sono i residtimes
cos = 'resid_times'
for mod in model_names_all:
    rmsall = []
    rmsall_90 = []
    for ke in results.keys():
        if mod == ke.split('_')[0]:
            rmsall.append([np.mean(results[ke][cos][reg]) for reg in range(4)])
            rmsall_90.append([np.percentile(results[ke][cos][reg], 90) for reg in range(4)])

    if len(rmsall) == 0:
        print('aFJAFJAJFAFJAJFAJFJFJAFJFAJ!!!!')
        continue

    nuvars[mod][cos+'_mean'] = np.mean(rmsall, axis = 0)
    nuvars[mod][cos+'_90'] = np.mean(rmsall_90, axis = 0)

    nuvars[mod][cos+'_mean_err'] = np.std(rmsall, axis = 0)
    nuvars[mod][cos+'_90_err'] = np.std(rmsall_90, axis = 0)


for cos, tit in zip(allnams[5:], longnams[5:]):
    for tipres, resolution in zip(['atm', 'oce'], [atm_resolution, oce_resolution]):
        fig = plt.figure(figsize = (16, 12))
        ind = 0
        axes = []
        for reg in range(4):
            ind += 1
            ax = plt.subplot(2, 2, ind)

            allvals = np.array([np.mean(nuvars[mod][cos][reg]) for mod in model_names])
            allerr = np.array([nuvars[mod][cos+'_err'][reg] for mod in model_names])

            for res, err, val, col, mark in zip(resolution, allerr, allvals, colors, sym_all):
                ax.scatter(res, val, color = col, marker = mark, s = 100, linewidth = lw)
                if not np.isnan(err):
                    ax.errorbar(res, val, yerr = err, color = col, linewidth = lw, capsize = 5)

            ax.axhline(nuvars['ERA'][cos][reg], color = 'grey', alpha = 0.6)

            ax.set_ylabel(tit)
            ax.set_xlabel('Eff. atm. resolution (km)')

            ax.axhline(0., color = 'grey', alpha = 0.6)

            #ax.set_title('Reg {}'.format(reg))
            ax.set_title(regnam[reg])
            axes.append(ax)
        ctl.adjust_ax_scale(axes)

        fig.suptitle(tit)
        plt.subplots_adjust(top = 0.9)
        fig = ctl.custom_legend(fig, colors, model_names, ncol = 6)
        fig.savefig(cart_out + cos+'_{}_{}resolution.pdf'.format(vtag, tipres))
        figures.append(fig)


    fig = plt.figure(figsize = (16, 12))
    ind = 0
    axes = []
    for reg in range(4):
        ind += 1
        ax = plt.subplot(2, 2, ind)

        eraval = nuvars['ERA'][cos][reg]
        allvals = np.array([np.mean(nuvars[mod][cos][reg]) for mod in model_names])
        allerr = np.array([nuvars[mod][cos+'_err'][reg] for mod in model_names])

        allvalsdiff = []
        allerrdiff = []
        for modcou in model_coups:
            v_hr = [val for val, mod, vv in zip(allvals, model_names, vers) if modcou in mod and vv == 'HR'][0]
            v_lr = [val for val, mod, vv in zip(allvals, model_names, vers) if modcou in mod and vv == 'LR'][0]
            #vdiff = abs(v_lr-eraval)-abs(v_hr-eraval)
            vdiff = v_hr-v_lr
            allvalsdiff.append(vdiff)

            v_err_hr = [val for val, mod, vv in zip(allerr, model_names, vers) if modcou in mod and vv == 'HR'][0]
            v_err_lr = [val for val, mod, vv in zip(allerr, model_names, vers) if modcou in mod and vv == 'LR'][0]
            errdiff = np.nanmax([v_err_hr, v_err_lr])
            allerrdiff.append(errdiff)

        i = 0
        for val, err, col, modcou in zip(allvalsdiff, allerrdiff, colorscoup, model_coups):
            if np.isnan(err):
                ax.bar(i, val, color = col, width = wi)
            else:
                ax.bar(i, val, color = col, width = wi, yerr = err, capsize = 5)
            i+=0.7

        i+=0.3
        ax.bar(i, np.mean(allvalsdiff), color = 'grey', width = wi, yerr = np.std(allvalsdiff), capsize = 5)
        ax.set_xticks([])
        ax.set_ylabel(tit + ' (difference: HR - LR)')
        #ax.set_xlabel('Eff. atm. resolution (km)')
        ax.axhline(0., color = 'grey', alpha = 0.6)

        ax.set_title(regnam[reg])
        axes.append(ax)
    ctl.adjust_ax_scale(axes)

    fig.suptitle(tit)
    plt.subplots_adjust(top = 0.9)
    fig = ctl.custom_legend(fig, colorscoup+['grey'], model_coups+['avg'])

    fig.savefig(cart_out + cos+'_{}_1vs1.pdf'.format(vtag))
    figures.append(fig)

ctl.plot_pdfpages(cart_out + 'allfigs_primacoup_{}_wresolution.pdf'.format(vtag), figures)
