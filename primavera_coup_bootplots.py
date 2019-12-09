#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
from matplotlib import pyplot as plt
from matplotlib import cm
import pickle

import climtools_lib as ctl
#import climdiags as cd

from importlib import reload
from matplotlib.colors import LogNorm

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 15
#######################################
cart_in = '/home/fabiano/Research/articoli/Papers/primavera_regimes/figures/figures_bootstraps_v7/'

plot_sig = True
plot_mean = True

if plot_sig:
    cart_out = cart_in + 'v7_rmshi_500/'
    if not os.path.exists(cart_out): os.mkdir(cart_out)
    filon = open(cart_in + 'res_bootstrap_v7_rmshi_500.p', 'rb')
else:
    cart_out = cart_in + 'v7_final/'
    if not os.path.exists(cart_out): os.mkdir(cart_out)
    filon = open(cart_in + 'res_bootstrap_v7_500_relent2_nosig.p', 'rb')

vtag = 'v7'
cart = '/home/fabiano/Research/lavori/WeatherRegimes/prima_coup_{}/'.format(vtag)
filogen = cart + 'out_prima_coup_{}_DJF_EAT_4clus_4pcs_1957-2014_refEOF.p'.format(vtag)
results, results_ref = pickle.load(open(filogen, 'rb'))

results.pop('HadGEM3-GC31-LL_r1i1p2f1')
results.pop('EC-Earth3P_r1i1p1f1')
results.pop('EC-Earth3P-HR_r1i1p1f1')
allresmembers = list(results.keys())
del results

##############################################################

model_names = ['AWI-CM-1-1-LR', 'AWI-CM-1-1-HR', 'CMCC-CM2-HR4', 'CMCC-CM2-VHR4', 'CNRM-CM6-1', 'CNRM-CM6-1-HR', 'EC-Earth3P', 'EC-Earth3P-HR', 'ECMWF-IFS-LR', 'ECMWF-IFS-MR', 'ECMWF-IFS-HR', 'MPI-ESM1-2-HR', 'MPI-ESM1-2-XR', 'HadGEM3-GC31-LL', 'HadGEM3-GC31-MM', 'HadGEM3-GC31-HM', 'HadGEM3-GC31-HH']

vers = ['LR', 'HR', 'LR', 'HR', 'LR', 'HR', 'LR', 'HR', 'LR', 'MR', 'HR', 'LR', 'HR', 'LR', 'MR', 'MR', 'HR'] + ['OBS']
model_coups = ['AWI-CM-1-1', 'CMCC-CM2', 'CNRM-CM6-1', 'EC-Earth3P', 'ECMWF-IFS', 'MPI-ESM1-2', 'HadGEM3-GC31']

model_names_all = model_names + ['ERA']

colors = ctl.color_set(len(model_coups)*2, sns_palette = 'Paired')
colorscoup = [np.mean([col1, col2], axis = 0) for col1, col2 in zip(colors[:-1:2], colors[1::2])]
color_main = []
for col2 in colors[1::2]:
    color_main += [col2, col2]

mr_cos = np.where(np.array(vers) == 'MR')[0]
for gi in mr_cos:
    colors.insert(gi, np.mean(colors[gi-1:gi+1], axis = 0))
    color_main.insert(gi, colors[gi+1])

colors_wERA = colors + ['black']
color_main.append('black')

#############################################################################################
#############################################################################################

nmods = len(model_names_all)

allnams = ['significance', 'varopt', 'autocorr', 'freq', 'dist_cen', 'resid_times_av', 'resid_times_90', 'trans_matrix', 'centroids', 'relative_entropy', 'patcor', 'filt_dist_cen', 'filt_relative_entropy', 'filt_patcor']
alltits = ['Sharpness of regime structure', 'Optimal variance ratio', 'lag-1 autocorrelation', 'Regime frequency', 'Centroid distance to reference in phase space', 'Average residence time', '90th percentile of residence time', 'Transition matrix', 'Regime centroid', 'Relative entropy of cluster cloud in phase space', 'Pattern correlation of cluster centroid in phase space', 'Centroid distance to reference in phase space (5-day filter)', 'Relative entropy of cluster cloud in phase space (5-day filter)', 'Pattern correlation of cluster centroid in phase space (5-day filter)']
allylabels = ['sharpness', 'optimal ratio', 'lag-1 autocorrelation', 'frequency', 'centroid distance (m)', 'resid. time (days)', 'resid. time (days)', 'trans. probability', 'centroid (m)', 'rel. entropy', 'patt. corr.', 'centroid distance (m)', 'rel. entropy', 'patt. corr.']

alltits = dict(zip(allnams, alltits))
allylabels = dict(zip(allnams, allylabels))

bootstri = dict()
for mod in model_names_all:
    print(mod)
    allmems = []
    for ke in allresmembers:
        if mod == ke.split('_')[0]:
            allmems.append(ke.split('_')[1])

    print(allmems)
    if mod == 'ERA':
        allmems = ['0']

    bootstraps_all = pickle.load(filon)
    # bootstraps_all = dict()
    # for mem in allmems:
    #     bootstraps_all[mem] = pickle.load(filon)

    #bootstraps_all = pickle.load(filon) # new version
    #allmems = list(bootstraps_all.keys())
    allkeysss = list(bootstraps_all[allmems[0]].keys())

    # list of all boot means
    for cos in ['mean', 'std', 'p10', 'p25', 'p50', 'p75', 'p90']:
        bootstraps_all['boot_'+cos] = dict()
    # mean of all boot means
    for cos in ['mean', 'std', 'min', 'max']:
        bootstraps_all['ens_'+cos] = dict()

    for ke in allkeysss:
        for mem in allmems:
            bootstraps_all[mem][ke] = np.array(bootstraps_all[mem][ke])

        #print(mem, ke, bootstraps_all[mem][ke].shape, allmems, bootstraps_all['boot_mean'].keys())
        bootstraps_all['boot_mean'][ke] = np.array([np.mean(bootstraps_all[mem][ke], axis = 0) for mem in allmems])
        bootstraps_all['ens_mean'][ke] = np.mean(bootstraps_all['boot_mean'][ke], axis = 0)
        bootstraps_all['ens_min'][ke] = np.min(bootstraps_all['boot_mean'][ke], axis = 0)
        bootstraps_all['ens_max'][ke] = np.max(bootstraps_all['boot_mean'][ke], axis = 0)
        bootstraps_all['ens_std'][ke] = np.std(bootstraps_all['boot_mean'][ke], axis = 0)

        bootstraps_all['boot_mean'][ke] = np.mean(bootstraps_all['boot_mean'][ke], axis = 0)
        bootstraps_all['boot_std'][ke] = np.std(np.concatenate([bootstraps_all[mem][ke] for mem in allmems]), axis = 0)

        ndim = bootstraps_all[allmems[0]][ke].ndim
        if ndim == 1:
            for cos, num in zip(['p10', 'p25', 'p50', 'p75', 'p90'], [10, 25, 50, 75, 90]):
                bootstraps_all['boot_'+cos][ke] = np.percentile(np.concatenate([bootstraps_all[mem][ke] for mem in allmems]), num)

            # bootstraps_all['ens_min'][ke] = np.min([np.percentile(bootstraps_all[mem][ke], 50) for mem in allmems])
            # bootstraps_all['ens_max'][ke] = np.max([np.percentile(bootstraps_all[mem][ke], 50) for mem in allmems])
            #print(ndim, type(bootstraps_all['ens_min'][ke]))
        elif ndim == 2:
            for cos, num in zip(['p10', 'p25', 'p50', 'p75', 'p90'], [10, 25, 50, 75, 90]):
                bootstraps_all['boot_'+cos][ke] = np.array([np.percentile(np.concatenate([bootstraps_all[mem][ke][:, reg] for mem in allmems]), num) for reg in range(4)])
            # bootstraps_all['ens_min'][ke] = np.array([np.min([np.percentile(bootstraps_all[mem][ke][:,reg], 50) for mem in allmems]) for reg in range(4)])
            # bootstraps_all['ens_max'][ke] = np.array([np.max([np.percentile(bootstraps_all[mem][ke][:,reg], 50) for mem in allmems]) for reg in range(4)])
            #print(ndim, type(bootstraps_all['ens_min'][ke]))
        else:
            for cos, num in zip(['p10', 'p25', 'p50', 'p75', 'p90'], [10, 25, 50, 75, 90]):
                bootstraps_all['boot_'+cos][ke] = np.zeros((4,4))
                for i in range(4):
                    for j in range(4):
                        bootstraps_all['boot_'+cos][ke][i, j] = np.percentile(np.concatenate([bootstraps_all[mem][ke][:, i, j] for mem in allmems]), num)

                        # if cos == 'p50':
                        #     bootstraps_all['ens_min'][ke][i,j] = np.min([np.percentile(bootstraps_all[mem][ke][:, i, j], 50) for mem in allmems])
                        #     bootstraps_all['ens_max'][ke][i,j] = np.max([np.percentile(bootstraps_all[mem][ke][:, i, j], 50) for mem in allmems])
            #print(ndim, type(bootstraps_all['ens_min'][ke]))

    bootstri[mod] = bootstraps_all


# colors = ctl.color_set(len(model_names), sns_palette = 'Paired') + [cm.colors.ColorConverter.to_rgb('black')]
nam = 'significance'
regnam = ['NAO +', 'Sc. BL', 'AR', 'NAO -']

#plt.style.use('seaborn-whitegrid')
#plt.style.use('default')
plt.rc('axes', axisbelow=True)
#plt.ion()

ens_names = ['LR', 'HR']
ens_colors = ['teal', 'indianred']

allfigs = []
if plot_sig:
    for nam in ['significance', 'varopt', 'autocorr']:
        fig, ax = plt.subplots(figsize=(16,12))

        allpercs = dict()
        for cos in ['mean', 'std', 'p10', 'p25', 'p50', 'p75', 'p90']:
            allpercs[cos] = [bootstri[mod]['boot_'+cos][nam] for mod in model_names_all]
        for cos in ['ens_min', 'ens_max', 'ens_std']:
            allpercs[cos] = [bootstri[mod][cos][nam] for mod in model_names_all]

        ctl.primavera_boxplot_on_ax(ax, allpercs, model_names_all, colors_wERA, color_main, vers, ens_names, ens_colors, plot_mean = plot_mean)
        ax.set_ylabel(allylabels[nam])

        #plt.title(alltits[nam])
        #fig.suptitle(alltits[nam], fontsize = titlefont)
        ctl.custom_legend(fig, colors_wERA + ens_colors, model_names_all + ens_names, ncol = 6)
        fig.savefig(cart_out + '{}_bootstraps_v7.pdf'.format(nam))
        allfigs.append(fig)


for nam in ['freq', 'dist_cen', 'resid_times_av', 'resid_times_90', 'relative_entropy', 'patcor', 'filt_dist_cen', 'filt_relative_entropy', 'filt_patcor']:
    fig = plt.figure(figsize=(16,12))
    axes = []
    for ii, reg in enumerate(regnam):
        ax = plt.subplot(2, 2, ii+1)

        allpercs = dict()
        for cos in ['mean', 'std', 'p10', 'p25', 'p50', 'p75', 'p90']:
            allpercs[cos] = [bootstri[mod]['boot_'+cos][nam][ii] for mod in model_names_all]
        for cos in ['ens_min', 'ens_max', 'ens_std']:
            allpercs[cos] = [bootstri[mod][cos][nam][ii] for mod in model_names_all]

        ctl.primavera_boxplot_on_ax(ax, allpercs, model_names_all, colors_wERA, color_main, vers, ens_names, ens_colors, wi = 0.5, plot_mean = plot_mean)

        ax.set_title(reg)
        ax.set_ylabel(allylabels[nam])
        axes.append(ax)

    ctl.adjust_ax_scale(axes)

    #fig.suptitle(alltits[nam], fontsize = titlefont)
    plt.subplots_adjust(top = 0.9)

    ctl.custom_legend(fig, colors_wERA + ens_colors, model_names_all + ens_names, ncol = 6)
    fig.savefig(cart_out + '{}_bootstraps_v7.pdf'.format(nam))
    allfigs.append(fig)


ax_pers_2 = []
fig_pers = plt.figure(figsize=(16,12))
for i in range(4):
    ax_pers_2.append(plt.subplot(2,2,i+1))

nam = 'trans_matrix'
figs = []
axes_diff = []
axes_pers = []
for ireg, reg in enumerate(regnam):
    fig = plt.figure(figsize=(16,12))
    for ii in range(4):
        ax = plt.subplot(2, 2, ii+1)

        allpercs = dict()
        for cos in ['mean', 'std', 'p10', 'p25', 'p50', 'p75', 'p90']:
            allpercs[cos] = [bootstri[mod]['boot_'+cos][nam][ireg, ii] for mod in model_names_all]
        for cos in ['ens_min', 'ens_max', 'ens_std']:
            allpercs[cos] = [bootstri[mod][cos][nam][ireg, ii] for mod in model_names_all]

        ctl.primavera_boxplot_on_ax(ax, allpercs, model_names_all, colors_wERA, color_main, vers, ens_names, ens_colors, wi = 0.5, plot_mean = plot_mean)
        if ireg == ii:
            ctl.primavera_boxplot_on_ax(ax_pers_2[ii], allpercs, model_names_all, colors_wERA, color_main, vers, ens_names, ens_colors, wi = 0.5, plot_mean = plot_mean)

        ax.set_title('trans {} -> {}'.format(regnam[ireg], regnam[ii]))
        ax.set_ylabel(allylabels[nam])
        if ii != ireg:
            axes_diff.append(ax)
        else:
            axes_pers.append(ax)
            axes_pers.append(ax_pers_2[ii])

    #fig.suptitle(alltits[nam], fontsize = titlefont)
    plt.subplots_adjust(top = 0.9)

    ctl.custom_legend(fig, colors_wERA + ens_colors, model_names_all + ens_names, ncol = 6)
    figs.append(fig)
    allfigs.append(fig)

ctl.adjust_ax_scale(axes_diff)
ctl.adjust_ax_scale(axes_pers)

for fig, ireg in zip(figs, range(4)):
    fig.savefig(cart_out + '{}_reg{}_bootstraps_v7.pdf'.format(nam, ireg))


#fig_pers.suptitle('Regime persistence probability', fontsize = titlefont)
plt.subplots_adjust(top = 0.9)

ctl.custom_legend(fig_pers, colors_wERA + ens_colors, model_names_all + ens_names, ncol = 6)
fig_pers.savefig(cart_out + 'Persistence_bootstraps_v7.pdf')
allfigs.append(fig_pers)

#
#
# nam = 'trans_matrix'
# figs = []
# axes_diff = []
# axes_pers = []
# for ireg, reg in enumerate(regnam):
#     fig = plt.figure(figsize=(24,8))
#     ind = 0
#     normpers = [np.array([(1-cos[ireg, ireg]) for cos in bootstri[mod][nam]]) for mod in model_names_all]
#     normpers = dict(zip(model_names_all, normpers))
#     for ii in range(4):
#         if ii == ireg: continue
#         ind += 1
#         ax = plt.subplot(1, 3, ind)
#         allsigs = [np.mean(np.array([cos[ireg, ii] for cos in bootstri[mod][nam]])/normpers[mod]) for mod in model_names_all]
#         allerrs = [np.std(np.array([cos[ireg, ii] for cos in bootstri[mod][nam]])/normpers[mod]) for mod in model_names_all]
#         allp90 = [(np.percentile(np.array([cos[ireg, ii] for cos in bootstri[mod][nam]])/normpers[mod], 10), np.percentile(np.array([cos[ireg, ii] for cos in bootstri[mod][nam]])/normpers[mod], 90)) for mod in model_names_all]
#         plt.scatter(list(range(nmods)), allsigs, c = colors, s = 20)
#         for imod, sig, err, errp9, col in zip(range(nmods), allsigs, allerrs, allp90, colors):
#             #plt.errorbar(imod, sig, yerr = err, ecolor = col, linestyle = 'None', elinewidth = 2, capsize = 10)
#             plt.errorbar([imod], [sig], yerr = [[sig - errp9[0]], [errp9[1]-sig]], ecolor = col, linestyle = 'None', elinewidth = 2, capsize = 5)
#         ax.set_title('trans {} -> {}'.format(regnam[ireg], regnam[ii]))
#         axes_diff.append(ax)
#
#         ax.axhline(allsigs[-1], color = 'grey', alpha = 0.6)
#         ax.axhline(allp90[-1][0], color = 'grey', alpha = 0.6, ls = '--')
#         ax.axhline(allp90[-1][1], color = 'grey', alpha = 0.6, ls = '--')
#
#     fig.suptitle(alltits[nam], fontsize = titlefont)
#     plt.subplots_adjust(top = 0.9)
#
#     ctl.custom_legend(fig, colors, model_names_all)
#     figs.append(fig)
#     allfigs.append(fig)
#
# ctl.adjust_ax_scale(axes_diff)
#
# for fig, ireg in zip(figs, range(4)):
#     fig.savefig(cart_out + 'transonly_reg{}_bootstraps_v7.pdf'.format(ireg))


ctl.plot_pdfpages(cart_out + 'all_bootplots.pdf', allfigs, save_single_figs = False)
