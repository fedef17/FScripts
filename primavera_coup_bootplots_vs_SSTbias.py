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
from scipy import stats

from matplotlib.colors import LogNorm
from copy import deepcopy as cp

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['figure.max_open_warning'] = 0
#######################################
cart_out = '/home/fabiano/Research/articoli/Papers/primavera_regimes/figures/figures_bootstraps_v7_vs_biases/'
if not os.path.exists(cart_out): os.mkdir(cart_out)

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

colors_wERA = colors + [cm.colors.ColorConverter.to_rgb('black')]
colors_wERA = [tuple(col) for col in colors_wERA]
color_main.append('black')

cart_boot = '/home/fabiano/Research/articoli/Papers/primavera_regimes/figures/figures_bootstraps_v7/'
filon = open(cart_boot + 'res_bootstrap_v7.p', 'rb')

vtag = 'v7'
cart = '/home/fabiano/Research/lavori/WeatherRegimes/prima_coup_{}/'.format(vtag)
filogen = cart + 'out_prima_coup_{}_DJF_EAT_4clus_4pcs_1957-2014_refEOF.p'.format(vtag)
results, results_ref = pickle.load(open(filogen, 'rb'))

results.pop('HadGEM3-GC31-LL_r1i1p2f1')
results.pop('EC-Earth3P_r1i1p1f1')
results.pop('EC-Earth3P-HR_r1i1p1f1')
allresmembers = results.keys()
del results

cart_sstbias = '/home/fabiano/Research/articoli/Papers/primavera_regimes/figures/SST_biases/'
allrms, allpatcor = pickle.load(open(cart_sstbias + 'sst_bias_rms_djf_eat.p', 'rb'))
for mod in model_names:
    if mod == 'AWI-CM-1-1-LR':
        modok = 'AWI-CM-1-0-LR'
    elif mod == 'AWI-CM-1-1-HR':
        modok == 'AWI-CM-1-0-HR'
    else:
        modok = mod

    rmsall = []
    patcorall = []
    for ke in allpatcor.keys():
        if modok in ke:
            rmsall.append(allrms[ke])
            patcorall.append(allpatcor[ke])

    allrms[mod] = np.mean(rmsall)
    allpatcor[mod] = np.mean(patcorall)

sstbias_rms = [allrms[mod] for mod in model_names]
sstbias_pat = [allpatcor[mod] for mod in model_names]

cart_jet = '/home/fabiano/Research/articoli/Papers/primavera_regimes/figures/jet_and_blocking/'
relent = pickle.load(open(cart_jet + 'jetdist_relent.p', 'rb'))
jetbias = []
for mod in model_names:
    if (mod, 'all', '') in relent.keys():
        jetbias.append(relent[(mod, 'all', '')])
    else:
        jetbias.append(np.nan)


jetlat_q = pickle.load(open(cart_jet + 'jetlat_q_all.p', 'rb'))
jetmean_bias = []
for mod in model_names:
    if (mod, 'mean') in jetlat_q.keys():
        jetmean_bias.append(abs(jetlat_q[(mod, 'mean')] - jetlat_q[('ERA', 'mean')]))
    else:
        jetmean_bias.append(np.nan)

jetmeansign_bias = []
for mod in model_names:
    if (mod, 'mean') in jetlat_q.keys():
        jetmeansign_bias.append(jetlat_q[(mod, 'mean')] - jetlat_q[('ERA', 'mean')])
    else:
        jetmeansign_bias.append(np.nan)

jetmedian_bias = []
for mod in model_names:
    if (mod, 'median') in jetlat_q.keys():
        jetmedian_bias.append(abs(jetlat_q[(mod, 'median')] - jetlat_q[('ERA', 'median')]))
    else:
        jetmedian_bias.append(np.nan)

jetvariab = []
for mod in model_names:
    if (mod, 'q90') in jetlat_q.keys():
        jetvariab.append(abs(jetlat_q[(mod, 'q90')] - jetlat_q[(mod, 'q10')]))
    else:
        jetvariab.append(np.nan)

cart_block = '/home/fabiano/Research/articoli/Papers/primavera_regimes/figures/Reinhard_blocking/'
allmaps, allrms, allpatcor, allsums = pickle.load(open(cart_block + 'out_bloccomp_v7.p', 'rb'))
res_boot, res_boot_ERA, EAT_mean = pickle.load(open(cart_block + 'out_blocboot_v7.p', 'rb'))

for mod in model_names:
    rmsall = []
    patcorall = []
    meanall = []
    for ke in allpatcor.keys():
        if mod in ke and 'full' in ke:
            rmsall.append(allrms[ke])
            patcorall.append(allpatcor[ke])

    for ke in EAT_mean.keys():
        if mod in ke:
            meanall.append(EAT_mean[ke])

    allrms[mod] = np.mean(rmsall)
    allpatcor[mod] = np.mean(patcorall)
    EAT_mean[mod] = np.mean(meanall)

    if mod == 'ECMWF-IFS-MR':
        allrms[mod] = np.nan
        allpatcor[mod] = np.nan
        EAT_mean[mod] = np.nan

blockbias_rms = [allrms[mod] for mod in model_names]
blockbias_pat = [allpatcor[mod] for mod in model_names]
blockfreq = [EAT_mean[mod] for mod in model_names]

#### MEAN FIELD, LOW FR VAR, HIGH FR VAR, ECCECC

cart_lohi = '/home/fabiano/Research/lavori/PRIMAVERA_bias/hist_1950/'
mean_field_all, lowfrvar, highfrvar, stat_eddy_all = pickle.load(open(cart_lohi + 'out_lowhighstat_var.p', 'rb'))

mf_bias = []
lf_bias = []
hf_bias = []
se_bias = []

print(mean_field_all['ERA'].shape)
n_latsss = mean_field_all['ERA'].shape[0]
latsss = np.arange(30., 88, 2.5)

i = 0
for cos, cosbia in zip([mean_field_all, lowfrvar, highfrvar, stat_eddy_all], [mf_bias, lf_bias, hf_bias, se_bias]):
    i+=1
    for mod in model_names:
        rmsall = []
        for ke in cos.keys():
            if mod in ke:
                rmsall.append(cos[ke])

        allrms = np.array([ctl.E_rms(gigi, cos['ERA'], latitude = latsss) for gigi in rmsall])
        allpatcor = np.array([ctl.Rcorr(gigi, cos['ERA'], latitude = latsss) for gigi in rmsall])

        cos[mod] = np.mean(rmsall, axis = 0)
        rmsmed = ctl.E_rms(cos[mod], cos['ERA'], latitude = latsss)
        patmed = ctl.Rcorr(cos[mod], cos['ERA'], latitude = latsss)
        #cosbia.append(allrms.mean())
        cosbia.append(allpatcor.mean())

        #if abs(rmsmed-allrms.mean())/rmsmed > 0.1:
        #    print(i, mod, rmsmed, np.mean(allrms), np.min(allrms), np.max(allrms))
        #    print(i, mod, patmed, np.mean(allpatcor), np.min(allpatcor), np.max(allpatcor))


# atm_resolution = np.array([245, 94, 97, 23, 250, 44, 100, 46, 48, 50, 27, 103, 52, 255, 106, 54, 56])
# oce_resolution = np.array([50, 29, 21, 22, 96, 23, 98, 24, 102, 25, 26, 38, 42, 104, 27, 28, 8])
atm_resolution = np.array([250, 100, 100, 25, 250, 50, 100, 50, 50, 50, 25, 100, 50, 250, 100, 50, 50])
oce_resolution = np.array([50, 25, 25, 25, 100, 25, 100, 25, 100, 25, 25, 40, 40, 100, 25, 25, 8])

# modelspec = np.array([0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 6, 6, 6, 6])


allbiasnam = ['_sstrms', '_sstpat', '_jet', '_blockrms', '_blockpat', '_blockfreq', '_atmres', '_oceres', '_meanfield', '_lowfrvar', '_hifrvar', '_stateddy', '_jetmean', '_jetmedian', '_jetmeansign', '_jetvariab']
allbias = [sstbias_rms, sstbias_pat, jetbias, blockbias_rms, blockbias_pat, blockfreq, atm_resolution, oce_resolution, mf_bias, lf_bias, hf_bias, se_bias, jetmean_bias, jetmedian_bias, jetmeansign_bias, jetvariab]


#############################################################################################
#############################################################################################

nmods = len(model_names_all)

allnams = ['significance', 'varopt', 'autocorr', 'freq', 'dist_cen', 'resid_times_av', 'resid_times_90', 'trans_matrix', 'centroids', 'relative_entropy', 'patcor', 'filt_dist_cen', 'filt_relative_entropy', 'filt_patcor']
alltits = ['Sharpness of regime structure', 'Optimal variance ratio', 'lag-1 autocorrelation', 'Regime frequency', 'Centroid distance to reference in phase space', 'Average residence time', '90th percentile of residence time', 'Transition matrix', 'Regime centroid', 'Relative entropy of cluster cloud in phase space', 'Pattern correlation of cluster centroid in phase space', 'Centroid distance to reference in phase space (5-day filter)', 'Relative entropy of cluster cloud in phase space (5-day filter)', 'Pattern correlation of cluster centroid in phase space (5-day filter)']
alltits = dict(zip(allnams, alltits))

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
    allkeysss = bootstraps_all[allmems[0]].keys()

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

        if bootstraps_all[allmems[0]][ke].ndim == 1:
            for cos, num in zip(['p10', 'p25', 'p50', 'p75', 'p90'], [10, 25, 50, 75, 90]):
                bootstraps_all['boot_'+cos][ke] = np.percentile(np.concatenate([bootstraps_all[mem][ke] for mem in allmems]), num)
        elif bootstraps_all[allmems[0]][ke].ndim == 2:
            for cos, num in zip(['p10', 'p25', 'p50', 'p75', 'p90'], [10, 25, 50, 75, 90]):
                bootstraps_all['boot_'+cos][ke] = np.array([np.percentile(np.concatenate([bootstraps_all[mem][ke][:, reg] for mem in allmems]), num) for reg in range(4)])
        else:
            for cos, num in zip(['p10', 'p25', 'p50', 'p75', 'p90'], [10, 25, 50, 75, 90]):
                bootstraps_all['boot_'+cos][ke] = np.zeros((4,4))
                for i in range(4):
                    for j in range(4):
                        bootstraps_all['boot_'+cos][ke][i, j] = np.percentile(np.concatenate([bootstraps_all[mem][ke][:, i, j] for mem in allmems]), num)

    bootstri[mod] = bootstraps_all

for mod in model_names_all:
    for cos in ['boot_mean', 'boot_std', 'ens_mean', 'ens_min', 'ens_max', 'ens_std', 'boot_p10', 'boot_p25', 'boot_p50', 'boot_p75', 'boot_p90']:
        bootstri[mod][cos]['mean_dist_cen'] = np.mean(bootstri[mod][cos]['dist_cen'])
        bootstri[mod][cos]['mean_patcor'] = np.mean(bootstri[mod][cos]['patcor'])
        bootstri[mod][cos]['mean_freq_bias'] = np.sqrt(np.sum(np.array([bootstri[mod][cos]['freq'][i]-bootstri['ERA'][cos]['freq'][i] for i in range(4)])**2)/4)


newnam = ['mean_dist_cen', 'mean_patcor', 'mean_freq_bias']
allnams += newnam
newtit = ['Mean distance to reference centroids', 'Mean regime pattern correlation', 'Mean frequency bias']
gigi = dict(zip(newnam, newtit))
alltits.update(gigi)

# colors = ctl.color_set(len(model_names), sns_palette = 'Paired') + [cm.colors.ColorConverter.to_rgb('black')]
nam = 'significance'
regnam = ['NAO +', 'Sc. BL', 'AR', 'NAO -']

#plt.style.use('seaborn-whitegrid')
#plt.style.use('default')
plt.rc('axes', axisbelow=True)
#plt.ion()

ens_names = ['LR', 'HR']
ens_colors = ['teal', 'indianred']

rcorrsall = dict()
pearsall = dict()
pears_nocmcc = dict()
pears_nooutli = dict()

biasdict = dict()
biasdict['model_names'] = model_names

for biasnam, bias in zip(allbiasnam, allbias):
    print(biasnam)
    allfigs = []
    bias = np.array(bias)
    biasdict[biasnam] = bias

    for nam in ['significance', 'varopt', 'autocorr', 'mean_dist_cen', 'mean_patcor', 'mean_freq_bias']:
        fig, ax = plt.subplots(figsize=(16,12))

        allpercs = dict()
        for cos in ['mean', 'std', 'p10', 'p25', 'p50', 'p75', 'p90']:
            allpercs[cos] = np.array([bootstri[mod]['boot_'+cos][nam] for mod in model_names_all])
        for cos in ['ens_min', 'ens_max', 'ens_std']:
            allpercs[cos] = np.array([bootstri[mod][cos][nam] for mod in model_names_all])

        #ctl.primavera_boxplot_on_ax(ax, allpercs, model_names_all, colors_wERA, color_main, vers, ens_names, ens_colors, x = bias)
        ax.scatter(bias, allpercs['mean'][:-1], c = colors_wERA[:-1], s = 20)
        yerr1 = np.array(allpercs['mean']) - np.array(allpercs['p25'])
        yerr2 = -np.array(allpercs['mean']) + np.array(allpercs['p75'])
        yerr = np.stack([yerr1, yerr2])[:, :-1]
        for bi, mea, ye, col in zip(bias, allpercs['mean'][:-1], yerr.T, colors_wERA[:-1]):
            ax.errorbar(bi, mea, yerr = np.vstack([ye]).T, c = col, capsize = 5, linewidth = 2)

        #filt = (bias != 0.)
        filt = ~np.isnan(bias)
        rcorr = ctl.Rcorr(bias[filt], allpercs['mean'][:-1][filt])
        pears, pval = stats.pearsonr(bias[filt], allpercs['mean'][:-1][filt])
        rcorrsall[(biasnam, nam)] = rcorr
        pearsall[(biasnam, nam)] = (pears, pval)

        filt2 = cp(filt)
        filt2[2] = False
        filt2[3] = False
        pears, pval = stats.pearsonr(bias[filt2], allpercs['mean'][:-1][filt2])
        pears_nocmcc[(biasnam, nam)] = (pears, pval)

        biasoutli = abs(bias-np.nanmean(bias)) > 2*np.nanstd(bias)
        metroutli = abs(allpercs['mean'][:-1]-np.nanmean(allpercs['mean'][:-1])) > 2*np.nanstd(allpercs['mean'][:-1])
        if 'res' not in biasnam:
            filt3 = (filt) & ~(biasoutli) & ~(metroutli)
        else:
            filt3 = (filt) & ~(metroutli)
        print(nam, np.sum(filt3), np.sum(filt), filt3[2:4])
        pears, pval = stats.pearsonr(bias[filt3], allpercs['mean'][:-1][filt3])
        pears_nooutli[(biasnam, nam)] = (pears, pval)

        plt.title(alltits[nam])
        ctl.custom_legend(fig, colors_wERA + ens_colors, model_names_all + ens_names, ncol = 6)
        fig.savefig(cart_out + '{}_bootstraps_v7{}.pdf'.format(nam, biasnam))
        allfigs.append(fig)


    for nam in ['freq', 'dist_cen', 'resid_times_av', 'resid_times_90', 'relative_entropy', 'patcor', 'filt_relative_entropy', 'filt_patcor']:
        fig = plt.figure(figsize=(16,12))
        axes = []
        for ii, reg in enumerate(regnam):
            ax = plt.subplot(2, 2, ii+1)

            allpercs = dict()
            for cos in ['mean', 'std', 'p10', 'p25', 'p50', 'p75', 'p90']:
                allpercs[cos] = np.array([bootstri[mod]['boot_'+cos][nam][ii] for mod in model_names_all])
            for cos in ['ens_min', 'ens_max', 'ens_std']:
                allpercs[cos] = np.array([bootstri[mod][cos][nam][ii] for mod in model_names_all])

            #ctl.primavera_boxplot_on_ax(ax, allpercs, model_names_all, colors_wERA, color_main, vers, ens_names, ens_colors, wi = 0.5, x = bias)
            ax.scatter(bias, allpercs['mean'][:-1], c = colors_wERA[:-1], s = 20)
            yerr1 = np.array(allpercs['mean']) - np.array(allpercs['p25'])
            yerr2 = -np.array(allpercs['mean']) + np.array(allpercs['p75'])
            yerr = np.stack([yerr1, yerr2])[:, :-1]
            for bi, mea, ye, col in zip(bias, allpercs['mean'][:-1], yerr.T, colors_wERA[:-1]):
                ax.errorbar(bi, mea, yerr = np.vstack([ye]).T, c = col, capsize = 5)
            #ax.errorbar(bias, allpercs['mean'][:-1], yerr = yerr, c = colors_wERA[:-1], capsize = 5)

            ax.set_title(reg)
            axes.append(ax)

            #filt = (bias != 0.)
            filt = ~np.isnan(bias)
            rcorr = ctl.Rcorr(bias[filt], allpercs['mean'][:-1][filt])
            rcorrsall[(biasnam, nam, reg)] = rcorr
            pears, pval = stats.pearsonr(bias[filt], allpercs['mean'][:-1][filt])
            pearsall[(biasnam, nam, reg)] = (pears, pval)

            filt2 = filt
            filt[2] = False
            filt[3] = False
            pears, pval = stats.pearsonr(bias[filt], allpercs['mean'][:-1][filt])
            pears_nocmcc[(biasnam, nam, reg)] = (pears, pval)

        ctl.adjust_ax_scale(axes)

        fig.suptitle(alltits[nam])
        plt.subplots_adjust(top = 0.9)

        ctl.custom_legend(fig, colors_wERA + ens_colors, model_names_all + ens_names, ncol = 6)
        fig.savefig(cart_out + '{}_bootstraps_v7{}.pdf'.format(nam, biasnam))
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
                allpercs[cos] = np.array([bootstri[mod]['boot_'+cos][nam][ireg, ii] for mod in model_names_all])
            for cos in ['ens_min', 'ens_max', 'ens_std']:
                allpercs[cos] = np.array([bootstri[mod][cos][nam][ireg, ii] for mod in model_names_all])

            #ctl.primavera_boxplot_on_ax(ax, allpercs, model_names_all, colors_wERA, color_main, vers, ens_names, ens_colors, wi = 0.5)
            ax.scatter(bias, allpercs['mean'][:-1], c = colors_wERA[:-1], s = 20)
            yerr1 = np.array(allpercs['mean']) - np.array(allpercs['p25'])
            yerr2 = -np.array(allpercs['mean']) + np.array(allpercs['p75'])
            yerr = np.stack([yerr1, yerr2])[:, :-1]
            for bi, mea, ye, col in zip(bias, allpercs['mean'][:-1], yerr.T, colors_wERA[:-1]):
                ax.errorbar(bi, mea, yerr = np.vstack([ye]).T, c = col, capsize = 5)
            #ax.errorbar(bias, allpercs['mean'][:-1], yerr = yerr, c = colors_wERA[:-1], capsize = 5)
            if ireg == ii:
                #ctl.primavera_boxplot_on_ax(ax_pers_2[ii], allpercs, model_names_all, colors_wERA, color_main, vers, ens_names, ens_colors, wi = 0.5, x = bias)
                ax_pers_2[ii].scatter(bias, allpercs['mean'][:-1], c = colors_wERA[:-1], s = 20)
                yerr1 = np.array(allpercs['mean']) - np.array(allpercs['p25'])
                yerr2 = -np.array(allpercs['mean']) + np.array(allpercs['p75'])
                yerr = np.stack([yerr1, yerr2])[:, :-1]
                for bi, mea, ye, col in zip(bias, allpercs['mean'][:-1], yerr.T, colors_wERA[:-1]):
                    ax_pers_2[ii].errorbar(bi, mea, yerr = np.vstack([ye]).T, c = col, capsize = 5)
                #ax_pers_2[ii].errorbar(bias, allpercs['mean'][:-1], yerr = yerr, c = colors_wERA[:-1], capsize = 5)

            ax.set_title('trans {} -> {}'.format(regnam[ireg], regnam[ii]))
            if ii != ireg:
                axes_diff.append(ax)
            else:
                axes_pers.append(ax)
                axes_pers.append(ax_pers_2[ii])

        fig.suptitle(alltits[nam])
        plt.subplots_adjust(top = 0.9)

        ctl.custom_legend(fig, colors_wERA + ens_colors, model_names_all + ens_names, ncol = 6)
        figs.append(fig)
        allfigs.append(fig)

    ctl.adjust_ax_scale(axes_diff)
    ctl.adjust_ax_scale(axes_pers)

    for fig, ireg in zip(figs, range(4)):
        fig.savefig(cart_out + '{}_reg{}_bootstraps_v7{}.pdf'.format(nam, ireg, biasnam))


    fig_pers.suptitle('Regime persistence probability')
    plt.subplots_adjust(top = 0.9)

    ctl.custom_legend(fig_pers, colors_wERA + ens_colors, model_names_all + ens_names, ncol = 6)
    fig_pers.savefig(cart_out + 'Persistence_bootstraps_v7{}.pdf'.format(biasnam))
    allfigs.append(fig_pers)


    ctl.plot_pdfpages(cart_out + 'all_bootplots{}.pdf'.format(biasnam), allfigs, save_single_figs = False)


pickle.dump(biasdict, open(cart_out + 'allbiases.p', 'wb'))

print('Higher correlations:\n')
for ke in rcorrsall:
    if abs(rcorrsall[ke]) > 0.6:
        print(ke, rcorrsall[ke])

pickle.dump(rcorrsall, open(cart_out + 'rcorrsall.p', 'wb'))

fi = open(cart_out + 'all_correlations.txt', 'w')
fi.write('All correlations: \n')
for ke in rcorrsall:
    fi.write('{:50s}  {:7.3f} \n'.format(str(ke), rcorrsall[ke]))
fi.write('\n Higher correlations: \n')
for ke in rcorrsall:
    if abs(rcorrsall[ke]) > 0.4:
        fi.write('{:50s}  {:7.3f} \n'.format(str(ke), rcorrsall[ke]))
fi.close()


fi = open(cart_out + 'all_pearsonr.txt', 'w')
fi.write('All correlations: \n')
for ke in pearsall:
    fi.write('{:50s}  {:7.3f}  {:10.3e} \n'.format(str(ke), pearsall[ke][0], pearsall[ke][1]))
fi.write('\n Higher correlations: \n')
for ke in pearsall:
    if abs(pearsall[ke][0]) > 0.4:
        fi.write('{:50s}  {:7.3f}  {:10.3e} \n'.format(str(ke), pearsall[ke][0], pearsall[ke][1]))
fi.close()


fi = open(cart_out + 'all_pearsonr_nocmcc.txt', 'w')
fi.write('All correlations: \n')
for ke in pears_nocmcc:
    fi.write('{:50s}  {:7.3f}  {:10.3e} \n'.format(str(ke), pears_nocmcc[ke][0], pears_nocmcc[ke][1]))
fi.write('\n Higher correlations: \n')
for ke in pears_nocmcc:
    if abs(pears_nocmcc[ke][0]) > 0.4:
        fi.write('{:50s}  {:7.3f}  {:10.3e} \n'.format(str(ke), pears_nocmcc[ke][0], pears_nocmcc[ke][1]))
fi.close()


fi = open(cart_out + 'best_pearsonr.txt', 'w')
fi.write('Best correlations: \n')
for ke in pearsall:
    if np.any([cos in ke for cos in ['varopt', 'significance', 'autocorr', 'mean_dist_cen', 'mean_patcor', 'mean_freq_bias']]):
        fi.write('{:50s}  {:7.3f}  {:10.3e}   {:7.3f}  {:10.3e} \n'.format(str(ke), pearsall[ke][0], pearsall[ke][1], pears_nocmcc[ke][0], pears_nocmcc[ke][1]))
fi.close()


fi = open(cart_out + 'best_pearsonr_nooutli.txt', 'w')
fi.write('Best correlations: \n')
for ke in pearsall:
    if np.any([cos in ke for cos in ['varopt', 'significance', 'autocorr', 'mean_dist_cen', 'mean_patcor', 'mean_freq_bias']]):
        fi.write('{:50s}  {:7.3f}  {:10.3e}   {:7.3f}  {:10.3e}   {:7.3f}  {:10.3e}\n'.format(str(ke), pearsall[ke][0], pearsall[ke][1], pears_nocmcc[ke][0], pears_nocmcc[ke][1], pears_nooutli[ke][0], pears_nooutli[ke][1]))
fi.close()
