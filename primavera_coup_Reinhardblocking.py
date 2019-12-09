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
import pandas as pd

plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 14
#######################################
cart_out = '/home/fabiano/Research/articoli/Papers/primavera_regimes/figures/Reinhard_blocking/'
if not os.path.exists(cart_out): os.mkdir(cart_out)

cart_out_maps = cart_out + 'maps/'
if not os.path.exists(cart_out_maps): os.mkdir(cart_out_maps)

cart = '/home/fabiano/Research/lavori/WeatherRegimes/prima_coup_v6/'
filogen = cart + 'out_prima_coup_v6_DJF_EAT_4clus_4pcs_1957-2014_refEOF.p'

cart_bloc = '/data-hobbes/fabiano/PRIMAVERA/Reinhard_blocking/all_hist-1950/'
all_fils_bloc = os.listdir(cart_bloc)
#fil_bloc = 'bf.5day.daily.daily_mean.{}.{}.hist-1950.r1i1p1f1.v20170915.{}-{}.nc'

model_names = ['AWI-CM-1-1-LR', 'AWI-CM-1-1-HR', 'CMCC-CM2-HR4', 'CMCC-CM2-VHR4', 'CNRM-CM6-1', 'CNRM-CM6-1-HR', 'EC-Earth3P', 'EC-Earth3P-HR', 'ECMWF-IFS-LR', 'ECMWF-IFS-MR', 'ECMWF-IFS-HR', 'MPI-ESM1-2-HR', 'MPI-ESM1-2-XR', 'HadGEM3-GC31-LL', 'HadGEM3-GC31-MM', 'HadGEM3-GC31-HM', 'HadGEM3-GC31-HH']
#model_names = ['CMCC-CM2-HR4', 'CMCC-CM2-VHR4', 'CNRM-CM6-1', 'CNRM-CM6-1-HR', 'EC-Earth3P', 'EC-Earth3P-HR', 'ECMWF-IFS-LR', 'ECMWF-IFS-HR', 'MPI-ESM1-2-HR', 'MPI-ESM1-2-XR', 'HadGEM3-GC31-LL', 'HadGEM3-GC31-MM']#, 'HadGEM3-GC31-HM']
#ens_mems = ['r1i1p1f002', 'r1i1p1f002', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f2', 'r1i1p1f2', 'r1i1p2f1', 'r1i1p2f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1']
#ens_mems = ['r1i1p1f1', 'r1i1p1f1', 'r1i1p1f2', 'r1i1p1f2', 'r3i1p2f1', 'r1i1p2f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1']

vers = ['LR', 'HR', 'LR', 'HR', 'LR', 'HR', 'LR', 'HR', 'LR', 'MR', 'HR', 'LR', 'HR', 'LR', 'MR', 'MR', 'HR'] + ['OBS']
#vers = 6*['LR', 'HR'] + ['OBS']
model_coups = ['AWI-CM-1-1', 'CMCC-CM2', 'CNRM-CM6-1', 'EC-Earth3P', 'ECMWF-IFS', 'MPI-ESM1-2', 'HadGEM3-GC31']
#model_coups = ['CMCC-CM2', 'CNRM-CM6-1', 'EC-Earth3P', 'ECMWF-IFS', 'MPI-ESM1-2', 'HadGEM3-GC31']
model_names_all = model_names + ['ERA']

# #colors = ctl.color_set(len(model_names), sns_palette = 'Paired')
# colors = ctl.color_set(len(model_coups)*2, sns_palette = 'Paired')
# mr_cos = np.where(np.array(vers) == 'MR')[0]
# for gi in mr_cos:
#     colors.insert(gi, np.mean(colors[gi-1:gi+1], axis = 0))
#
# colors_wERA = colors + ['black']
# # colors_wERA2 = np.array(colors_wERA)
# # colors_wERA2[0:-1:2] = colors_wERA2[1::2]


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


ens_names = ['LR', 'HR']
ens_colors = ['teal', 'indianred']
col_LR = ens_colors[0]
col_HR = ens_colors[1]

regnam = ['NAO +', 'Sc. BL', 'AR', 'NAO -']

################################################################################

# results_refEOF, results_ref = pickle.load(open(filogen, 'rb'))
# results_refEOF['ERA'] = results_ref

# cart = '/home/fabiano/Research/lavori/WeatherRegimes/panosjet/'
# filogen = cart + 'out_panosjet_DJF_EAT_4clus_4pcs_1957-2014_refCLUS.p'
# results_refCLUS, results_ref = pickle.load(open(filogen, 'rb'))
# results_refCLUS['ERA'] = results_ref
#
# cart = '/home/fabiano/Research/lavori/WeatherRegimes/panosjet/'
# filogen = cart + 'out_panosjet_DJF_EAT_4clus_4pcs_1957-2014_refEOF.p'
# results_refEOF, results_ref = pickle.load(open(filogen, 'rb'))
# results_refEOF['ERA'] = results_ref

cart = '/home/fabiano/Research/lavori/WeatherRegimes/prima_coup_v7/'
filogen = cart + 'out_prima_coup_v7_DJF_EAT_4clus_4pcs_1957-2014_refCLUS.p'
results_refCLUS, results_ref = pickle.load(open(filogen, 'rb'))
results_refCLUS['ERA'] = results_ref

cart = '/home/fabiano/Research/lavori/WeatherRegimes/prima_coup_v7/'
filogen = cart + 'out_prima_coup_v7_DJF_EAT_4clus_4pcs_1957-2014_refEOF.p'
results_refEOF, results_ref = pickle.load(open(filogen, 'rb'))
results_refEOF['ERA'] = results_ref

results_refEOF.pop('HadGEM3-GC31-LL_r1i1p2f1')
results_refEOF.pop('EC-Earth3P_r1i1p1f1')
results_refEOF.pop('EC-Earth3P-HR_r1i1p1f1')
allresmembers = list(results_refEOF.keys())

# mod = model_names[5]
# mem = ens_mems[5]
# print(mod, mem)
mod = 'ERA'
filok = 'bf.5day.daily.daily_mean.ERA.ERA.1958-2017.nc'

blocked_days, datacoords, aux_info = ctl.readxDncfield(cart_bloc + filok)
dates_block = datacoords['dates']
lat = datacoords['lat']
lon = datacoords['lon']

WR_index = results_refEOF[mod]['labels']
WR_dates = results_refEOF[mod]['dates']

blok, wri, datcom = ctl.extract_common_dates(dates_block, WR_dates, blocked_days, WR_index)
ERA_map_full = np.mean(blok, axis = 0)
ctl.plot_map_contour(ERA_map_full, lat, lon, visualization = 'Nstereo', plot_anomalies = False, filename = cart_out_maps + 'map_full_ERA.pdf')

def func_sum(blok_ind, reg_ind):
    alsums = []
    for reg in range(4):
        blokok = blok_ind[reg_ind == reg]
        okmap = np.mean(blokok, axis = 0)
        alsums.append(np.mean(okmap))

    return np.array(alsums)

blok_anom = blok - ERA_map_full
blok_anom_area, _, _ = ctl.sel_area(lat, lon, blok_anom, 'EAT')
res_boot_ERA = ctl.bootstrap(blok_anom_area, datcom, 'DJF', y = wri, apply_func = func_sum, n_choice = 30, n_bootstrap = 500)

EAT_mean = dict()
ERA_map_area, _, _ =  ctl.sel_area(lat, lon, ERA_map_full, 'EAT')
EAT_mean['ERA_0'] = np.mean(ERA_map_area)

# allmaps[('ERA', 'full')], lato, lono = ctl.sel_area(lat, lon, ERA_map_full, 'EAT')
#
# ERA_allmaps = []
# for reg in range(4):
#     okind = wri == reg
#     okblo_map = np.mean(blok[okind, ...], axis = 0) - ERA_map_full
#     ERA_allmaps.append(okblo_map)
#     allmaps[('ERA', reg)], lato, lono = ctl.sel_area(lat, lon, okblo_map, 'EAT')
#
# ERA_allmaps = np.stack(ERA_allmaps)
# ctl.plot_multimap_contour(ERA_allmaps, lat, lon, filename = cart_out_maps + 'map_regs_ERA.pdf', visualization = 'Nstereo')

#for results, tag in zip([results_refEOF, results_refCLUS], ['', '_refCLUS']):
for results, tag in zip([results_refEOF], ['']):
    allmaps = dict()
    allrms = dict()
    allpatcor = dict()
    allsums = dict()
    res_boot = dict()

    allmaps[('ERA', 'full')], lato, lono = ctl.sel_area(lat, lon, ERA_map_full, 'EAT')

    ERA_allmaps = []
    for reg in range(4):
        okind = wri == reg
        okblo_map = np.mean(blok[okind, ...], axis = 0) - ERA_map_full
        ERA_allmaps.append(okblo_map)
        allmaps[('ERA', reg)], lato, lono = ctl.sel_area(lat, lon, okblo_map, 'EAT')
        allsums[('ERA', reg)] = np.sum(allmaps[('ERA', reg)])

    ERA_allmaps = np.stack(ERA_allmaps)
    ctl.plot_multimap_contour(ERA_allmaps, lat, lon, filename = cart_out_maps + 'map_regs_ERA.pdf', visualization = 'Nstereo')

    for mod in model_names:
        filok = [fi for fi in all_fils_bloc if mod == fi.split('.')[5]]
        all_mems = [fi.split('.')[7] for fi in filok]
        print(filok)

        for filo, mem in zip(filok, all_mems):
            if 'EC-Earth' in mod and mem == 'r1i1p1f1':
                continue
            if mod == 'HadGEM3-GC31-LL' and mem == 'r1i1p2f1':
                continue
            blocked_days, datacoords, aux_info = ctl.readxDncfield(cart_bloc + filo)
            dates_block = datacoords['dates']
            lat = datacoords['lat']
            lon = datacoords['lon']

            modmem = mod + '_' + mem
            WR_index = results[modmem]['labels']
            WR_dates = results[modmem]['dates']

            blok, wri, datcom = ctl.extract_common_dates(dates_block, WR_dates, blocked_days, WR_index)
            map_full = np.mean(blok, axis = 0)
            ctl.plot_map_contour(map_full, lat, lon, visualization = 'Nstereo', plot_anomalies = False, filename = cart_out_maps + 'map_full_{}_{}.pdf'.format(mod,mem), cbar_range = (0., 0.1))

            allmaps[(mod, mem, 'full')], lato, lono = ctl.sel_area(lat, lon, map_full, 'EAT')
            allrms[(mod, mem, 'full')] = ctl.E_rms(allmaps[(mod, mem, 'full')], allmaps[('ERA', 'full')], lato)
            allpatcor[(mod, mem, 'full')] = ctl.Rcorr(allmaps[(mod, mem, 'full')], allmaps[('ERA', 'full')], lato)

            allma = []
            for reg in range(4):
                okind = wri == reg
                okblo_map = np.mean(blok[okind, ...], axis = 0) - map_full
                allma.append(okblo_map)
                allmaps[(mod, mem, reg)], lato, lono = ctl.sel_area(lat, lon, okblo_map, 'EAT')
                allrms[(mod, mem, reg)] = ctl.E_rms(allmaps[(mod, mem, reg)], allmaps[('ERA', reg)], lato)
                allpatcor[(mod, mem, reg)] = ctl.Rcorr(allmaps[(mod, mem, reg)], allmaps[('ERA', reg)], lato)
                allsums[(mod, mem, reg)] = np.sum(allmaps[(mod, mem, reg)])

            ctl.plot_multimap_contour(allma, lat, lon, filename = cart_out_maps + 'map_regs_{}_{}{}.pdf'.format(mod,mem,tag), visualization = 'Nstereo', cbar_range = (-0.2, 0.2))

            blok_anom = blok - map_full
            blok_anom_area, _, _ = ctl.sel_area(lat, lon, blok_anom, 'EAT')
            res_boot[modmem] = ctl.bootstrap(blok_anom_area, datcom, 'DJF', y = wri, apply_func = func_sum, n_choice = 30, n_bootstrap = 500)

            map_area, _, _ =  ctl.sel_area(lat, lon, map_full, 'EAT')
            EAT_mean[modmem] = np.mean(map_area)

    pickle.dump([allmaps, allrms, allpatcor, allsums], open(cart_out + 'out_bloccomp_v7{}.p'.format(tag), 'wb'))
    pickle.dump([res_boot, res_boot_ERA, EAT_mean], open(cart_out + 'out_blocboot_v7{}.p'.format(tag), 'wb'))

    res_boot, res_boot_ERA, EAT_mean = pickle.load(open(cart_out + 'out_blocboot_v7{}.p'.format(tag), 'rb'))
    allmaps, allrms, allpatcor, allsums = pickle.load(open(cart_out + 'out_bloccomp_v7{}.p'.format(tag), 'rb'))

    for ke in allsums:
        allsums[ke] = allsums[ke]/allmaps[ke].size

    #### PLOTTTTTS
    patnames = ['NAO +', 'Sc. Blocking', 'Atl. Ridge', 'NAO -']

    fig = plt.figure(figsize = (16, 12))
    axes = []
    for reg in range(4):
        ax = plt.subplot(2, 2, reg+1)
        axes.append(ax)
        #ax.set_xlabel('Latitude')
        ax.set_ylabel('RMS')
        ax.set_title(patnames[reg])
        ax.set_xticks([])

        i = 0
        wi = 0.6

        for mod, col, vv in zip(model_names, colors, vers):
            print('aaaaaaaaaaaaaaaaaaa', mod)
            rmsall = []
            for ke in allrms.keys():
                if mod in ke and reg in ke:
                    print(ke)
                    rmsall.append(allrms[ke])

            if len(rmsall) == 0: continue

            allrms[(mod, reg)] = np.mean(rmsall)

            if len(rmsall) == 1:
                ax.bar(i, allrms[(mod, reg)], width = wi, color = col)
            else:
                print('{} has {} ens members'.format(mod, len(rmsall)))
                okrms = np.mean(rmsall)
                rmserr_lo = np.min(rmsall)
                rmserr_hi = np.max(rmsall)
                print(okrms, rmserr_lo, rmserr_hi)
                yerr = np.array([[okrms-rmserr_lo, rmserr_hi-okrms]]).T

                ax.bar(i, okrms, width = wi, color = col, yerr = yerr, capsize = 5)

            i += 0.7
            if vv == 'HR': i += 0.2

        all_LR = np.mean([allrms[(mod, reg)] for mod, vv in zip(model_names_all, vers) if vv == 'LR'])
        all_HR = np.mean([allrms[(mod, reg)] for mod, vv in zip(model_names_all, vers) if vv == 'HR'])
        all_LR_std = np.std([allrms[(mod, reg)] for mod, vv in zip(model_names_all, vers) if vv == 'LR'])
        all_HR_std = np.std([allrms[(mod, reg)] for mod, vv in zip(model_names_all, vers) if vv == 'HR'])

        i+=0.4

        ax.bar(i, all_LR, width = wi, color = col_LR, yerr = all_LR_std, capsize = 5)
        i += 0.7
        ax.bar(i, all_HR, width = wi, color = col_HR, yerr = all_HR_std, capsize = 5)

    ctl.custom_legend(fig, colors + ens_colors, model_names + ens_names)

    #fig.suptitle('RMS of blocking anomalies')
    fig.savefig(cart_out + 'block_allregs_rms{}.pdf'.format(tag))



    fig = plt.figure(figsize = (16, 12))
    axes = []
    for reg in range(4):
        ax = plt.subplot(2, 2, reg+1)
        axes.append(ax)
        #ax.set_xlabel('Latitude')
        ax.set_ylabel('patcor')
        ax.set_title(patnames[reg])
        ax.set_xticks([])

        i = 0
        wi = 0.6

        for mod, col, vv in zip(model_names, colors, vers):
            rmsall = []
            for ke in allpatcor.keys():
                if mod in ke and reg in ke:
                    rmsall.append(allpatcor[ke])

            if len(rmsall) == 0: continue

            allpatcor[(mod, reg)] = np.mean(rmsall)

            if len(rmsall) == 1:
                ax.bar(i, allpatcor[(mod, reg)], width = wi, color = col)
            else:
                print('{} has {} ens members'.format(mod, len(rmsall)))
                okrms = np.mean(rmsall)
                rmserr_lo = np.min(rmsall)
                rmserr_hi = np.max(rmsall)
                print(okrms, rmserr_lo, rmserr_hi)
                yerr = np.array([[okrms-rmserr_lo, rmserr_hi-okrms]]).T

                ax.bar(i, okrms, width = wi, color = col, yerr = yerr, capsize = 5)

            i += 0.7
            if vv == 'HR': i += 0.2

        all_LR = np.mean([allpatcor[(mod, reg)] for mod, vv in zip(model_names_all, vers) if vv == 'LR'])
        all_HR = np.mean([allpatcor[(mod, reg)] for mod, vv in zip(model_names_all, vers) if vv == 'HR'])
        all_LR_std = np.std([allpatcor[(mod, reg)] for mod, vv in zip(model_names_all, vers) if vv == 'LR'])
        all_HR_std = np.std([allpatcor[(mod, reg)] for mod, vv in zip(model_names_all, vers) if vv == 'HR'])

        i+=0.4

        ax.bar(i, all_LR, width = wi, color = col_LR, yerr = all_LR_std, capsize = 5)
        i += 0.7
        ax.bar(i, all_HR, width = wi, color = col_HR, yerr = all_HR_std, capsize = 5)

    ctl.custom_legend(fig, colors + ens_colors, model_names + ens_names)

    #fig.suptitle('Pattern correlation of blocking anomalies')
    fig.savefig(cart_out + 'block_allregs_patcor{}.pdf'.format(tag))


    fig = plt.figure(figsize = (16, 12))
    axes = []
    for reg in range(4):
        ax = plt.subplot(2, 2, reg+1)
        axes.append(ax)
        #ax.set_xlabel('Latitude')
        ax.set_ylabel('Freq. anomaly')
        ax.set_title(patnames[reg])
        ax.set_xticks([])

        i = 0
        wi = 0.6

        for mod, col, vv in zip(model_names, colors, vers):
            rmsall = []
            for ke in allsums.keys():
                if mod in ke and reg in ke:
                    rmsall.append(allsums[ke])

            if len(rmsall) == 0: continue

            allsums[(mod, reg)] = np.mean(rmsall)

            if len(rmsall) == 1:
                ax.bar(i, allsums[(mod, reg)], width = wi, color = col)
            else:
                print('{} has {} ens members'.format(mod, len(rmsall)))
                okrms = np.mean(rmsall)
                rmserr_lo = np.min(rmsall)
                rmserr_hi = np.max(rmsall)
                print(okrms, rmserr_lo, rmserr_hi)
                yerr = np.array([[okrms-rmserr_lo, rmserr_hi-okrms]]).T

                ax.bar(i, okrms, width = wi, color = col, yerr = yerr, capsize = 5)

            i += 0.7
            if vv == 'HR': i += 0.2

        all_LR = np.mean([allsums[(mod, reg)] for mod, vv in zip(model_names_all, vers) if vv == 'LR'])
        all_HR = np.mean([allsums[(mod, reg)] for mod, vv in zip(model_names_all, vers) if vv == 'HR'])
        all_LR_std = np.std([allsums[(mod, reg)] for mod, vv in zip(model_names_all, vers) if vv == 'LR'])
        all_HR_std = np.std([allsums[(mod, reg)] for mod, vv in zip(model_names_all, vers) if vv == 'HR'])

        i += 0.4
        ax.bar(i, allsums[('ERA', reg)], width = wi, color = 'black', capsize = 5)
        i += 0.7

        ax.bar(i, all_LR, width = wi, color = col_LR, yerr = all_LR_std, capsize = 5)
        i += 0.7
        ax.bar(i, all_HR, width = wi, color = col_HR, yerr = all_HR_std, capsize = 5)


    ctl.adjust_ax_scale(axes)
    ctl.custom_legend(fig, colors_wERA + ens_colors, model_names_all + ens_names)

    #fig.suptitle('Blocking frequency anomaly in the EAT sector')
    fig.savefig(cart_out + 'block_allregs_sumeve{}.pdf'.format(tag))

    res_boot['ERA_0'] = res_boot_ERA

    bootstri = dict()
    models_ok = []
    for mod in model_names_all:
        print(mod)
        allmems = []
        #for ke in allresmembers:
        for ke in res_boot.keys():
            if mod == ke.split('_')[0]:
                allmems.append(ke.split('_')[1])

        if len(allmems) == 0:
            models_ok.append(False)
            continue
        models_ok.append(True)

        print(allmems)
        if mod == 'ERA':
            allmems = ['0']

        res_boot[mod] = np.concatenate([res_boot[mod+'_'+mem] for mem in allmems], axis = 0)
        bootstri[mod] = dict()

        # list of all boot means
        bootstri[mod]['boot_mean'] = np.mean(res_boot[mod], axis = 0)
        bootstri[mod]['boot_std'] = np.std(res_boot[mod], axis = 0)
        for cos, num in zip(['p10', 'p25', 'p50', 'p75', 'p90'], [10, 25, 50, 75, 90]):
            bootstri[mod]['boot_'+cos] = np.array([np.percentile(res_boot[mod][:, reg], num) for reg in range(4)])


    mind = model_names_all.index('ECMWF-IFS-MR')
    model_names_ok = model_names_all[0:mind] + model_names_all[mind+1:]
    colors_wERA_ok = colors_wERA[0:mind] + colors_wERA[mind+1:]
    color_main_ok = color_main[0:mind] + color_main[mind+1:]
    vers_ok = vers[0:mind] + vers[mind+1:]

    fig = plt.figure(figsize=(16,12))
    axes = []
    for ii, reg in enumerate(regnam):
        ax = plt.subplot(2, 2, ii+1)

        allpercs = dict()
        for cos in ['mean', 'std', 'p10', 'p25', 'p50', 'p75', 'p90']:
            allpercs[cos] = [bootstri[mod]['boot_'+cos][ii] for mod in model_names_ok]
        for cos in ['ens_min', 'ens_max']:
            allpercs[cos] = [bootstri[mod]['boot_mean'][ii] for mod in model_names_ok]

        allpercs['ens_std'] = np.zeros(len(model_names_ok))

        ctl.primavera_boxplot_on_ax(ax, allpercs, model_names_ok, colors_wERA_ok, color_main_ok, vers_ok, ens_names, ens_colors, wi = 0.5)

        ax.set_title(reg)
        ax.set_ylabel('Freq. anomaly')
        ax.axhline(0., color = 'grey', linewidth = 0.5, linestyle = ':')
        axes.append(ax)

    ctl.adjust_ax_scale(axes)

    #fig.suptitle(alltits[nam], fontsize = titlefont)
    plt.subplots_adjust(top = 0.9)

    ctl.custom_legend(fig, colors_wERA_ok + ens_colors, model_names_ok + ens_names, ncol = 6)
    fig.savefig(cart_out + 'blockave_bootplot_v7.pdf')

    for mod in model_names_ok:
        rmsall = []
        for ke in EAT_mean.keys():
            if mod in ke:
                rmsall.append(EAT_mean[ke])

        EAT_mean[mod] = np.mean(rmsall)


    fig = plt.figure(figsize=(16,12))
    axes = []
    for ii, reg in enumerate(regnam):
        ax = plt.subplot(2, 2, ii+1)

        allpercs = dict()
        for cos in ['mean', 'std', 'p10', 'p25', 'p50', 'p75', 'p90']:
            allpercs[cos] = [bootstri[mod]['boot_'+cos][ii]+EAT_mean[mod] for mod in model_names_ok]
        for cos in ['ens_min', 'ens_max']:
            allpercs[cos] = [bootstri[mod]['boot_mean'][ii]+EAT_mean[mod] for mod in model_names_ok]

        allpercs['ens_std'] = np.zeros(len(model_names_ok))

        ctl.primavera_boxplot_on_ax(ax, allpercs, model_names_ok, colors_wERA_ok, color_main_ok, vers_ok, ens_names, ens_colors, wi = 0.5)

        ax.set_title(reg)
        ax.set_ylabel('Block. freq.')
        ax.axhline(0., color = 'grey', linewidth = 0.5, linestyle = ':')
        axes.append(ax)

    ctl.adjust_ax_scale(axes)

    #fig.suptitle(alltits[nam], fontsize = titlefont)
    plt.subplots_adjust(top = 0.9)

    ctl.custom_legend(fig, colors_wERA_ok + ens_colors, model_names_ok + ens_names, ncol = 6)
    fig.savefig(cart_out + 'blockave_ABS_bootplot_v7.pdf')


##### Bootstrapping della mean anomaly
