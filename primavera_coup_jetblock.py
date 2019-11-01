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

#######################################
cart_out = '/home/fabiano/Research/articoli/Papers/primavera_regimes/figures/jet_and_blocking/'
if not os.path.exists(cart_out): os.mkdir(cart_out)

cart = '/home/fabiano/Research/lavori/WeatherRegimes/prima_coup_v6/'
filogen = cart + 'out_prima_coup_v6_DJF_EAT_4clus_4pcs_1957-2014_refEOF.p'

cart_jet = '/data-hobbes/fabiano/PRIMAVERA/jetlat_panos/hist-1950/'
fil_jet = 'JetLatDist_850hPa_{}_hist-1950_DJF.npy'

cart_bloc = '/data-hobbes/fabiano/PRIMAVERA/Reinhard_blocking/all_hist-1950/'
fil_bloc = 'bf.5day.daily.daily_mean.{}.{}.hist-1950.r1i1p1f1.v20170915.{}-{}.nc'

# model_names = ['AWI-CM-1-1-LR', 'AWI-CM-1-1-HR', 'CMCC-CM2-HR4', 'CMCC-CM2-VHR4', 'CNRM-CM6-1', 'CNRM-CM6-1-HR', 'EC-Earth3P', 'EC-Earth3P-HR', 'ECMWF-IFS-LR', 'ECMWF-IFS-HR', 'MPI-ESM1-2-HR', 'MPI-ESM1-2-XR', 'HadGEM3-GC31-LL', 'HadGEM3-GC31-MM']#'HadGEM3-GC31-HM']
model_names = ['CMCC-CM2-HR4', 'CMCC-CM2-VHR4', 'CNRM-CM6-1', 'CNRM-CM6-1-HR', 'EC-Earth3P', 'EC-Earth3P-HR', 'ECMWF-IFS-LR', 'ECMWF-IFS-HR', 'MPI-ESM1-2-HR', 'MPI-ESM1-2-XR', 'HadGEM3-GC31-LL', 'HadGEM3-GC31-MM']
model_mems = ['r1i1p1f1', 'r1i1p1f1', 'r1i1p1f2', 'r1i1p1f2', 'r3i1p2f1', 'r1i1p2f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', '']
vers = 6*['LR', 'HR'] + ['OBS']
# model_coups = ['AWI-CM-1-1', 'CMCC-CM2', 'CNRM-CM6-1', 'EC-Earth3P', 'ECMWF-IFS', 'MPI-ESM1-2', 'HadGEM3-GC31']
model_coups = ['CMCC-CM2', 'CNRM-CM6-1', 'EC-Earth3P', 'ECMWF-IFS', 'MPI-ESM1-2', 'HadGEM3-GC31']
model_names_all = model_names + ['ERA']

colors2 = ctl.color_set(len(model_coups))
col22 = list(np.concatenate([[co, co] for co in colors2]))
colors_wERA2 = col22 + ['black']
# styall = np.empty(len(colors_wERA2), dtype = str)
# styall[0::2] = ':'
# styall[1::2] = '-'
# styall[-1] = '-'

colors = ctl.color_set(len(model_names), sns_palette = 'Paired')
colors_wERA = colors + ['black']
# colors_wERA2 = np.array(colors_wERA)
# colors_wERA2[0:-1:2] = colors_wERA2[1::2]

regnam = ['NAO +', 'Sc. BL', 'AR', 'NAO -']

################################################################################

# results_refEOF, results_ref = pickle.load(open(filogen, 'rb'))
# results_refEOF['ERA'] = results_ref
#
# fil_ECE = '/home/fabiano/Research/lavori/WeatherRegimes/prima_ECE_only_allens/out_prima_ECE_only_allens_DJF_EAT_4clus_4pcs_1957-2014_refEOF.p'
# result_ECE, _ = pickle.load(open(fil_ECE, 'rb'))
# # Panos uses r1i1p2f1 (HR) and r3i1p2f1 (LR)
# # results['EC-Earth3P'] = result_ECE['ECE-Stream2-LR']
# results_refEOF['EC-Earth3P-HR'] = result_ECE['ECE-r1i1p2f1-HR']
# results_refEOF['EC-Earth3P'] = result_ECE['ECE-r3i1p2f1-LR']
#
# fil_Had = '/home/fabiano/Research/lavori/WeatherRegimes/prima_Had_only/out_prima_Had_only_DJF_EAT_4clus_4pcs_1957-2014_refEOF.p'
# result_Had, _ = pickle.load(open(fil_Had, 'rb'))
# results_refEOF['HadGEM3-GC31-MM'] = result_Had['HadGEM-MM']
# results_refEOF['HadGEM3-GC31-LL'] = result_Had['HadGEM-LL']
# results_refEOF['HadGEM3-GC31-HM'] = result_Had['HadGEM-HM']


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

cart = '/home/fabiano/Research/lavori/WeatherRegimes/prima_coup_v7_rmshi/'
filogen = cart + 'out_prima_coup_v7_rmshi_DJF_EAT_4clus_4pcs_1957-2014_refEOF.p'
results_refEOF_2, results_ref = pickle.load(open(filogen, 'rb'))
results_refEOF_2['ERA'] = results_ref

reg_change = dict()
reg_change['CMCC-CM2-HR4'] = [0, 3, 1, 2]
reg_change['EC-Earth3P-HR'] = [2, 1, 0, 3]

for results, tag in zip([results_refEOF, results_refCLUS, results_refEOF_2], ['', '_refCLUS', '_v7_optmatch']):
    latgri = np.arange(20,71,0.5)

    patnames = ['NAO +', 'Sc. Blocking', 'Atl. Ridge', 'NAO -']

    jetlat_comp = dict()
    # first for the jet
    fig_all = plt.figure(figsize = (16, 12))
    ax_all = plt.subplot()

    fig = plt.figure(figsize = (16, 12))
    axes = []
    for reg in range(4):
        ax = plt.subplot(2, 2, reg+1)
        axes.append(ax)
        ax.set_xlabel('Latitude')
        ax.set_ylabel('Rel. frequency')
        ax.set_title(patnames[reg])
        ax.grid()

    fig_all_HR = plt.figure(figsize = (16, 12))
    ax_all_HR = plt.subplot()

    fig_HR = plt.figure(figsize = (16, 12))
    axes_HR = []
    for reg in range(4):
        ax = plt.subplot(2, 2, reg+1)
        axes_HR.append(ax)
        ax.set_xlabel('Latitude')
        ax.set_ylabel('Rel. frequency')
        ax.set_title(patnames[reg])
        ax.grid()

    sty = '-'
    for mod, col, vv, mem in zip(model_names_all, colors_wERA2, vers, model_mems):
        print(mod, mem, sty)
        fil_ok = fil_jet.format(mod)
        # if mod == 'HadGEM3-GC31-LL-stoc':
        #     fil_ok = fil_jet.format('HadGEM3-GC31-LL')
        if mod == 'ERA':
            fil_ok = 'JetLatDist_850hPa_NCEP_DJF.npy'
        print(fil_ok)
        # elif mod == 'EC-Earth3P':
        #     print('Missing 1991')
        #     continue

        if not os.path.exists(cart_jet + fil_ok):
            print('waiting for panos..\n\n')
            continue

        jetind = np.load(cart_jet + fil_ok)[0]
        jetlat_comp[(mod, 'all', tag)] = jetind

        modmem = mod + '_' + mem

        if mod == 'ERA':
            regind = results[mod]['labels']
            dates = results[mod]['dates']
        else:
            regind = results[modmem]['labels']
            dates = results[modmem]['dates']

        dates_pdh = pd.to_datetime(dates)
        okdat = ~((dates_pdh.month == 2) & (dates_pdh.day == 29))
        print('Num leap days : {}'.format(len(regind) - np.sum(okdat)))
        regind = regind[okdat]
        dates = dates[okdat]

        if 'HadGEM' in mod:
            # 31 31 28 (panos); 30 31 28 (io)
            # splitto panos in cosi da 90; tolgo primo giorno e sto a posto
            nch = len(jetind)//90
            print(nch)
            jetind_sp = np.split(jetind, nch)
            jetind = np.concatenate([ji[1:] for ji in jetind_sp])
            print(len(jetind))
            jetind = jetind[7*89:]

            data1 = pd.to_datetime('{}1201'.format(1957), format='%Y%m%d')
            data2 = pd.to_datetime('{}0228'.format(2014), format='%Y%m%d')
        elif mod == 'EC-Earth3P':
            jetind = jetind[1260:] # the first day is 01-12-1964, the last is 28/02/2014
            data1 = pd.to_datetime('{}1201'.format(1964), format='%Y%m%d')
            data2 = pd.to_datetime('{}0228'.format(2014), format='%Y%m%d')
        else:
            # Skip the first 7 seasons (panos data start from 1950): 7*90 = 630
            # Panos uses 28 day february
            jetind = jetind[630:] # the first day is 01-12-1957, the last is 28/02/2014
            data1 = pd.to_datetime('{}1201'.format(1957), format='%Y%m%d')
            data2 = pd.to_datetime('{}0228'.format(2014), format='%Y%m%d')

        print(len(regind), len(dates), data1, data2)
        regok, datok = ctl.sel_time_range(regind, dates, (data1, data2))
        print(len(regok), len(jetind))

        #if len(regok) != 5130:
        if len(regok) != len(jetind):
            print('SKIPPING!!!!\n')
            continue
            #raise ValueError('OPS!')

        linewidth = 1.5
        if mod == 'ERA': linewidth = 2.5

        for reg in range(4):
            okin = regok == reg
            okjet = jetind[okin]
            if tag == '_v7_optmatch':
                if mod == 'CMCC-CM2-HR4' or mod == 'EC-Earth3P-HR':
                    okin = regok == reg_change[mod][reg]
                    okjet = jetind[okin]
            jetlat_comp[(mod, reg, tag)] = okjet

            pdf = ctl.calc_pdf(okjet)
            if vv in ['LR', 'OBS']:
                axes[reg].plot(latgri, pdf(latgri), color = col, label = mod, linewidth = linewidth, linestyle = sty)
            if vv in ['HR', 'OBS']:
                axes_HR[reg].plot(latgri, pdf(latgri), color = col, label = mod, linewidth = linewidth, linestyle = sty)

        pdf = ctl.calc_pdf(jetind)
        if vv in ['LR', 'OBS']:
            ax_all.plot(latgri, pdf(latgri), color = col, label = mod, linewidth = linewidth, linestyle = sty)
        if vv in ['HR', 'OBS']:
            ax_all_HR.plot(latgri, pdf(latgri), color = col, label = mod, linewidth = linewidth, linestyle = sty)

    ctl.adjust_ax_scale(axes+axes_HR)
    # ctl.custom_legend(fig, colors_wERA, model_names_all)
    ctl.custom_legend(fig, colors_wERA2[0::2], model_coups + ['OBS'])
    ctl.custom_legend(fig_HR, colors_wERA2[0::2], model_coups + ['OBS'])
    fig.suptitle('Jet Latitude and Weather Regimes (SR)')
    fig.savefig(cart_out + 'jetlat_compos_LR{}.pdf'.format(tag))
    fig_HR.suptitle('Jet Latitude and Weather Regimes (HR)')
    fig_HR.savefig(cart_out + 'jetlat_compos_HR{}.pdf'.format(tag))
    # fig.savefig(cart_out + 'jetlat_compos_wECEstream2.pdf')

    ctl.adjust_ax_scale([ax_all, ax_all_HR])
    ctl.custom_legend(fig_all, colors_wERA2[0::2], model_coups + ['OBS'])
    #ctl.custom_legend(fig_all, colors_wERA, model_names_all)
    fig_all.savefig(cart_out + 'jetlat_allregs_LR{}.pdf'.format(tag))
    ctl.custom_legend(fig_all_HR, colors_wERA2[0::2], model_coups + ['OBS'])
    fig_all_HR.savefig(cart_out + 'jetlat_allregs_HR{}.pdf'.format(tag))

    pickle.dump(jetlat_comp, open(cart_out + 'jetlat_compos_all{}.p'.format(tag), 'wb'))

    # Calculate relative entropy

    fig = plt.figure(figsize = (16, 12))
    axes = []
    for reg in range(4):
        ax = plt.subplot(2, 2, reg+1)
        axes.append(ax)
        #ax.set_xlabel('Latitude')
        ax.set_ylabel('Rel. entropy')
        ax.set_title(patnames[reg])
        ax.set_xticks([])
        #ax.grid()

    sty = '-'
    relent_all = dict()
    i = 0
    wi = 0.6
    for mod, col, vv in zip(model_names_all, colors_wERA, vers):
        for reg in range(4):
            obsjet = jetlat_comp[('ERA', reg, tag)]
            pdfob = ctl.calc_pdf(obsjet)
            okjet = jetlat_comp[(mod, reg, tag)]
            pdf = ctl.calc_pdf(okjet)

            relent = stats.entropy(pdf(latgri), pdfob(latgri))
            relent_all[(mod, reg, tag)] = relent
            axes[reg].bar(i, relent, width = wi, color = col)
        i+=0.7

    ctl.adjust_ax_scale(axes)
    ctl.custom_legend(fig, colors_wERA, model_names_all)

    fig.suptitle('Relative entropy of jet lat composites')
    fig.savefig(cart_out + 'jetlat_compos_relent{}.pdf'.format(tag))


    fig_all = plt.figure(figsize = (16, 12))
    ax = fig_all.add_axes([0.1, 0.55, 0.8, 0.4])
    ax2 = fig_all.add_axes([0.1, 0.1, 0.8, 0.45], sharex = ax)

    i = 0
    wi = 0.6
    for mod, col, vv in zip(model_names_all, colors_wERA, vers):
        obsjet = jetlat_comp[('ERA', 'all', tag)]
        pdfob = ctl.calc_pdf(obsjet)
        okjet = jetlat_comp[(mod, 'all', tag)]
        pdf = ctl.calc_pdf(okjet)

        relent = stats.entropy(pdf(latgri), pdfob(latgri))
        relent_all[(mod, 'all', tag)] = relent
        ax.bar(i, relent, width = wi, color = col)
        #print(i)

        if vv == 'HR':
            #print('a2 ',np.mean([i, lasti]))
            ax2.bar(np.mean([i, lasti]), relent - relent_all[(lastmod, 'all', tag)], width = wi, color = np.mean([lastcol, col], axis = 0))

        lastmod = mod
        lasti = i
        lastcol = col

        i += 0.7
        if vv == 'HR': i += 0.2

    all_LR = np.mean([relent_all[(mod, 'all', tag)] for mod, vv in zip(model_names_all, vers) if vv == 'LR'])
    all_LR_std = np.std([relent_all[(mod, 'all', tag)] for mod, vv in zip(model_names_all, vers) if vv == 'LR'])
    all_HR = np.mean([relent_all[(mod, 'all', tag)] for mod, vv in zip(model_names_all, vers) if vv == 'HR'])
    all_HR_std = np.std([relent_all[(mod, 'all', tag)] for mod, vv in zip(model_names_all, vers) if vv == 'HR'])

    i+=0.2
    col_LR = 'grey'
    col_HR = 'black'
    ax.bar(i, all_LR, width = wi, color = col_LR)
    i += 0.7
    ax.bar(i, all_HR, width = wi, color = col_HR)
    i -= 0.35
    ax2.bar(i, all_HR-all_LR, width = wi, color = 'darkslategray')

    ctl.custom_legend(fig_all, colors_wERA, model_names_all)
    ax2.axhline(0.0, color = 'grey')

    fig_all.suptitle('Relative entropy of jet lat')
    fig_all.savefig(cart_out + 'jetlat_allregs_relent{}.pdf'.format(tag))

    pickle.dump(relent_all, open(cart_out + 'jetdist_relent{}.p'.format(tag), 'wb'))
