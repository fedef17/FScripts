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

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 15
###############################################################################

def read_file_index(filename, var_name = None):
    fh = nc.Dataset(filename,'r')
    ### TIME ###
    time = fh.variables['time'][:]
    time_units = fh.variables['time'].units # Reading in the time units
    time_cal = fh.variables['time'].calendar # Calendar in use (proleptic_gregorian))
    dates = nc.num2date(time, time_units, time_cal)

    if time_cal == '365_day' or time_cal == 'noleap':
        dates = ctl.adjust_noleap_dates(dates)
    elif time_cal == '360_day':
        dates = ctl.adjust_360day_dates(dates)

    if var_name is not None:
        var = fh.variables[var_name][:]
    else:
        #print(fh.variables.keys())
        #print('Returning the last variable')
        var_name = list(fh.variables.keys())[-1]
        var = fh.variables[var_name][:]

    fh.close()
    return var, dates

#####################################################################################

indexname = 'AMV'

cart_out_all = '/home/fabiano/Research/lavori/WR_and_modes/'

file_AMV_mod = '/nas/PRIMAVERA/indices/Stream1/hist-1950/AMV/AMVindex-monthly_{}_{}_hist-1950_{}.nc'
file_ENSO_mod = '/nas/PRIMAVERA/indices/Stream1/hist-1950/ENSO/N3.4Index_{}_{}_hist-1950_{}.nc'
sedat = '1950-2014'

file_amv_ref = '/nas/reference/indices/AMV/AMVindex-monthly_ERA40+Interim_1957-2018.nc'
file_enso_ref = '/nas/reference/indices/ENSO/N3.4Index.nc-monthly_ERA40+Interim_1957-2018.nc'

#cart = '/home/fabiano/Research/lavori/WeatherRegimes/prima_coup_v7/'
cart = '/nas/PRIMAVERA/indices/Stream1/hist-1950/WeatherRegimes/'
filogen_EAT = cart + 'DJFM_EAT/out_prima_coup_v7_DJFM_EAT_4clus_4pcs_1957-2014_refEOF.p'
filogen_PNA = cart + 'DJFM_PNA/out_prima_coup_v7_DJFM_PNA_4clus_4pcs_1957-2014_refEOF.p'
results, results_ref = pickle.load(open(filogen_EAT, 'rb'))

#####################################################################################
indexes = ['AMV', 'ENSO']
file_refs = [file_amv_ref, file_enso_ref]
file_mods = [file_AMV_mod, file_ENSO_mod]
filogens = [filogen_EAT, filogen_PNA]
n_yrs = [5, 1]
reg_names_all = [['NAO+', 'SBL', 'AR', 'NAO-'], ['PT', 'PNA-', 'PNA+', 'AR']]

seasons = ['DJFM', 'DJ', 'FM', 'DJF']

allcorrs = dict()
allfreqs_posneg = dict()
allfreqs_terc = dict()

for indexname, file_ref, file_mod, filogen, reg_names, n_yr in zip(indexes, file_refs, file_mods, filogens, reg_names_all, n_yrs):
    cart_out_ind = cart_out_all + indexname + '/all_models/'
    if not os.path.exists(cart_out_ind): os.mkdir(cart_out_ind)

    results, results_ref = pickle.load(open(filogen, 'rb'))
    #results['ERA'] = results_ref

    for modmem in results.keys():
        plt.close('all')
        print(modmem)
        cart_mod = cart_out_ind + modmem + '/'
        if not os.path.exists(cart_mod): os.mkdir(cart_mod)

        mod, mem = modmem.split('_')

        try:
            if 'AWI' in mod:
                sedat = '1950-2010'
                y2 = 2010
                if 'LR' in mod:
                    modok = 'AWI-CM-1-0-LR'
                else:
                    modok = 'AWI-CM-1-0-HR'
                amv_ref, dates = read_file_index(file_mod.format(modok, mem, sedat))
            else:
                sedat = '1950-2014'
                y2 = 2014
                amv_ref, dates = read_file_index(file_mod.format(mod, mem, sedat))
        except FileNotFoundError:
            print('NOT FOUND!!!')
            continue

        amv_ref, dates = ctl.sel_time_range(amv_ref, dates, ctl.range_years(1957, y2-1))
        # Yearly amv index
        amv_ref_yr, yrdates = ctl.yearly_average(amv_ref, dates)
        amv_ref_yr = np.squeeze(amv_ref_yr)
        if len(amv_ref_yr) < y2-1957:
            for ye in range(1957, y2):
                if ye not in [da.year for da in yrdates]: print('Missing {}'.format(ye))
            continue
        # amv_ref_djf, dates_djf = ctl.sel_season(amv_ref, dates, seas)
        # amv_ref_djf = np.squeeze(amv_ref_djf)
        amv_nyr = np.array(ctl.running_mean(amv_ref_yr, n_yr))
        amvc = amv_nyr[~np.isnan(amv_nyr)]

        fig_all = plt.figure(figsize = (16,12))
        ax_all = []
        for i in range(4): ax_all.append(fig_all.add_subplot(2, 2, i+1))

        fig_all_ul = plt.figure(figsize = (16,12))
        ax_all_ul = []
        for i in range(4): ax_all_ul.append(fig_all_ul.add_subplot(2, 2, i+1))

        for ise, seas in enumerate(seasons):
            print(seas)
            y1 = 1957
            if seas == 'FM': y1 = 1958
            print(indexname, seas, y1)
            cart_out = cart_out_ind + seas + '/'
            if not os.path.exists(cart_out): os.mkdir(cart_out)

            # plt.ion()
            fig = plt.figure(figsize = (16,12))
            for reg in range(4):
                ax = fig.add_subplot(2, 2, reg+1)

                #freq_ok, dates_ok = ctl.sel_season(results_ref['monthly_freq']['freq'][reg], results_ref['monthly_freq']['dates'], seas)
                freq_seas, dates_seas = ctl.seasonal_set(results[modmem]['monthly_freq']['freq'][reg], results[modmem]['monthly_freq']['dates'], seas, dates_range = ctl.range_years(y1, y2))
                freq_seas = np.mean(freq_seas, axis = 1)

                if len(freq_seas) < len(amv_ref_yr):
                    print('aaaaaaaaaaaa', len(freq_seas), len(amv_ref_yr))
                    continue

                years = np.array([da.year for da in yrdates])
                yealen = np.arange(len(years))

                # freq = np.array(ctl.running_mean(freq_ok, 15))
                freq = np.array(ctl.running_mean(freq_seas, n_yr))
                oks = ~np.isnan(freq)
                freq = freq[oks]
                #amvc = np.array(ctl.running_mean(amv_ref_djf, 15))

                #print(len(freq), len(years))
                years = years[oks]
                yealen = yealen[oks]

                rco = ctl.Rcorr(amvc, freq)
                print('AAAAAAAAAAA', rco)
                allcorrs[(indexname, modmem, reg, seas)] = rco
                print(allcorrs.keys())

                ax.set_title('Corr {}: {:5.2f}'.format(reg, rco))
                ax.plot(yealen, freq, color = 'steelblue')

                ax2 = ax.twinx()
                ax2.plot(yealen, amvc, color = 'indianred')

                ax.set_xticks(yealen[2::15])
                ax.set_xticklabels(years[2::15])
                ax.set_xlabel('Years')
                ax.set_ylabel('WR frequency')

            fig.suptitle('Correlation of WR frequency with {} index'.format(indexname))
            fig.savefig(cart_mod + '{}_{}_vs_WRfreq_{}.pdf'.format(modmem, indexname, seas))

            if len(freq_seas) < len(amv_ref_yr):
                print('aaaaaaaaaaaa', len(freq_seas), len(amv_ref_yr))
                continue


            hi_amv = amvc > 0.
            lo_amv = amvc <= 0.

            fig = plt.figure(figsize = (16,12))
            ax = fig.add_subplot(111)

            allfr = []
            for reg in range(4):
                #freq_ok, dates_ok = ctl.sel_season(results_ref['monthly_freq']['freq'][reg], results_ref['monthly_freq']['dates'], seas)
                freq_seas, dates_seas = ctl.seasonal_set(results[modmem]['monthly_freq']['freq'][reg], results[modmem]['monthly_freq']['dates'], seas, dates_range = ctl.range_years(y1, y2))
                freq_seas = np.mean(freq_seas, axis = 1)

                freq = np.array(ctl.running_mean(freq_seas, n_yr))
                oks = ~np.isnan(freq)
                freq = freq[oks]

                # freqs_hi = freq_seas[hi_amv]
                # freqs_lo = freq_seas[lo_amv]
                freqs_hi = freq[hi_amv]
                freqs_lo = freq[lo_amv]

                allfreqs_posneg[(indexname, modmem, reg, seas)] = np.array([np.mean(freqs_hi), np.mean(freqs_lo)])

                allfr.append(freqs_hi)
                allfr.append(freqs_lo)

            pos = [1,2,4,5,7,8,10,11]
            posthicks = [1.5, 4.5, 7.5, 10.5]

            ax.boxplot(allfr, positions = np.array(pos))
            ax.set_xticks(posthicks)
            ax.set_xticklabels(reg_names)
            ax.set_ylabel('WR frequency')

            ax_all[ise].boxplot(allfr, positions = np.array(pos))
            ax_all[ise].set_xticks(posthicks)
            ax_all[ise].set_xticklabels(reg_names)
            ax_all[ise].set_ylabel('WR frequency')
            ax_all[ise].set_title('season: {}'.format(seas))


            fig.suptitle('Positive vs negative {}'.format(indexname))

            fig.savefig(cart_mod + '{}_{}_vs_WRfreq_posneg_{}.pdf'.format(modmem, indexname, seas))


            # hi_amv = amv_ref_yr > np.percentile(amv_ref_yr, 67)
            # lo_amv = amv_ref_yr <= np.percentile(amv_ref_yr, 33)
            hi_amv = amvc > np.percentile(amvc, 67)
            lo_amv = amvc <= np.percentile(amvc, 33)

            fig = plt.figure(figsize = (16,12))
            ax = fig.add_subplot(111)

            allfr = []
            for reg in range(4):
                #freq_ok, dates_ok = ctl.sel_season(results_ref['monthly_freq']['freq'][reg], results_ref['monthly_freq']['dates'], seas)
                freq_seas, dates_seas = ctl.seasonal_set(results[modmem]['monthly_freq']['freq'][reg], results[modmem]['monthly_freq']['dates'], seas, dates_range = ctl.range_years(y1, y2))
                freq_seas = np.mean(freq_seas, axis = 1)
                freq = np.array(ctl.running_mean(freq_seas, n_yr))
                oks = ~np.isnan(freq)
                freq = freq[oks]

                # freqs_hi = freq_seas[hi_amv]
                # freqs_lo = freq_seas[lo_amv]
                freqs_hi = freq[hi_amv]
                freqs_lo = freq[lo_amv]

                allfreqs_terc[(indexname, modmem, reg, seas)] = np.array([np.mean(freqs_hi), np.mean(freqs_lo)])

                allfr.append(freqs_hi)
                allfr.append(freqs_lo)

            pos = [1,2,4,5,7,8,10,11]

            posthicks = [1.5, 4.5, 7.5, 10.5]
            ax.boxplot(allfr, positions = np.array(pos))
            ax.set_xticks(posthicks)
            ax.set_xticklabels(reg_names)
            ax.set_ylabel('WR frequency')

            ax_all_ul[ise].boxplot(allfr, positions = np.array(pos))
            ax_all_ul[ise].set_xticks(posthicks)
            ax_all_ul[ise].set_xticklabels(reg_names)
            ax_all_ul[ise].set_ylabel('WR frequency')
            ax_all_ul[ise].set_title('season: {}'.format(seas))

            fig.suptitle('Upper vs lower tercile of {} index'.format(indexname))

            fig.savefig(cart_mod + '{}_{}_vs_WRfreq_upperlowertercile_{}.pdf'.format(modmem, indexname, seas))

        ctl.adjust_ax_scale(ax_all+ax_all_ul)
        fig_all.suptitle('Positive vs negative {} index'.format(indexname))
        fig_all.savefig(cart_mod + '{}_{}_vs_WRfreq_posneg_allseas.pdf'.format(modmem, indexname))

        fig_all_ul.suptitle('Upper vs lower tercile of {} index'.format(indexname))
        fig_all_ul.savefig(cart_mod + '{}_{}_vs_WRfreq_upperlowertercile_allseas.pdf'.format(modmem, indexname))

pickle.dump([allcorrs, allfreqs_posneg, allfreqs_terc], open(cart_out_all + 'corrfreq_AMV_ENSO.p', 'wb'))

[allcorrs, allfreqs_posneg, allfreqs_terc] = pickle.load(open(cart_out_all + 'corrfreq_AMV_ENSO.p', 'rb'))
[ref_corrs, ref_freqs_posneg, ref_freqs_terc] = pickle.load(open(cart_out_all + 'ref_corrfreq_AMV_ENSO.p', 'rb'))

sea_ok = 'DJFM'

for ind in indexes:
    print('Correlations with {}:\n'.format(ind))
    print(('{:25s}:'+4*' {:7.2f}').format('ERA', *[ref_corrs[(ind, reg, sea_ok)] for reg in range(4)]))
    for modmem in results.keys():
        if (ind, modmem, 0) in allcorrs.keys():
            print(('{:25s}:'+4*' {:7.2f}').format(modmem, *[allcorrs[(ind, modmem, reg, sea_ok)] for reg in range(4)]))

    print('\n')
    print('freq diff with {} pos/neg:\n'.format(ind))
    print(('{:25s}:'+4*' {:7.2f}').format('ERA', *[ref_freqs_posneg[(ind, reg, sea_ok)][0]-ref_freqs_posneg[(ind, reg, sea_ok)][1] for reg in range(4)]))
    for modmem in results.keys():
        if (ind, modmem, 0, sea_ok) in allcorrs.keys():
            print(('{:25s}:'+4*' {:7.2f}').format(modmem, *[allfreqs_posneg[(ind, modmem, reg, sea_ok)][0]-allfreqs_posneg[(ind, modmem, reg, sea_ok)][1] for reg in range(4)]))
    print('\n')

    print('freq diff with {} upper/lower tercile:\n'.format(ind))
    print(('{:25s}:'+4*' {:7.2f}').format('ERA', *[ref_freqs_terc[(ind, reg, sea_ok)][0]-ref_freqs_terc[(ind, reg, sea_ok)][1] for reg in range(4)]))
    for modmem in results.keys():
        if (ind, modmem, 0) in allcorrs.keys():
            print(('{:25s}:'+4*' {:7.2f}').format(modmem, *[allfreqs_terc[(ind, modmem, reg, sea_ok)][0]-allfreqs_terc[(ind, modmem, reg, sea_ok)][1] for reg in range(4)]))
    print('\n')
    print('\n')

########################################################################################
### PLOTS!

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

col_LR, col_HR = ['teal', 'indianred']
########################################################################################

for seas in seasons:
    for ind, reg_names in zip(indexes, reg_names_all):
        cou = cart_out_all + ind + '/' + seas + '/'
    # PLOT 1: correlazioni con gli indici
        fig = plt.figure(figsize = (16,12))
        axes = []
        for reg in range(4):
            modoks = []
            coloks = []

            ax = fig.add_subplot(2, 2, reg + 1)
            axes.append(ax)

            HRs = []
            LRs = []
            i = 0
            for mod, col, ver in zip(model_names, colors, vers):
                keoks = [cos for cos in allcorrs.keys() if cos[2] == reg and mod == cos[1].split('_')[0] and cos[0] == ind and cos[3] == seas]
                allvals = [allcorrs[ke] for ke in keoks]

                if len(allvals) == 0:
                    continue

                modoks.append(mod)
                coloks.append(col)
                meanval = np.mean(allvals)
                ax.scatter(i, meanval, marker = 'D', color = col, s = 25)
                if len(allvals) > 1:
                    ax.scatter([i]*len(allvals), allvals, marker = 'o', color = col, s = 5)

                i += 1
                if ver == 'LR':
                    LRs.append(meanval)
                elif ver == 'HR':
                    HRs.append(meanval)
                    i += 0.5

            i += 0.3
            ax.axvline(i, color = 'gray', linewidth = 0.5)
            ax.axhline(0, color = 'gray', linewidth = 0.5)
            i += 0.7

            ax.scatter(i, ref_corrs[(ind, reg, seas)], marker = 'D', color = 'black', s = 25)
            i+=1
            ax.scatter(i, np.mean(LRs), marker = 'D', color = col_LR, s = 25)
            i+=1
            ax.scatter(i, np.mean(HRs), marker = 'D', color = col_HR, s = 25)
            ax.set_xticks([])
            ax.set_title(reg_names[reg])


        ctl.adjust_ax_scale(axes)
        ctl.custom_legend(fig, coloks + ['black', col_LR, col_HR], modoks + ['ERA', 'LR', 'HR'])
        fig.suptitle('Correlation of WR freq. with {} index - {}'.format(ind, seas))

        fig.savefig(cou + 'corr_models_WRfreq_vs_{}_{}.pdf'.format(ind, seas))


    # PLOT 2: freq diff per ogni regime pos/neg
        fig = plt.figure(figsize = (16,12))
        axes = []
        for reg in range(4):
            modoks = []
            coloks = []

            ax = fig.add_subplot(2, 2, reg + 1)
            axes.append(ax)

            HRs = []
            LRs = []
            i = 0
            for mod, col, ver in zip(model_names, colors, vers):
                keoks = [cos for cos in allfreqs_posneg.keys() if cos[2] == reg and mod == cos[1].split('_')[0] and cos[0] == ind and cos[3] == seas]
                allvals = [allfreqs_posneg[ke][0]-allfreqs_posneg[ke][1] for ke in keoks]

                if len(allvals) == 0:
                    continue

                modoks.append(mod)
                coloks.append(col)
                meanval = np.mean(allvals)
                ax.scatter(i, meanval, marker = 'D', color = col, s = 25)
                if len(allvals) > 1:
                    ax.scatter([i]*len(allvals), allvals, marker = 'o', color = col, s = 5)

                i += 1
                if ver == 'LR':
                    LRs.append(meanval)
                elif ver == 'HR':
                    HRs.append(meanval)
                    i += 0.5

            i += 0.3
            ax.axvline(i, color = 'gray', linewidth = 0.5)
            ax.axhline(0, color = 'gray', linewidth = 0.5)
            i += 0.7

            ax.scatter(i, ref_freqs_posneg[(ind, reg, seas)][0]-ref_freqs_posneg[(ind, reg, seas)][1], marker = 'D', color = 'black', s = 25)
            i+=1
            ax.scatter(i, np.mean(LRs), marker = 'D', color = col_LR, s = 25)
            i+=1
            ax.scatter(i, np.mean(HRs), marker = 'D', color = col_HR, s = 25)
            ax.set_xticks([])
            ax.set_title(reg_names[reg])


        ctl.adjust_ax_scale(axes)
        ctl.custom_legend(fig, coloks + ['black', col_LR, col_HR], modoks + ['ERA', 'LR', 'HR'])
        fig.suptitle('WR freq. difference between positive and negative phase of {} - {}'.format(ind, seas))

        fig.savefig(cou + 'WRfreq_posneg_models_{}_{}.pdf'.format(ind, seas))
