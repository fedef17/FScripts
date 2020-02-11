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

    if var_name is not None:
        var = fh.variables[var_name][:]
    else:
        print(fh.variables.keys())
        print('Returning the last variable')
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

#####################################################################################
indexes = ['AMV', 'ENSO']
file_refs = [file_amv_ref, file_enso_ref]
file_mods = [file_AMV_mod, file_ENSO_mod]
filogens = [filogen_EAT, filogen_PNA]
n_yrs = [5, 1]

reg_names_all = [['NAO+', 'SBL', 'AR', 'NAO-'], ['PT', 'PNA-', 'PNA+', 'AR']]

seasons = ['DJFM', 'DJ', 'FM', 'DJF']

ref_corrs = dict()
ref_freqs_posneg = dict()
ref_freqs_terc = dict()

for indexname, file_ref, file_mod, filogen, reg_names, n_yr in zip(indexes, file_refs, file_mods, filogens, reg_names_all, n_yrs):
    cart_out_ind = cart_out_all + indexname + '/'
    if not os.path.exists(cart_out_ind): os.mkdir(cart_out_ind)

    results, results_ref = pickle.load(open(filogen, 'rb'))

    amv_ref, dates = read_file_index(file_ref)
    amv_ref, dates = ctl.sel_time_range(amv_ref, dates, ctl.range_years(1957, 2013))
    # Yearly amv index
    amv_ref_yr, yrdates = ctl.yearly_average(amv_ref, dates)
    amv_ref_yr = np.squeeze(amv_ref_yr)
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
        y1 = 1958
        if seas == 'FM': y1 = 1959
        print(indexname, seas, y1)
        cart_out = cart_out_ind + seas + '/'
        if not os.path.exists(cart_out): os.mkdir(cart_out)

        # plt.ion()
        fig = plt.figure(figsize = (16,12))
        for reg in range(4):
            ax = fig.add_subplot(2, 2, reg+1)

            #freq_ok, dates_ok = ctl.sel_season(results_ref['monthly_freq']['freq'][reg], results_ref['monthly_freq']['dates'], seas)
            freq_seas, dates_seas = ctl.seasonal_set(results_ref['monthly_freq']['freq'][reg], results_ref['monthly_freq']['dates'], seas, dates_range = ctl.range_years(y1, 2014))
            freq_seas = np.mean(freq_seas, axis = 1)

            years = np.array([da.year for da in yrdates])
            yealen = np.arange(len(years))

            # freq = np.array(ctl.running_mean(freq_ok, 15))
            freq = np.array(ctl.running_mean(freq_seas, n_yr))
            oks = ~np.isnan(freq)
            freq = freq[oks]
            #amvc = np.array(ctl.running_mean(amv_ref_djf, 15))

            print(len(freq), len(years))
            years = years[oks]
            yealen = yealen[oks]

            rco = ctl.Rcorr(amvc, freq)
            ref_corrs[(indexname, reg, seas)] = rco

            ax.set_title('Corr {}: {:5.2f}'.format(reg, rco))
            ax.plot(yealen, freq, color = 'steelblue')

            ax2 = ax.twinx()
            ax2.plot(yealen, amvc, color = 'indianred')

            ax.set_xticks(yealen[2::15])
            ax.set_xticklabels(years[2::15])
            ax.set_xlabel('Years')
            ax.set_ylabel('WR frequency')

        fig.suptitle('Correlation of WR frequency with {} index'.format(indexname))
        fig.savefig(cart_out + 'ERA_{}_vs_WRfreq_{}.pdf'.format(indexname, seas))


        # hi_amv = amv_ref_yr > np.median(amv_ref_yr)
        # lo_amv = amv_ref_yr <= np.median(amv_ref_yr)
        # hi_amv = amvc > np.median(amvc)
        # lo_amv = amvc <= np.median(amvc)
        #
        # fig = plt.figure(figsize = (16,12))
        # ax = fig.add_subplot(111)
        #
        # allfr = []
        # for reg in range(4):
        #     #freq_ok, dates_ok = ctl.sel_season(results_ref['monthly_freq']['freq'][reg], results_ref['monthly_freq']['dates'], seas)
        #     freq_seas, dates_seas = ctl.seasonal_set(results_ref['monthly_freq']['freq'][reg], results_ref['monthly_freq']['dates'], seas, dates_range = ctl.range_years(y1, 2014))
        #     freq_seas = np.mean(freq_seas, axis = 1)
        #     freq = np.array(ctl.running_mean(freq_seas, n_yr))
        #     oks = ~np.isnan(freq)
        #     freq = freq[oks]
        #
        #     # freqs_hi = freq_seas[hi_amv]
        #     # freqs_lo = freq_seas[lo_amv]
        #     freqs_hi = freq[hi_amv]
        #     freqs_lo = freq[lo_amv]
        #
        #     allfr.append(freqs_hi)
        #     allfr.append(freqs_lo)
        #
        # pos = [1,2,4,5,7,8,10,11]
        #
        # ax.boxplot(allfr, positions = np.array(pos))
        # #ax.set_title('Corr: {:5.2f}'.format(rco))
        #
        # posthicks = [1.5, 4.5, 7.5, 10.5]
        # ax.set_xticks(posthicks)
        # ax.set_xticklabels(reg_names)
        # ax.set_ylabel('WR frequency')
        #
        # fig.suptitle('High vs low {}'.format(indexname))
        #
        # fig.savefig(cart_out + 'ERA_{}_vs_WRfreq_hilo_{}.pdf'.format(indexname, seas))


        # hi_amv = amv_ref_yr > 0.
        # lo_amv = amv_ref_yr <= 0.
        hi_amv = amvc > 0.
        lo_amv = amvc <= 0.

        fig = plt.figure(figsize = (16,12))
        ax = fig.add_subplot(111)

        allfr = []
        for reg in range(4):
            #freq_ok, dates_ok = ctl.sel_season(results_ref['monthly_freq']['freq'][reg], results_ref['monthly_freq']['dates'], seas)
            freq_seas, dates_seas = ctl.seasonal_set(results_ref['monthly_freq']['freq'][reg], results_ref['monthly_freq']['dates'], seas, dates_range = ctl.range_years(y1, 2014))
            freq_seas = np.mean(freq_seas, axis = 1)
            freq = np.array(ctl.running_mean(freq_seas, n_yr))
            oks = ~np.isnan(freq)
            freq = freq[oks]

            # freqs_hi = freq_seas[hi_amv]
            # freqs_lo = freq_seas[lo_amv]
            freqs_hi = freq[hi_amv]
            freqs_lo = freq[lo_amv]
            ref_freqs_posneg[(indexname, reg, seas)] = np.array([np.mean(freqs_hi), np.mean(freqs_lo)])

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

        fig.savefig(cart_out + 'ERA_{}_vs_WRfreq_posneg_{}.pdf'.format(indexname, seas))


        # hi_amv = amv_ref_yr > np.percentile(amv_ref_yr, 67)
        # lo_amv = amv_ref_yr <= np.percentile(amv_ref_yr, 33)
        hi_amv = amvc > np.percentile(amvc, 67)
        lo_amv = amvc <= np.percentile(amvc, 33)

        fig = plt.figure(figsize = (16,12))
        ax = fig.add_subplot(111)

        allfr = []
        for reg in range(4):
            #freq_ok, dates_ok = ctl.sel_season(results_ref['monthly_freq']['freq'][reg], results_ref['monthly_freq']['dates'], seas)
            freq_seas, dates_seas = ctl.seasonal_set(results_ref['monthly_freq']['freq'][reg], results_ref['monthly_freq']['dates'], seas, dates_range = ctl.range_years(y1, 2014))
            freq_seas = np.mean(freq_seas, axis = 1)
            freq = np.array(ctl.running_mean(freq_seas, n_yr))
            oks = ~np.isnan(freq)
            freq = freq[oks]

            # freqs_hi = freq_seas[hi_amv]
            # freqs_lo = freq_seas[lo_amv]
            freqs_hi = freq[hi_amv]
            freqs_lo = freq[lo_amv]
            ref_freqs_terc[(indexname, reg, seas)] = np.array([np.mean(freqs_hi), np.mean(freqs_lo)])

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

        fig.savefig(cart_out + 'ERA_{}_vs_WRfreq_upperlowertercile_{}.pdf'.format(indexname, seas))

    ctl.adjust_ax_scale(ax_all+ax_all_ul)
    fig_all.suptitle('Positive vs negative {} index'.format(indexname))
    fig_all.savefig(cart_out_ind + 'ERA_{}_vs_WRfreq_posneg_allseas.pdf'.format(indexname))

    fig_all_ul.suptitle('Upper vs lower tercile of {} index'.format(indexname))
    fig_all_ul.savefig(cart_out_ind + 'ERA_{}_vs_WRfreq_upperlowertercile_allseas.pdf'.format(indexname))

pickle.dump([ref_corrs, ref_freqs_posneg, ref_freqs_terc], open(cart_out_all + 'ref_corrfreq_AMV_ENSO.p', 'wb'))
