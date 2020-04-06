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
if os.uname()[1] == 'hobbes':
    cart_in = '/home/fabiano/Research/lavori/CMIP6/'
elif os.uname()[1] == 'ff-clevo':
    cart_in = '/home/fedefab/Scrivania/Research/Post-doc/lavori/CMIP6/'

cart_out_orig = cart_in + 'Results_cmip5/'
ctl.mkdir(cart_out_orig)

#file_hist_refEOF = cart_in + 'cmip6_hist/out_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refEOF_dtr.p'
file_hist = cart_in + 'cmip5_hist/out_cmip5_hist_NDJFM_{}_4clus_4pcs_allyrs_refCLUS_dtr.p'
#file_hist_refEOF = cart_in + 'cmip6_hist/out_cmip6_hist_NDJFM_EAT_4clus_4pcs_1964-2014_refEOF.p'
gen_file_ssp = cart_in + 'cmip5_{}/out_cmip5_{}_NDJFM_{}_4clus_4pcs_2005-2100_refCLUS_dtr.p'

numclus = 4
reg_names_area = dict()
reg_names_area['EAT'] = ['NAO+', 'SBL', 'AR', 'NAO-']
reg_names_area['PNA'] = ['PT', 'PNA+', 'PNA-', 'AR']


#allssps = 'ssp119 ssp126 ssp245 ssp370 ssp585'.split()
#allssps = 'ssp126 ssp245 ssp370 ssp585'.split()
allssps = ['rcp85']

for area in ['EAT']:#, 'PNA']:
    cart_out = cart_out_orig + '{}_NDJFM/'.format(area)
    ctl.mkdir(cart_out)

    results = pickle.load(open(file_hist.format(area), 'rb'))
    results_hist = results['models']
    results_ref = results['reference']

    yr = np.arange(1950, 2005)

    # Erasing incomplete runs
    avlen = np.max([len(results_hist[ke]['labels']) for ke in results_hist.keys()])
    for ke in tuple(results_hist.keys()):
        if len(results_hist[ke]['labels']) < avlen-50:
            del results_hist[ke]

    okmods = [cos for cos in results_hist.keys()]
    print(okmods)
    print(len(okmods))

    runfreq = dict()
    seasfreq = dict()

    fig = plt.figure(figsize = (16,12))
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg+1)
        cosi = []
        for mem in okmods:
            seasfr, yr = ctl.calc_seasonal_clus_freq(results_hist[mem]['labels'], results_hist[mem]['dates'], numclus)
            seasfreq[('hist', mem, reg)] = seasfr[reg, :]
            seas20 = np.array(ctl.running_mean(seasfr[reg, :], 20))
            ax.plot(yr, seas20)
            if len(seas20) == len(yr):
                cosi.append(seas20)
            else:
                print(mem, 'too short')
        coso = np.mean(cosi, axis = 0)
        runfreq[('hist', reg)] = coso
        ax.plot(yr, coso, color = 'black', linewidth = 3)
        ax.set_title(reg_names_area[area][reg])

    fig.savefig(cart_out + 'models_run20_{}_hist.pdf'.format(area))

    reg_names = reg_names_area[area]

    seasfr, yr_ref = ctl.calc_seasonal_clus_freq(results_ref['labels'], results_ref['dates'], numclus)
    for reg in range(4):
        seasfreq[('hist', 'ref', reg)] = seasfr[reg, :]

    fig = plt.figure(figsize = (16,12))
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg+1)
        cosi = []
        for mem in okmods:
            seas20 = np.array(ctl.running_mean(seasfreq[('hist', mem, reg)], 20))
            if len(seas20) == len(yr):
                cosi.append(seas20)
            else:
                print(mem, 'too short')
            cosi.append(seas20)
        coso = np.mean(cosi, axis = 0)
        runfreq[('run20', reg)] = coso
        coserr = np.std(cosi, axis = 0)
        runfreq[('run20err', reg)] = coserr/np.sqrt(len(okmods)-1)
        ax.fill_between(yr, coso-coserr, coso+coserr, color = 'steelblue', alpha = 0.3)
        ax.plot(yr, coso, color = 'black', linewidth = 3)

        seas20ref = np.array(ctl.running_mean(seasfreq[('hist', 'ref', reg)], 20))
        ax.plot(yr_ref, seas20ref, color = 'red', linewidth = 2, linestyle = '--')

        ax.set_title(reg_names_area[area][reg])

    fig.savefig(cart_out + 'long_run20_{}_hist.pdf'.format(area))

    # now for rcp85
    print('rcp85')
    results = pickle.load(open(gen_file_ssp.format('rcp85', 'rcp85', area), 'rb'))
    results_ssp = results['models']

    okmods = [cos for cos in results_ssp.keys()]
    yr = np.arange(2005, 2100)

    fig = plt.figure(figsize = (16,12))
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg+1)
        cosi = []
        for mem in okmods:
            seasfr, yr = ctl.calc_seasonal_clus_freq(results_ssp[mem]['labels'], results_ssp[mem]['dates'], numclus)
            seasfreq[('rcp85', mem, reg)] = seasfr[reg, :]
            seas20 = np.array(ctl.running_mean(seasfr[reg, :], 20))
            ax.plot(yr, seas20)
            if len(seas20) == len(yr):
                cosi.append(seas20)
            else:
                print(mem, 'too short')
            cosi.append(seas20)
        coso = np.mean(cosi, axis = 0)
        runfreq[('rcp85', reg)] = coso
        ax.plot(yr, coso, color = 'black', linewidth = 3)
        ax.set_title(reg_names_area[area][reg])

    fig.savefig(cart_out + 'models_run20_{}_rcp85.pdf'.format(area))

    reg_names = reg_names_area[area]

    fig = plt.figure(figsize = (16,12))
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg+1)
        cosi = []
        for mem in okmods:
            seas20 = np.array(ctl.running_mean(seasfreq[('rcp85', mem, reg)], 20))
            if len(seas20) == len(yr):
                cosi.append(seas20)
            else:
                print(mem, 'too short')
            cosi.append(seas20)
        coso = np.mean(cosi, axis = 0)
        runfreq[('run20', reg)] = coso
        coserr = np.std(cosi, axis = 0)
        runfreq[('run20err', reg)] = coserr/np.sqrt(len(okmods)-1)
        ax.fill_between(yr, coso-coserr, coso+coserr, color = 'steelblue', alpha = 0.3)
        ax.plot(yr, coso, color = 'black', linewidth = 3)

        ax.set_title(reg_names_area[area][reg])

    fig.savefig(cart_out + 'long_run20_{}_rcp85.pdf'.format(area))
