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

cart_out_orig = cart_in + 'Results_hist1850/'
ctl.mkdir(cart_out_orig)

file_hist = cart_in + 'cmip6_hist/out_cmip6_hist_NDJFM_{}_4clus_4pcs_1850-2014_refCLUS_dtr.p'

numclus = 4
reg_names_area = dict()
reg_names_area['EAT'] = ['NAO+', 'SBL', 'AR', 'NAO-']
reg_names_area['PNA'] = ['PT', 'PNA+', 'PNA-', 'AR']


for area in ['EAT', 'PNA']:
    cart_out = cart_out_orig + '{}_NDJFM/'.format(area)
    ctl.mkdir(cart_out)

    results = pickle.load(open(file_hist.format(area), 'rb'))
    results_hist = results['models']
    results_ref = results['reference']

    # Erasing incomplete runs
    avlen = np.max([len(results_hist[ke]['labels']) for ke in results_hist.keys()])
    for ke in tuple(results_hist.keys()):
        if len(results_hist[ke]['labels']) < 0.9*avlen:
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
            cosi.append(seas20)
        coso = np.mean(cosi, axis = 0)
        runfreq[('hist', reg)] = coso
        ax.plot(yr, coso, color = 'black', linewidth = 3)
        ax.set_title(reg_names_area[area][reg])

    fig.savefig(cart_out + 'models_run20_{}_hist.pdf'.format(area))

    reg_names = reg_names_area[area]


    seasfr, yr_ref = ctl.calc_seasonal_clus_freq(results_ref['labels'], results_ref['dates'], numclus)
    for reg in range(4):
        seasfreq[('hist', 'ref', reg)] = seasfr[reg, :]

    yr = np.arange(1850, 2014)

    fig = plt.figure(figsize = (16,12))
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg+1)
        cosi = []
        for mem in okmods:
            seas20 = np.array(ctl.running_mean(seasfreq[('hist', mem, reg)], 20))
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


    fig = plt.figure(figsize = (16,12))
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg+1)
        cosi = []
        for mem in okmods:
            seas20 = np.array(ctl.running_mean(seasfreq[('hist', mem, reg)], 20))
            metot = np.mean(seasfreq[('hist', mem, reg)])
            cosi.append(seas20-metot)
        coso = np.mean(cosi, axis = 0)
        runfreq[('run20anom', reg)] = coso
        coserr = np.std(cosi, axis = 0)
        runfreq[('run20erranom', reg)] = coserr/np.sqrt(len(okmods)-1)
        ax.fill_between(yr, coso-coserr, coso+coserr, color = 'steelblue', alpha = 0.3)
        ax.plot(yr, coso, color = 'black', linewidth = 3)

        seas20ref = np.array(ctl.running_mean(seasfreq[('hist', 'ref', reg)], 20))-np.mean(seasfreq[('hist', 'ref', reg)])
        ax.plot(yr_ref, seas20ref, color = 'red', linewidth = 2, linestyle = '--')

        ax.set_title(reg_names_area[area][reg])

    fig.savefig(cart_out + 'long_run20anom_{}_hist.pdf'.format(area))
