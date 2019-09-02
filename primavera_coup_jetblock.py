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

cart = '/home/fabiano/Research/lavori/WeatherRegimes/prima_coup_v4/'
filogen = cart + 'out_prima_coup_v4_DJF_EAT_4clus_4pcs_1957-2014_refEOF.p'

cart_jet = '/data-hobbes/fabiano/PRIMAVERA/jetlat_panos/hist-1950/'
fil_jet = 'JetLatDist_850hPa_{}_hist-1950_DJF.npy'

cart_bloc = '/data-hobbes/fabiano/PRIMAVERA/Reinhard_blocking/all_hist-1950/'
fil_bloc = 'bf.5day.daily.daily_mean.{}.{}.hist-1950.r1i1p1f1.v20170915.{}-{}.nc'

model_names = ['CMCC-CM2-HR4', 'CMCC-CM2-VHR4', 'CNRM-CM6-1', 'CNRM-CM6-1-HR', 'EC-Earth3P', 'EC-Earth3P-HR', 'ECMWF-IFS-LR', 'ECMWF-IFS-HR', 'MPI-ESM1-2-HR', 'MPI-ESM1-2-XR', 'HadGEM3-GC31-MM', 'HadGEM3-GC31-HM', 'HadGEM3-GC31-LL-det', 'HadGEM3-GC31-LL-stoc', 'EC-Earth3P-det', 'EC-Earth3P-stoc']
model_coups = ['CMCC-CM2', 'CNRM-CM6-1', 'EC-Earth3P', 'ECMWF-IFS', 'MPI-ESM1-2', 'HadGEM3-GC31', 'HadGEM3-GC31 (det vs stoc)', 'EC-Earth3P (det vs stoc)']
model_names_all = model_names + ['ERA']

colors = ctl.color_set(len(model_names), sns_palette = 'Paired')
colors_wERA = colors + ['black']

regnam = ['NAO +', 'Sc. BL', 'AR', 'NAO -']

################################################################################

results, results_ref = pickle.load(open(filogen, 'rb'))
results['ERA'] = results_ref

latgri = np.arange(30,70,1)

jetlat_comp = dict()
# first for the jet
fig = plt.figure(figsize = (16, 12))

for mod, col in zip(model_names_all, colors_wERA):
    print(mod)
    fil_ok = fil_jet.format(mod)
    if mod == 'HadGEM3-GC31-LL-stoc':
        fil_ok = fil_jet.format('HadGEM3-GC31-LL')
    elif mod == 'ERA':
        fil_ok = 'JetLatDist_850hPa_NCEP_DJF.npy'
    elif mod == 'EC-Earth3P':
        print('Missing 1991')
        continue

    if not os.path.exists(cart_jet + fil_ok):
        print('waiting for panos..')
        continue

    jetind = np.load(cart_jet + fil_ok)[0]
    jetlat_comp[(mod, 'all')] = jetind

    regind = results[mod]['labels']
    dates = results[mod]['dates']
    dates_pdh = pd.to_datetime(dates)
    okdat = ~((dates_pdh.month == 2) & (dates_pdh.day == 29))
    print(len(regind) - np.sum(okdat))
    regind = regind[okdat]
    dates = dates[okdat]

    data1 = pd.to_datetime('{}1201'.format(1957), format='%Y%m%d')
    data2 = pd.to_datetime('{}0228'.format(2014), format='%Y%m%d')
    regok, datok = ctl.sel_time_range(regind, dates, (data1, data2))

    print(len(regok), len(jetind))
    if len(regok) != 5130:
        print('SKIPPING')
        continue
        #raise ValueError('OPS!')

    jetind = jetind[630:]
    print(len(regok), len(jetind))

    for reg in range(4):
        okin = regok == reg
        okjet = jetind[okin]
        jetlat_comp[(mod, reg)] = okjet

        pdf = ctl.calc_pdf(okjet)

        ax = plt.subplot(2, 2, reg+1)
        ax.plot(latgri, pdf(latgri), color = col)

fig.savefig(cart_out + 'jetlat_compos.pdf')

pickle.dump(jetlat_comp, open(cart_out + 'jetlat_compos_all.p', 'wb'))
