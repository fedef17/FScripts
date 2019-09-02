#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
from matplotlib import pyplot as plt
import pickle

import climtools_lib as ctl
import climdiags as cd

#######################################

cart_out = '/home/fabiano/Research/lavori/WeatherRegimes/primavera_coupled_1957-2014/'
[tutti, ref] = pickle.load(open(cart_out + 'out_primavera_coupled_1957-2014_DJF_EAT_4clus_4pcs_1957-2014_dtr.p'))
groups = dict()
groups['HR'] = 'CMCC-CM2-VHR4,EC-Earth-3-HR,ECMWF-IFS-HR,HadGEM3-GC31-HM,MPI-ESM1-2-XR'.split(',')
groups['LR'] = 'CMCC-CM2-HR4,EC-Earth-3-LR,ECMWF-IFS-LR,HadGEM3-GC31-MM,MPI-ESM1-2-HR'.split(',')

nsqr = np.sqrt(ref['cluspattern_area'].size)

data = dict()
for cos in ['RMS', 'patcor']:
    fac = 1.
    if cos == 'RMS': fac = nsqr
    data[cos] = []
    data[cos+'_err'] = []
    for grp in ['LR', 'HR']:
        data[cos].append(np.mean([tutti[mod][cos] for mod in groups[grp]], axis = 0)/fac)
        data[cos+'_err'].append(np.std([tutti[mod][cos] for mod in groups[grp]], axis = 0)/(fac * np.sqrt(len(groups[grp])-1)))
    data[cos] = np.stack(data[cos])
    data[cos+'_err'] = np.stack(data[cos+'_err'])

# rms_HR = [13.17, 12.02, 12.75, 13.31]
# patcor_HR = [0.79, 0.84, 0.82, 0.89]
# err_rms_HR = [4.79, 5.33, 6.88, 7.64]
# err_patcor_HR = [0.14, 0.14, 0.19, 0.15]
#
# rms_LR = [16.11, 14.05, 14.33, 14.20]
# err_rms_LR = [6.80, 4.04, 5.99, 7.69]
# patcor_LR = [0.69, 0.79, 0.79, 0.87]
# err_patcor_LR = [0.23, 0.13, 0.21, 0.16]

patt = ['NAO +', 'Sc. Blocking', 'Atl. Ridge', 'NAO -']

fig = plt.figure(figsize = (16,12))

for i in range(4):
    ax = fig.add_subplot(2, 2, i+1)
    ax.set_title(patt[i], fontsize = 18, fontweight = 'bold')

    # x = [patcor_LR[i], patcor_HR[i]]
    # y = [rms_LR[i], rms_HR[i]]
    # errx = [err_patcor_LR[i], err_patcor_HR[i]]
    # erry = [err_rms_LR[i], err_rms_HR[i]]
    ctl.ellipse_plot(data['patcor'][:,i], data['RMS'][:,i], data['patcor_err'][:,i], data['RMS_err'][:,i], labels = ['LR', 'HR'], ax = ax, colors = ['indianred', 'steelblue'], alpha = 0.7)

    for col, grp, sim in zip(['indianred', 'steelblue'], ['LR', 'HR'], ['$L$', '$H$']):
        pats = [tutti[mod]['patcor'][i] for mod in groups[grp]]
        rmss = [tutti[mod]['RMS'][i]/nsqr for mod in groups[grp]]
        ax.scatter(pats, rmss, color = col, s = 25, marker = sim)

    ax.set_xlim(0.35, 1.0)
    ax.set_ylim(0., 27.0)
    ax.tick_params(labelsize=14)
    plt.gca().invert_xaxis()
    ax.set_xlabel('Pattern correlation', fontsize = 18)
    ax.set_ylabel('RMS (m)', fontsize = 18)
    ax.grid()

plt.tight_layout()
plt.subplots_adjust(top = 0.9)
plt.suptitle('Average performance of PRIMAVERA stream1 coupled models', fontsize = 28)
fig.savefig(cart_out + 'ellipse_plot.pdf')
