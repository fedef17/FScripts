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

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 18

#############################################################################
cart_in = '/data-hobbes/fabiano/WR_CMIP6/'

yr10 = 10 # length of running mean
#dtrtyp = 'light'
dtrtyp = 'histrebase'
cart_out_orig = '/home/fabiano/Research/lavori/CMIP6/cmip6_vs_cmip5_refEOF/'
ctl.mkdir(cart_out_orig)

file_hist_refEOF = cart_in + 'out_NEW_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refEOF_dtr.p'
file_hist = cart_in + 'out_NEW_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refCLUS_dtr_light.p'

file_hist_refEOF_cmip5 = '/home/fabiano/Research/lavori/CMIP6/cmip5_hist_reb/out_cmip5_hist_reb_NDJFM_{}_4clus_4pcs_allyrs_refEOF_dtr.p'
file_hist_cmip5 = '/home/fabiano/Research/lavori/CMIP6/cmip5_hist_reb/out_cmip5_hist_reb_NDJFM_{}_4clus_4pcs_allyrs_refCLUS_dtr.p'

numclus = 4
reg_names_area = dict()
reg_names_area['EAT'] = ['NAO+', 'SBL', 'AR', 'NAO-']
reg_names_area['PNA'] = ['PT', 'PNA+', 'PNA-', 'AR']

clatlo = dict()
clatlo['EAT'] = (70., -20.)
clatlo['PNA'] = (70., -120.)

colormip = dict()
colormip['cmip5'] = ctl.color_set(6)[4]
colormip['cmip6'] = ctl.color_set(6)[0]

area = 'EAT'
for area in ['EAT']:#, 'PNA']:
    cart_out = cart_out_orig + '{}_NDJFM/'.format(area)
    ctl.mkdir(cart_out)

    results_hist, results_ref = ctl.load_wrtool(file_hist.format(area))
    results_hist_refEOF, _ = ctl.load_wrtool(file_hist_refEOF.format(area))
    results_hist_cmip5, _ = ctl.load_wrtool(file_hist_cmip5.format(area))
    results_hist_refEOF_cmip5, _ = ctl.load_wrtool(file_hist_refEOF_cmip5.format(area))

    resdict = dict()
    resdict['cmip6'] = results_hist
    resdict['cmip5'] = results_hist_cmip5
    resdict['cmip6_refEOF'] = results_hist_refEOF
    resdict['cmip5_refEOF'] = results_hist_refEOF_cmip5

    var_ratio = dict()
    freqbias = dict()
    for ke in resdict:
        var_ratio[ke] = np.array([resdict[ke][mod]['var_ratio'] for mod in resdict[ke].keys()])
        freqbias[ke] = np.array([np.sum(np.abs(resdict[ke][mod]['freq_clus']-results_ref['freq_clus'])) for mod in resdict[ke].keys()])

    for tip in ['', '_refEOF']:
        fig = plt.figure(figsize = (16,12))
        ax = fig.add_subplot(111)
        for cos in ['cmip5', 'cmip6']:
            ax.scatter(var_ratio[cos+tip], freqbias[cos+tip], label = cos, color = colormip[cos])

            varme = np.mean(var_ratio[cos+tip])
            varstd = np.std(var_ratio[cos+tip])
            frme = np.mean(freqbias[cos+tip])
            frstd = np.std(freqbias[cos+tip])
            ctl.ellipse_plot(varme, frme, varstd, frstd, colors = colormip[cos], ax = ax, alpha = 0.5)

        plt.legend()
        plt.xlabel('Variance ratio')
        plt.ylabel('Frequency bias')
        plt.title('cmip6 vs cmip6 var_fr '+tip)
        fig.savefig(cart_out + 'var_fr_plot_{}{}.pdf'.format(area, tip))


    # taylor with clouds
    patnames = reg_names_area[area]
    for tip in ['', '_refEOF']:
        fig = plt.figure(figsize=(16,12))
        for num, patt in enumerate(patnames):
            ax = plt.subplot(2, 2, num+1, polar = True)

            obs = results_ref['cluspattern_area'][num, ...]
            meapats = dict()
            for cos in ['cmip5', 'cmip6']:
                models = resdict[cos+tip].keys()
                modpats = [resdict[cos+tip][mod]['cluspattern_area'][num, ...] for mod in models]
                colors = [colormip[cos]]*len(modpats)
                ctl.Taylor_plot(modpats, obs, ax = ax, title = None, colors = colors, only_first_quarter = True, plot_ellipse = True, ellipse_color = colors[0])#, mod_points_size = taylor_mark_dim, obs_points_size = int(1.1*taylor_mark_dim), max_val_sd = max_val_sd)

        fig.savefig(cart_out + 'taylor_{}{}.pdf'.format(area, tip))
