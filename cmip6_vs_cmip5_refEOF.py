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
plt.rcParams['axes.titlesize'] = 28
plt.rcParams['axes.labelsize'] = 18

#############################################################################
cart_in = '/data-hobbes/fabiano/WR_CMIP6/'

yr10 = 10 # length of running mean
#dtrtyp = 'light'
dtrtyp = 'histrebase'
cart_out_orig = '/home/fabiano/Research/lavori/CMIP6/cmip6_vs_cmip5_refEOF/'
ctl.mkdir(cart_out_orig)

file_hist_refEOF = cart_in + 'out_NEW_cmip6_hist_NDJFM_{}_4clus_4pcs_1957-2005_refEOF_dtr.p'
file_hist = cart_in + 'out_NEW_cmip6_hist_NDJFM_{}_4clus_4pcs_1957-2005_refCLUS_dtr.p'

file_hist_refEOF_cmip5 = '/home/fabiano/Research/lavori/CMIP6/cmip5_hist_reb/out_cmip5_hist_reb_NDJFM_{}_4clus_4pcs_1957-2005_refEOF_dtr.p'
file_hist_cmip5 = '/home/fabiano/Research/lavori/CMIP6/cmip5_hist_reb/out_cmip5_hist_reb_NDJFM_{}_4clus_4pcs_1957-2005_refCLUS_dtr.p'

numclus = 4
reg_names_area = dict()
reg_names_area['EAT'] = ['NAO+', 'SBL', 'AR', 'NAO-']
reg_names_area['PNA'] = ['PT', 'PNA+', 'PNA-', 'BR']

clatlo = dict()
clatlo['EAT'] = (70., -20.)
clatlo['PNA'] = (70., -120.)

colormip = dict()
colormip[('cmip5', 'EAT')] = ctl.color_set(7)[4]
colormip[('cmip6', 'EAT')] = ctl.color_set(7)[0]
colormip[('cmip5', 'PNA')] = ctl.color_set(7)[2]
colormip[('cmip6', 'PNA')] = ctl.color_set(7)[6]

area = 'EAT'
plocos = dict()
for area in ['EAT', 'PNA']:
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
        plocos[('var_ratio', ke, area)] = var_ratio[ke]
        freqbias[ke] = np.array([np.mean(np.abs(resdict[ke][mod]['freq_clus']-results_ref['freq_clus'])) for mod in resdict[ke].keys()])
        plocos[('freqbias', ke, area)] = freqbias[ke]

    for tip in ['', '_refEOF']:
        fig = plt.figure(figsize = (16,12))
        ax = fig.add_subplot(111)
        for cos in ['cmip5', 'cmip6']:
            ax.scatter(var_ratio[cos+tip], freqbias[cos+tip], label = cos, color = colormip[(cos, area)], s = 50)

            varme = np.mean(var_ratio[cos+tip])
            varstd = np.std(var_ratio[cos+tip])
            frme = np.mean(freqbias[cos+tip])
            frstd = np.std(freqbias[cos+tip])
            ctl.ellipse_plot(varme, frme, varstd, frstd, colors = colormip[(cos, area)], ax = ax, alpha = 0.5)

        plocos[('var_ratio', 'ref', area)] = results_ref['var_ratio']
        ax.scatter(results_ref['var_ratio'], 0, s = 250, marker = '*', color = 'black', label = 'ERA')
        plt.legend(fontsize = 20)
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
            #obs = results_ref['centroids'][num, ...]
            plocos[('patterns', 'ref', area, num)] = obs
            meapats = dict()
            for cos in ['cmip5', 'cmip6']:
                models = resdict[cos+tip].keys()
                print(cos, models)
                modpats = [resdict[cos+tip][mod]['cluspattern_area'][num, ...] for mod in models]
                #modpats = [resdict[cos+tip][mod]['eff_centroids'][num, ...] for mod in models]
                plocos[('patterns', cos+tip, area, num)] = modpats
                colors = [colormip[(cos, area)]]*len(modpats)
                ctl.Taylor_plot(modpats, obs, ax = ax, title = patt, colors = colors, only_first_quarter = True, plot_ellipse = True, ellipse_color = colors[0])#, mod_points_size = taylor_mark_dim, obs_points_size = int(1.1*taylor_mark_dim), max_val_sd = max_val_sd)

        fig.savefig(cart_out + 'taylor_{}{}.pdf'.format(area, tip))

with open(cart_out_orig + 'histbiases.p', 'wb') as filo:
    pickle.dump(plocos, filo)

cart_out = cart_out_orig + 'EAT_and_PNA/'
ctl.mkdir(cart_out)

for tip in ['', '_refEOF']:
    fig = plt.figure(figsize = (24,12))
    for i, area in enumerate(['EAT', 'PNA']):
        ax = plt.subplot(1, 2, i+1)
        for cos in ['cmip5', 'cmip6']:
            var_ratio = plocos[('var_ratio', cos+tip, area)]
            freqbias = plocos[('freqbias', cos+tip, area)]
            ax.scatter(var_ratio, freqbias, label = cos, color = colormip[(cos, area)], s = 50)

            varme = np.mean(var_ratio)
            varstd = np.std(var_ratio)
            frme = np.mean(freqbias)
            frstd = np.std(freqbias)
            ctl.ellipse_plot(varme, frme, varstd, frstd, colors = colormip[(cos, area)], ax = ax, alpha = 0.5)

        ax.scatter(plocos[('var_ratio', 'ref', area)], 0, s = 250, marker = '*', color = 'black', label = 'ERA')
        ax.legend(fontsize = 20)
        ax.set_xlabel('Variance ratio')
        ax.set_ylabel('Frequency bias')
        ax.set_title(area)
    fig.savefig(cart_out + 'var_fr_plot_{}.pdf'.format(tip))


# taylor with clouds
for tip in ['', '_refEOF']:
    fig = plt.figure(figsize=(24,12))
    for i, area in enumerate(['EAT', 'PNA']):
        patnames = reg_names_area[area]
        for num, patt in enumerate(patnames):
            ax = plt.subplot(2, 4, num+1+4*i, polar = True)

            obs = plocos[('patterns', 'ref', area, num)]
            meapats = dict()
            for cos in ['cmip5', 'cmip6']:
                modpats = plocos[('patterns', cos+tip, area, num)]

                colors = [colormip[(cos, area)]]*len(modpats)
                ctl.Taylor_plot(modpats, obs, ax = ax, title = patt, colors = colors, only_first_quarter = True, plot_ellipse = True, ellipse_color = colors[0], max_val_sd = 1.6)

    ax.text(0.05, 0.75, 'EAT', horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 35)
    ax.text(0.05, 0.25, 'PNA', horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 35)
    fig.savefig(cart_out + 'taylor_{}.pdf'.format(tip))
