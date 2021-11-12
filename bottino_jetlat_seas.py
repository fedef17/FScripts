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
import xarray as xr
import glob
import xclim

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 28
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.axisbelow'] = True

#############################################################################

cart_out = '/home/fabiano/Research/lavori/BOTTINO/jetlat/'
ctl.mkdir(cart_out)

filna = '/nas/BOTTINO/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/{}/{}/day_r25/ua/ua*nc'

allru = ['pi', 'b025', 'b050', 'b100']
allnams = ['piControl', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']

colors = ['black', 'forestgreen', 'orange', 'violet']

allru2 = allru + ['c1950']
allnams2 = allnams + ['control-1950']
colors2 = colors + ['steelblue']

#############################################################################
areasel = dict()
areasel['NATL'] = [-60., 0., 20., 70.]
areasel['NEPAC'] = [200., 240., 20., 65.]
areasel['NCPAC'] = [150., 200., 20., 65.]
areasel['SPAC'] = [180., 240., -70., -20.]
areasel['SATL'] = [-45., 10., -70., -20.]
areasel['SIND'] = [50., 110., -70., -20.]

allseas = ['DJF', 'MAM', 'JJA', 'SON']

with open(cart_out + '../analisi/res_jli200_v2.p', 'rb') as filox:
    resdict = pickle.load(filox)

for season, lensea in zip(['DJF', 'JJA'], [90, 92]):
    for area in ['NATL', 'NEPAC', 'SPAC']:
        xose = dict()
        for ru in allru:
            gigi = resdict[(ru, 'jli', area, season)]
            seasme = gigi.reshape((-1, lensea)).mean(axis = 1)
            xose[ru] = seasme

        allpercs = dict()
        for nu in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(nu)] = [np.percentile(xose[ru], nu) for ru in allru]
        allpercs['mean'] = [np.mean(xose[ru]) for ru in allru]
        allpercs['min'] = [np.min(xose[ru]) for ru in allru]
        allpercs['max'] = [np.max(xose[ru]) for ru in allru]

        fig, ax = plt.subplots(figsize = (12,8))
        ctl.boxplot_on_ax(ax, allpercs, allru, colors, plot_ensmeans=False)
        ax.set_xticklabels(allru)
        ax.set_ylabel('{} {} mean jet lat'.format(area, season))
        fig.savefig(cart_out + 'seamean_jetlat_{}_{}.pdf'.format(area, season))

        xose2 = dict()
        for ru in allru:
            gigi = resdict[(ru, 'jli', area, season)]
            seasstd = gigi.reshape((-1, 90)).std(axis = 1)
            xose2[ru] = seasstd

        fig, ax = plt.subplots(figsize = (12,8))
        allpercs = dict()
        for nu in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(nu)] = [np.percentile(xose2[ru], nu) for ru in allru]
        allpercs['mean'] = [np.mean(xose2[ru]) for ru in allru]
        allpercs['min'] = [np.min(xose2[ru]) for ru in allru]
        allpercs['max'] = [np.max(xose2[ru]) for ru in allru]
        ctl.boxplot_on_ax(ax, allpercs, allru, colors, plot_ensmeans=False)
        ax.set_xticklabels(allru)
        ax.set_ylabel('{} {} intra-seasonal stddev of jet lat'.format(area, season))
        fig.savefig(cart_out + 'seastd_jetlat_{}_{}.pdf'.format(area, season))


        xose3 = dict()
        for ru in allru:
            gigi = resdict[(ru, 'jspeed', area, season)]
            seas = gigi.reshape((-1, 90)).mean(axis = 1)
            xose3[ru] = seas

        allpercs = dict()
        for nu in [10, 25, 50, 75, 90]:
            allpercs['p{}'.format(nu)] = [np.percentile(xose3[ru], nu) for ru in allru]
        allpercs['mean'] = [np.mean(xose3[ru]) for ru in allru]
        allpercs['min'] = [np.min(xose3[ru]) for ru in allru]
        allpercs['max'] = [np.max(xose3[ru]) for ru in allru]

        fig, ax = plt.subplots(figsize = (12,8))
        ctl.boxplot_on_ax(ax, allpercs, allru, colors, plot_ensmeans=False)
        ax.set_xticklabels(allru)
        ax.set_ylabel('{} {} mean jet speed'.format(area, season))
        fig.savefig(cart_out + 'seamean_jetspeed_{}_{}.pdf'.format(area, season))
