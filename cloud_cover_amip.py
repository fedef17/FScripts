#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os

from matplotlib import pyplot as plt
import netCDF4 as nc
import pandas as pd
from numpy import linalg as LA
from scipy import stats
from scipy import interpolate as itrp
import itertools as itt

from datetime import datetime
import pickle

import climtools_lib as ctl
import climdiags as cd

from copy import deepcopy as cp
###############################################

# Calculates radiation balances
cart_in = '/data-hobbes/fabiano/SPHINX/AMIP/cloud_cover/'
cart_out = '/home/fabiano/Research/lavori/SPHINX_for_lisboa/cloud_cover_AMIP/'
if not os.path.exists(cart_out): os.mkdir(cart_out)

varnames = ['hcc', 'lcc', 'mcc', 'tcw']
namefi = '{}_mon_{}_{}.nc'

ann = np.arange(1979,2008,6)
annme = [(a1+a2)//2 for a1,a2 in zip(ann[:-1], ann[1:])]
print(annme)

radclim = dict()
ensmems = ['labm', 'lasm']
for exp in ensmems:
    for a1, a2 in zip(ann[:-1], ann[1:]):
        am = (a2+a1)/2
        vardict = dict()
        for year in range(a1+1,a2+1):
            for varna in varnames:
                varniuu, lat, lon, dates, time_units, var_units = ctl.read3Dncfield(cart_in+namefi.format(exp,year,varna))
                if year == a1+1:
                    vardict[varna] = [np.mean(varniuu, axis = 0)]
                else:
                    vardict[varna].append(np.mean(varniuu, axis = 0))

        for varna in varnames:
            vardict[varna] = np.stack(vardict[varna])

        for key in vardict:
            print(vardict[key].shape)
            radclim[(exp, 'map', key, am)] = np.mean(vardict[key], axis = 0)
            radclim[(exp, 'map_std', key, am)] = np.std(vardict[key], axis = 0)
            radclim[(exp, 'zonal', key, am)] = np.mean(radclim[(exp, 'map', key, am)], axis = -1)
            radclim[(exp, 'zonal_std', key, am)] = np.mean(radclim[(exp, 'map_std', key, am)], axis = -1)
            radclim[(exp, 'global', key, am)] = ctl.global_mean(radclim[(exp, 'map', key, am)], lat)
            radclim[(exp, 'global_std', key, am)] = ctl.global_mean(radclim[(exp, 'map_std', key, am)], lat)


pickle.dump(radclim, open(cart_out+'cloudcover_allens.p', 'wb'))
# radclim = pickle.load(open(cart_out+'cloudcover_allens.p'))
# varniuu, lat, lon, dates, time_units, var_units = ctl.read3Dncfield(cart_in+namefi.format('lcb0',1988,'hcc'))
# del varniuu

# figure
#voglio figura con globalmean base e stoc anno per anno (con err da std? forse)
titlevar = dict()
titlevar['hcc'] = 'High clouds cover'
titlevar['lcc'] = 'Low clouds cover'
titlevar['mcc'] = 'Mid clouds cover'
titlevar['tcw'] = 'Total column water'

#figure_file = cart_out+'rad_forcing_TOA_lcb0_vs_lcs0.pdf'
figure_file = cart_out+'cloud_cover_TOA_base_vs_stoc.pdf'
for varna in titlevar:
    figures = []
    figure_file = cart_out+'cloud_{}_base_vs_stoc.pdf'.format(varna)

    stoc_glob = np.stack([radclim[('lasm', 'global', varna, ye)] for ye in annme])
    stoc_glob_err = [radclim[('lasm', 'global_std', varna, ye)] for ye in annme]
    base_glob = np.stack([radclim[('labm', 'global', varna, ye)] for ye in annme])
    base_glob_err = [radclim[('labm', 'global_std', varna, ye)] for ye in annme]

    fig = plt.figure()
    plt.title(titlevar[varna])
    for ens in ensmems[:3]:
        coso = np.stack([radclim[(ens, 'global', varna, ye)] for ye in annme])
        plt.plot(annme, coso, label = None, color = 'grey', linestyle = '-', linewidth = 0.7)
    for ens in ensmems[3:]:
        coso = np.stack([radclim[(ens, 'global', varna, ye)] for ye in annme])
        plt.plot(annme, coso, label = None, color = 'grey', linestyle = '--', linewidth = 0.7)
    plt.plot(annme, base_glob, label = 'base', linewidth = 2.0)
    plt.plot(annme, stoc_glob, label = 'stoc', linewidth = 2.0)
    plt.xlabel('Year')
    plt.ylabel('Cloud cover')
    plt.legend()
    figures.append(fig)

    fig = plt.figure()
    plt.title(titlevar[varna])
    plt.plot(annme, stoc_glob-base_glob, label = 'diff stoc-base')
    plt.xlabel('Year')
    plt.ylabel('Cloud cover')
    plt.grid()
    plt.legend()
    figures.append(fig)

    for ann in annme:
        print(varna, ann)
        base = radclim[('labm', 'map', varna, ann)]
        stoc = radclim[('lasm', 'map', varna, ann)]
        print(base.shape, stoc.shape)
        figures.append(ctl.plot_map_contour(stoc-base, lat, lon, title = titlevar[varna]+' (stoc-base diff): {}-{}'.format(ann-5, ann+5), cb_label = 'cloud cover'))#, cbar_range = (-20.,20.)))

    mino = np.min([radclim[('lasm', 'zonal', varna, ann)] for ann in annme])
    maxo = np.max([radclim[('labm', 'zonal', varna, ann)] for ann in annme])
    mino = mino - 0.1*abs(maxo-mino)
    maxo = maxo + 0.1*abs(maxo-mino)

    mino_diff = 10.1*np.min([radclim[('lasm', 'zonal', varna, ann)]-radclim[('labm', 'zonal', varna, ann)] for ann in annme])
    maxo_diff = 10.1*np.max([radclim[('lasm', 'zonal', varna, ann)]-radclim[('labm', 'zonal', varna, ann)] for ann in annme])
    mino = np.min([mino, mino_diff])
    maxo = np.max([maxo, maxo_diff])

    for ann in annme:
        print(varna, ann)
        fig = plt.figure()
        stoc = radclim[('lasm', 'zonal', varna, ann)]
        base = radclim[('labm', 'zonal', varna, ann)]
        plt.title('Zonal '+titlevar[varna]+': {}-{}'.format(ann-5, ann+5))
        plt.plot(lat, base, label = 'base', linewidth = 2)
        plt.plot(lat, stoc, label = 'stoc', linewidth = 2)
        plt.plot(lat, 10*(stoc-base), label = 'diff x 10', linestyle = '--', linewidth = 1)
        plt.ylim(mino, maxo)
        plt.grid()
        plt.xlabel('Latitude')
        plt.ylabel('cloud cover')
        plt.legend()
        figures.append(fig)

    ctl.plot_pdfpages(figure_file, figures)
    plt.close('all')
