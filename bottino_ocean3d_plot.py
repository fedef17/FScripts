#!/usr/bin/python
# -*- coding: utf-8 -*-


import numpy as np
import sys
import os
from matplotlib import pyplot as plt
from matplotlib import cm

import pickle
#import netCDF4 as nc

import climtools_lib as ctl
#import climdiags as cd
#from tunlib import gregplot_on_ax

#from matplotlib.colors import LogNorm
#from datetime import datetime

#from scipy import stats
import xarray as xr
import glob
#import xclim

import multiprocessing as mp
import psutil

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 22
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['axes.axisbelow'] = True

#############################################################################

if os.uname()[1] == 'hobbes':
    cart_out = '/home/fabiano/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'xaru':
    cart_out = '/home/fedef/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'tintin':
    cart_out = '/home/fabiano/work/lavori/BOTTINO/'

#cart_out = '/g100_work/IscrB_QUECLIM/BOTTINO/bottino_an/'
cart_out = cart_out + 'ocean3d/'
ctl.mkdir(cart_out)


allru = ['b025', 'b050', 'b100']
allsyear = [2030, 2050, 2100]
allnams = ['stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']

colors = ['forestgreen', 'orange', 'violet']


ru = 'b100'
for ru in allru:
    syear = allsyear[allru.index(ru)]

    filo = open(cart_out + 'thetao_{}.p'.format(ru), 'rb')
    zuk = []
    for ii in range(500):
        try:
            pitu = pickle.load(filo)
        except EOFError:
            print('Only {} entries found'.format(ii))
            break

        zuk.append(pitu.assign_coords(year = syear + ii))

    thetao = xr.concat(zuk, 'year')

    trends = dict()
    for nam in thetao.data_vars:
        print(nam)
        coso = thetao[nam]

        trendmat = np.empty_like(coso.values[0])
        errtrendmat = np.empty_like(coso.values[0])

        for i in np.arange(trendmat.shape[0]):
            for j in np.arange(trendmat.shape[1]):
                m, c, err_m, err_c = ctl.linear_regre_witherr(coso.year.values, coso.values[:,i,j])
                #coeffs, covmat = np.polyfit(years, var_set[i,j], deg = deg, cov = True)
                trendmat[i,j] = m
                errtrendmat[i,j] = err_m

        trends[(ru, nam)] = 100*trendmat

        coso = thetao[nam].sel(year = slice(syear + 300, syear + 500))

        trendmat = np.empty_like(coso[0].values)
        errtrendmat = np.empty_like(coso[0].values)

        for i in np.arange(trendmat.shape[0]):
            for j in np.arange(trendmat.shape[1]):
                m, c, err_m, err_c = ctl.linear_regre_witherr(coso.year.values, coso.values[:,i,j])
                #coeffs, covmat = np.polyfit(years, var_set[i,j], deg = deg, cov = True)
                trendmat[i,j] = m
                errtrendmat[i,j] = err_m

        trends[(ru, nam, 'last200')] = 100*trendmat

        coso = thetao[nam].sel(year = slice(syear, syear + 200))

        trendmat = np.empty_like(coso[0].values)
        errtrendmat = np.empty_like(coso[0].values)

        for i in np.arange(trendmat.shape[0]):
            for j in np.arange(trendmat.shape[1]):
                m, c, err_m, err_c = ctl.linear_regre_witherr(coso.year.values, coso.values[:,i,j])
                #coeffs, covmat = np.polyfit(years, var_set[i,j], deg = deg, cov = True)
                trendmat[i,j] = m
                errtrendmat[i,j] = err_m

        trends[(ru, nam, 'first200')] = 100*trendmat

    ############################################################
    #### Plots ####

    cbran = (0., 1)
    top    = 0.88  # the top of the subplots
    bottom = 0.20    # the bottom of the subplots
    left   = 0.1    # the left side
    right  = 0.98  # the right side
    hspace = 0.20   # height reserved for white space
    wspace = 0.05    # width reserved for blank space

    nams = ['thetao_atl', 'thetao_ind', 'thetao_pac']
    tits = ['Atlantic', 'Indian', 'Pacific']

    fig, axes = plt.subplots(1,3, figsize = (24,9))
    for nam, tit, ax in zip(nams, tits, axes.flatten()):
        map_plot = ctl.plot_lat_crosssection(trends[(ru, nam)], thetao.lat, thetao.lev, ax = ax, set_logscale_levels=False, cmap = 'viridis', ylim = (None, 10), ylabel='Depth (m)', xlabel = 'Latitude', cb_label='Trend (K/year)', cbar_range = cbran)
        ax.set_title(tit)
        if tits.index(tit) > 0:
            ax.set_yticklabels([])
            ax.set_ylabel('')

    cax = plt.axes([0.1, 0.11, 0.8, 0.05]) #horizontal
    cb = plt.colorbar(map_plot, cax=cax, orientation='horizontal')#, labelsize=18)
    cb.ax.tick_params(labelsize=14)
    cb.set_label('Trend (K/cent)', fontsize=16)
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

    fig.savefig(cart_out + 'oce_trends_{}.pdf'.format(ru))

    fig, axes = plt.subplots(1,3, figsize = (24,9))
    for nam, tit, ax in zip(nams, tits, axes.flatten()):
        map_plot = ctl.plot_lat_crosssection(trends[(ru, nam, 'last200')], thetao.lat, thetao.lev, ax = ax, set_logscale_levels=False, cmap = 'viridis', ylim = (None, 10), ylabel='Depth (m)', xlabel = 'Latitude', cb_label='Trend (K/year)', cbar_range = cbran)
        ax.set_title(tit)
        if tits.index(tit) > 0:
            ax.set_yticklabels([])
            ax.set_ylabel('')

    cax = plt.axes([0.1, 0.11, 0.8, 0.05]) #horizontal
    cb = plt.colorbar(map_plot, cax=cax, orientation='horizontal')#,
    cb.ax.tick_params(labelsize=14)
    cb.set_label('Trend (K/cent)', fontsize=16)
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

    fig.savefig(cart_out + 'oce_trends_{}_last200.pdf'.format(ru))

    cbran = (0., 2)

    fig, axes = plt.subplots(1,3, figsize = (24,9))
    for nam, tit, ax in zip(nams, tits, axes.flatten()):
        map_plot = ctl.plot_lat_crosssection(trends[(ru, nam, 'first200')], thetao.lat, thetao.lev, ax = ax, set_logscale_levels=False, cmap = 'viridis', ylim = (None, 10), ylabel='Depth (m)', xlabel = 'Latitude', cb_label='Trend (K/year)', cbar_range = cbran)
        ax.set_title(tit)
        if tits.index(tit) > 0:
            ax.set_yticklabels([])
            ax.set_ylabel('')

    cax = plt.axes([0.1, 0.11, 0.8, 0.05]) #horizontal
    cb = plt.colorbar(map_plot, cax=cax, orientation='horizontal')#,
    cb.ax.tick_params(labelsize=14)
    cb.set_label('Trend (K/cent)', fontsize=16)
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

    fig.savefig(cart_out + 'oce_trends_{}_first200.pdf'.format(ru))

    cbran = (-0.3, 0.3)

    fig, axes = plt.subplots(1,3, figsize = (24,9))
    for nam, tit, ax in zip(nams, tits, axes.flatten()):
        map_plot = ctl.plot_lat_crosssection((trends[(ru, nam, 'last200')]-trends[(ru, nam, 'first200')]), thetao.lat, thetao.lev, ax = ax, cmap = 'RdBu_r', ylabel='Depth (m)', xlabel = 'Latitude', cb_label='Trend diff (K/year)', cbar_range = cbran)
        ax.set_title(tit)
        if tits.index(tit) > 0:
            ax.set_yticklabels([])
            ax.set_ylabel('')

    cax = plt.axes([0.1, 0.11, 0.8, 0.05]) #horizontal
    cb = plt.colorbar(map_plot, cax=cax, orientation='horizontal')#,
    cb.ax.tick_params(labelsize=14)
    cb.set_label('Trend diff (K/cent)', fontsize=16)
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

    fig.savefig(cart_out + 'oce_trends_{}_diff.pdf'.format(ru))


for ke in trends:
    print(ke, 'points with negative trend: ', np.sum(trends[ke] < 0.))


# AAAAAAAAAAAAAAAAAAAAA
# fig, axes = plt.subplots(2, 3, figsize = (16,9))
# for nam, tit, ax in zip(nams, tits, axes.flatten()):
#     map_plot = ctl.plot_lat_crosssection((trends[(ru, nam, 'last200')]-trends[(ru, nam, 'first200')]), thetao.lat, thetao.lev, ax = ax, cmap = 'RdBu_r', ylabel='Depth (m)', xlabel = 'Latitude', cb_label='Trend diff (K/year)', cbar_range = cbran)
#     ax.set_title(tit)
#     if tits.index(tit) > 0:
#         ax.set_yticklabels([])
#         ax.set_ylabel('')
#
# cax = plt.axes([0.1, 0.11, 0.8, 0.05]) #horizontal
# cb = plt.colorbar(map_plot, cax=cax, orientation='horizontal')#,
# cb.ax.tick_params(labelsize=14)
# cb.set_label('Trend diff (K/cent)', fontsize=16)
# plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
#
# fig.savefig(cart_out + 'oce_trends_{}_diff.pdf'.format(ru))
