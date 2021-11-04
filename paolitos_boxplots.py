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

import xarray as xr

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 28
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.axisbelow'] = True

#############################################################################

# leggi cose
cart_in = '/home/fedef/Research/lavori/paolitos_LWA/LWA_timeseries/'
cart_out = '/home/fedef/Research/lavori/paolitos_LWA/'

filin = cart_in + 'timeseries_{}_{}.p'

areas = ['NH', 'EAT']#, 'PNA']
versions = ['LR', 'MR', 'HR']

colors = ctl.color_set(12, sns_palette = 'Paired')

cose = dict()
for area in areas:
    cose[area] = dict()
    for vers in versions:
        print(area, vers)
        gigi = pickle.load(open(filin.format(area, vers), 'rb'))
        print(gigi['domain'], gigi['resolution'])
        #cose[(area, vers)] = gigi['timeseries']
        mod_names = gigi['labels']
        if 'ERA5' in mod_names:
            cose[area]['obs'] = gigi['timeseries'][0]
            mod_names = mod_names[1:]
            ts = gigi['timeseries'][1:]
        for mod, tise in zip(mod_names, ts):
            cose[area][(mod, vers)] = tise

mod_list = []
mod_list2 = []
allvers = []
for mod in mod_names:
    for vers in versions:
        if (mod, vers) in cose['NH']:
            mod_list.append(mod + '-' + vers)
            mod_list2.append((mod, vers))
            allvers.append(vers)

positions, posticks = ctl.positions(mod_list)

### Adding colors for mid res
in1 = mod_list2.index(('HadGem', 'LR'))
colors.insert(in1+1, np.mean([colors[in1], colors[in1+1]], axis = 0))
in1 = mod_list2.index(('ECMWF', 'LR'))
colors.insert(in1+1, np.mean([colors[in1], colors[in1+1]], axis = 0))

# positions = [0., 0.7]
# posticks = [0.35]
# for i in range(len(mod_names)-2):
#     positions.append(positions[-1]+0.7+0.4)
#     positions.append(positions[-1]+0.7)
#     posticks.append(np.mean(positions[-2:]))

positions.append(positions[-1]+0.6+0.7)
positions.append(positions[-1]+0.7+0.4)
positions.append(positions[-1]+0.7)

posticks += positions[-3:]

for area in areas:
    fig, ax = plt.subplots(figsize = (12,8))

    allpercs = dict()
    for nu in [10, 25, 50, 75, 90]:
        allpercs['p{}'.format(nu)] = [np.percentile(cose[area][(modvers[0], modvers[1])], nu) for modvers in mod_list2]
    allpercs['mean'] = [np.mean(cose[area][(modvers[0], modvers[1])]) for modvers in mod_list2]

    obsperc = dict()
    for nu in [10, 25, 50, 75, 90]:
        obsperc['p{}'.format(nu)] = np.percentile(cose[area]['obs'], nu)
    obsperc['mean'] = np.mean(cose[area]['obs'])

    #nomi = [mod + '_' + vers for mod in mod_names for vers in versions]
    # allvers = [vers for mod in mod_names for vers in versions]
    nomi = mod_list

    ctl.boxplot_on_ax(ax, allpercs, nomi, colors, versions = allvers, plot_mean = True, plot_ensmeans = True, ens_colors = ['steelblue', 'indianred'], ens_names = ['LR', 'HR'], obsperc = obsperc, obs_color = 'black', obs_name = 'ERA5', plot_minmax = False, positions = positions)
    # ax.axhline(0, color = 'gray', linewidth = 0.5)
    ax.set_xticks(posticks)
    ax.set_xticklabels(mod_names + ['ERA5', 'LR', 'HR'])
    #ax.set_title(tit)

    #ctl.custom_legend(fig, colors, nomi, ncol = 2, add_space_below = 0.1)
    ax.set_ylabel('LWA (m/s)')

    #plt.title('cmip6 vs cmip6 var_fr '+tip)
    ax.set_title(area)
    fig.savefig(cart_out + 'LWA_boxplot_{}.pdf'.format(area))


sys.exit()
###########################################################################################

cart = '/home/fedef/Research/lavori/paolitos_LWA/LWA_patterns/'
# che devo fa: leggere i pattern di ERA, fare patcor e RMS o il taylor?
# farei il taylor e poi uno scatter tra montg e transLWA patcors

resu = dict()
for fil in os.listdir(cart):
    if 'WR_' in fil:
        mod = fil[3:-2]
        resu[mod] = pickle.load(open(cart+fil, 'rb'))
        resu[mod]['Montg streamf'] = np.stack(resu[mod]['Montg streamf'])[:, ::-1, :] # il montg è invertito

mods_all = []
for mod in resu:
    print(mod, resu[mod]['WR labels'])
    if mod != 'ERA5':
        mods_all.append(mod)

mods_all.sort()

mods_all = ['CMCC-CM2-HR4', 'CMCC-CM2-VHR4', 'CNRM-CM6-1', 'CNRM-CM6-1-HR', 'EC-Earth3', 'EC-Earth3-HR', 'ECMWF-IFS-LR', 'ECMWF-IFS-MR', 'ECMWF-IFS-HR', 'HadGEM3-GC31-LL', 'HadGEM3-GC31-MM', 'HadGEM3-GC31-HH', 'MPI-ESM1-2-HR', 'MPI-ESM1-2-XR']

positions, posticks = ctl.positions(mods_all)

modstip = dict()
modstip['HR'] = ['CMCC-CM2-VHR4', 'CNRM-CM6-1-HR', 'EC-Earth3-HR', 'ECMWF-IFS-HR', 'HadGEM3-GC31-HH', 'MPI-ESM1-2-XR']
modstip['LR'] = ['CMCC-CM2-HR4', 'CNRM-CM6-1', 'EC-Earth3', 'ECMWF-IFS-LR', 'HadGEM3-GC31-LL', 'MPI-ESM1-2-HR']
modstip['MR'] = ['ECMWF-IFS-MR', 'HadGEM3-GC31-MM']

patnames = resu[mod]['WR labels']

modshort = ['CMCC', 'CNRM', 'EC-Earth', 'ECMWF', 'HadGEM', 'MPI']
markers = ['o', '^', 's', '*', 'X', 'd']

colors = ['steelblue', 'forestgreen', 'indianred']
cols = []
for mod in mods_all:
    if mod in modstip['LR']:
        cols.append(colors[0])
    elif mod in modstip['HR']:
        cols.append(colors[2])
    else:
        cols.append(colors[1])

marks = []
for mod in mods_all:
    for modsh, ma in zip(modshort, markers):
        if modsh in mod:
            marks.append(ma)

coldict = dict(zip(mods_all, cols))
markdict = dict(zip(mods_all, marks))

colors = ctl.color_set(12, sns_palette = 'Paired')
mr_cos = np.where(np.array([mod in modstip['MR'] for mod in mods_all]))[0]
for gi in mr_cos:
    colors.insert(gi, np.mean(colors[gi-1:gi+1], axis = 0))

gigi = resu['ERA5']
gogo = dict()
gogo['lat'] = resu['EC-Earth3']['lat']
gogo['lon'] = gigi['lon']
for nom in ['tot LWA', 'trans LWA', 'Montg streamf']:
    gigixr = xr.DataArray(np.stack(gigi[nom]).squeeze(), dims = ('reg', 'lat', 'lon'), coords = {'reg': [0, 1, 2, 3], 'lat': gigi['lat'], 'lon': gigi['lon']})
    gigixr_int = gigixr.interp(lat=gogo['lat'], lon = gogo['lon'])
    gogo[nom] = gigixr_int.values
    #gigixr_int.sel(reg = 0).plot()

resu['ERA5'] = gogo


areas = dict() # lonW, lonE, latS, latN
areas['NML'] = (-180., 180., 30., 80.)
areas['EAT'] = (-80., 40., 30., 80.)
areas['EAText'] = (-80., 100., 30., 80.)

for aaa in ['NML', 'EAT', 'EAText']:
    # Taylor semplice
    cart_out = cart + aaa + '/'
    ctl.mkdir(cart_out)

    for tip in ['tot LWA', 'trans LWA', 'Montg streamf']:
        patt_ref = resu['ERA5'][tip]
        olat = resu['ERA5']['lat']
        olon = resu['ERA5']['lon']

        fig = plt.figure(figsize=(16,12))
        for num, patt in enumerate(patnames):
            ax = plt.subplot(2, 2, num+1, polar = True)

            obs = patt_ref[num]
            obs, lat_area, lon_area = ctl.sel_area(olat, olon, obs, areas[aaa])
            modpats = [ctl.sel_area(resu[mod]['lat'], resu[mod]['lon'], resu[mod][tip][num], areas[aaa])[0] for mod in mods_all]

            ctl.Taylor_plot(modpats, obs, latitude = lat_area, ax = ax, colors = cols, markers = marks, only_first_quarter = True, plot_ellipse = False, ellipse_color = colors[0], max_val_sd = 1.6, title = patt, alpha_markers = 0.8)

            #ax.text(0.95, 0.95, patt, horizontalalignment='center', verticalalignment='center', rotation='horizontal',transform=ax.transAxes, fontsize = 18)

        ctl.custom_legend(fig, cols, mods_all, markers = marks, ncol = 4, add_space_below = 0.05, loc = 'right', fontsize = 14)
        fig.savefig(cart_out + 'taylor_{}_{}.pdf'.format(tip.split()[0], aaa))


    # Taylor con ellipse
    for tip in ['tot LWA', 'trans LWA', 'Montg streamf']:
        patt_ref = resu['ERA5'][tip]

        fig = plt.figure(figsize=(16,12))
        for num, patt in enumerate(patnames):
            ax = plt.subplot(2, 2, num+1, polar = True)

            obs = patt_ref[num]
            obs, lat_area, lon_area = ctl.sel_area(olat, olon, obs, areas[aaa])
            for vers in ['LR', 'HR']:
                modpats = [ctl.sel_area(resu[mod]['lat'], resu[mod]['lon'], resu[mod][tip][num], areas[aaa])[0] for mod in modstip[vers]]
                cols_ok = [coldict[mod] for mod in modstip[vers]]
                marks_ok = [markdict[mod] for mod in modstip[vers]]
                ctl.Taylor_plot(modpats, obs, latitude = lat_area, ax = ax, colors = cols_ok, markers = marks_ok, only_first_quarter = True, plot_ellipse = True, ellipse_color = cols_ok[0], max_val_sd = 1.6, alpha_markers = 0.8)#, title = patt)

            vers = 'MR'
            modpats = [ctl.sel_area(resu[mod]['lat'], resu[mod]['lon'], resu[mod][tip][num], areas[aaa])[0] for mod in modstip[vers]]
            cols_ok = [coldict[mod] for mod in modstip[vers]]
            marks_ok = [markdict[mod] for mod in modstip[vers]]
            ctl.Taylor_plot(modpats, obs, latitude = lat_area, ax = ax, colors = cols_ok, markers = marks_ok, only_first_quarter = True, plot_ellipse = False, ellipse_color = cols_ok[0], max_val_sd = 1.6, title = patt, alpha_markers = 0.8)

            #ax.text(1, 1, patt, horizontalalignment='center', verticalalignment='center', rotation='horizontal',transform=ax.transAxes, fontsize = 16)

        # ax.text(0.05, 0.75, 'EAT', horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 35)
        ctl.custom_legend(fig, cols, mods_all, markers = marks, ncol = 4, add_space_below = 0.05, loc = 'right', fontsize = 14)
        fig.savefig(cart_out + 'taylor_{}_well_{}.pdf'.format(tip.split()[0], aaa))


    for tip in ['tot LWA', 'trans LWA', 'Montg streamf']:
        # scatter/bar plot montg/trans
        fig = plt.figure(figsize=(16,12))
        axes = []
        for num, patt in enumerate(patnames):
            ax = plt.subplot(2, 2, num+1)

            obs = resu['ERA5'][tip][num]
            obs, lat_area, lon_area = ctl.sel_area(olat, olon, obs, areas[aaa])
            modpats = [ctl.sel_area(resu[mod]['lat'], resu[mod]['lon'], resu[mod][tip][num], areas[aaa])[0] for mod in mods_all]

            patcors = [ctl.Rcorr(obs, patt, lat_area) for patt in modpats]

            for pos, tra, col in zip(positions, patcors, colors):
                ax.bar(pos, tra, color = col, width = 0.4)

            ax.set_xticks(posticks)
            ax.set_xticklabels([])
            if num in [2,3]: ax.set_xticklabels(modshort, rotation = 45.)

            ax.set_title(patt, fontsize = 16)
            axes.append(ax)
            ax.grid(axis = 'y')
            #if num in [2,3]: ax.set_xlabel('regime streamf. pattern correlation')
            if num in [0,2]: ax.set_ylabel('{} pattern correlation'.format(tip))

        ctl.adjust_ax_scale(axes)
        # ax.text(0.05, 0.75, 'EAT', horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 35)
        #ctl.custom_legend(fig, colors, mods_all, ncol = 4, fontsize = 14)
        fig.savefig(cart_out + '{}_barplot_{}.pdf'.format(tip.split()[0], aaa))


    # scatter tra montg e trans
    fig = plt.figure(figsize=(16,12))
    axes = []
    for num, patt in enumerate(patnames):
        ax = plt.subplot(2, 2, num+1)

        tip = 'trans LWA'
        obs = resu['ERA5'][tip][num]
        obs, lat_area, lon_area = ctl.sel_area(olat, olon, obs, areas[aaa])
        modpats = [ctl.sel_area(resu[mod]['lat'], resu[mod]['lon'], resu[mod][tip][num], areas[aaa])[0] for mod in mods_all]

        patcors_trans = [ctl.Rcorr(obs, patt, lat_area) for patt in modpats]

        tip = 'Montg streamf'
        obs = resu['ERA5'][tip][num]
        obs, lat_area, lon_area = ctl.sel_area(olat, olon, obs, areas[aaa])
        modpats = [ctl.sel_area(resu[mod]['lat'], resu[mod]['lon'], resu[mod][tip][num], areas[aaa])[0] for mod in mods_all]
        patcors_montg = [ctl.Rcorr(obs, patt, lat_area) for patt in modpats]

        for mon, tra, col, mar in zip(patcors_montg, patcors_trans, cols, marks):
            ax.scatter(mon, tra, color = col, marker = mar, s = 30)
        ax.set_title(patt)
        axes.append(ax)
        if num in [2,3]: ax.set_xlabel('regime streamf. pattern correlation')
        if num in [0,2]: ax.set_ylabel('trans LWA pattern correlation')

    ctl.adjust_ax_scale(axes)
    # ax.text(0.05, 0.75, 'EAT', horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 35)
    ctl.custom_legend(fig, cols, mods_all, markers = marks, ncol = 4, add_space_below = 0.05, loc = 'right', fontsize = 14)
    fig.savefig(cart_out + 'scatt_montg_vs_trans_{}.pdf'.format(aaa))
