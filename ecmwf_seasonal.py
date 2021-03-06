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

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 15
###############################################################################

cart_in = '/home/fedefab/Scrivania/Research/Post-doc/lavori/ecmwf_seasonal/WRtool_output/'
cart_out = '/home/fedefab/Scrivania/Research/Post-doc/lavori/ecmwf_seasonal/results/'
ctl.mkdir(cart_out)

exps = 'efpo3 efxtv egt3s egv90 egwve egwvf esys3 esys4 esys5'.split()
finam = cart_in + '{}/out_{}_DJF_EAT_4clus_4pcs_allyrs_refEOF.p'

## Load all
results = dict()
for ex in exps:
    fiok = finam.format(ex, ex)
    gigi = pickle.load(open(fiok, 'rb'))
    modnam = list(gigi['models'].keys())[0]
    results[ex] = gigi['models'][modnam]

results_ref = gigi['reference']

gigi = pickle.load(open('/home/fedefab/Scrivania/Research/Post-doc/lavori/ecmwf_seasonal/ERA_test/out_ERA_test_DJF_EAT_4clus_4pcs_1981-2010_refEOF.p', 'rb'))
results['ERA'] = gigi['models']['ERA']

colors = ctl.color_set(len(exps))
colors_all = colors + ['black']
exps_all = exps + ['ERA']

patnames = ['NAO+', 'SBL', 'AR', 'NAO-']

## Diagnostics
## Plots: variance ratio, patcor & phase space, persistence, frequency
## Two versions: simple and with ensemble member spread
# Variance ratio

figs = []

fig = plt.figure(figsize = (16,12))
#ax = fig.add_subplot(111)

allye = results['ERA']['freq_clus_seasonal_years']
obsfr = results['ERA']['freq_clus_seasonal'].T

axes = []
for reg in range(4):
    ax = fig.add_subplot(2, 2, reg+1)

    modskill = []
    modskill_p50 = []
    obserie = np.array([ob[reg] for ob in obsfr])
    for mod in exps:
        modserie = []
        modserie_p50 = []
        for ye in allye:
            cose = [results[mod]['freq_clus_seasonal']['{:02d}_{}1101'.format(nu, ye)][reg] for nu in range(25)]
            modserie.append(np.mean(cose))
            modserie_p50.append(np.median(cose))

        modskill.append(ctl.Rcorr(obserie, modserie))
        modskill_p50.append(ctl.Rcorr(obserie, modserie_p50))

    xs = np.arange(len(modskill))
    ax.scatter(xs, modskill, c = colors, s = 100)
    ax.scatter(xs, modskill_p50, c = colors, s = 100, marker = 'x')
    ax.grid()

    ax.set_xticks([])
    ax.set_title(patnames[reg])
    axes.append(ax)

ctl.adjust_ax_scale(axes)

fig.suptitle('Seasonal skill')
ctl.custom_legend(fig, colors, exps)
fig.savefig(cart_out + 'seas_skill.pdf')
figs.append(fig)


fig = plt.figure(figsize = (16,12))
ax = fig.add_subplot(111)

allye = results['ERA']['freq_clus_seasonal_years']
obsfr = results['ERA']['freq_clus_seasonal'].T

modskill = []
modskill_p50 = []
modskill_p75 = []
modskill_p25 = []
for mod in exps:
    modserie = []
    modserie_p50 = []
    modserie_p25 = []
    modserie_p75 = []
    for ye, obf in zip(allye, obsfr):
        cose = [ctl.Rcorr(results[mod]['freq_clus_seasonal']['{:02d}_{}1101'.format(nu, ye)], obf) for nu in range(25)]
        modserie.append(np.mean(cose))
        modserie_p50.append(np.median(cose))
        modserie_p25.append(np.percentile(cose, 25))
        modserie_p75.append(np.percentile(cose, 75))

    modskill.append(np.mean(modserie))
    modskill_p50.append(np.mean(modserie_p50))
    modskill_p25.append(np.mean(modserie_p25))
    modskill_p75.append(np.mean(modserie_p75))

xs = np.arange(len(modskill))
ax.scatter(xs, modskill, c = colors, s = 100)
ax.scatter(xs, modskill_p50, c = colors, s = 100, marker = 'x')
ax.scatter(xs, modskill_p25, c = colors, s = 100, marker = '<')
ax.scatter(xs, modskill_p75, c = colors, s = 100, marker = '>')
ax.grid()

fig.suptitle('Seasonal skill allregs')
ctl.custom_legend(fig, colors, exps)
fig.savefig(cart_out + 'seas_skill_allregs.pdf')
figs.append(fig)


for month in range(3):
    fig = plt.figure(figsize = (16,12))
    #ax = fig.add_subplot(111)

    allye = results['ERA']['freq_clus_seasonal_years']
    obsfr = results['ERA']['freq_clus_monthly'].T[month::3]

    axes = []
    for reg in range(4):
        ax = fig.add_subplot(2, 2, reg+1)

        modskill = []
        modskill_p50 = []
        obserie = np.array([ob[reg] for ob in obsfr])
        for mod in exps:
            modserie = []
            modserie_p50 = []
            for ye in allye:
                cose = [results[mod]['freq_clus_monthly']['{:02d}_{}1101'.format(nu, ye)][reg, month] for nu in range(25)]
                modserie.append(np.mean(cose))
                modserie_p50.append(np.median(cose))

            modskill.append(ctl.Rcorr(obserie, modserie))
            modskill_p50.append(ctl.Rcorr(obserie, modserie_p50))

        xs = np.arange(len(modskill))
        ax.scatter(xs, modskill, c = colors, s = 100)
        ax.scatter(xs, modskill_p50, c = colors, s = 100, marker = 'x')
        ax.grid()

        ax.set_xticks([])
        ax.set_title(patnames[reg])
        axes.append(ax)

    ctl.adjust_ax_scale(axes)

    fig.suptitle('Month {} skill'.format(month))
    ctl.custom_legend(fig, colors, exps)
    fig.savefig(cart_out + 'month{}_skill.pdf'.format(month))
    figs.append(fig)


    fig = plt.figure(figsize = (16,12))
    ax = fig.add_subplot(111)

    allye = results['ERA']['freq_clus_seasonal_years']
    obsfr = results['ERA']['freq_clus_monthly'].T[month::3]

    modskill = []
    modskill_p50 = []
    modskill_p75 = []
    modskill_p25 = []
    for mod in exps:
        modserie = []
        modserie_p50 = []
        modserie_p25 = []
        modserie_p75 = []
        for ye, obf in zip(allye, obsfr):
            cose = [ctl.Rcorr(results[mod]['freq_clus_seasonal']['{:02d}_{}1101'.format(nu, ye)], obf) for nu in range(25)]
            cose = [ctl.Rcorr(results[mod]['freq_clus_monthly']['{:02d}_{}1101'.format(nu, ye)][:, month], obf) for nu in range(25)]
            modserie.append(np.mean(cose))
            modserie_p50.append(np.median(cose))
            modserie_p25.append(np.percentile(cose, 25))
            modserie_p75.append(np.percentile(cose, 75))

        modskill.append(np.mean(modserie))
        modskill_p50.append(np.mean(modserie_p50))
        modskill_p25.append(np.mean(modserie_p25))
        modskill_p75.append(np.mean(modserie_p75))

    xs = np.arange(len(modskill))
    ax.scatter(xs, modskill, c = colors, s = 100)
    ax.scatter(xs, modskill_p50, c = colors, s = 100, marker = 'x')
    ax.scatter(xs, modskill_p25, c = colors, s = 100, marker = '<')
    ax.scatter(xs, modskill_p75, c = colors, s = 100, marker = '>')
    ax.grid()

    fig.suptitle('Month {} skill allregs'.format(month))
    ctl.custom_legend(fig, colors, exps)
    fig.savefig(cart_out + 'month{}_skill_allregs.pdf'.format(month))
    figs.append(fig)

sys.exit()

fig = plt.figure(figsize = (16,12))
ax = fig.add_subplot(111)

varrats = np.array([results[ex]['var_ratio'] for ex in exps_all])
xs = np.arange(len(varrats))
ax.scatter(xs, varrats, c = colors_all, s = 100)

ax.set_xticks([])
fig.suptitle('Variance ratio')
ctl.custom_legend(fig, colors_all, exps_all)
fig.savefig(cart_out + 'var_ratio.pdf')
figs.append(fig)

## Patcor
fig = plt.figure(figsize = (16,12))
axes = []
for reg in range(4):
    ax = fig.add_subplot(2, 2, reg+1)

    patcor = np.array([results[ex]['patcor'][reg] for ex in exps_all])
    xs = np.arange(len(patcor))

    ax.scatter(xs, patcor, c = colors_all, s = 100)

    ax.set_xticks([])
    ax.set_title(patnames[reg])
    axes.append(ax)
ctl.adjust_ax_scale(axes)
fig.suptitle('Pattern correlation')
ctl.custom_legend(fig, colors_all, exps_all)
fig.savefig(cart_out + 'patcor.pdf')
figs.append(fig)


## RMS
gg = 9.81
for ex in exps:
    results[ex]['dist_cen'] = np.array([ctl.distance(results[ex]['centroids'][reg]/gg, results_ref['centroids'][reg]) for reg in range(4)])
results['ERA']['dist_cen'] = np.array([ctl.distance(results['ERA']['centroids'][reg], results_ref['centroids'][reg]) for reg in range(4)])

fig = plt.figure(figsize = (16,12))
axes = []
for reg in range(4):
    ax = fig.add_subplot(2, 2, reg+1)

    patcor = np.array([results[ex]['dist_cen'][reg] for ex in exps_all])
    xs = np.arange(len(patcor))

    ax.scatter(xs, patcor, c = colors_all, s = 100)

    ax.set_title(patnames[reg])
    ax.set_xticks([])
    axes.append(ax)
ctl.adjust_ax_scale(axes)
fig.suptitle('Centroid distance')
ctl.custom_legend(fig, colors_all, exps)
fig.savefig(cart_out + 'distcen.pdf')
figs.append(fig)


## Frequency
fig = plt.figure(figsize = (16,12))
axes = []
for reg in range(4):
    ax = fig.add_subplot(2, 2, reg+1)

    patcor = np.array([results[ex]['freq_clus'][reg] for ex in exps_all])
    xs = np.arange(len(patcor))

    ax.scatter(xs, patcor, c = colors_all, s = 100)

    ax.set_title(patnames[reg])
    ax.set_xticks([])
    axes.append(ax)
ctl.adjust_ax_scale(axes)
fig.suptitle('Frequency')
ctl.custom_legend(fig, colors_all, exps_all)
fig.savefig(cart_out + 'freq_clus.pdf')
figs.append(fig)


## Persistence
for ex in exps_all:
    results[ex]['av_res'] = np.array([np.mean(results[ex]['resid_times'][reg]) for reg in range(4)])
results_ref['av_res'] = np.array([np.mean(results_ref['resid_times'][reg]) for reg in range(4)])

fig = plt.figure(figsize = (16,12))
axes = []
for reg in range(4):
    ax = fig.add_subplot(2, 2, reg+1)

    patcor = np.array([results[ex]['av_res'][reg] for ex in exps_all])
    xs = np.arange(len(patcor))

    ax.scatter(xs, patcor, c = colors_all, s = 100)

    ax.set_title(patnames[reg])
    ax.set_xticks([])
    axes.append(ax)
ctl.adjust_ax_scale(axes)
fig.suptitle('Av. persistence')
ctl.custom_legend(fig, colors_all, exps_all)
fig.savefig(cart_out + 'persistence.pdf')
figs.append(fig)


ctl.plot_pdfpages(cart_out + 'all_metrics.pdf', figs)


all_figures = []
for ex in exps_all:
    patt = results[ex]['cluspattern_area']
    if ex != 'ERA': patt = patt/gg
    patt_ref = results_ref['cluspattern_area']
    lat_ref = results[ex]['lat_area']
    lon_ref = results[ex]['lon_area']
    figs = ctl.plot_multimap_contour(patt, lat_ref, lon_ref, None, visualization = 'nearside', central_lat_lon = (70, -20), cmap = 'RdBu_r', title = ex, subtitles = patnames, cb_label = 'Geop. height anomaly (m)', cbar_range = (-150, 150), number_subplots = False, bounding_lat = 0, plot_margins = None, draw_grid = True, plot_type = 'pcolormesh')
    all_figures += figs

figs = ctl.plot_multimap_contour(patt_ref, lat_ref, lon_ref, None, visualization = 'nearside', central_lat_lon = (70, -20), cmap = 'RdBu_r', title = 'ERA (1957-2018)', subtitles = patnames, cb_label = 'Geop. height anomaly (m)', cbar_range = (-150, 150), number_subplots = False, bounding_lat = 0, plot_margins = None, draw_grid = True, plot_type = 'pcolormesh')
all_figures += figs

ctl.plot_pdfpages(cart_out + 'all_patterns.pdf', all_figures)


all_figures = []
for ex in exps_all:
    patt = results[ex]['cluspattern_area']
    if ex != 'ERA': patt = patt/gg
    patt_ref = results_ref['cluspattern_area']
    lat_ref = results[ex]['lat_area']
    lon_ref = results[ex]['lon_area']
    figs = ctl.plot_multimap_contour(patt-patt_ref, lat_ref, lon_ref, None, visualization = 'nearside', central_lat_lon = (70, -20), cmap = 'RdBu_r', title = ex, subtitles = patnames, cb_label = 'Diff to reference (m)', cbar_range = (-60, 60), number_subplots = False, bounding_lat = 0, plot_margins = None, draw_grid = True, plot_type = 'pcolormesh')
    all_figures += figs

ctl.plot_pdfpages(cart_out + 'all_patterns_diff.pdf', all_figures)

## Predictability
## Is there a correlation of seasonal frequencies with obs?
## On a daily basis which is the number of predictable days?
