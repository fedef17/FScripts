#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os

from matplotlib import pyplot as plt
import matplotlib.cm as cm
import matplotlib.patheffects as PathEffects

import netCDF4 as nc
import cartopy.crs as ccrs
import pandas as pd

from numpy import linalg as LA
from eofs.standard import Eof
from scipy import stats
import itertools as itt

from sklearn.cluster import KMeans



from datetime import datetime
import pickle

import climtools_lib as ctl
import climdiags as cd

from mpl_toolkits.mplot3d import Axes3D

#######################################

cart_in = '/data-hobbes/fabiano/SPHINX/zg_daily/'
cart_out = '/home/fabiano/Research/lavori/SPHINX_for_lisboa/WRtool/test_trend_traj_lcb0/'

if not os.path.exists(cart_out): os.mkdir(cart_out)

fil = cart_in + 'lcb0-1850-2100-NDJFM_zg500_NH_14473.nc'

var, level, lat, lon, dates, time_units, var_units, time_cal = ctl.read4Dncfield(fil, extract_level = 50000)

#climat_mean, dates_climate_mean = ctl.trend_daily_climat(var, dates)

#var_anom = ctl.anomalies_daily_detrended(var, dates, climat_mean = climat_mean, dates_climate_mean = dates_climate_mean)
var_season, dates_season = ctl.sel_season(var, dates, 'DJF')

erafile = '/data-hobbes/fabiano/OBS/ERA/ERAInterim/zg500/zg500_Aday_ERAInterim_2deg_1979-2014.nc'
ERA_ref_EAT = cd.WRtool_from_file(erafile, 'DJF', 'EAT', extract_level_4D = 50000., numclus = 4, heavy_output = True, run_significance_calc = False)

area = 'EAT'
ref_solver = ERA_ref_EAT['solver']
ref_patterns_area = ERA_ref_EAT['cluspattern_area']

#resu_nodet = cd.WRtool_core(var_season, lat, lon, dates_season, area, run_significance_calc = False, ref_solver = ref_solver, ref_patterns_area = ref_patterns_area, detrended_eof_calculation = False, detrended_anom_for_clustering = False, heavy_output = True)

resu1 = cd.WRtool_core(var_season, lat, lon, dates_season, area, run_significance_calc = False, ref_solver = ref_solver, ref_patterns_area = ref_patterns_area, detrended_eof_calculation = True, detrended_anom_for_clustering = True, heavy_output = False)

resu2 = cd.WRtool_core(var_season, lat, lon, dates_season, area, run_significance_calc = False, ref_solver = ref_solver, ref_patterns_area = ref_patterns_area, detrended_eof_calculation = True, detrended_anom_for_clustering = False, heavy_output = False)

freq1 = ctl.calc_seasonal_clus_freq(resu1['labels'], resu1['dates'])
freq2 = ctl.calc_seasonal_clus_freq(resu2['labels'], resu2['dates'])

years = np.unique(pd.to_datetime(resu1['dates']).year)[:-1]
patnames = ['NAO +', 'Blocking', 'NAO -', 'Atl. Ridge']
fig = plt.figure()
plt.ylim(15.,35.)
plt.title('detr.EOFs + detrended anomalies')
for clu, clunam in enumerate(patnames):
    smut = ctl.running_mean(freq1[:,clu], wnd = 30)
    plt.plot(years, smut, label = clunam)
fig = plt.figure()
plt.ylim(20.,30.)
plt.title('detr.EOFs + non-detrended anomalies')
for clu, clunam in enumerate(patnames):
    smut = ctl.running_mean(freq2[:,clu], wnd = 30)
    plt.plot(years, smut, label = clunam)


trans_matrix = ctl.calc_regime_transmatrix(1, resu1['labels'], dates_season)
print(trans_matrix)

filter = resu1['dist_centroid'] < np.percentile(resu1['dist_centroid'], 70)
trans_pcs = ctl.find_transition_pcs(1, resu1['labels'][filter], dates_season[filter], resu1['pcs'][filter], filter_longer_than = 3, max_days_between = 1)

ngp = 100
(x0, x1) = (np.percentile(resu1['pcs'][:,0], 1), np.percentile(resu1['pcs'][:,0], 99))
(y0, y1) = (np.percentile(resu1['pcs'][:,1], 1), np.percentile(resu1['pcs'][:,1], 99))
xss = np.linspace(x0,x1,ngp)
yss = np.linspace(y0,y1,ngp)
xi2, yi2 = np.meshgrid(xss, yss)

(z0, z1) = (np.percentile(resu1['pcs'][:,2], 1), np.percentile(resu1['pcs'][:,2], 99))
zss = np.linspace(z0,z1,ngp)
xi3, yi3, zi3 = np.meshgrid(xss, yss, zss)

xib, zib = np.meshgrid(xss, zss)
yic, zic = np.meshgrid(yss, zss)

cordss = [xss, yss, zss]
namz = ['x', 'y', 'z']
coups = [(0,1), (0,2), (1,2)]

for cou in coups:
    regime_pdf_2D = []
    regime_pdf_2D_func = []
    xi, yi = np.meshgrid(cordss[cou[0]], cordss[cou[1]])

    for clus in range(4):
        print(clus)
        okclus = resu1['labels'] == clus
        okpc = resu1['pcs'][okclus, :]
        kufu = ctl.calc_pdf(okpc[:,cou].T)
        print('fatto kufu\n')
        zi = kufu(np.vstack([xi.flatten(), yi.flatten()]))
        regime_pdf_2D.append(zi)
        regime_pdf_2D_func.append(kufu)

    fig = plt.figure(figsize=(8, 6), dpi=150)
    for li in trans_pcs[0,1]:
        li3 = li[:,cou].T
        plt.plot(li3[0], li3[1], color = 'grey')

    for clus, namcm in enumerate(['Purples','Blues','Greens','Oranges']):
        plt.contour(xi, yi, regime_pdf_2D[clus].reshape(xi.shape), cmap = cm.get_cmap(namcm))

    plt.title('Projection on {} and {} eofs'.format(cou[0], cou[1]))
    plt.xlabel('EOF {}'.format(cou[0]))
    plt.ylabel('EOF {}'.format(cou[1]))

for i in [1., 10., 100., 1000.]:
    kufu = ctl.calc_pdf(okpc[:, cou].T, bnd_width = i)
    print('fatto kufu\n')
    zi = kufu(np.vstack([xi.flatten(), yi.flatten()]))
    fig = plt.figure(figsize=(8, 6), dpi=150)
    # plt.title('Projection on {} and {} eofs'.format(cou[0], cou[1]))
    plt.title('BAND {}'.format(i))
    ax = fig.add_subplot(111, projection='3d')

    ax.contour(xi, yi, zi.reshape(xi.shape), cmap = cm.get_cmap('Oranges'))
    ax.plot_surface(xi, yi, zi.reshape(xi.shape), cmap = cm.get_cmap('Oranges'))

    # fig = plt.figure(figsize=(8, 6), dpi=150)
    # ax = fig.add_subplot(111, projection='3d')
    # for li in trans_pcs[0,1]:
    #     li3 = li[:,cou].T
    #     zi1 = regime_pdf_2D_func[0](li3)
    #     ax.plot(li3[0], li3[1], zi1, color = 'grey')
    #
    # for clus, namcm in enumerate(['Purples','Blues','Greens','Oranges']):
    #     ax.plot_surface(xi, yi, regime_pdf_2D[clus].reshape(xi.shape), cmap = cm.get_cmap(namcm), alpha = 0.8)
    #
    # plt.title('Projection on {} and {} eofs'.format(cou[0], cou[1]))
    # plt.xlabel('EOF {}'.format(cou[0]))
    # plt.ylabel('EOF {}'.format(cou[1]))
    #
    # okpc = resu1['pcs']
    # kufu = ctl.calc_pdf(okpc[:, cou].T)
    # print('fatto kufu\n')
    # zi = kufu(np.vstack([xi.flatten(), yi.flatten()]))
    # fig = plt.figure(figsize=(8, 6), dpi=150)
    # plt.title('Projection on {} and {} eofs'.format(cou[0], cou[1]))
    # ax = fig.add_subplot(111, projection='3d')
    #
    # ax.contour(xi, yi, zi.reshape(xi.shape), cmap = cm.get_cmap('Oranges'))
    # ax.plot_surface(xi, yi, zi.reshape(xi.shape), cmap = cm.get_cmap('Oranges'), alpha = 0.8)

cou = (1,2)
okpc = resu1['pcs']
for i in [1., 10., 100., 1000.]:
    kufu = ctl.calc_pdf(okpc[:, cou].T, bnd_width = i)
    print('fatto kufu\n')
    zi = kufu(np.vstack([xi.flatten(), yi.flatten()]))
    fig = plt.figure(figsize=(8, 6), dpi=150)
    # plt.title('Projection on {} and {} eofs'.format(cou[0], cou[1]))
    plt.title('BAND {}'.format(i))
    ax = fig.add_subplot(111, projection='3d')

    ax.contour(xi, yi, zi.reshape(xi.shape), cmap = cm.get_cmap('Oranges'))
    ax.plot_surface(xi, yi, zi.reshape(xi.shape), cmap = cm.get_cmap('Oranges'))

sys.exit()
#3d

ngp = 20
(x0, x1) = (np.percentile(resu1['pcs'][:,0], 1), np.percentile(resu1['pcs'][:,0], 99))
(y0, y1) = (np.percentile(resu1['pcs'][:,1], 1), np.percentile(resu1['pcs'][:,1], 99))
xss = np.linspace(x0,x1,ngp)
yss = np.linspace(y0,y1,ngp)
xi2, yi2 = np.meshgrid(xss, yss)

(z0, z1) = (np.percentile(resu1['pcs'][:,2], 1), np.percentile(resu1['pcs'][:,2], 99))
zss = np.linspace(z0,z1,ngp)
xi3, yi3, zi3 = np.meshgrid(xss, yss, zss)

regime_pdf_3D = []
for clus in range(4):
    print(clus)
    okclus = resu1['labels'] == clus
    okpc = resu1['pcs'][okclus, :]

    kufu = ctl.calc_pdf(okpc[:,:3].T)
    print('fatto kufu\n')
    zi = kufu(np.vstack([okpc[:,0].flatten(), okpc[:,1].flatten(), okpc[:,2].flatten()]))
    regime_pdf_3D.append(zi)

fig = plt.figure(figsize=(8, 6), dpi=150)
ax = fig.add_subplot(111, projection='3d')

for clus, namcm in enumerate(['Purples','Blues','Greens','Oranges']):
    okclus = resu1['labels'] == clus
    okpc = resu1['pcs'][okclus, :3]
    #ax.scatter(okpc[:,0], okpc[:,1], okpc[:,2], color = namcm[:-1], s = 2, alpha = 0.5)
    ax.scatter(okpc[:,0], okpc[:,1], okpc[:,2], c = regime_pdf_3D[clus], cmap = cm.get_cmap(namcm), s = 2, alpha = 0.5)

for li in trans_pcs[0,1]:
    li3 = li[:,:3].T
    ax.plot(li3[0], li3[1], li3[2], color = 'grey')

# years = [dtclm[0].year for dtclm in dates_climate_mean]
# gigi = np.stack([np.mean(ctl.sel_season(va, dat, 'DJF', cut = False)[0], axis = 0) for va, dat in zip(climat_mean, dates_climate_mean)])
#
# #gigi = gigi - gigi[0]
#
# plt.ion()
# ctl.plot_animation_map(gigi, lat, lon, fps_anim = 5, labels = years, filename = cart_out + 'mean_geopot.gif', visualization = 'standard', cmap = 'RdBu_r', title = 'Mean geopotential field trend', cb_label = 'Geopotential height anomaly (m)', color_percentiles = (0.5,99.5), plot_anomalies = False, figsize = (16,8))
