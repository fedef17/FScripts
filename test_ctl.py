#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import netCDF4 as nc
import climtools_lib as ctl
import pandas as pd
import numpy as np
from numpy import linalg as LA
from matplotlib import pyplot as plt
import pickle

cart = '/home/fabiano/Research/lavori/WP2_deliverable_Oct2018/Results_WP2/regime_indices/'
indxfi = 'regime_indices_ECE_LR.txt'

labels = np.loadtxt(cart+indxfi)
var, level, lat, lon, dates, time_units, var_units, time_cal = ctl.read4Dncfield('/data-hobbes/fabiano/PRIMAVERA_Stream1_Z500remap/zg500_Aday_EC-Earth3_T255_regrid25_1979-2014.nc', extract_level = 50000.)

resid_times, dates_init = ctl.calc_regime_residtimes(labels, dates = dates)
resid_times_noskip, dt = ctl.calc_regime_residtimes(labels, dates = dates, skip_singleday_pause = False)

plt.ion()
patnames = ['NAO +', 'Blocking', 'NAO -', 'Atl. Ridge']
patnames_short = ['NP', 'BL', 'NN', 'AR']
binzzz = np.arange(0,37,2)
for clu, clunam in zip(list(range(4)), patnames):
    pts = patnames_short[clu]
    fig = plt.figure()
    plt.title(clunam)
    n, bins, patches = plt.hist(resid_times[clu], bins = binzzz, alpha = 0.5, density = True, label = 'skip 1 day')
    n2, bins2, patches2 = plt.hist(resid_times_noskip[clu], bins = binzzz, alpha = 0.5, density = True, label = 'noskip')
    plt.legend()
    plt.xlabel('days')
    plt.ylabel('freq')
    fig.savefig(cart+'persistence_{}_wskip.pdf'.format(pts))

sys.exit()
#plt.ion()

# filename = '/data-hobbes/fabiano/SPHINX/zg_daily/ridottone.nc'
# var_day, lat_day, lon_day, dates_day, time_units, var_units = ctl.read3Dncfield(filename)
#
# dates_pdh_day = pd.to_datetime(dates_day)
#
# # filename='/data-hobbes/fabiano/SPHINX/tas_mon/ridotto.nc'
# # var_mon, lat_mon, lon_mon, dates_mon, time_units, var_units = ctl.read3Dncfield(filename)
# #
# # dates_pdh_mon = pd.to_datetime(dates_mon)
#
# anom = ctl.anomalies_daily(var_day, dates_day)
#
# anom_area, lat_area, lon_area = ctl.sel_area(lat_day, lon_day, anom, 'EAT')
# anom_areasea, dates_sea = ctl.sel_season(anom_area, dates_day, 'DJF')
#
# indicini = np.random.choice(anom_areasea.shape[0], 500, replace=False)
# anok = anom_areasea[indicini]
#
# pickle.dump([anok, lat_area, lon_area, dates_sea], open('anok_for_testing.p','wb'))
anok, lat_area, lon_area, dates_sea = pickle.load(open('anok_for_testing.p'))
#del var_day, anom, anom_area, anom_areasea, dates_pdh_day, dates_day

solver = ctl.eof_computation(anok, lat_area)

allclus_sk = []
allclus_ml = []
for i in range(50):
    cluster_sk, labels_sk = ctl.Kmeans_clustering_from_solver(solver, 4, 5)

    cluspatt_sk = ctl.compute_clusterpatterns(anok, labels_sk)

    cluster_ml, labels_ml = ctl.Kmeans_clustering_from_solver(solver, 4, 5, algorithm = 'molteni')

    cluspatt_ml = ctl.compute_clusterpatterns(anok, labels_ml)

    cluspatt_proj_ml = []
    for el in cluspatt_ml:
        cluspatt_proj_ml.append(solver.projectField(el)[:5])
    cluspatt_proj_ml = np.stack(cluspatt_proj_ml)

    cluspatt_proj_sk = []
    for el in cluspatt_sk:
        cluspatt_proj_sk.append(solver.projectField(el)[:5])
    cluspatt_proj_sk = np.stack(cluspatt_proj_sk)

    diff_mlsk = cluster_ml-cluster_sk
    diff_mlsk_proj = cluspatt_proj_ml -cluspatt_proj_sk

    diff_ml = cluspatt_proj_ml - cluster_ml
    diff_sk = cluspatt_proj_sk - cluster_sk

    for i, (el,elpr) in enumerate(zip(diff_mlsk, diff_mlsk_proj)): print(i, LA.norm(el), LA.norm(elpr))

    allclus_ml.append(cluster_ml)
    allclus_sk.append(cluster_sk)

allclus_ml = np.stack(allclus_ml)
allclus_sk = np.stack(allclus_sk)

mean_ml = np.mean(allclus_ml, axis = 0)
mean_sk = np.mean(allclus_sk, axis = 0)
std_ml = np.std(allclus_ml, axis = 0)
std_sk = np.std(allclus_sk, axis = 0)

print('Standard deviation of the different clusters:\n')
for i, (el,elpr) in enumerate(zip(std_ml, std_sk)): print(i, LA.norm(el), LA.norm(elpr))
print('Mean of the different clusters:\n')
for i, (el,elpr) in enumerate(zip(mean_ml, mean_sk)): print(i, LA.norm(el-elpr))

#plt.ion()
ctl.plot_multimap_contour(cluspatt_ml, lat_area, lon_area, filename = cart+'clusmap_ml.pdf')

ctl.plot_multimap_contour(cluspatt_sk, lat_area, lon_area, filename = cart+'clusmap_sk.pdf')
