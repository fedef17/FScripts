#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib import colors

import pickle

import climtools_lib as ctl
import climdiags as cd

from matplotlib.colors import LogNorm
from datetime import datetime

from scipy import stats
import seaborn as sns

#######################################
colo = '#d73027 #f46d43 #fdae61 #fee090 #ffffff #e0f3f8 #abd9e9 #74add1 #4575b4'
#colo = '#a50026 #d73027 #f46d43 #fdae61 #fee090 #e0f3f8 #abd9e9 #74add1 #4575b4 #313695'
colo = colo.split()
colo = colo[::-1]
# sns.palplot(colo)
cmappa = colors.ListedColormap(colo)
cmappa.set_over('#800026') #662506
cmappa.set_under('#023858') #542788

cart = '/home/fabiano/Research/lavori/WeatherRegimes/ERA_ref_r25_v4/'
fil = 'out_ERA_DJF_EAT_4clus_4pcs_1957-2014.p'

result_obs = pickle.load(open(cart+fil, 'rb'))

patt_ref = result_obs['cluspattern']
freq_ref = result_obs['freq_clus']
lat = result_obs['lat']
lon = result_obs['lon']
patnames = ['NAO +', 'Sc. Blocking', 'Atl. Ridge', 'NAO -']

clatlo = (70, -20)

clevels = np.arange(-135., 136., 30.)

# plt.rcParams['lines.dashed_pattern'] = [5, 5]
# for proj, blat in zip(['nearside', 'Npolar', 'Nstereo'], [0, 5, 25]):
#     print(proj)
#     filename = cart+'Allclus_OBSERVED_{}.pdf'.format(proj)
#     figs = ctl.plot_multimap_contour(patt_ref, lat, lon, filename, visualization = proj, central_lat_lon = clatlo, cmap = cmappa, title = 'Observed weather regimes', subtitles = patnames, cb_label = 'Geopotential height anomaly (m)', color_percentiles = (0.5,99.5), number_subplots = False, bounding_lat = blat, draw_grid = True, n_color_levels = 10, draw_contour_lines = True, clevels = clevels, lw_contour = 0.7)
#
#
cart = '/home/fabiano/Research/lavori/WeatherRegimes/prima_coup_v7/'
filogen = cart + 'out_prima_coup_v7_DJF_EAT_4clus_4pcs_1957-2014_refEOF.p'
results, results_ref = pickle.load(open(filogen, 'rb'))
results['ERA_0'] = results_ref
all_mods = np.array([ke.split('_')[0] for ke in results.keys()])
all_mems = np.array([ke.split('_')[1] for ke in results.keys()])

proj = 'nearside'
blat = 0
cart_out = '/home/fabiano/Research/articoli/Papers/EC-Earth-paper_Rein/'
#
# figs_all = []
# filename = cart_out + 'Allclus_{}.pdf'
# meanpatts = dict()
# meanfreqs = dict()
#
# for mod in ['EC-Earth3P', 'EC-Earth3P-HR']:
#     modmems = []
#     for mem in ['r1i1p2f1', 'r2i1p2f1', 'r3i1p2f1']:
#         modmem = mod+'_'+mem
#         modmems.append(modmem)
#         patt = results[modmem]['cluspattern']
#         figs = ctl.plot_multimap_contour(patt, lat, lon, filename.format(modmem), visualization = proj, central_lat_lon = clatlo, cmap = cmappa, title = modmem, subtitles = patnames, cb_label = 'Geopotential height anomaly (m)', color_percentiles = (0.5,99.5), number_subplots = False, bounding_lat = blat, draw_grid = True, n_color_levels = 10, draw_contour_lines = True, clevels = clevels, lw_contour = 0.7)
#         figs_all += figs
#     meanpatts[(mod, 'mean')] = np.mean([results[modmem]['cluspattern'] for modmem in modmems], axis = 0)
#     meanpatts[(mod, 'std')] = np.std([results[modmem]['cluspattern'] for modmem in modmems], axis = 0)
#     meanfreqs[(mod, 'mean')] = np.mean([results[modmem]['freq_clus'] for modmem in modmems], axis = 0)
#     meanfreqs[(mod, 'std')] = np.std([results[modmem]['freq_clus'] for modmem in modmems], axis = 0)
#
#
# filename = cart_out + 'Allclus_mean_compare.pdf'
# pattall = np.concatenate([patt_ref, meanpatts[('EC-Earth3P', 'mean')], meanpatts[('EC-Earth3P-HR', 'mean')]])
#
# subnames = 4*['NAO + ({:4.1f} %)', 'Sc. Blocking ({:4.1f} %)', 'Atl. Ridge ({:4.1f} %)', 'NAO - ({:4.1f} %)']
# allfreqs = np.concatenate([freq_ref, meanfreqs[('EC-Earth3P', 'mean')], meanfreqs[('EC-Earth3P-HR', 'mean')]])
#
# subtits = [nam.format(freq) for nam, freq in zip(subnames, allfreqs)]
#
# figs = ctl.plot_multimap_contour(pattall, lat, lon, filename, visualization = proj, central_lat_lon = clatlo, cmap = cmappa, title = None, subtitles = subtits, cb_label = 'Geopotential height anomaly (m)', color_percentiles = (0.5,99.5), number_subplots = False, bounding_lat = blat, draw_grid = True, n_color_levels = 10, draw_contour_lines = True, clevels = clevels, lw_contour = 0.5, fix_subplots_shape = (3, 4))
#
# filename = cart_out + 'Allmems_compare.pdf'
# ctl.plot_pdfpages(filename, figs_all, save_single_figs = False)


allmods = ['EC-Earth3P', 'EC-Earth3P-HR']
allmems = ['r1i1p2f1', 'r2i1p2f1', 'r3i1p2f1']
for mod in allmods:
    for mem in allmems:
        modmem = mod+'_'+mem
        results[modmem]['varopt'] = ctl.calc_varopt_molt(results[modmem]['pcs'], results[modmem]['centroids'], results[modmem]['labels'])
        results[modmem]['av_res_time'] = np.array([np.mean(results[modmem]['resid_times'][reg]) for reg in range(4)])

result_obs['varopt'] = ctl.calc_varopt_molt(result_obs['pcs'], result_obs['centroids'], result_obs['labels'])

cosee = ['significance', 'varopt', 'freq', 'RMS', 'patcor', 'trans_matrix', 'av_res_time']
allcose = dict()
for cos in cosee:
    try:
        for mod in allmods:
            modmems = [mod+'_'+mem for mem in allmems]
            allcose[(mod, cos, 'mean')] = np.mean([results[modmem][cos] for modmem in modmems], axis = 0)
            allcose[(mod, cos, 'std')] = np.std([results[modmem][cos] for modmem in modmems], axis = 0)
            print(cos, mod, allcose[(mod, cos, 'mean')], allcose[(mod, cos, 'std')])
        if cos in result_obs.keys():
            print(cos, result_obs[cos])
    except:
        continue

for cos in ['RMS', 'patcor']:
    for mod in allmods:
        allcose[(mod, cos, 'allregs')] = np.mean(allcose[(mod, cos, 'mean')])
        allcose[(mod, cos, 'allregs_std')] = np.mean(allcose[(mod, cos, 'std')])
        print(cos, mod, allcose[(mod, cos, 'allregs')], allcose[(mod, cos, 'allregs_std')])

cos = 'pers_prob'
for mod in allmods:
    allcose[(mod, cos, 'mean')] = np.array([allcose[(mod, 'trans_matrix', 'mean')][j,j] for j in range(4)])
    allcose[(mod, cos, 'std')] = np.array([allcose[(mod, 'trans_matrix', 'std')][j,j] for j in range(4)])
    print(cos, mod, allcose[(mod, cos, 'mean')], allcose[(mod, cos, 'std')])

result_obs[cos] = np.array([result_obs['trans_matrix'][j,j] for j in range(4)])
print(cos, result_obs[cos])



# import matplotlib as mpl
# from importlib import reload
# #gui_env = ['TKAgg','GTKAgg','Qt4Agg','WXAgg', 'Agg']
# for gui in mpl.rcsetup.all_backends:
#     print("testing", gui)
#     try:
#         mpl.use(gui,warn=False, force=True)
#         reload(plt)
#         print("Using:", mpl.get_backend())
#         filename = cart+'Allclus_OBSERVED_{}.pdf'.format(gui)
#         figs = ctl.plot_multimap_contour(patt, lat, lon, filename, visualization = 'Npolar', central_lat_lon = clatlo, cmap = cmappa, title = 'Observed weather regimes', subtitles = patnames, cb_label = 'Geopotential height anomaly (m)', color_percentiles = (0.5,99.5), number_subplots = False, bounding_lat = 5, add_rectangles = [ctl.sel_area_translate('EAT')], draw_grid = True, n_color_levels = 10, draw_contour_lines = True)
#     except Exception as expt:
#         print(expt)
#         continue
