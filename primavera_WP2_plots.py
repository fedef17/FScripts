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

##############################################
############ INPUTS #########

cart = '/home/fabiano/Research/lavori/WP2_deliverable_Oct2018/Results_WP2/'
results = pickle.load(open(cart+'res_primavera.p','rb'))

cart_ecmwf = '/home/fabiano/Research/lavori/WP2_deliverable_Oct2018/Results_ECMWF/'
results_ecmwf = pickle.load(open(cart_ecmwf+'res_primavera.p','rb'))

cart_ecmwf_full = '/home/fabiano/Research/lavori/WP2_deliverable_Oct2018/Results_ECMWF_full/'
results_ecmwf_full = pickle.load(open(cart_ecmwf_full+'res_primavera.p','rb'))

res_tags = ['labels', 'freq_clus', 'patcor', 'cluspattern', 'significance', 'et', 'cluspattern_area']

mod_tags = ['ECE_HR', 'HadGEM_HR', 'MPI_LR', 'ECE_LR', 'CMCC_LR','CNRM_LR', 'ECMWF_LR', 'CMCC_HR', 'CNRM_HR','ECMWF_HR', 'NCEP', 'MPI_HR', 'HadGEM_LR', 'ERA']
models = ['ECE', 'HadGEM', 'CMCC', 'CNRM', 'MPI', 'ECMWF']
# Significance plot
sig = results['significance']
sig_ecmwf = results_ecmwf['significance']


fig = plt.figure()
lab1 = 'LR'
lab2 = 'HR'
key_HR = ['HR_{}'.format(ke) for ke in range(1,5)]
key_LR = ['LR_{}'.format(ke) for ke in range(1,7)]
sig_HR = [sig_ecmwf[ke] for ke in key_HR]
sig_LR = [sig_ecmwf[ke] for ke in key_LR]
plt.scatter(list(range(len(sig_LR))), sig_LR, color = 'green', s=30, marker = '$L$', label = 'LR')
plt.scatter(6+np.arange(len(sig_HR)), sig_HR, color = 'orange', s=30, marker = '$H$', label = 'HR')

mean_HR = np.mean(sig_HR)
std_HR = np.std(sig_HR)
mean_LR = np.mean(sig_LR)
std_LR = np.std(sig_LR)
plt.errorbar(6+(len(key_HR)-1)/2., mean_HR, std_HR, color = 'orange', capsize = 5)
plt.scatter(6+(len(key_HR)-1)/2., mean_HR, color = 'orange', s = 20)
plt.errorbar((len(key_LR)-1)/2., mean_LR, std_LR, color = 'green', capsize = 5)
plt.scatter((len(key_LR)-1)/2., mean_LR, color = 'green', s = 20)

plt.legend(fontsize = 'small', loc = 1)
plt.title('Significance of regime structure - ECMWF ensemble')
keytot = key_LR+key_HR
plt.xticks(list(range(len(keytot))), keytot, size='small')
plt.ylabel('Significance')
fig.savefig(cart+'Significance_ECMWF.pdf')

# print('\n \n \n SOSTITUISCO LA MEDIAAA DI ECMWFFFFFFF\n \n \n')
# sig['ECMWF_LR'] = mean_LR
# sig['ECMWF_HR'] = mean_HR

syms = ['$H$']*len(models) + ['$L$']*len(models) + ['D']
labels = ['ECE', 'HadGEM', 'CMCC', 'CNRM', 'MPI', 'ECMWF']+6*[None]+['NCEP']
colors = ctl.color_set(len(models)+1)

wi = 0.3
#plt.ion()
fig = plt.figure()
lab1 = 'LR'
lab2 = 'HR'
for i, mod in enumerate(models):
    # plt.scatter(i, sig[mod+'_LR'], color = 'green', s=20, marker = '$L$')
    # plt.scatter(i, sig[mod+'_HR'], color = 'orange', s=20, marker = '$H$')
    if i > 0:
        lab1 = None
        lab2 = None
    plt.bar(i-wi/2, sig[mod+'_LR'], width = wi, color = 'green', label = lab1)
    plt.bar(i+wi/2,sig[mod+'_HR'], width = wi,  color = 'orange', label = lab2)
plt.bar(i+1-wi/2, sig['ERA'],width = wi,  color = 'black', label = 'ERA')
plt.bar(i+1+wi/2, sig['NCEP'],width = wi,  color = colors[-1], label = 'NCEP')
# plt.scatter(i+1, sig['ERA'], color = 'black', s=20, marker = 'D')
# plt.scatter(i+1, sig['NCEP'], color = colors[-1], s=15, marker = 'D')
models2 = models+['Obs']
plt.legend(fontsize = 'small', loc = 4)
plt.title('Significance of regime structure - Stream 1')
plt.xticks(list(range(len(models2))), models2, size='small')
plt.ylabel('Significance')
fig.savefig(cart+'Significance.pdf')


hrsigs = sig_HR + [sig[mod] for mod in sig.keys() if 'HR' in mod]
lrsigs = sig_LR + [sig[mod] for mod in sig.keys() if 'LR' in mod]

np.savetxt(cart+'HRsigs.dat', hrsigs, fmt = '%7.3f')
np.savetxt(cart+'LRsigs.dat', lrsigs, fmt = '%7.3f')

patt_ref = results['cluspattern']['ERA']
lat = np.arange(-90., 91., 2.5)
lon = np.arange(0., 360., 2.5)

patnames = ['NAO +', 'Blocking', 'NAO -', 'Atl. Ridge']
patnames_short = ['NP', 'BL', 'NN', 'AR']
for tag in mod_tags:
    print('adesso {}\n'.format(tag))
    patt = results['cluspattern'][tag]
    if np.any(np.isnan(patt)):
        print('There are {} NaNs in this patt.. replacing with zeros\n'.format(np.sum(np.isnan(patt))))
        patt[np.isnan(patt)] = 0.0
    cartout = cart+'Model_{}/'.format(tag)
    if not os.path.exists(cartout): os.mkdir(cartout)
    filename = cartout+'Allclus_'+tag+'.pdf'
    ctl.plot_multimap_contour(patt, lat, lon, filename, visualization = 'polar', central_lat_lon = (50.,0.), cmap = 'RdBu_r', title = 'North-Atlantic weather regimes - {}'.format(tag), subtitles = patnames, cb_label = 'Geopotential height anomaly (m)', color_percentiles = (0.5,99.5), fix_subplots_shape = (2,2), number_subplots = False)
    for patuno, patuno_ref, pp, pps in zip(patt, patt_ref, patnames, patnames_short):
        nunam = cartout+'clus_'+pps+'_'+tag+'.pdf'
        ctl.plot_double_sidebyside(patuno, patuno_ref, lat, lon, filename = nunam, visualization = 'polar', central_lat_lon = (50., 0.), title = pp, cb_label = 'Geopotential height anomaly (m)', stitle_1 = tag, stitle_2 = 'ERA', color_percentiles = (0.5,99.5))

sys.exit()
filename = cart+'HadGEM_blocking_example.pdf'
subtitles = ['ERA-Interim', 'HadGEM LR', 'HadGEM HR']
patts = [results['cluspattern']['ERA'][1], results['cluspattern']['HadGEM_LR'][1], results['cluspattern']['HadGEM_HR'][1]]
ctl.plot_multimap_contour(patts, lat, lon, filename, visualization = 'polar', central_lat_lon = (50.,0.), cmap = 'RdBu_r', title = 'Blocking regime', subtitles = subtitles, cb_label = 'Geopotential height anomaly (m)', color_percentiles = (0.5,99.5), fix_subplots_shape = (1,3), number_subplots = False, figsize = (16,8))

# Taylor plots
fig = plt.figure(figsize=(16,12))

for num, patt in enumerate(patnames):
    ax = plt.subplot(2, 2, num+1, polar = True)

    obs = results['cluspattern_area']['ERA'][num, ...]
    modpats_HR = [results['cluspattern_area'][tag+'_HR'][num, ...] for tag in models]
    tags = [tag+'_HR' for tag in models]
    modpats_LR = [results['cluspattern_area'][tag+'_LR'][num, ...] for tag in models]
    tags += [tag+'_LR' for tag in models]
    modpats = modpats_HR+modpats_LR
    if num == 2:
        modpats += [results['cluspattern_area']['NCEP'][3, ...]]
    elif num == 3:
        modpats += [results['cluspattern_area']['NCEP'][2, ...]]
    else:
        modpats += [results['cluspattern_area']['NCEP'][num, ...]]

    tags += ['NCEP']

    colors = ctl.color_set(len(modpats_HR)+1, bright_thres = 0.3)
    colors = colors[:-1] + colors[:-1] + [colors[-1]]
    syms = ['$H$']*len(modpats_HR) + ['$L$']*len(modpats_HR) + ['D']
    labels = ['ECE', 'HadGEM', 'CMCC', 'CNRM', 'MPI', 'ECMWF']+6*[None]+['NCEP']

    filename = cart+'TaylorPlot_{}.pdf'.format(patnames_short[num])
    label_ERMS_axis = 'Total RMS error (m)'
    label_bias_axis = 'Pattern mean (m)'
    #ctl.Taylor_plot(modpats, obs, filename, title = patt, label_bias_axis = label_bias_axis, label_ERMS_axis = label_ERMS_axis, colors = colors, markers = syms, only_first_quarter = True, legend = True, marker_edge = None, labels = labels, obs_label = 'ERA')
    legok = False
    ctl.Taylor_plot(modpats, obs, ax = ax, title = None, colors = colors, markers = syms, only_first_quarter = True, legend = legok, labels = labels, obs_label = 'ERA', mod_points_size = 50, obs_points_size = 70)

#Custom legend
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

legend_elements = []
for col, lab in zip(colors, labels):
    if lab is None:
        break
    legend_elements.append(Line2D([0], [0], marker='o', color=col, label=lab, linestyle = ''))

legend_elements.append(Line2D([0], [0], marker='D', color=colors[-1], label=labels[-1], linestyle = ''))
legend_elements.append(Line2D([0], [0], marker='D', color='black', label='ERA', linestyle = ''))
fig.legend(handles=legend_elements, loc=1, fontsize = 'large')

fig.tight_layout()
n1 = plt.text(0.15,0.6,'NAO +',transform=fig.transFigure, fontsize = 20)
n2 = plt.text(0.15,0.1,'NAO -',transform=fig.transFigure, fontsize = 20)
n3 = plt.text(0.6,0.6,'Blocking',transform=fig.transFigure, fontsize = 20)
n4 = plt.text(0.6,0.1,'Atl.Ridge',transform=fig.transFigure, fontsize = 20)
bbox=dict(facecolor = 'lightsteelblue', alpha = 0.7, edgecolor='black', boxstyle='round,pad=0.2')
n1.set_bbox(bbox)
n2.set_bbox(bbox)
n3.set_bbox(bbox)
n4.set_bbox(bbox)

fig.savefig(cart+'TaylorPlot_Stream1.pdf')
sys.exit()

colors = 4*['orange']+6*['green']+[colors[-1]]
syms = 4*['$H$'] + 6*['$L$'] + ['D']
labels = ['HR'] + 3*[None] + ['LR'] + 5*[None] + ['NCEP']

for num, patt in enumerate(patnames):
    obs = results['cluspattern_area']['ERA'][num, ...]
    modpats_HR = [results_ecmwf['cluspattern_area'][tag][num, ...] for tag in key_HR]
    modpats_LR = [results_ecmwf['cluspattern_area'][tag][num, ...] for tag in key_LR]
    tags = key_HR+key_LR+['NCEP']
    modpats = modpats_HR+modpats_LR
    if num == 2:
        modpats += [results['cluspattern_area']['NCEP'][3, ...]]
    elif num == 3:
        modpats += [results['cluspattern_area']['NCEP'][2, ...]]
    else:
        modpats += [results['cluspattern_area']['NCEP'][num, ...]]

    filename = cart_ecmwf+'TaylorPlot_{}_ECMWF.pdf'.format(patnames_short[num])
    label_ERMS_axis = 'Total RMS error (m)'
    label_bias_axis = 'Pattern mean (m)'
    ctl.Taylor_plot(modpats, obs, filename, title = patt, label_bias_axis = label_bias_axis, label_ERMS_axis = label_ERMS_axis, colors = colors, markers = syms, only_first_quarter = True, legend = True, marker_edge = None, labels = labels, obs_label = 'ERA')

key_HR = ['HR_{}_full'.format(ke) for ke in range(1,5)]
key_LR = ['LR_{}_full'.format(ke) for ke in range(1,7)]
for num, patt in enumerate(patnames):
    obs = results['cluspattern_area']['ERA'][num, ...]
    modpats_HR = [results_ecmwf_full['cluspattern_area'][tag][num, ...] for tag in key_HR]
    modpats_LR = [results_ecmwf_full['cluspattern_area'][tag][num, ...] for tag in key_LR]
    tags = key_HR+key_LR+['NCEP']
    modpats = modpats_HR+modpats_LR
    if num == 2:
        modpats += [results['cluspattern_area']['NCEP'][3, ...]]
    elif num == 3:
        modpats += [results['cluspattern_area']['NCEP'][2, ...]]
    else:
        modpats += [results['cluspattern_area']['NCEP'][num, ...]]

    filename = cart_ecmwf_full+'TaylorPlot_{}_ECMWF_full.pdf'.format(patnames_short[num])
    label_ERMS_axis = 'Total RMS error (m)'
    label_bias_axis = 'Pattern mean (m)'
    ctl.Taylor_plot(modpats, obs, filename, title = patt, label_bias_axis = label_bias_axis, label_ERMS_axis = label_ERMS_axis, colors = colors, markers = syms, only_first_quarter = True, legend = True, marker_edge = None, labels = labels, obs_label = 'ERA')
