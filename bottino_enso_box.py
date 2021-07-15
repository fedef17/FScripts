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
#from tunlib import gregplot_on_ax

from matplotlib.colors import LogNorm
from datetime import datetime

from scipy import stats
import xarray as xr
import glob
import xclim

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

cart_out = cart_out + 'indices/'
ctl.mkdir(cart_out)

allru = ['pi', 'b025', 'b050', 'b100']
allnams = ['piControl', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']

colors = ['black', 'forestgreen', 'orange', 'violet']
colors_vtr = ['black', 'lightgreen', 'forestgreen', 'moccasin', 'orange', 'thistle', 'violet']

####################################################################################################

enso = dict()
cartind = '/nas/BOTTINO/indices/enso500_xr/'
for ru in allru:
    if ru == 'pi':
        enso[ru] = xr.load_dataset(cartind + '{}_enso_360day.nc'.format(ru), use_cftime = True)
    else:
        enso[ru] = xr.load_dataset(cartind + '{}_enso_360day.nc'.format(ru), use_cftime = True)

        firstye = enso[ru].time.values[0].year
        lastye = enso[ru].time.values[-1].year
        enso[ru+'_st'] = enso[ru].sel(time = slice('{}-01-01'.format(lastye-200), '{}-12-30'.format(lastye)))
        enso[ru+'_tr'] = enso[ru].sel(time = slice('{}-01-01'.format(firstye), '{}-12-30'.format(firstye+50)))


fig, ax = plt.subplots(figsize = (12,8))

allpercs = dict()
for nu in [10, 25, 50, 75, 90]:
    allpercs['p{}'.format(nu)] = [np.percentile(enso[ru+'_'+vers]['tos'], nu) for ru in allru[1:] for vers in ['tr', 'st']]
allpercs['mean'] = [np.mean(enso[ru]['tos']).values for ru in allru[1:] for vers in ['tr', 'st']]
allpercs['min'] = [np.min(enso[ru]['tos']).values for ru in allru[1:] for vers in ['tr', 'st']]
allpercs['max'] = [np.max(enso[ru]['tos']).values for ru in allru[1:] for vers in ['tr', 'st']]

ru = 'pi'
obsperc = dict()
for nu in [10, 25, 50, 75, 90]:
    obsperc['p{}'.format(nu)] = np.percentile(enso[ru]['tos'], nu)
obsperc['mean'] = np.mean(enso[ru]['tos']).values
obsperc['min'] = np.min(enso[ru]['tos']).values
obsperc['max'] = np.max(enso[ru]['tos']).values

nams = [ru+'_'+vers for ru in allru[1:] for vers in ['tr', 'st']]
edgecol = np.concatenate([(col, col) for col in colors[1:]])

ctl.boxplot_on_ax(ax, allpercs, nams, colors_vtr[1:], edge_colors = edgecol, plot_mean = False, plot_minmax = True, plot_ensmeans = False, obsperc = obsperc, obs_color = 'black', obs_name = 'pi')
# ax.axhline(0, color = 'gray', linewidth = 0.5)
#ax.set_xticks(np.arange(4))
ax.set_xticklabels(nams + ['pi'])
#ax.set_title(tit)

#ctl.custom_legend(fig, colors_vtr, ['pi'] + nams, ncol = 4, add_space_below = 0.1)
ax.set_ylabel('ENSO index (K)')

fig.savefig(cart_out + 'enso_boxplot.pdf')
