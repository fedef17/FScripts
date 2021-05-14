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

cart_out = '/home/fabiano/Research/lavori/BOTTINO/analisi/'

filna = '/nas/BOTTINO/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/{}/r1i1p1f1/day_r25/ua/ua*nc'

allru = ['pi', 'b025', 'b050', 'b100']
allnams = ['piControl', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']

colors = ['black', 'forestgreen', 'orange', 'violet']

#############################################################################
# resdict = dict()
#
# for na, ru, col in zip(allnams, allru, colors):
#     print(ru)
#     filist = glob.glob(filna.format(na))
#     filist.sort()
#
#     jli, jspeed, dates = cd.jli_from_files(filist)
#
#     resdict[(ru, 'jli')] = jli
#     resdict[(ru, 'jspeed')] = jspeed
#     resdict[(ru, 'dates')] = dates
#
# with open(cart_out + 'res_jli200.p', 'wb') as filox:
#     pickle.dump(resdict, filox)

with open(cart_out + 'res_jli200.p', 'rb') as filox:
    resdict = pickle.load(filox)

fig = plt.figure(figsize = (24,12))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

latsel = np.arange(30, 71, 0.5)
kos = np.concatenate([resdict[(ru, 'jspeed')] for ru in allru])
vmin, vmax = (np.min(kos), np.max(kos))
print(vmin, vmax)
vsel = np.linspace(vmin, vmax, 100)

def dopdf(var, xi, bnd_width):
    pdf = ctl.calc_pdf(var, bnd_width = bnd_width)
    pdfok = pdf(xi)
    pdfok /= np.sum(pdfok)
    return pdfok

for ru, col in zip(allru, colors):
    jli = resdict[(ru, 'jli')]
    jspeed = resdict[(ru, 'jspeed')]

    jliserie = ctl.bootstrap(jli, resdict[(ru, 'dates')], None, apply_func = dopdf, func_args = [latsel, 0.22], n_choice = 50, n_bootstrap = 1000)
    jspedserie = ctl.bootstrap(jspeed, resdict[(ru, 'dates')], None, apply_func = dopdf, func_args = [vsel, None], n_choice = 50, n_bootstrap = 1000)

    pdf = ctl.calc_pdf(jli, bnd_width = 0.22)
    pdfok = pdf(latsel)
    pdfok /= np.sum(pdfok)

    jlimin = np.percentile(jliserie, 10, axis = 0)
    jlimax = np.percentile(jliserie, 90, axis = 0)
    ax1.fill_between(latsel, jlimin, jlimax, color = col, alpha = 0.3)
    ax1.plot(latsel, pdfok, label = ru, color = col, linewidth = 3)

    pdf = ctl.calc_pdf(jspeed)
    pdfok = pdf(vsel)
    pdfok /= np.sum(pdfok)

    jlimin = np.percentile(jspedserie, 10, axis = 0)
    jlimax = np.percentile(jspedserie, 90, axis = 0)
    ax2.fill_between(vsel, jlimin, jlimax, color = col, alpha = 0.3)
    ax2.plot(vsel, pdfok, label = ru, color = col, linewidth = 3)

ax1.grid()
ax2.grid()
ax1.set_xlabel('Latitude')
ax2.set_xlabel('u wind (m/s)')
ax1.set_title('Jet latitude index')
ax2.set_title('Jet speed')
ax1.legend()
ax2.legend()
fig.savefig(cart_out + 'jlinspeed_bottino.pdf')
