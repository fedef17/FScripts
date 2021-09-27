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
    cart = '/home/fabiano/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'xaru':
    cart = '/home/fedef/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'tintin':
    cart = '/home/fabiano/work/lavori/BOTTINO/'

cart_in = cart + 'seasmean/'
cart_out = cart + 'internal_var/'
ctl.mkdir(cart_out)

allru = ['pi', 'b025', 'b050', 'b100']
allnams = ['piControl', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']

colors = ['black', 'forestgreen', 'orange', 'violet']

allnams2 = allnams + ['ssp585', 'historical']
allru2 = allru + ['ssp585', 'hist']
colors2 = colors + ['indianred', 'steelblue']
####################################################################################################

glomeans, pimean, yeamean, mapmean = pickle.load(open(cart_in + 'bottino_seasmean_2D.p', 'rb'))

for ru in allru:
    print(ru)
    years, gtas = glomeans[(ru, 'tas')]

    if ru != 'pi':
        #coef2, covmat2 = np.polyfit(coso[0], coso[1], deg = 2, cov = True)
        coef3, covmat3 = np.polyfit(years, gtas, deg = 3, cov = True)

        #fitco2 = np.polyval(coef2, years)
        fitco3 = np.polyval(coef3, years)
        gtas3 = gtas-fitco3
    else:
        gtas3 = gtas

    g50_3 = ctl.butter_filter(gtas3, 50)

    allfigs = []
    fig = plt.figure()
    plt.plot(years, gtas3)
    plt.plot(years, g50_3, linewidth = 3)
    allfigs.append(fig)

    fig, ax = plt.subplots(figsize = (16,12))

    ps = np.abs(np.fft.rfft(gtas3))**2
    frq = np.fft.rfftfreq(gtas3.size, 1)
    invfr = 1/frq

    frbins = [5, 10, 20, 50, 100]

    barz = []
    xba = []
    for bi0, bi1 in zip(frbins[:-1], frbins[1:]):
        xba.append('{} - {}'.format(bi0, bi1))
        okke = (bi0 <= invfr) & (invfr < bi1)
        gig = np.sum(ps[okke])
        barz.append(gig)

    ax.bar(np.arange(len(barz)), barz)
    ax.set_xticks(np.arange(len(barz)))
    ax.set_xticklabels(xba, rotation = 30)
    ax.set_xlabel('Period (yr)')
    ax.set_ylabel(r'Integrated spectral power ($K^2$)')

    allfigs.append(fig)

    for var in ['tas', 'pr', 'clt', 'rlut']:
        tama = yeamean[(ru, var)][50:]
        if ru != 'pi':
            coef3_var, covmat3_var = np.polyfit(years, glomeans[(ru, var)][1], deg = 3, cov = True)
            fitco3_var = np.polyval(coef3_var, years)
            tama3 = tama - fitco3_var[50:, np.newaxis, np.newaxis]
        else:
            tama3 = tama

        # plt.figure()
        # plt.plot(years, np.gradient(g50_3))
        # plt.grid()
        incr = np.gradient(g50_3[50:]) > 0
        decr = np.gradient(g50_3[50:]) < 0

        tasdecr = tama3[decr].mean('year') - tama3.mean('year')
        tasincr = tama3[incr].mean('year') - tama3.mean('year')

        fig = ctl.plot_multimap_contour([tasincr, tasdecr], figsize = (16,9), plot_anomalies=True, subtitles= ['gtas increasing', 'gtas decreasing'], color_percentiles = (5,95), title = ru+' - '+var)
        allfigs.append(fig[0])

        var_trend, var_intercept, var_trend_err, var_intercept_err, var_pval = ctl.calc_trend_climatevar(g50_3[50:], tama3)

        fig = ctl.plot_map_contour(var_trend, tama3.lat, tama3.lon, figsize = (16,9), plot_anomalies=True, color_percentiles = (5,95), title = 'regr with gtas: '+ru+' - '+var, add_hatching = var_pval > 0.05)

        allfigs.append(fig)

    ctl.plot_pdfpages(cart_out + 'gtas_oscillations_{}.pdf'.format(ru), allfigs)
