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

filna = '/nas/BOTTINO/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/{}/{}/day_r25/ua/ua*nc'

allru = ['pi', 'b025', 'b050', 'b100']
allnams = ['piControl', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']

colors = ['black', 'forestgreen', 'orange', 'violet']

allru2 = allru + ['c1950']
allnams2 = allnams + ['control-1950']
colors2 = colors + ['steelblue']

#############################################################################
areasel = dict()
areasel['NATL'] = [-60., 0., 20., 70.]
areasel['NEPAC'] = [200., 240., 20., 65.]
areasel['NCPAC'] = [150., 200., 20., 65.]
areasel['SPAC'] = [180., 240., -70., -20.]
areasel['SATL'] = [-45., 10., -70., -20.]
areasel['SIND'] = [50., 110., -70., -20.]

allseas = ['DJF', 'MAM', 'JJA', 'SON']

resdict = dict()
for na, ru, col in zip(allnams2, allru2, colors2):
    print(ru)
    mem = 'r1i1p1f1'
    if ru == 'c1950':
        mem = 'r1i1p2f1'

    filist = glob.glob(filna.format(na, mem))
    filist.sort()

    # for area in areasel.keys():
    #     jli, jspeed, dates = cd.jli_from_files(filist, area = areasel[area])
    #
    #     resdict[(ru, 'jli', area)] = jli
    #     resdict[(ru, 'jspeed', area)] = jspeed
    #     resdict[(ru, 'dates')] = dates
    for season in allseas:
        jli, jspeed, dates = cd.jli_from_files(filist, allareadict = areasel, season = season)

        for area in areasel.keys():
            resdict[(ru, 'jli', area, season)] = jli[area]
            resdict[(ru, 'jspeed', area, season)] = jspeed[area]
        resdict[(ru, 'dates', season)] = dates

with open(cart_out + 'res_jli200_v2.p', 'wb') as filox:
    pickle.dump(resdict, filox)

with open(cart_out + 'res_jli200_v2.p', 'rb') as filox:
    resdict = pickle.load(filox)

for tip, aruok, colok in zip(['', '_wc1950'], [allru, allru2], [colors, colors2]):
    figN, axN = plt.subplots(1,3,figsize = (24,8))
    figS, axS = plt.subplots(1,3,figsize = (24,8))
    iS = 0
    iN = 0
    for area in areasel.keys():
        for season in allseas:
            fig = plt.figure(figsize = (24,12))
            ax1 = fig.add_subplot(121)
            ax2 = fig.add_subplot(122)

            if area[0] == 'N' and season == 'DJF':
                axlu = axN[iN]
                iN += 1
            elif area[0] == 'S' and season == 'JJA':
                axlu = axS[iS]
                iS += 1
            else:
                axlu = None

            latsel = np.arange(areasel[area][-2], areasel[area][-1]+0.1, 0.5)
            kos = np.concatenate([resdict[(ru, 'jspeed', area, season)] for ru in aruok])
            vmin, vmax = (np.min(kos), np.max(kos))
            print(vmin, vmax)
            vsel = np.linspace(vmin, vmax, 100)

            def dopdf(var, xi, bnd_width):
                pdf = ctl.calc_pdf(var, bnd_width = bnd_width)
                pdfok = pdf(xi)
                pdfok /= np.sum(pdfok)
                return pdfok

            for ru, col in zip(aruok, colok):
                jli = resdict[(ru, 'jli', area, season)]
                jspeed = resdict[(ru, 'jspeed', area, season)]

                jliserie = ctl.bootstrap(jli, resdict[(ru, 'dates', season)], None, apply_func = dopdf, func_args = [latsel, 0.22], n_choice = 50, n_bootstrap = 1000)
                jspedserie = ctl.bootstrap(jspeed, resdict[(ru, 'dates', season)], None, apply_func = dopdf, func_args = [vsel, None], n_choice = 50, n_bootstrap = 1000)

                pdf = ctl.calc_pdf(jli, bnd_width = 0.22)
                pdfok = pdf(latsel)
                pdfok /= np.sum(pdfok)

                jlimin = np.percentile(jliserie, 10, axis = 0)
                jlimax = np.percentile(jliserie, 90, axis = 0)
                ax1.fill_between(latsel, jlimin, jlimax, color = col, alpha = 0.3)
                ax1.plot(latsel, pdfok, label = ru, color = col, linewidth = 3)

                if axlu is not None:
                    axlu.fill_between(latsel, jlimin, jlimax, color = col, alpha = 0.3)
                    axlu.plot(latsel, pdfok, label = ru, color = col, linewidth = 3)

                pdf = ctl.calc_pdf(jspeed)
                pdfok = pdf(vsel)
                pdfok /= np.sum(pdfok)

                jlimin = np.percentile(jspedserie, 10, axis = 0)
                jlimax = np.percentile(jspedserie, 90, axis = 0)
                ax2.fill_between(vsel, jlimin, jlimax, color = col, alpha = 0.2)
                ax2.plot(vsel, pdfok, label = ru, color = col, linewidth = 3)

        axlu.grid()
        axlu.set_xlabel('Latitude')
        axlu.set_title(area)
        axlu.legend()

        ax1.grid()
        ax2.grid()
        ax1.set_xlabel('Latitude')
        ax2.set_xlabel('u wind (m/s)')
        ax1.set_title('Jet latitude index')
        ax2.set_title('Jet speed')
        ax1.legend()
        ax2.legend()
        fig.suptitle('{} - {}'.format(area, season))
        #fig.savefig(cart_out + 'jlinspeed_bottino_{}{}.pdf'.format(area, tip))
        figs.append(fig)

    ctl.plot_pdfpages(cart_out + 'jlinspeed_bottinoall{}.pdf'.format(tip), figs)

    figN.savefig(cart_out + 'jlibott_allareas_N{}.pdf'.format(tip))
    figS.savefig(cart_out + 'jlibott_allareas_S{}.pdf'.format(tip))
