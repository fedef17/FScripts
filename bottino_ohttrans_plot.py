#!/usr/bin/python
# -*- coding: utf-8 -*-


import numpy as np
import sys
import os
from matplotlib import pyplot as plt
from matplotlib import cm

import pickle
#import netCDF4 as nc

import climtools_lib as ctl
#import climdiags as cd
#from tunlib import gregplot_on_ax

#from matplotlib.colors import LogNorm
#from datetime import datetime

#from scipy import stats
import xarray as xr
import glob
#import xclim

import multiprocessing as mp
import psutil

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

#cart_out = '/g100_work/IscrB_QUECLIM/BOTTINO/bottino_an/'
cart_out = cart_out + 'ocean3d/'
ctl.mkdir(cart_out)

cart_in = cart_out + '../seasmean/'
gogo = pickle.load(open(cart_in + 'bottino_glomeans_1000.p', 'rb'))
glomeans, pimean = gogo

# allru = ['b990', 'b025', 'b050', 'b100']
# allsyear = [1990, 2025, 2050, 2100]
# allnams = ['stabilization-hist-1990', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']
# colors = ['teal', 'forestgreen', 'orange', 'violet']

allru = ['pi', 'b990', 'b025', 'b050', 'b065', 'b080', 'b100']
colors = ['black', 'lightslategray', 'forestgreen', 'orange', 'chocolate', 'maroon', 'violet']

cp0 = 3989.245 # J/kg/K
rho = 1025
oce_area = 3.6481962e+14 # m2

oht_all = dict()

Adiff = 1.2e-5

read_ts = False

oht_all = pickle.load(open(cart_out + 'oht_ts_deep.p', 'rb'))

fig, axs = plt.subplots(1, 3, figsize = (18,6))
for ru, col in zip(allru, colors):
    if not read_ts:
        oht_lev = []
        dtdz_lev = []

        if ru not in ['pi', 'b065']:
            filo = open(cart_out + 'ohtrans_{}.p'.format(ru), 'rb')
            for i in range(1000):
                try:
                    gigi = pickle.load(filo)
                except:
                    break
                oht_lev.append(gigi[0])
                dtdz_lev.append(gigi[1])

            filo.close()
        elif ru == 'b065':
            filo = open(cart_out + 'ohtrans_{}_f800.p'.format(ru), 'rb')
            for i in range(735):
                try:
                    gigi = pickle.load(filo)
                except:
                    break
                oht_lev.append(gigi[0])
                dtdz_lev.append(gigi[1])
                
            filo.close()

            filo = open(cart_out + 'ohtrans_{}_l200.p'.format(ru), 'rb')
            for i in range(265):
                try:
                    gigi = pickle.load(filo)
                except:
                    break
                oht_lev.append(gigi[0])
                dtdz_lev.append(gigi[1])
                
            filo.close()

            print(len(oht_lev))
            if len(oht_lev) < 1000:
                print('AAAAA MISSING 1 YEAR FOR b065!! copying last year twice for now')
                oht_lev.append(gigi[0])
                dtdz_lev.append(gigi[1])
        elif ru == 'pi':
            filo = open(cart_out + 'ohtrans_{}.p'.format(ru), 'rb')
            for i in range(501):
                try:
                    gigi = pickle.load(filo)
                except:
                    break
                oht_lev.append(gigi[0])
                dtdz_lev.append(gigi[1])

            filo.close()


        oht_lev = xr.concat(oht_lev, dim = 'year')*cp0*rho*86400*365*100 # J/cent
        dtdz_lev = xr.concat(dtdz_lev, dim = 'year')*Adiff*cp0*rho*86400*365*100 # J

        oht100 = oht_lev.sel(lev = slice(96., 98.)).mean('lev')
        oht700 = oht_lev.sel(lev = slice(696., 698.)).mean('lev')
        oht2000 = oht_lev.sel(lev = slice(1944., 1946.)).mean('lev')

        diff100 = dtdz_lev.sel(lev = slice(96., 98.)).mean('lev')
        diff700 = dtdz_lev.sel(lev = slice(696., 698.)).mean('lev')
        diff2000 = dtdz_lev.sel(lev = slice(1944., 1946.)).mean('lev')

        oht_all[(ru, 'moc', 100)] = oht100
        oht_all[(ru, 'moc', 700)] = oht700
        oht_all[(ru, 'moc', 2000)] = oht2000

        oht_all[(ru, 'diff', 100)] = diff100
        oht_all[(ru, 'diff', 700)] = diff700
        oht_all[(ru, 'diff', 2000)] = diff2000

    gtas = glomeans[(ru, 'tas')][1] - pimean['tas']
    yeas = np.arange(500)

    grun = ctl.running_mean(gtas, 20, remove_nans = True)
    for cosu, ax in zip([oht100, oht700, oht2000], axs.flatten()):
        #ax.scatter(gtas, cosu, s = 5, color = col)
        larun = ctl.running_mean(cosu, 20, remove_nans = True)
        ax.plot(grun, -larun, color = col, label = ru + ' - moc', lw = 2)

    for cosu, ax in zip([diff100, diff700, diff2000], axs.flatten()):
        #ax.scatter(gtas, cosu, s = 5, color = col)
        larun = ctl.running_mean(cosu, 20, remove_nans = True)
        ax.plot(grun, -larun, color = col, label = ru + ' - diff', lw = 2, ls = ':')

for ax, tit in zip(axs.flatten(), ['100 m', '700 m', '2000 m']):
    ax.grid()
    ax.set_title(tit)
    ax.set_xlabel('GTAS (K)')

#axs[2].legend()
axs[0].set_ylabel('Downward Heat transport (J/cent)')

fig.savefig(cart_out + 'ohtrans_vs_gtas.pdf')

fig, ax = plt.subplots(1, 1, figsize = (16,9))
for ru, col in zip(allru, colors):
    gtas = glomeans[(ru, 'tas')][1] - pimean['tas']
    grun = ctl.running_mean(gtas, 20, remove_nans = True)

    tot_trans = 100*oht_all[(ru, 'deep')].differentiate('year')
    larun = ctl.running_mean(tot_trans, 20, remove_nans = True)
    if len(larun) == len(grun):
        ax.plot(grun, larun, color = col, label = ru + ' - TOT', lw = 2)
    else:
        ax.plot(grun[:-1], larun, color = col, label = ru + ' - TOT', lw = 2)

    larun1 = ctl.running_mean(oht_all[(ru, 'moc', 2000)], 20, remove_nans = True)
    ax.plot(grun, -larun1, color = col, label = ru + ' - moc', ls = '--', lw = 2)

    larun2 = ctl.running_mean(oht_all[(ru, 'diff', 2000)], 20, remove_nans = True)
    ax.plot(grun, -larun2, color = col, label = ru + ' - diff', lw = 2, ls = ':')

    #ax.plot(grun, -(larun1+larun2), color = col, label = ru + ' - diff+moc', lw = 2, ls = '-.')

    #print(ru, ctl.Rcorr(-tot_trans.values, oht_all[(ru, 'moc', 2000)].values+oht_all[(ru, 'diff', 2000)].values))

ax.grid()
ax.set_xlabel('GTAS (K)')
ax.set_ylabel('Downward Heat transport (J/cent)')

fig.savefig(cart_out + 'ohtrans_vs_gtas_deep.pdf')



fig, ax = plt.subplots(1, 1, figsize = (16,9))

ru = 'b100'

gtas = glomeans[(ru, 'tas')][1] - pimean['tas']
grun = ctl.running_mean(gtas, 20, remove_nans = True)

tot_trans = 100*oht_all[(ru, 'deep')].differentiate('year')
larun = ctl.running_mean(tot_trans, 20, remove_nans = True)

ax.plot(grun, larun, label = 'TOT', lw = 2)

larun1 = ctl.running_mean(oht_all[(ru, 'moc', 2000)], 20, remove_nans = True)
ax.plot(grun, -larun1, label = 'ADV', lw = 2)

larun2 = ctl.running_mean(oht_all[(ru, 'diff', 2000)], 20, remove_nans = True)
ax.plot(grun, -larun2, label = 'dia. diff')

ax.grid()
ax.set_xlabel('GTAS (K)')
ax.set_ylabel('Downward Heat transport (J/cent)')

fig.savefig(cart_out + 'ohtrans_vs_gtas_deep_b100.pdf')


pickle.dump(oht_all, open(cart_out + 'ohtrans_mean.p', 'wb'))


tot_trans_pi = 100*oht_all[('pi', 'deep')].differentiate('year')
larun_pi = np.mean(tot_trans_pi).values

larun1_pi = np.mean(oht_all[('pi', 'moc', 2000)]).values

larun2_pi = np.mean(oht_all[('pi', 'diff', 2000)]).values

#resid = tot_trans_pi - (oht_all[('pi', 'diff', 2000)]+oht_all[('pi', 'moc', 2000)])
resid_pi = larun_pi + larun1_pi + larun2_pi


fig, ax = plt.subplots(1, 1, figsize = (16,9))
for ru, col in zip(allru[1:], colors[1:]):
    gtas = glomeans[(ru, 'tas')][1] - pimean['tas']
    grun = ctl.running_mean(gtas, 20, remove_nans = True)

    tot_trans = 100*oht_all[(ru, 'deep')].differentiate('year')
    larun = ctl.running_mean(tot_trans, 20, remove_nans = True)
    if len(larun) == len(grun):
        ax.plot(grun, larun-larun_pi, color = col, label = ru + ' - TOT', lw = 2)
    else:
        ax.plot(grun[:-1], larun-larun_pi, color = col, label = ru + ' - TOT', lw = 2)

    larun1 = ctl.running_mean(oht_all[(ru, 'moc', 2000)], 20, remove_nans = True)
    ax.plot(grun, -larun1+larun1_pi, color = col, label = ru + ' - moc', ls = '--', lw = 2)

    larun2 = ctl.running_mean(oht_all[(ru, 'diff', 2000)], 20, remove_nans = True)
    ax.plot(grun, -larun2+larun2_pi, color = col, label = ru + ' - diff', lw = 2, ls = ':')


    resid = tot_trans + (oht_all[(ru, 'diff', 2000)]+oht_all[(ru, 'moc', 2000)])
    resid = ctl.running_mean(resid, 20, remove_nans = True)
    ax.plot(grun, resid - resid_pi, color = col, label = ru + ' - resid', lw = 2, ls = '-.')

    #ax.plot(grun, -(larun1+larun2), color = col, label = ru + ' - diff+moc', lw = 2, ls = '-.')

    #print(ru, ctl.Rcorr(-tot_trans.values, oht_all[(ru, 'moc', 2000)].values+oht_all[(ru, 'diff', 2000)].values))

ax.grid()
ax.set_xlabel('GTAS (K)')
ax.set_ylabel('Downward Heat transport (J/cent)')

fig.savefig(cart_out + 'ohtrans_vs_gtas_deep_anom.pdf')


sys.exit()
#############################################################

#############################################################
### maps of OHT trends
lats = np.linspace(-89.5, 89.5, 180)
lons = np.linspace(0, 359, 360)

oht_patt = dict()

for ru, col in zip(allru, colors):
    print(ru)
    filo = open(cart_out + 'ohtrans_{}.p'.format(ru), 'rb')

    oht100 = []
    oht700 = []
    oht2000 = []
    for i in range(500):
        try:
            oht_lev, oht100_i, oht700_i, oht2000_i = pickle.load(filo)
        except:
            break
        oht100.append(oht100_i)
        oht700.append(oht700_i)
        oht2000.append(oht2000_i)

    # if ru == 'b065':
    #     print('AAAAA MISSING 1 YEAR FOR b065!! copying last year twice for now')
    #     oht100.append(oht100_i)
    #     oht700.append(oht700_i)
    #     oht2000.append(oht2000_i)

    filo.close()

    oht100 = np.stack(oht100)
    oht700 = np.stack(oht700)
    oht2000 = np.stack(oht2000)

    for var, lab in zip([oht100, oht700, oht2000], [100, 700, 2000]):
        var[var == 0.0] = np.nan

        oht_patt[(ru, lab, 'ini')] = var[:50].mean(axis = 0)*cp0*rho
        oht_patt[(ru, lab, 'fin')] = var[-50:].mean(axis = 0)*cp0*rho
        oht_patt[(ru, lab, 'change')] = oht_patt[(ru, lab, 'fin')] - oht_patt[(ru, lab, 'ini')]

pickle.dump(oht_patt, open(cart_out + 'ohtrans_patt.p', 'wb'))
oht_patt = pickle.load(open(cart_out + 'ohtrans_patt.p', 'rb'))

#for lev, tit in zip([700, 2000, 'deep'], ['0-700m', '700-2000m', '> 2000 m']):
for lev, tit in zip([100, 700, 2000], ['100 m', '700 m', '2000 m']):
    plpa = []
    subt = []
    hatch = []
    for ru in allru:
        plpa.append(oht_patt[(ru, lev, 'fin')])
        subt.append(ru + ': ' + tit)

    [fig] = ctl.plot_multimap_contour(plpa, lats, lons, visualization = 'Robinson', central_lat_lon = (0., -120.), filename = cart_out + 'ohtrans_patt_fin_{}.pdf'.format(lev), subtitles = subt, plot_anomalies = True, cbar_range = (-1, 1), cmap = ctl.heatmap(), figsize = (16,9), fix_subplots_shape = (2,3), cb_label = 'Downward heat trasport (W/m2)')#, add_hatching = hatch, hatch_styles = ['///', '', ''])

    for ax in fig.axes[:-1]:
        ax.set_facecolor('gainsboro')
    fig.savefig(cart_out + 'ohtrans_patt_fin_{}.pdf'.format(lev))
