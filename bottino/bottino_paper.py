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
import pymannkendall as mk

plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
titlefont = 22
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.axisbelow'] = True
plt.rcParams['legend.fontsize'] = 18

#############################################################################

if os.uname()[1] == 'hobbes':
    cart_out = '/home/fabiano/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'xaru':
    cart_out = '/home/fedef/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'tintin':
    cart_out = '/home/fabiano/work/lavori/BOTTINO/'

cart_out += 'nonlin_evol_1000/'
ctl.mkdir(cart_out)

colors = ['black', 'royalblue', 'lightslategray', 'forestgreen', 'orange', 'chocolate', 'maroon', 'violet', 'crimson']
allru = ['pi', 'hist', 'b990', 'b025', 'b050', 'b065', 'b080', 'b100', 'ssp585']

####################################################################################################

cart_in = cart_out + '../seasmean/'
# gogo = pickle.load(open(cart_in + 'bottino_seasmean_2D.p', 'rb'))
# glomeans, pimean, yeamean, _ = gogo
gogo = pickle.load(open(cart_in + 'bottino_glomeans_1000.p', 'rb'))
glomeans, pimean = gogo


## calc tas trend at the end of the sim
filo = open(cart_out + 'trend_gtas.txt', 'w')
filo.write('Trend of GTAS in last 100 years:\n')
for ru in allru:
    yea, tas = glomeans[(ru, 'tas')]

    res = stats.linregress(yea[-100:], tas[-100:]) # trend of last 100 year
    rel_trend = 100*res.slope
    rel_trend_err = 100*res.stderr

    print('{} - {:6.4f} +/- {:6.4f} K/cent\n'.format(ru, rel_trend, rel_trend_err))
    filo.write('{} - {:6.4f} +/- {:6.4f} K/cent\n'.format(ru, rel_trend, rel_trend_err))

filo.close()

print('200 years')
for ru in allru:
    yea, tas = glomeans[(ru, 'tas')]

    res = stats.linregress(yea[-200:], tas[-200:]) # trend of last 200 year
    rel_trend = 100*res.slope
    rel_trend_err = 100*res.stderr

    print('{} - {:6.4f} +/- {:6.4f} K/cent\n'.format(ru, rel_trend, rel_trend_err))
    #filo.write('{} - {:6.4f} +/- {:6.4f} K/cent\n'.format(ru, rel_trend, rel_trend_err))
print('500 years')
for ru in allru:
    yea, tas = glomeans[(ru, 'tas')]

    res = stats.linregress(yea[-500:], tas[-500:]) # trend of last 500 year
    rel_trend = 100*res.slope
    rel_trend_err = 100*res.stderr

    print('{} - {:6.4f} +/- {:6.4f} K/cent\n'.format(ru, rel_trend, rel_trend_err))

filo = open(cart_out + 'final_gtas.txt', 'w')
filo.write('Final GTAS increase (last 30 years):\n')
for ru in allru:
    yea, tas = glomeans[(ru, 'tas')]

    deltatas = np.mean(tas[-30:])-pimean['tas']
    print('{} - {:6.4f} K\n'.format(ru, deltatas))
    filo.write('{} - {:6.4f} K\n'.format(ru, deltatas))

filo.close()

print(pimean['net_toa'])

filo = open(cart_out + 'final_netTOA.txt', 'w')
filo.write('Final net TOA (last 30 years):\n')
for ru in allru:
    yea, tas = glomeans[(ru, 'net_toa')]

    deltatas = np.mean(tas[-30:])
    print('{} - {:6.4f} W/m2\n'.format(ru, deltatas))
    filo.write('{} - {:6.4f} W/m2\n'.format(ru, deltatas))

filo.close()

filo = open(cart_out + 'final_netSRF.txt', 'w')
filo.write('Final net SFC (last 30 years):\n')
for ru in allru[2:-1]:
    yea, tas = glomeans[(ru, 'net_srf')]

    deltatas = np.mean(tas[-30:])
    print('{} - {:6.4f} W/m2\n'.format(ru, deltatas))
    filo.write('{} - {:6.4f} W/m2\n'.format(ru, deltatas))

filo.close()

### Determine up to which point the trend is significant
for ru in allru:
    print(ru)
    tas = glomeans[(ru, 'tas')][1]
    for n in np.arange(0, 951, 10):
        if n > len(tas)-10: continue
        pio = mk.original_test(tas[n:])
        print(n, pio.trend, pio.p)

fact = 60*60*24
do_pr = True
if do_pr:
    plt.rcParams['xtick.labelsize'] = 28
    plt.rcParams['ytick.labelsize'] = 28
    titlefont = 28
    plt.rcParams['figure.titlesize'] = titlefont
    plt.rcParams['axes.titlesize'] = 28
    plt.rcParams['axes.labelsize'] = 28
    plt.rcParams['axes.axisbelow'] = True
    plt.rcParams['legend.fontsize'] = 28


    # A coefficient for clausius-clapeyron applied to global pr
    exptas = lambda t : np.exp(17.625*(t-273.15)/(t-30.11))
    A = pimean['pr']/exptas(pimean['tas'])

    fig = plt.figure(figsize = (16,12))
    for ru, col in zip(allru, colors):
        if ru[0] == 'b':
            s = 1
        else:
            s = 10
        plt.scatter(glomeans[(ru, 'tas')][1]-pimean['tas'], fact*(glomeans[(ru, 'pr')][1]-pimean['pr']), color = col, s = s, label = ru)

        if ru == 'ssp585':
            pr = glomeans[(ru, 'pr')][1]
            pr_anom = fact*(pr-pimean['pr'])
            tas = glomeans[(ru, 'tas')][1]
            tas_anom = tas-pimean['tas']

            coeffs, covmat = np.polyfit(tas_anom, pr_anom, deg = 1, cov = True)
            x_nu = np.arange(-0.5, 10.1, 0.5)
            fitted = np.polyval(coeffs, x_nu)
            plt.plot(x_nu, fitted, color = col, ls = '--')
            x_nu = np.arange(1, 6.6, 0.5)
            fitted = np.polyval(coeffs, x_nu)
            plt.plot(x_nu, fitted, color = col)

            #coeffs, covmat = np.polyfit(tas_anom, np.log(pr_anom), deg = 1, cov = True)
            #coeffs, covmat = np.polyfit(tas_anom, np.log(pr), deg = 2, cov = True)
            #x_nu = np.arange(-0.5, 9.6, 0.5)
            #fitted = np.exp(np.polyval(coeffs, x_nu))
            #plt.plot(x_nu, fitted-pimean['pr'], color = col)

            # fitted2 = A*exptas(x_nu+pimean['tas'])
            # plt.plot(x_nu, fitted2-pimean['pr'], color = col, linestyle = '--')

    plt.grid()
    #plt.legend()
    ctl.custom_legend(fig, colors, allru, fontsize = 24, add_space_below = 0.03)
    plt.subplots_adjust(bottom = 0.2)
    plt.xlabel('GTAS anomaly (K)')
    plt.ylabel(r'Prec. anomaly (mm/day)')
    fig.savefig(cart_out + 'pr_gtas_scatter_linear.pdf')


    fig = plt.figure(figsize = (16,12))
    for ru, col in zip(allru, colors):
        if ru[0] == 'b':
            s = 1
        else:
            s = 10
        plt.scatter(glomeans[(ru, 'tas')][1]-pimean['tas'], glomeans[(ru, 'pr')][1]-pimean['pr'], color = col, s = s, label = ru)

        if ru == 'ssp585':
            pr = glomeans[(ru, 'pr')][1]
            pr_anom = fact*(pr-pimean['pr'])
            tas = glomeans[(ru, 'tas')][1]
            tas_anom = tas-pimean['tas']

            coeffs, covmat = np.polyfit(tas_anom, pr_anom, deg = 2, cov = True)
            x_nu = np.arange(-0.5, 10.1, 0.5)
            fitted = np.polyval(coeffs, x_nu)
            plt.plot(x_nu, fitted, color = col)

            #coeffs, covmat = np.polyfit(tas_anom, np.log(pr_anom), deg = 1, cov = True)
            #coeffs, covmat = np.polyfit(tas_anom, np.log(pr), deg = 2, cov = True)
            #x_nu = np.arange(-0.5, 9.6, 0.5)
            #fitted = np.exp(np.polyval(coeffs, x_nu))
            #plt.plot(x_nu, fitted-pimean['pr'], color = col)

            # fitted2 = A*exptas(x_nu+pimean['tas'])
            # plt.plot(x_nu, fitted2-pimean['pr'], color = col, linestyle = '--')

    plt.grid()
    #plt.legend()
    ctl.custom_legend(fig, colors, allru, fontsize = 24, add_space_below = 0.03)
    plt.subplots_adjust(bottom = 0.2)
    plt.xlabel('GTAS anomaly (K)')
    plt.ylabel(r'Prec. anomaly (mm/day)')
    fig.savefig(cart_out + 'pr_gtas_scatter_quadratic.pdf')

    fig, ax = plt.subplots(1, 1, figsize = (16,9))
    allruok = ['hist', 'ssp585', 'b990', 'b025', 'b050', 'b065', 'b080', 'b100']
    colok = [colors[allru.index(ru)] for ru in allruok]
    for i, (ru, col) in enumerate(zip(allruok, colok)):
        res = stats.linregress(glomeans[(ru, 'tas')][1], glomeans[(ru, 'pr')][1]/pimean['pr'])
        rel_trend = 100*res.slope
        rel_trend_err = 100*res.stderr
        ax.errorbar(i, rel_trend, yerr = rel_trend_err, color = col, capsize = 5, zorder = 6, lw = 4)
        ax.scatter(i, rel_trend, color = col, s = 200, label = ru)

        deltaT = np.mean(glomeans[(ru, 'tas')][1][-10:])-np.mean(glomeans[(ru, 'tas')][1][:10])
        print(ru, rel_trend, np.log(1 + rel_trend * deltaT)/deltaT)

    ax.grid(axis = 'y')
    #plt.legend()
    ax.set_xticks(range(len(allruok)))
    ax.set_xticklabels(allruok)

    ax.set_ylabel(r'Prec. trend (%/K)')
    fig.savefig(cart_out + 'pr_trend_scatter.pdf')

    fig = plt.figure(figsize = (16,9))
    for ru, col in zip(allru, colors):
        if ru in ['pi', 'hist', 'ssp585']: continue
        pr = glomeans[(ru, 'pr')][1]
        pr_anom = fact*(pr-pimean['pr'])
        tas_anom = glomeans[(ru, 'tas')][1]-pimean['tas']

        pr_anom_nl = pr - np.exp(np.polyval(coeffs, tas_anom))

        coeffs_ru, covmat_ru = np.polyfit(tas_anom, pr_anom_nl, deg = 2, cov = True)
        x_nu = np.arange(tas_anom.min()-0.1, tas_anom.max()+0.2, 0.1)
        fitted_ru = np.polyval(coeffs_ru, x_nu)
        #plt.plot(x_nu, fitted_ru, color = col, ls = '--')

        # coeffs_ru, covmat_ru = np.polyfit(tas_anom, np.log(pr_anom_nl), deg = 1, cov = True)
        # x_nu = np.arange(tas_anom.min()-0.1, tas_anom.max()+0.2, 0.1)
        # fitted_ru = np.exp(np.polyval(coeffs_ru, x_nu))
        # plt.plot(x_nu, fitted_ru, color = col, ls = '--')

        plt.scatter(tas_anom, pr_anom_nl, color = col, s = 5, label = ru)
        #plt.plot(x_nu, fitted_ru, color = col, lw = 2)

        plt.scatter(tas_anom[:50].mean(), pr_anom_nl[:50].mean(), facecolor = None, edgecolor = col, s = 50)
        plt.scatter(tas_anom[-50:].mean(), pr_anom_nl[-50:].mean(), facecolor = col, edgecolor = col, s = 50)

    plt.grid()
    #plt.legend()
    ctl.custom_legend(fig, colors[2:-1], allru[2:-1])
    plt.subplots_adjust(bottom = 0.2)
    plt.xlabel('GTAS anomaly (K)')
    plt.ylabel(r'Stabilization prec. anomaly (mm/day)')
    fig.savefig(cart_out + 'pr_gtas_nonlin_scatter.pdf')

    fig = plt.figure(figsize = (16,9))
    for ru, col in zip(allru, colors):
        if ru in ['pi', 'hist', 'ssp585']: continue
        pr = glomeans[(ru, 'pr')][1]
        pr_anom = fact*(pr-pimean['pr'])
        tas_anom = glomeans[(ru, 'tas')][1]-pimean['tas']

        # pr_anom_nl = pr_anom - np.polyval(coeffs, tas_anom)
        # coeffs_ru, covmat_ru = np.polyfit(tas_anom, pr_anom_nl, deg = 2, cov = True)
        # fitted_ru = np.polyval(coeffs_ru, tas_anom)

        # pr_anom_nl = pr_anom - np.exp(np.polyval(coeffs, tas_anom))
        pr_anom_nl = pr - np.exp(np.polyval(coeffs, tas_anom))

        pr_std = ctl.running_std(pr_anom_nl, 50)/np.sqrt(50-1)
        pr_me = ctl.running_mean(pr_anom_nl, 50)

        plt.fill_between(np.arange(len(pr_anom_nl)), pr_me-pr_std, pr_me+pr_std, color = col, alpha = 0.3)
        plt.plot(np.arange(len(pr_anom_nl)), pr_me, color = col, lw = 2)
        # plt.plot(np.arange(len(pr_anom_nl)), fitted_ru, color = col, lw = 0.5, ls = ':')
        # plt.plot(np.arange(len(pr_anom_nl)), ctl.running_mean(fitted_ru, 20), color = col, lw = 2)

    plt.grid()
    #plt.legend()
    ctl.custom_legend(fig, colors[2:-1], allru[2:-1])
    plt.subplots_adjust(bottom = 0.2)
    plt.xlabel('Years after stabilization')
    plt.ylabel(r'Stabilization prec. anomaly (mm/day)')
    fig.savefig(cart_out + 'pr_gtas_nonlin_scatter_years.pdf')


    #############################################
    # precipitation over Mediterranean
    med_area = [-8.0, 38.0, 27.0, 48.0]
    tropics = [0, 0, -20, 20]
    midlat_nh = [0, 0, 30, 60]
    midlat_sh = [0, 0, -60, -30]
    #med_area = [-5.0, 35.0, 30.0, 45.0]

    areas = [med_area, tropics, midlat_nh, midlat_sh]
    anams = ['Med', 'Trop', 'NML', 'SML']

    # for area, anam in zip(areas, anams):
    #     pino = yeamean[('pi', 'pr')]
    #     pino_med_pi = ctl.global_mean(ctl.sel_area_xr(pino, area)).mean()
    #
    #     fig = plt.figure(figsize = (16,9))
    #     for ru, col in zip(allru, colors):
    #         tas_anom = glomeans[(ru, 'tas')][1]-pimean['tas']
    #
    #         pino = yeamean[(ru, 'pr')]
    #         pino_med = ctl.global_mean(ctl.sel_area_xr(pino, area))
    #
    #         plt.scatter(tas_anom, pino_med - pino_med_pi, color = col, s = 10, label = ru)
    #
    #         if ru == 'ssp585':
    #             coeffs_ru, covmat_ru = np.polyfit(tas_anom, pino_med - pino_med_pi, deg = 1, cov = True)
    #             x_nu = np.arange(tas_anom.min()-0.1, tas_anom.max()+0.2, 0.1)
    #             fitted_ru = np.polyval(coeffs_ru, x_nu)
    #             plt.plot(x_nu, fitted_ru, color = col, ls = '-', lw = 2)
    #             # rume = ctl.running_mean(pino_med-pino_med_pi, 3s0)
    #             # teme = ctl.running_mean(tas_anom, 30)
    #             # plt.plot(teme, rume, color = col)
    #     plt.grid()
    #     #plt.legend()
    #     ctl.custom_legend(fig, colors, allru)
    #     plt.subplots_adjust(bottom = 0.2)
    #     plt.xlabel('GTAS anomaly (K)')
    #     plt.ylabel(r'Prec. anomaly (kg m$^{-2}$ s$^{-1}$)')
    #     plt.title('Area: {}'.format(anam))
    #
    #     fig.savefig(cart_out + 'pr_{}_gtas_scatter.pdf'.format(anam))



    for var in ['rlut', 'rsut', 'clt']:
        fig = plt.figure(figsize = (16,9))
        for ru, col in zip(allru, colors):
            plt.scatter(glomeans[(ru, 'tas')][1]-pimean['tas'], glomeans[(ru, var)][1]-pimean[var], color = col, s = 10, label = ru)

            if ru == 'ssp585':
                pr = glomeans[(ru, var)][1]
                pr_anom = pr-pimean[var]
                tas_anom = glomeans[(ru, 'tas')][1]-pimean['tas']

                # coeffs, covmat = np.polyfit(tas_anom, pr_anom, deg = 2, cov = True)
                # x_nu = np.arange(-0.5, 9.6, 0.5)
                # fitted = np.polyval(coeffs, x_nu)
                # plt.plot(x_nu, fitted, color = col)

                #coeffs, covmat = np.polyfit(tas_anom, np.log(pr_anom), deg = 1, cov = True)
                coeffs, covmat = np.polyfit(tas_anom, np.log(pr), deg = 1, cov = True)
                x_nu = np.arange(-0.5, 9.6, 0.5)
                fitted = np.exp(np.polyval(coeffs, x_nu))
                plt.plot(x_nu, fitted-pimean[var], color = col)

        plt.grid()
        #plt.legend()
        ctl.custom_legend(fig, colors, allru)
        plt.subplots_adjust(bottom = 0.2)
        plt.xlabel('GTAS anomaly (K)')
        if var == 'clt':
            plt.ylabel('Cloud cover anomaly')
        elif var == 'rsut':
            plt.ylabel(r'SW outgoing radiation at TOA ($W/m^2$)')
        elif var == 'rlut':
            plt.ylabel(r'LW outgoing radiation at TOA ($W/m^2$)')
        fig.savefig(cart_out + '{}_gtas_scatter.pdf'.format(var))


do_fb = True
if do_fb:
    var = 'net_toa'
    fig, ax = plt.subplots(figsize = (16,9))
    #fig2, ax2 = plt.subplots(figsize = (16,9))
    for ru, col in zip(allru[2:-1], colors[2:-1]):
        pr = glomeans[(ru, var)][1]

        pr_anom = pr-np.mean(pr[:10])
        #pr_anom = ctl.running_mean(pr_anom, 5, remove_nans = True)[::5]

        tas = glomeans[(ru, 'tas')][1]
        tas_anom = tas - tas[:10].mean()
        #tas_anom = ctl.running_mean(tas_anom, 5, remove_nans = True)[::5]

        # ordering for temp # NO! this is different than removing the fast manifold
        #gino = np.argsort(tas_anom)
        #tas_anom = tas_anom[gino]
        #pr_anom = pr_anom[gino]

        pr_anom = np.mean(np.split(pr_anom, int(len(pr_anom)/5)), axis = 1)
        tas_anom = np.mean(np.split(tas_anom, int(len(tas_anom)/5)), axis = 1)
        ax.scatter(tas_anom, pr_anom, color = col, s = 5)

        x_nu = np.arange(tas_anom.min()-0.2, tas_anom.max()+0.2, 0.1)
        # coeffs2, covmat = np.polyfit(tas_anom, pr_anom, deg = 2, cov = True)
        # print(coeffs2)
        # fitted = np.polyval(coeffs2, x_nu)
        # ax.plot(x_nu, fitted, color = col, ls = ':')

        coeffs, covmat = np.polyfit(tas_anom, pr_anom, deg = 1, cov = True)
        fitted = np.polyval(coeffs, x_nu)
        ax.plot(x_nu, fitted, color = col, lw = 2, label = ru)

        x_nu_ok = x_nu + tas[:10].mean() - pimean['tas']
        # ax2.plot(x_nu_ok, coeffs2[1]+2*coeffs2[0]*x_nu, color = col, lw = 2, label = ru)

        m, c, err_m, err_c = ctl.linear_regre_witherr(tas_anom, pr_anom)
        stron = r'{}: $\lambda = {:5.3f} \pm {:5.3f} \ W/m^2/K$'.format(ru, m, err_m)
        print(stron)
        ax.text(-0.3, -1.5-(allru.index(ru)-2)*0.2, stron, fontsize = 14)

    ax.grid()
    ax.legend()
    #ctl.custom_legend(fig, colors[2:-1], allru[2:-1])
    #ax.subplots_adjust(bottom = 0.2)
    ax.set_xlabel('GTAS change (K)')
    ax.set_ylabel(r'Change in net incoming TOA ($W/m^2$)')
    fig.savefig(cart_out + 'feedback_evolution.pdf'.format(var))

    # ax2.grid()
    # ax2.legend()
    # ax2.set_xlabel('GTAS (K)')
    # ax2.set_ylabel(r'feedback parameter ($W m^{-2} K^{-1}$)')
    # fig2.savefig(cart_out + 'feedback_evolution_poly2.pdf'.format(var))

    fig, ax = plt.subplots(figsize = (16,9))
    for ru, col in zip(allru[2:-1], colors[2:-1]):
        pr = glomeans[(ru, var)][1]

        pr_anom = pr-np.mean(pr[:10])
        #pr_anom = ctl.running_mean(pr_anom, 5, remove_nans = True)[::5]
        #pr_anom = np.mean(np.split(pr_anom, int(len(pr_anom)/5)), axis = 1)

        tas = glomeans[(ru, 'tas')][1]
        tas_anom = tas - tas[:10].mean()
        #tas_anom = np.mean(np.split(tas_anom, int(len(tas_anom)/5)), axis = 1)
        #tas_anom = ctl.running_mean(tas_anom, 5, remove_nans = True)[::5]

        x_nu = np.arange(tas_anom.min()-0.2, tas_anom.max()+0.2, 0.1)
        coeffs2, covmat = np.polyfit(tas_anom, pr_anom, deg = 2, cov = True)
        print(coeffs2)
        fitted = np.polyval(coeffs2, x_nu)
        ax.scatter(tas_anom, pr_anom, color = col, s = 5)
        ax.plot(x_nu, fitted, color = col, ls = ':')

        coeffs, covmat = np.polyfit(tas_anom, pr_anom, deg = 1, cov = True)
        fitted = np.polyval(coeffs, x_nu)
        ax.plot(x_nu, fitted, color = col, lw = 2, label = ru)

        x_nu_ok = x_nu + tas[:10].mean() - pimean['tas']
        #ax2.plot(x_nu_ok, coeffs2[1]+2*coeffs2[0]*x_nu, color = col, lw = 2, label = ru)

        m, c, err_m, err_c = ctl.linear_regre_witherr(tas_anom, pr_anom)
        stron = r'{}: $\lambda = {:5.3f} \pm {:5.3f} \ W/m^2/K$'.format(ru, m, err_m)
        print(stron)
        ax.text(-0.3, -2-(allru.index(ru)-2)*0.2, stron, fontsize = 14)

    ax.grid()
    ax.legend()
    #ctl.custom_legend(fig, colors[2:-1], allru[2:-1])
    #ax.subplots_adjust(bottom = 0.2)
    ax.set_xlabel('GTAS change (K)')
    ax.set_ylabel(r'Change in net incoming TOA ($W/m^2$)')
    fig.savefig(cart_out + 'feedback_evolution_yearly.pdf'.format(var))


    #
    # fig, ax = plt.subplots(figsize = (16,9))
    # for ru, col in zip(allru[2:-1], colors[2:-1]):
    #     pr = glomeans[(ru, var)][1]
    #
    #     pr_anom = pr-np.mean(pr[:10])
    #     #pr_anom = ctl.running_mean(pr_anom, 5, remove_nans = True)[::5]
    #     #pr_anom = np.mean(np.split(pr_anom, int(len(pr_anom)/5)), axis = 1)
    #
    #     tas = glomeans[(ru, 'tas')][1]
    #     tas_anom = tas - tas[:10].mean()
    #     #tas_anom = np.mean(np.split(tas_anom, int(len(tas_anom)/5)), axis = 1)
    #     #tas_anom = ctl.running_mean(tas_anom, 5, remove_nans = True)[::5]
    #
    #     x_nu = np.arange(tas_anom.min(), tas_anom.max(), 0.1)
    #
    #     pr_ok = []
    #     tas_ok = []
    #     for xi in x_nu:
    #         indi = (tas_anom >= xi) & (tas_anom < xi + 0.1)
    #         pr_ok.append(np.mean(pr_anom[indi]))
    #         tas_ok.append(np.mean(tas_anom[indi]))
    #
    #     pr_ok = np.array(pr_ok)
    #     tas_ok = np.array(tas_ok)
    #
    #     ax.scatter(tas_anom, pr_anom, color = col, s = 1)
    #     ax.scatter(tas_ok, pr_ok, color = col, s = 5, facecolor = 'none', marker = 'o')
    #
    #     coeffs, covmat = np.polyfit(tas_ok, pr_ok, deg = 1, cov = True)
    #     fitted = np.polyval(coeffs, x_nu)
    #     ax.plot(x_nu, fitted, color = col, lw = 2, label = ru)
    #
    #     coeffs, covmat = np.polyfit(tas_anom, pr_anom, deg = 1, cov = True)
    #     fitted = np.polyval(coeffs, x_nu)
    #     ax.plot(x_nu, fitted, color = col, lw = 1, label = ru, ls = '--')
    #
    #     m, c, err_m, err_c = ctl.linear_regre_witherr(tas_ok, pr_ok)
    #     stron = r'{}: $\lambda = {:5.3f} \pm {:5.3f} \ W/m^2/K$'.format(ru, m, err_m)
    #     print(stron)
    #     ax.text(-0.3, -2-(allru.index(ru)-2)*0.2, stron, fontsize = 14)
    #
    # ax.grid()
    # ax.legend()
    # #ctl.custom_legend(fig, colors[2:-1], allru[2:-1])
    # #ax.subplots_adjust(bottom = 0.2)
    # ax.set_xlabel('GTAS change (K)')
    # ax.set_ylabel(r'Change in net incoming TOA ($W/m^2$)')
    # fig.savefig(cart_out + 'feedback_evolution_tasbin.pdf'.format(var))
    ##################################################################

    fig_greg, ax_greg = plt.subplots(figsize = (16,9))

    for ru, col in zip(allru, colors):
        print(ru)
        checkdi = True
        if ru in ['hist', 'pi']: checkdi = False
        ctl.gregplot_on_ax(ax_greg, glomeans[(ru, 'tas')][1]-pimean['tas'], glomeans[(ru, 'net_toa')][1], color = col, label = ru, calc_ERF = False, calc_ECS = False, mean_nyr = False, check_diff=checkdi)

    #ax_greg.legend()
    ax_greg.grid()

    ax_greg.set_xlabel('Global mean tas (K)')
    ax_greg.set_ylabel(r'Global net incoming TOA flux ($W/m^2$)')
    fig_greg.savefig(cart_out + 'bottino_gregory_ok_yr.pdf')

    #################################################################

    fig, ax = plt.subplots(figsize = (16,9))
    for ru, col in zip(allru[2:-1], colors[2:-1]):
        pr = glomeans[(ru, var)][1]

        pr_anom = pr-np.mean(pr[:10])
        #pr_anom = ctl.running_mean(pr_anom, 5, remove_nans = True)[::5]
        tas_anom = glomeans[(ru, 'tas')][1]-pimean['tas']
        # tas_anom = tas - tas[:10].mean()
        #tas_anom = ctl.running_mean(tas_anom, 5, remove_nans = True)[::5]
        xs = []
        ys = []
        yerrs = []

        nspli = 9
        wind = 100
        shift = 25
        ini = 0

        # gino = np.argsort(tas_anom)
        # pr_anom = pr_anom[gino]
        # tas_anom = tas_anom[gino]
        #for spli in range(nspli):
        while ini <= 1000-wind:
            fin = ini + wind
            print(fin)

            in0 = max(0, ini-5)

            rutas = ctl.running_mean(tas_anom[in0:fin+5], 5, remove_nans = True)
            rupr = ctl.running_mean(pr_anom[in0:fin+5], 5, remove_nans = True)

            #m, c, err_m, err_c = ctl.linear_regre_witherr(tas_anom[ini:fin], pr_anom[ini:fin])
            m, c, err_m, err_c = ctl.linear_regre_witherr(rutas, rupr)
            stron = r'{}: $\lambda = {:5.3f} \pm {:5.3f} \ W/m^2/K$'.format(ru, m, err_m)
            print(stron)
            #ax.text(-0.3, -1.5-(allru.index(ru)-2)*0.2, stron, fontsize = 14)
            if ini == 0:
                lab = ru
            else:
                lab = None
            #plt.scatter(np.mean(tas_anom[ini:fin]), m, s = 10, color = col, label = lab)
            #xs.append(np.mean(tas_anom[ini:fin]))
            xs.append(np.mean([fin, ini]))
            ys.append(m)
            yerrs.append(err_m)
            ini = int(ini + shift)

        plt.errorbar(xs, ys, yerr = yerrs, color = col, label = lab, capsize = 2, lw = 2)
        plt.scatter(xs, ys, s = 20, color = col)#, label = lab)


    plt.grid()
    plt.legend()
    #ctl.custom_legend(fig, colors[2:-1], allru[2:-1])
    plt.subplots_adjust(bottom = 0.2)
    plt.xlabel('GTAS (K)')
    plt.ylabel(r'feedback parameter ($W m^{-2} K^{-1}$)')
    fig.savefig(cart_out + 'feedback_evolution_2.pdf'.format(var))


#########
## EOFs of the various variables
calc_eofs = False

if calc_eofs:
    lat = yeamean[('b100', 'tas')].lat.values
    lon = yeamean[('b100', 'tas')].lon.values
    first_eof = dict()
    for var in ['tas', 'pr', 'rlut', 'rsut', 'clt']:
        var_low_ssp = ctl.running_mean(yeamean[('ssp585', var)], 10, remove_nans=True)
        var_low = ctl.running_mean(yeamean[('b100', var)], 10, remove_nans=True)
        solver_ssp = ctl.eof_computation(var_low_ssp, latitude = lat)

        #for ru in ['b025', 'b050', 'b100']:
        solver = ctl.eof_computation(var_low, latitude = lat)
        first_eof[(var, 'ssp585')] = solver_ssp.eofs()[0]
        first_eof[(var, 'b100')] = solver.eofs()[0]

        if var == 'tas':
            cma = 'viridis'
            plotanom = False
        else:
            cma = 'RdBu_r'
            plotanom = True
        ctl.plot_multimap_contour([solver_ssp.eofs()[0], solver.eofs()[0]], lat, lon, cmap = cma, plot_anomalies = plotanom, figsize = (16,9), filename = cart_out + '{}_eof_sspvsb100.pdf'.format(var), color_percentiles = (2,98), subtitles = ['ssp585', 'b100'], cb_label = 'First EOF of {} (normalized)'.format(var))

plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
titlefont = 22
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.axisbelow'] = True
plt.rcParams['legend.fontsize'] = 18

###################################################
###########
carto = cart_out + '../ocean3d/'

masscell = xr.load_dataset(carto + 'masscello_Omon_EC-Earth3_stabilization-ssp585-2050_r1i1p1f1_gn_222201-222212.nc', use_cftime = True)['masscello']
areacell = xr.load_dataset(carto + 'areacello_Ofx_EC-Earth3_stabilization-ssp585-2050_r1i1p1f1_gn.nc', use_cftime = True)['areacello']
massvol = masscell*areacell
mavo = massvol.mean('time')
mavo_lev = mavo.sum(['i','j'])

oce_mass = 1.381107e+21 # global and vertical sum of masscello*areacello (year 2222)
ml_mass = 3.7352024e+19 # first 100 m globally
bulk_mass = 1.3436264e+21 # below 100 m globally

oce_area = 3.6481962e+14 # m2

cp0 = 3989.245 # J/kg/K

ml_depth = 150. # as in Gregory, 2000.
#max_depth = 2500. # as well
max_depth = 6000.
ml_mass = mavo_lev.sel(lev = slice(0.,ml_depth)).sum(['lev'])
bulk_mass = mavo_lev.sel(lev = slice(ml_depth, max_depth)).sum(['lev'])

t_oceall = dict()

read_ts = False
if read_ts:
    oht_all = pickle.load(open(carto + 'oht_ts_deep.p', 'rb'))
else:
    oht_all = dict()

if not read_ts:
    oht_lev = []
    filo = open(carto + 'oht_piControl.p', 'rb')
    for i in range(500):
        try:
            gigi = pickle.load(filo)
        except:
            break
        oht_lev.append(gigi[0])
    filo.close()

    oht_lev = xr.concat(oht_lev, dim = 'year')
    # oht_lev_pi = oht_lev.mean('year')*cp0
    # oht_lev_pi_std = oht_lev.std('year')*cp0
    oht_all[('pi', 700)] = cp0*oht_lev.sel(lev = slice(0., 700.)).sum('lev')
    oht_all[('pi', 2000)] = cp0*oht_lev.sel(lev = slice(700., 2000.)).sum('lev')
    oht_all[('pi', 'deep')] = cp0*oht_lev.sel(lev = slice(2000., 6000.)).sum('lev')

    oht_lev = oht_lev.sel(year = slice(70, 120))  ### CONSIDERING LAST 50 YEARS BEFORE BRANCHING OF HISTORICAL r4
    oht_lev_pi = oht_lev.mean('year')*cp0
    oht_lev_pi_std = oht_lev.std('year')*cp0
    oht_tot_pi = oht_lev_pi.sum('lev')
    t_deep_pi = 273.15 + oht_tot_pi/oce_mass/cp0

    oht1_pi = oht_lev_pi.sel(lev = slice(0., 700.)).sum('lev')
    oht2_pi = oht_lev_pi.sel(lev = slice(700., 2000.)).sum('lev')
    oht3_pi = oht_lev_pi.sel(lev = slice(2000., 6000.)).sum('lev')

    oht_ml_pi = oht_lev_pi.sel(lev = slice(0., ml_depth)).sum('lev')
    oht_bulk_pi = oht_lev_pi.sel(lev = slice(ml_depth, max_depth)).sum('lev')
    t_ml_pi = 273.15 + oht_ml_pi/ml_mass/cp0
    t_bulk_pi = 273.15 + oht_bulk_pi/bulk_mass/cp0
    t_oceall[('pi', 'ml')] = t_ml_pi
    t_oceall[('pi', 'bulk')] = t_bulk_pi

    oht_tot_pi = oht_lev_pi.sum('lev')

    pickle.dump([oht1_pi, oht2_pi, oht3_pi, oht_tot_pi, oht_ml_pi, oht_bulk_pi], open(carto + 'oht_ts_picontrol.p', 'wb'))
else:
    oht1_pi, oht2_pi, oht3_pi, oht_tot_pi, oht_ml_pi, oht_bulk_pi = pickle.load(open(carto + 'oht_ts_picontrol.p', 'rb'))
    t_deep_pi = 273.15 + oht_tot_pi/oce_mass/cp0
    t_ml_pi = 273.15 + oht_ml_pi/ml_mass/cp0
    t_bulk_pi = 273.15 + oht_bulk_pi/bulk_mass/cp0

oht_all[('pi', 't_prof')] = oht_lev_pi/mavo_lev/cp0
oht_all[('pi', 't_prof_std')] = oht_lev_pi_std/mavo_lev/cp0

fig, ax = plt.subplots(figsize = (16,9))
fig2, ax2 = plt.subplots(figsize = (16,9))
fig3, ax3 = plt.subplots(figsize = (16,9))
for ru, col in zip(allru[2:-1], colors[2:-1]):
    if not read_ts:
        oht_lev = []
        filo = open(carto + 'oht_{}_1000.p'.format(ru), 'rb')
        for i in range(1000):
            try:
                gigi = pickle.load(filo)
            except:
                break
            oht_lev.append(gigi[0])

        if ru == 'b065':
            print('AAAAA MISSING 1 YEAR FOR b065!! copying last year twice for now')
            oht_lev.append(gigi[0])

        filo.close()

        oht_lev = xr.concat(oht_lev, dim = 'year')*cp0
        oht_ml = oht_lev.sel(lev = slice(0., ml_depth)).sum('lev')
        oht_bulk = oht_lev.sel(lev = slice(ml_depth, max_depth)).sum('lev')

        oht_lev = oht_lev - oht_lev_pi # removing pi base level
        oht1 = oht_lev.sel(lev = slice(0., 700.)).sum('lev')
        oht2 = oht_lev.sel(lev = slice(700., 2000.)).sum('lev')
        oht3 = oht_lev.sel(lev = slice(2000., 6000.)).sum('lev')
        oht_tot = oht_lev.sum('lev')

        oht_all[(ru, 't_ini')] = oht_lev[:30].mean('year')/mavo_lev/cp0
        oht_all[(ru, 't_final')] = oht_lev[-30:].mean('year')/mavo_lev/cp0

        oht_all[(ru, 700)] = oht1
        oht_all[(ru, 2000)] = oht2
        oht_all[(ru, 'deep')] = oht3
        oht_all[(ru, 'ml')] = oht_ml
        oht_all[(ru, 'bulk')] = oht_bulk
        oht_all[(ru, 'tot')] = oht_tot

    else:
        oht1 = oht_all[(ru, 700)]
        oht2 = oht_all[(ru, 2000)]
        oht3 = oht_all[(ru, 'deep')]
        oht_tot = oht_all[(ru, 'tot')]
        oht_ml = oht_all[(ru, 'ml')]
        oht_bulk = oht_all[(ru, 'bulk')]

    t_deep = 273.15 + oht_tot/oce_mass/cp0

    t_ml = 273.15 + oht_ml/ml_mass/cp0
    t_bulk = 273.15 + oht_bulk/bulk_mass/cp0
    t_oceall[(ru, 'ml')] = t_ml
    t_oceall[(ru, 'bulk')] = t_bulk

    gtas = glomeans[(ru, 'tas')][1]
    yeas = np.arange(1000)
    # if ru == 'b025':
    #     gtas = gtas[5:]
    #     yeas = yeas[5:]

    ax.scatter(gtas-273.15, t_deep-273.15, s = 5, color = col, label = ru)
    grun = ctl.running_mean(gtas, 5, remove_nans = True)
    larun = ctl.running_mean(t_deep, 5, remove_nans = True)
    # ax.plot(grun, larun, color = col, label = ru)
    # x_nu = np.arange(gtas.min(), gtas.max(), 0.1)
    # coeffs, covmat = np.polyfit(grun, larun, deg = 2, cov = True)
    # fitted = np.polyval(coeffs, x_nu)
    # ax.plot(x_nu, fitted, color = col, label = ru, lw = 2)

    ax2.plot(oht1, color = col, ls = '-', label = ru, lw = 2)
    ax2.plot(oht2, color = col, ls = '--', lw = 2)
    ax2.plot(oht3, color = col, ls = ':', lw = 2)

    grun = ctl.running_mean(gtas, 10, remove_nans = True)
    oht1l = ctl.running_mean(oht1, 10, remove_nans = True)
    oht2l = ctl.running_mean(oht2, 10, remove_nans = True)
    oht3l = ctl.running_mean(oht3, 10, remove_nans = True)

    # ax3.scatter(grun, oht1l, color = col, label = ru, s = 20)
    # ax3.scatter(grun, oht2l, edgecolor = col, s = 20, marker = 'd', facecolor = 'none')
    # ax3.scatter(grun, oht3l, color = col, s = 20, marker = '*')
    ax3.plot(grun, oht1l, color = col, label = ru, lw = 2)
    ax3.plot(grun, oht2l, color = col, ls = '--', lw = 2)
    ax3.plot(grun, oht3l, color = col, ls = ':', lw = 2)

ax.legend()
ax.grid()
#ax.set_ylabel(r'$(1-T_d/T_s)$')
ax.set_ylabel('T bulk ocean (K)')
ax.set_xlabel('GTAS (K)')

ax2.legend()
ax2.grid()
ax2.set_ylabel('OHC anomaly (J)')
ax2.set_xlabel('Years after stabilization')

ax3.legend()
ax3.grid()
ax3.set_ylabel('OHC anomaly (J)')
ax3.set_xlabel('GTAS (K)')

fig.savefig(carto + 'lambda_factor.pdf')
fig2.savefig(carto + 'oht_deep_time.pdf')
fig3.savefig(carto + 'oht_deep_gtas.pdf')

pickle.dump(oht_all, open(carto + 'oht_ts_deep.p', 'wb'))

####
fig, ax = plt.subplots(figsize = (16,9))
for ru, col in zip(allru[2:-1], colors[2:-1]):
    ax.plot(oht_all[(ru, 't_final')], oht_all[(ru, 't_final')].lev, color = col, label = ru, lw = 2)
ax.legend()
ax.grid()
ax.invert_yaxis()
ax.set_ylabel('Depth (m)')
#ax.set_xlabel('Heat content anomaly (J)')
ax.set_xlabel('Cons. temperature anomaly (K)')
fig.savefig(carto + 'otemp_anom_profile.pdf')

####
fig, ax = plt.subplots(figsize = (16,9))
piprof = oht_all[('pi', 't_prof')]
piprof_std = oht_all[('pi', 't_prof_std')]
ax.fill_betweenx(piprof.lev, piprof-piprof_std, piprof+piprof_std, color = 'grey', alpha = 0.3)
ax.plot(piprof, piprof.lev, color = 'black', label = 'pi', lw = 2)
for ru, col in zip(allru[2:-1], colors[2:-1]):
    ax.plot(piprof+oht_all[(ru, 't_ini')], oht_all[(ru, 't_final')].lev, color = col, lw = 2, ls = ':')
for ru, col in zip(allru[2:-1], colors[2:-1]):
    ax.plot(piprof+oht_all[(ru, 't_final')], oht_all[(ru, 't_final')].lev, color = col, label = ru, lw = 2)
ax.legend()
ax.grid()
ax.invert_yaxis()
ax.set_ylabel('Depth (m)')
#ax.set_xlabel('Heat content anomaly (J)')
ax.set_xlabel('Cons. temperature (K)')
fig.savefig(carto + 'otemp_abs_profile.pdf')

#######################################################
## Conto stima vertical heat diffusion
print('Conto stima heat diffusion')
fig, ax = plt.subplots(figsize = (16,9))
for ru, col in zip(allru[2:-1], colors[2:-1]):
    tgrad_ini = (piprof+oht_all[(ru, 't_ini')]).differentiate('lev')
    tgrad_fin = (piprof+oht_all[(ru, 't_final')]).differentiate('lev')
    ax.plot(tgrad_ini, oht_all[(ru, 't_final')].lev, color = col, lw = 2, ls = ':')
    ax.plot(tgrad_fin, oht_all[(ru, 't_final')].lev, color = col, lw = 2, label = ru)

ax.grid()
ax.invert_yaxis()
ax.set_ylim(900., 3000.)
ax.set_xlim(-0.01, 0.)
ax.invert_yaxis()
ax.set_ylabel('Depth (m)')
#ax.set_xlabel('Heat content anomaly (J)')
ax.set_xlabel('Cons. temperature gradient (K/m)')
fig.savefig(carto + 'otemp_gradient_profile.pdf')

#######################################################

fig, ax = plt.subplots(figsize = (12,8))
ax.plot(mavo_lev, mavo_lev.lev, lw = 2)
ax.grid()
ax.invert_yaxis()
ax.set_ylabel('Depth (m)')
ax.set_xlabel('Ocean mass (kg)')
fig.savefig(carto + 'ocemass_profile.pdf')

fig, ax = plt.subplots(figsize = (12,8))
for ru, col in zip(allru[2:-1], colors[2:-1]):
    ax.plot(oht_all[(ru, 't_ini')], oht_all[(ru, 't_ini')].lev, color = col, lw = 2, ls = ':')
for ru, col in zip(allru[2:-1], colors[2:-1]):
    ax.plot(oht_all[(ru, 't_final')], oht_all[(ru, 't_final')].lev, color = col, label = ru, lw = 2)
ax.legend()
ax.grid()
ax.invert_yaxis()
ax.set_ylabel('Depth (m)')
#ax.set_xlabel('Heat content anomaly (J)')
ax.set_xlabel('Cons. temperature anomaly (K)')
fig.savefig(carto + 'otemp_anom_vsini_profile.pdf')


sys.exit()

### Plot efficiency of heat transfer, as in Armour (2017)
nlow = 50
fig, ax = plt.subplots(figsize = (16,9))
for ru, col in zip(allru[2:-1], colors[2:-1]):
    oht_ml = oht_all[(ru, 'ml')] - oht_ml_pi
    oht_bulk = oht_all[(ru, 'bulk')] - oht_bulk_pi
    oht_ml = ctl.running_mean(oht_ml, nlow)
    oht_bulk = ctl.running_mean(oht_bulk, nlow)

    t_ml = t_oceall[(ru, 'ml')] - t_oceall[('pi', 'ml')]
    t_bulk = t_oceall[(ru, 'bulk')] - t_oceall[('pi', 'bulk')]

    t_ml = ctl.running_mean(t_ml, nlow)
    t_bulk = ctl.running_mean(t_bulk, nlow)

    gtas = glomeans[(ru, 'tas')][1]-pimean['tas']
    gtas = ctl.running_mean(gtas, nlow)
    #gtas = 273.15 + oht1/oce_mass/cp0
    yeas = np.arange(1000)
    dtbulk = np.gradient(oht_bulk)/(86400*365)
    #gama = dtbulk/(gtas-t_deep)
    gama = dtbulk/(t_ml-t_bulk)/oce_area
    ax.scatter(gtas, gama, label = ru, color = col)
    #plt.scatter(yeas, gama, label = ru, color = col, s = 2)

ax.grid()
ax.legend()
#ax.set_title('Heat uptake efficiency')
ax.set_ylabel(r'$\gamma$ ($W m^{-2} K^{-1}$)')
ax.set_xlabel('GTAS (K)')
fig.savefig(carto + 'gamma_vs_gtas.pdf')


### Plot efficiency of heat transfer, as in Armour (2017)
nlow = 50
fig, ax = plt.subplots(figsize = (16,9))
for ru, col in zip(allru[2:-1], colors[2:-1]):
    oht_ml = oht_all[(ru, 'ml')] - oht_ml_pi
    oht_bulk = oht_all[(ru, 'bulk')] - oht_bulk_pi
    oht_ml = ctl.running_mean(oht_ml, nlow)
    oht_bulk = ctl.running_mean(oht_bulk, nlow)

    t_ml = t_oceall[(ru, 'ml')] - t_oceall[('pi', 'ml')]
    t_bulk = t_oceall[(ru, 'bulk')] - t_oceall[('pi', 'bulk')]

    t_ml = ctl.running_mean(t_ml, nlow)
    t_bulk = ctl.running_mean(t_bulk, nlow)

    gtas = glomeans[(ru, 'tas')][1]-pimean['tas']
    gtas = ctl.running_mean(gtas, nlow)
    #gtas = 273.15 + oht1/oce_mass/cp0
    yeas = np.arange(1000)
    dtbulk = np.gradient(oht_bulk)/(86400*365)
    #gama = dtbulk/(gtas-t_deep)
    gama = dtbulk/(t_ml-t_bulk)/oce_area
    #ax.scatter(gtas, gama, label = ru, color = col)
    ax.scatter(t_ml-t_bulk, gama, label = ru, color = col)
    #plt.scatter(yeas, gama, label = ru, color = col, s = 2)

ax.grid()
ax.legend()
#ax.set_title('Heat uptake efficiency')
ax.set_ylabel(r'$\gamma$ ($W m^{-2} K^{-1}$)')
ax.set_xlabel(r'$T_s-T_d$ (K)')
fig.savefig(carto + 'gamma_vs_deltaT.pdf')

### Plot efficiency of heat transfer, as in Armour (2017)
nlow = 50
fig, ax = plt.subplots(figsize = (16,9))
for ru, col in zip(allru[2:-1], colors[2:-1]):
    oht_ml = oht_all[(ru, 'ml')] - oht_ml_pi
    oht_bulk = oht_all[(ru, 'bulk')] - oht_bulk_pi
    oht_ml = ctl.running_mean(oht_ml, nlow)
    oht_bulk = ctl.running_mean(oht_bulk, nlow)

    t_ml = t_oceall[(ru, 'ml')] - t_oceall[('pi', 'ml')]
    t_bulk = t_oceall[(ru, 'bulk')] - t_oceall[('pi', 'bulk')]

    t_ml = ctl.running_mean(t_ml, nlow)
    t_bulk = ctl.running_mean(t_bulk, nlow)

    gtas = glomeans[(ru, 'tas')][1]-pimean['tas']
    gtas = ctl.running_mean(gtas, nlow)
    #gtas = 273.15 + oht1/oce_mass/cp0
    yeas = np.arange(1000)
    dtbulk = np.gradient(oht_bulk)/(86400*365)
    #gama = dtbulk/(gtas-t_deep)
    gama = dtbulk/(t_ml-t_bulk)/oce_area
    #ax.scatter(gtas, gama, label = ru, color = col)
    ax.scatter(yeas, gama, label = ru, color = col)
    #plt.scatter(yeas, gama, label = ru, color = col, s = 2)

ax.grid()
ax.legend()
#ax.set_title('Heat uptake efficiency')
ax.set_ylabel(r'$\gamma$ ($W m^{-2} K^{-1}$)')
ax.set_xlabel('Years after stabilization')
fig.savefig(carto + 'gamma_vs_time.pdf')

### Plot efficiency of heat transfer, as in Armour (2017)
nlow = 50
fig, ax = plt.subplots(figsize = (16,9))
for ru, col in zip(allru[2:-1], colors[2:-1]):
    oht_ml = oht_all[(ru, 'ml')] - oht_ml_pi
    oht_bulk = oht_all[(ru, 'bulk')] - oht_bulk_pi
    oht_ml = ctl.running_mean(oht_ml, nlow)
    oht_bulk = ctl.running_mean(oht_bulk, nlow)

    t_ml = t_oceall[(ru, 'ml')] - t_oceall[('pi', 'ml')]
    t_bulk = t_oceall[(ru, 'bulk')] - t_oceall[('pi', 'bulk')]
    t_ml = ctl.running_mean(t_ml, nlow)
    t_bulk = ctl.running_mean(t_bulk, nlow)

    yeas = np.arange(1000)
    dtbulk = np.gradient(oht_bulk)/(86400*365)

    nonan = np.isnan(dtbulk)
    res = stats.linregress(t_ml[~nonan]-t_bulk[~nonan], dtbulk[~nonan])
    trend = res.slope
    trend_err = res.stderr
    print(ru + 'Corr: ', res.rvalue)

    gtas = glomeans[(ru, 'tas')][1]-pimean['tas']
    gtas = ctl.running_mean(gtas, nlow)
    #gtas = 273.15 + oht1/oce_mass/cp0
    #gama = dtbulk/(gtas-t_deep)
    gama = dtbulk/(t_ml-t_bulk)/oce_area
    #ax.scatter(gtas, gama, label = ru, color = col)
    #ax.scatter(yeas, gama, label = ru, color = col)
    ax.scatter(t_ml-t_bulk, dtbulk/oce_area, label = ru, color = col)

    #plt.scatter(yeas, gama, label = ru, color = col, s = 2)

ax.grid()
ax.legend()
#ax.set_title('Heat uptake efficiency')
#ax.set_ylabel(r'$\gamma/A_{oce}$ ($W m^{-2} K^{-1}$)')
ax.set_ylabel(r'$C_d \frac{dT_d}{dt}$ (W/$m^2$)')
ax.set_xlabel(r'$T_s-T_d$ (K)')
fig.savefig(carto + 'greg2000_dt_vs_tdiff.pdf')
#############################################################


fig, axs = plt.subplots(2, 3, figsize = (16,9))

for ru, col, ax in zip(allru[2:-1], colors[2:-1], axs.flatten()):
    oht1 = oht_all[(ru, 700)]
    oht2 = oht_all[(ru, 2000)]
    oht3 = oht_all[(ru, 'deep')]
    oht_tot = oht_all[(ru, 'tot')]

    ax.fill_between(np.arange(len(oht3)), np.zeros(len(oht3)), oht3, color = 'steelblue', alpha = 0.5)
    ax.fill_between(np.arange(len(oht3)), oht3, oht3+oht2, color = 'forestgreen', alpha = 0.5)
    ax.fill_between(np.arange(len(oht3)), oht3+oht2, oht_tot, color = 'gold', alpha = 0.5)

    ax.plot(oht3, color = 'steelblue', ls = '-', lw = 2)
    ax.plot(oht3+oht2, color = 'forestgreen', ls = '-', lw = 2)
    ax.plot(oht_tot, color = 'gold', ls = '-', lw = 2)

    ax.set_title(ru)
    ax.grid()

axs[0,0].set_ylabel('OHC anomaly (J)')
axs[1,0].set_ylabel('OHC anomaly (J)')

for j in range(3):
    axs[0,j].set_xticklabels([])
    axs[1,j].set_xlabel('Years after stabilization')

for j in range(1, 3):
    axs[0,j].set_yticklabels([])
    axs[1,j].set_yticklabels([])

ctl.adjust_ax_scale(axs.flatten())
ctl.custom_legend(fig, ['steelblue', 'forestgreen', 'gold'], ['> 2000m', '700-2000m', '0-700m'], ncol = 3)
plt.subplots_adjust(bottom = 0.25)
fig.savefig(carto + 'oht_deep_all_evol.pdf')


fig, axs = plt.subplots(2, 3, figsize = (16,9))

for ru, col, ax in zip(allru[2:-1], colors[2:-1], axs.flatten()):
    oht_tot = oht_all[(ru, 'tot')]
    oht1 = oht_all[(ru, 700)]/oht_tot
    oht2 = oht_all[(ru, 2000)]/oht_tot
    oht3 = oht_all[(ru, 'deep')]/oht_tot

    for ohcos in [oht1, oht2, oht3, oht_tot]:
        print(np.mean(ohcos).values)

    ax.fill_between(np.arange(len(oht3)), np.zeros(len(oht3)), oht3, color = 'steelblue', alpha = 0.5)
    ax.fill_between(np.arange(len(oht3)), oht3, oht3+oht2, color = 'forestgreen', alpha = 0.5)
    ax.fill_between(np.arange(len(oht3)), oht3+oht2, 1., color = 'gold', alpha = 0.5)

    ax.plot(oht3, color = 'steelblue', ls = '-', lw = 2)
    ax.plot(oht3+oht2, color = 'forestgreen', ls = '-', lw = 2)
    #ax.plot(oht_tot, color = 'gold', ls = '-', lw = 2)

    ax.set_title(ru)
    ax.grid()

ax.text(0.05, 0.55, 'Fraction of heat absorbed', horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 18)
#axs[0,0].set_ylabel('Fraction of heat absorbed')
#axs[1,0].set_ylabel('Fraction of heat absorbed')

for j in range(3):
    axs[0,j].set_xticklabels([])
    axs[1,j].set_xlabel('Years after stabilization')

for j in range(1, 3):
    axs[0,j].set_yticklabels([])
    axs[1,j].set_yticklabels([])

ctl.adjust_ax_scale(axs.flatten())
ctl.custom_legend(fig, ['steelblue', 'forestgreen', 'gold'], ['> 2000m', '700-2000m', '0-700m'], ncol = 3)
plt.subplots_adjust(bottom = 0.25)
fig.savefig(carto + 'oht_deep_all_evol_rel.pdf')

## Trends in ocean heat content, per layer
filo = open(cart_out + 'trend_ohc.txt', 'w')
filo.write('Trend of OHC in last 100 years:\n')
for cos in [700, 2000, 'deep', 'ml', 'bulk']:
    filo.write('Layer: '+str(cos)+'\n')
    for ru in allru:
        if 'b' in ru:
            yea = np.arange(len(oht_all[(ru, 700)]))
            res = stats.linregress(yea[-100:], oht_all[(ru, cos)][-100:]) # trend of last 100 year
            rel_trend = 100*res.slope
            rel_trend_err = 100*res.stderr

            print('{} {} - {:10.3e} +/- {:10.3e} J/cent\n'.format(ru, cos, rel_trend, rel_trend_err))
            filo.write('  {} - {:10.3e} +/- {:10.3e} J/cent\n'.format(ru, rel_trend, rel_trend_err))

filo.close()

for nom, lev in zip(['Upper (< 700 m)', 'Mid (700-2000 m)', 'Deep (> 2000 m)', 'Mixed-layer (<100 m)', 'Bulk'], [700, 2000, 'deep', 'ml', 'bulk']):
    cose = [100*np.mean(oht_all[(ru, lev)][-30:]).values/np.mean(oht_all[(ru, 'tot')][-30:]).values for ru in allru if 'b' in ru]
    strin = ' '+6*'& {:5.0f}\%'+'\\\\'
    print(nom+strin.format(*cose))

for nom, lev in zip(['Upper (< 700 m)', 'Mid (700-2000 m)', 'Deep (> 2000 m)', 'Mixed-layer (<100 m)', 'Bulk', 'Tot'], [700, 2000, 'deep', 'ml', 'bulk', 'tot']):
    cose = [np.mean(oht_all[(ru, lev)][-30:]).values/1.e24 for ru in allru if 'b' in ru]
    strin = ' '+6*'& {:5.2f}'+'\\\\'
    print(nom+strin.format(*cose))

#############################################################

## historical r4 starts at year 2380 of piControl. piControl (official r1) starts from 2259, so historical r4 starts at year 121.


## Check deep trend pi
filo = open(carto + 'oht_piControl.p', 'rb')
oht_lev = []
for i in range(500):
    try:
        gigi = pickle.load(filo)
    except:
        break
    oht_lev.append(gigi[0])
filo.close()

oht_lev = xr.concat(oht_lev, dim = 'year')*cp0
oht3 = oht_lev.sel(lev = slice(2000., 6000.)).sum('lev')
oht3_pi_mean = oht3.mean('year')

for lay, tit, levok in zip([700, 2000, 'deep'], ['0-700m', '700-2000m', '> 2000 m'], [(0, 700), (700, 2000), (2000, 6000)]):
    fig = plt.figure(figsize = (16,9))
    plt.plot(oht3-oht3[:30].mean(), color = 'black', label = 'pi', lw = 2)
    for ru, col in zip(allru, colors):
        if 'b' in ru:
            plt.plot(oht_all[(ru, lay)] + oht3_pi_mean - oht3[:30].mean(), color = col, label = ru, lw = 2)
    
    if lay =='deep':
        plt.ylabel('Deep ocean heat content anomaly (J)')
    else:
        plt.ylabel('Ocean heat content anomaly (J) in: '+tit)

    plt.xlabel('years')
    plt.legend()
    fig.savefig(carto + 'drift_{}_ocean_vs_pi.pdf'.format(lay))

    deep_mass = mavo_lev.sel(lev = slice(levok[0], levok[1])).sum('lev')

    fig, ax = plt.subplots(figsize = (16,9))

    i=0
    res = stats.linregress(np.arange(500), oht3/deep_mass/cp0)
    ax.scatter(i, 100*res.slope, marker = 'D', color = 'black', label = 'pi', s = 100)
    for ru, col in zip(allru, colors):
        if 'b' in ru:
            i += 1
            res = stats.linregress(np.arange(1000), oht_all[(ru, lay)]/deep_mass/cp0)
            ax.scatter(i, 100*res.slope, marker = 'D', color = col, label = ru, s = 100)


    if lay =='deep':
        plt.ylabel('Deep ocean temperature trend (K/cent)')
    else:
        plt.ylabel('Ocean temperature trend (K/cent) in: '+tit)
    
    ax.set_ylim([0., 0.12])
    ax.set_xticks(np.arange(i+1))
    ax.set_xticklabels(['pi'] + [ru for ru in allru if 'b' in ru])
    ax.grid(axis = 'y')

    plt.legend()
    fig.savefig(carto + 'drift_{}_ocean_vs_pi_trend.pdf'.format(lay))


sys.exit()
#############################################################
### maps of OHT trends
lats = np.linspace(-89.5, 89.5, 180)
lons = np.linspace(0, 359, 360)

refoht = dict() #### THIS IS now assigned to b025, but should be pi
oht_patt = dict()

areacello = xr.load_dataset(carto + 'areacello.nc')['areacello']
regr_acello = ctl.regrid_dataset(areacello, lats, lons)

for ru, col in zip(allru[2:-1], colors[2:-1]):
    print(ru)
    filo = open(carto + 'oht_{}_1000.p'.format(ru), 'rb')

    oht700 = []
    oht2000 = []
    ohtdeep = []
    for i in range(1000):
        try:
            oht_lev_i, oht700_i, oht2000_i, ohtdeep_i = pickle.load(filo)
        except:
            break
        oht700.append(oht700_i)
        oht2000.append(oht2000_i)
        ohtdeep.append(ohtdeep_i)

    if ru == 'b065':
        print('AAAAA MISSING 1 YEAR FOR b065!! copying last year twice for now')
        oht700.append(oht700_i)
        oht2000.append(oht2000_i)
        ohtdeep.append(ohtdeep_i)

    filo.close()

    oht700 = np.stack(oht700)
    oht2000 = np.stack(oht2000)
    ohtdeep = np.stack(ohtdeep)

    # oht700 = oht700/regr_acello.values
    # oht2000 = oht2000/regr_acello.values
    # ohtdeep = ohtdeep/regr_acello.values

    for var, lab in zip([oht700, oht2000, ohtdeep], [700, 2000, 'deep']):
        var[var == 0.0] = np.nan

        var_trend, var_intercept, var_trend_err, var_intercept_err, var_pval = ctl.calc_trend_climatevar(np.arange(len(var)), var)

        oht_patt[(ru, lab)] = var_trend
        oht_patt[(ru, lab, 'err')] = var_trend_err
        oht_patt[(ru, lab, 'pval')] = var_pval

        nonan = ~np.isnan(var[0])
        gmea = ctl.global_mean(var, lats, mask = nonan)

        var_trend, var_intercept, var_trend_err, var_intercept_err, var_pval = ctl.calc_trend_climatevar(gmea, var)

        oht_patt[(ru, lab, 'rel')] = var_trend
        oht_patt[(ru, lab, 'rel_err')] = var_trend_err
        oht_patt[(ru, lab, 'rel_pval')] = var_pval

        var_trend, var_intercept, var_trend_err, var_intercept_err, var_pval = ctl.calc_trend_climatevar(np.arange(100), var[:100])

        oht_patt[(ru, lab, '100')] = var_trend
        oht_patt[(ru, lab, '100_err')] = var_trend_err
        oht_patt[(ru, lab, '100_pval')] = var_pval

        var_trend, var_intercept, var_trend_err, var_intercept_err, var_pval = ctl.calc_trend_climatevar(np.arange(100), var[-100:])

        oht_patt[(ru, lab, '1000')] = var_trend
        oht_patt[(ru, lab, '1000_err')] = var_trend_err
        oht_patt[(ru, lab, '1000_pval')] = var_pval


        if lab not in refoht:
            refoht[lab] = var[:20].mean(axis = 0) # taking b025 as reference

        oht_patt[(ru, lab, 'change')] = var[-20:].mean(axis = 0) - refoht[lab]

pickle.dump(oht_patt, open(carto + 'temp_patt_deep.p', 'wb'))
oht_patt = pickle.load(open(carto + 'temp_patt_deep.p', 'rb'))

#for lev, tit in zip([700, 2000, 'deep'], ['0-700m', '700-2000m', '> 2000 m']):
for lev, tit in zip([700, 2000, 'deep'], ['700 m', '2000 m', '4000 m']):
    plpa = []
    subt = []
    hatch = []
    for ru, col in zip(allru[2:-1], colors[2:-1]):
        plpa.append(oht_patt[(ru, lev)])
        hatch.append(oht_patt[(ru, lev, 'pval')] < 0.05)
        subt.append(ru + ': ' + tit)

    [fig] = ctl.plot_multimap_contour(plpa, lats, lons, visualization = 'Robinson', central_lat_lon = (0., -120.), filename = carto + 'temp_patt_ocean_{}.pdf'.format(lev), subtitles = subt, plot_anomalies = False, cmap = 'viridis', figsize = (16,9), fix_subplots_shape = (2,3), cb_label = 'Conservative temperature trend (K/yr)')#, cbar_range = (0., 0.008) , add_hatching = hatch, hatch_styles = ['///', '', ''])

    for ax in fig.axes[:-1]:
        ax.set_facecolor('gainsboro')

    fig.savefig(carto + 'temp_patt_ocean_{}.pdf'.format(lev))

    plpa = []
    subt = []
    hatch = []
    for ru, col in zip(allru[2:-1], colors[2:-1]):
        plpa.append(oht_patt[(ru, lev, 'rel')])
        hatch.append(oht_patt[(ru, lev, 'rel_pval')] < 0.05)
        subt.append(ru + ': ' + tit)

    #divnorm = colors.TwoSlopeNorm(vmin=-1., vcenter=1., vmax=3.5)
    [fig] = ctl.plot_multimap_contour(plpa, lats, lons, visualization = 'Robinson', central_lat_lon = (0., -120.), filename = carto + 'temp_patt_ocean_{}_rel.pdf'.format(lev), subtitles = subt, plot_anomalies = False, cmap = ctl.heatmap(), cbar_range = (-1, 3), figsize = (16,9), fix_subplots_shape = (2, 3), cb_label = 'Cons. temp. trend pattern', n_color_levels = 37)#, add_hatching = hatch, hatch_styles = ['///', '', ''])
    for ax in fig.axes[:-1]:
        ax.set_facecolor('gainsboro')
    #fig.savefig(carto + 'oht_patt_deep_rel.pdf')
    fig.savefig(carto + 'temp_patt_ocean_{}_rel.pdf'.format(lev))

    plpa = []
    subt = []
    hatch = []
    for ru, col in zip(allru[2:-1], colors[2:-1]):
        plpa.append(oht_patt[(ru, lev, 'change')])
        subt.append(ru + ': ' + tit)

    [fig] = ctl.plot_multimap_contour(plpa, lats, lons, visualization = 'Robinson', central_lat_lon = (0., -120.), filename = carto + 'temp_patt_ocean_{}_change.pdf'.format(lev), subtitles = subt, plot_anomalies = False, cmap = 'viridis', figsize = (16,9), fix_subplots_shape = (2, 3), cb_label = 'Conservative temperature change (K)') #, cbar_range = (0., 5.)

    for ax in fig.axes[:-1]:
        ax.set_facecolor('gainsboro')
    fig.savefig(carto + 'temp_patt_ocean_{}_change.pdf'.format(lev))


    plpa = []
    subt = []
    hatch = []
    for ru, col in zip(allru[2:-1], colors[2:-1]):
        plpa.append(oht_patt[(ru, lev, '1000')]-oht_patt[(ru, lev, '100')])
        subt.append(ru + ': ' + tit)

    [fig] = ctl.plot_multimap_contour(plpa, lats, lons, visualization = 'Robinson', central_lat_lon = (0., -120.), filename = carto + 'temp_patt_ocean_{}_1000-100.pdf'.format(lev), subtitles = subt, plot_anomalies = True, cmap = ctl.heatmap(), n_color_levels = 37, figsize = (16,9), fix_subplots_shape = (2,3), cb_label = 'Cons. temp. trend difference (K/yr) (last - first century)') #cbar_range = (-0.02, 0.02)

    for ax in fig.axes[:-1]:
        ax.set_facecolor('gainsboro')
    fig.savefig(carto + 'temp_patt_ocean_{}_1000-100.pdf'.format(lev))
