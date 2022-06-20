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

cart_out += 'nonlin_evol/'
ctl.mkdir(cart_out)

colors = ['black', 'steelblue', 'teal', 'forestgreen', 'orange', 'violet', 'indianred']
allru = ['pi', 'hist', 'b990', 'b025', 'b050', 'b100', 'ssp585']

####################################################################################################

cart_in = cart_out + '../seasmean/'
gogo = pickle.load(open(cart_in + 'bottino_seasmean_2D.p', 'rb'))
glomeans, pimean, yeamean, _ = gogo

do_pr = False

if do_pr:
    fig = plt.figure(figsize = (16,9))
    for ru, col in zip(allru, colors):
        plt.scatter(glomeans[(ru, 'tas')][1]-pimean['tas'], glomeans[(ru, 'pr')][1]-pimean['pr'], color = col, s = 10, label = ru)

        if ru == 'ssp585':
            pr = glomeans[(ru, 'pr')][1]
            pr_anom = pr-pimean['pr']
            tas_anom = glomeans[(ru, 'tas')][1]-pimean['tas']

            # coeffs, covmat = np.polyfit(tas_anom, pr_anom, deg = 2, cov = True)
            # x_nu = np.arange(-0.5, 9.6, 0.5)
            # fitted = np.polyval(coeffs, x_nu)
            # plt.plot(x_nu, fitted, color = col)

            #coeffs, covmat = np.polyfit(tas_anom, np.log(pr_anom), deg = 1, cov = True)
            coeffs, covmat = np.polyfit(tas_anom, np.log(pr), deg = 2, cov = True)
            x_nu = np.arange(-0.5, 9.6, 0.5)
            fitted = np.exp(np.polyval(coeffs, x_nu))
            plt.plot(x_nu, fitted-pimean['pr'], color = col)

    plt.grid()
    #plt.legend()
    ctl.custom_legend(fig, colors, allru)
    plt.subplots_adjust(bottom = 0.2)
    plt.xlabel('GTAS anomaly (K)')
    plt.ylabel(r'Prec. anomaly (kg m$^{-2}$ s$^{-1}$)')
    fig.savefig(cart_out + 'pr_gtas_scatter.pdf')


    fig = plt.figure(figsize = (16,9))
    for ru, col in zip(allru, colors):
        if ru in ['pi', 'hist', 'ssp585']: continue
        pr = glomeans[(ru, 'pr')][1]
        pr_anom = pr-pimean['pr']
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
    plt.ylabel(r'Stabilization prec. anomaly (kg m$^{-2}$ s$^{-1}$)')
    fig.savefig(cart_out + 'pr_gtas_nonlin_scatter.pdf')

    fig = plt.figure(figsize = (16,9))
    for ru, col in zip(allru, colors):
        if ru in ['pi', 'hist', 'ssp585']: continue
        pr = glomeans[(ru, 'pr')][1]
        pr_anom = pr-pimean['pr']
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
    plt.ylabel(r'Stabilization prec. anomaly (kg m$^{-2}$ s$^{-1}$)')
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

    for area, anam in zip(areas, anams):
        pino = yeamean[('pi', 'pr')]
        pino_med_pi = ctl.global_mean(ctl.sel_area_xr(pino, area)).mean()

        fig = plt.figure(figsize = (16,9))
        for ru, col in zip(allru, colors):
            tas_anom = glomeans[(ru, 'tas')][1]-pimean['tas']

            pino = yeamean[(ru, 'pr')]
            pino_med = ctl.global_mean(ctl.sel_area_xr(pino, area))

            plt.scatter(tas_anom, pino_med - pino_med_pi, color = col, s = 10, label = ru)

            if ru == 'ssp585':
                coeffs_ru, covmat_ru = np.polyfit(tas_anom, pino_med - pino_med_pi, deg = 1, cov = True)
                x_nu = np.arange(tas_anom.min()-0.1, tas_anom.max()+0.2, 0.1)
                fitted_ru = np.polyval(coeffs_ru, x_nu)
                plt.plot(x_nu, fitted_ru, color = col, ls = '-', lw = 2)
                # rume = ctl.running_mean(pino_med-pino_med_pi, 3s0)
                # teme = ctl.running_mean(tas_anom, 30)
                # plt.plot(teme, rume, color = col)
        plt.grid()
        #plt.legend()
        ctl.custom_legend(fig, colors, allru)
        plt.subplots_adjust(bottom = 0.2)
        plt.xlabel('GTAS anomaly (K)')
        plt.ylabel(r'Prec. anomaly (kg m$^{-2}$ s$^{-1}$)')
        plt.title('Area: {}'.format(anam))

        fig.savefig(cart_out + 'pr_{}_gtas_scatter.pdf'.format(anam))



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


do_fb = False
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

        # ordering for temp
        gino = np.argsort(tas_anom)
        tas_anom = tas_anom[gino]
        pr_anom = pr_anom[gino]

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

        x_nu = np.arange(tas_anom.min(), tas_anom.max(), 0.1)

        pr_ok = []
        tas_ok = []
        for xi in x_nu:
            indi = (tas_anom >= xi) & (tas_anom < xi + 0.1)
            pr_ok.append(np.mean(pr_anom[indi]))
            tas_ok.append(np.mean(tas_anom[indi]))

        pr_ok = np.array(pr_ok)
        tas_ok = np.array(tas_ok)

        ax.scatter(tas_anom, pr_anom, color = col, s = 1)
        ax.scatter(tas_ok, pr_ok, color = col, s = 5, facecolor = 'none', marker = 'o')

        coeffs, covmat = np.polyfit(tas_ok, pr_ok, deg = 1, cov = True)
        fitted = np.polyval(coeffs, x_nu)
        ax.plot(x_nu, fitted, color = col, lw = 2, label = ru)

        coeffs, covmat = np.polyfit(tas_anom, pr_anom, deg = 1, cov = True)
        fitted = np.polyval(coeffs, x_nu)
        ax.plot(x_nu, fitted, color = col, lw = 1, label = ru, ls = '--')

        m, c, err_m, err_c = ctl.linear_regre_witherr(tas_ok, pr_ok)
        stron = r'{}: $\lambda = {:5.3f} \pm {:5.3f} \ W/m^2/K$'.format(ru, m, err_m)
        print(stron)
        ax.text(-0.3, -2-(allru.index(ru)-2)*0.2, stron, fontsize = 14)

    ax.grid()
    ax.legend()
    #ctl.custom_legend(fig, colors[2:-1], allru[2:-1])
    #ax.subplots_adjust(bottom = 0.2)
    ax.set_xlabel('GTAS change (K)')
    ax.set_ylabel(r'Change in net incoming TOA ($W/m^2$)')
    fig.savefig(cart_out + 'feedback_evolution_tasbin.pdf'.format(var))
    ##################################################################

    fig_greg, ax_greg = plt.subplots(figsize = (16,9))

    for ru, col in zip(allru, colors):
        print(ru)
        checkdi = True
        if ru in ['hist', 'pi']: checkdi = False
        ctl.gregplot_on_ax(ax_greg, glomeans[(ru, 'tas')][1]-pimean['tas'], glomeans[(ru, 'net_toa')][1], color = col, label = ru, calc_ERF = False, calc_ECS = False, mean5yr = False, check_diff=checkdi)

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
        while ini <= 500-wind:
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
            if spli == 0:
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


###################################################
###########
carto = cart_out + '../ocean3d/'

oce_mass = 1.381107e+21 # global and vertical sum of masscello*areacello (year 2222)
cp0 = 3989.245 # J/kg/K

fig, ax = plt.subplots(figsize = (16,9))
fig2, ax2 = plt.subplots(figsize = (16,9))
fig3, ax3 = plt.subplots(figsize = (16,9))
for ru, col in zip(allru[3:-1], colors[3:-1]):
    oht_lev = []
    filo = open(carto + 'oht_{}.p'.format(ru), 'rb')
    for i in range(500):
        try:
            gigi = pickle.load(filo)
        except:
            break
        oht_lev.append(gigi[0])

    filo.close()

    oht_lev = xr.concat(oht_lev, dim = 'year')
    # filo = open(carto + 'oht_{}_global.p'.format(ru), 'rb')
    # gigi = pickle.load(filo)
    # filo.close()

    gtas = glomeans[(ru, 'tas')][1]
    yeas = np.arange(500)
    if ru == 'b025':
        gtas = gtas[5:]
        yeas = yeas[5:]

    oht_tot = oht_lev.sum('lev')*cp0
    t_deep = 273.15 + oht_tot/oce_mass/cp0

    #ax.scatter(gtas, (1-t_deep/gtas), s = 3, color = col)
    ax.scatter(gtas-273.15, t_deep-273.15, s = 5, color = col, label = ru)
    grun = ctl.running_mean(gtas, 5, remove_nans = True)
    #larun = ctl.running_mean(1-t_deep/gtas, 5, remove_nans = True)
    larun = ctl.running_mean(t_deep, 5, remove_nans = True)
    # ax.plot(grun, larun, color = col, label = ru)
    # x_nu = np.arange(gtas.min(), gtas.max(), 0.1)
    # coeffs, covmat = np.polyfit(grun, larun, deg = 2, cov = True)
    # fitted = np.polyval(coeffs, x_nu)
    # ax.plot(x_nu, fitted, color = col, label = ru, lw = 2)

    oht1 = oht_lev.sel(lev = slice(0., 700.)).sum('lev')
    oht2 = oht_lev.sel(lev = slice(700., 2000.)).sum('lev')
    oht3 = oht_lev.sel(lev = slice(2000., 6000.)).sum('lev')

    ax2.plot(oht1, color = col, ls = '-', label = ru, lw = 2)
    ax2.plot(oht2, color = col, ls = '--', lw = 2)
    ax2.plot(oht3, color = col, ls = '-.', lw = 2)

    grun = ctl.running_mean(gtas, 10, remove_nans = True)
    oht1l = ctl.running_mean(oht1, 10, remove_nans = True)
    oht2l = ctl.running_mean(oht2, 10, remove_nans = True)
    oht3l = ctl.running_mean(oht3, 10, remove_nans = True)

    ax3.scatter(grun, oht1l, color = col, label = ru, s = 20)
    ax3.scatter(grun, oht2l, edgecolor = col, s = 20, marker = 'd', facecolor = 'none')
    ax3.scatter(grun, oht3l, color = col, s = 20, marker = '*')

ax.legend()
ax.grid()
#ax.set_ylabel(r'$(1-T_d/T_s)$')
ax.set_ylabel('T bulk ocean (K)')
ax.set_xlabel('GTAS (K)')

ax2.legend()
ax2.grid()
ax2.set_ylabel('OHT (J)')
ax2.set_xlabel('Years after stabilization')

ax3.legend()
ax3.grid()
ax3.set_ylabel('OHT (J)')
ax3.set_xlabel('GTAS (K)')

fig.savefig(carto + 'lambda_factor.pdf')
fig2.savefig(carto + 'oht_deep_time.pdf')
fig3.savefig(carto + 'oht_deep_gtas.pdf')

#############################################################

fig, axs = plt.subplots(1, 3, figsize = (16,6))

oht_all = dict()

oht1ref = None
for ru, col, ax in zip(allru[3:-1], colors[3:-1], axs.flatten()):
    oht_lev = []
    filo = open(carto + 'oht_{}.p'.format(ru), 'rb')
    for i in range(500):
        try:
            gigi = pickle.load(filo)
        except:
            break
        oht_lev.append(gigi[0])

    filo.close()

    oht_lev = xr.concat(oht_lev, dim = 'year')

    gtas = glomeans[(ru, 'tas')][1]
    yeas = np.arange(500)
    if ru == 'b025':
        gtas = gtas[5:]
        yeas = yeas[5:]

    oht_tot = oht_lev.sum('lev')
    oht1 = oht_lev.sel(lev = slice(0., 700.)).sum('lev')
    oht2 = oht_lev.sel(lev = slice(700., 2000.)).sum('lev')
    oht3 = oht_lev.sel(lev = slice(2000., 6000.)).sum('lev')

    oht_all[(ru, 700)] = oht1
    oht_all[(ru, 2000)] = oht2
    oht_all[(ru, 'deep')] = oht3

    if oht1ref is None:
        oht1ref = oht1[:10].mean()
        oht2ref = oht2[:10].mean()
        oht3ref = oht3[:10].mean()
        ohttotref = oht_tot[:10].mean()

    # ax.fill_between(np.zeros(len(oht3)), (oht3-oht3ref)/(oht_tot-ohttotref), color = 'steelblue', alpha = 0.5)
    # ax.fill_between((oht3-oht3ref)/(oht_tot-ohttotref), ((oht3-oht3ref)+(oht2-oht2ref))/(oht_tot-ohttotref), color = 'forestgreen', alpha = 0.5)
    # ax.fill_between(((oht3-oht3ref)+(oht2-oht2ref))/(oht_tot-ohttotref), np.ones(len(oht1)), color = 'gold', alpha = 0.5)
    #
    # ax.plot((oht3-oht3ref)/(oht_tot-ohttotref), color = 'steelblue', ls = '-', lw = 2)
    # ax.plot(((oht3-oht3ref)+(oht2-oht2ref))/(oht_tot-ohttotref), color = 'forestgreen', ls = '-', lw = 2)
    # # ax.plot(oht3-oht3ref, color = col, ls = '-.', lw = 2)

    ax.fill_between(np.arange(len(oht3)), np.zeros(len(oht3)), (oht3-oht3ref), color = 'steelblue', alpha = 0.5)
    ax.fill_between(np.arange(len(oht3)), (oht3-oht3ref), (oht3-oht3ref)+(oht2-oht2ref), color = 'forestgreen', alpha = 0.5)
    ax.fill_between(np.arange(len(oht3)), (oht3-oht3ref)+(oht2-oht2ref), oht_tot-ohttotref, color = 'gold', alpha = 0.5)

    ax.plot((oht3-oht3ref), color = 'steelblue', ls = '-', lw = 2)
    ax.plot((oht3-oht3ref)+(oht2-oht2ref), color = 'forestgreen', ls = '-', lw = 2)
    ax.plot(oht_tot-ohttotref, color = 'gold', ls = '-', lw = 2)

    ax.set_title(ru)
    ax.grid()

axs[0].set_ylabel('OHT (J)')
axs[1].set_xlabel('Years after stabilization')

ctl.adjust_ax_scale(axs)
ctl.custom_legend(fig, ['steelblue', 'forestgreen', 'gold'], ['0-700m', '700-2000m', '> 2000m'], ncol = 3)
plt.subplots_adjust(bottom = 0.25)
fig.savefig(carto + 'oht_deep_all_evol.pdf')

pickle.dump(oht_all, open(carto + 'oht_ts_deep.p', 'wb'))
#############################################################
### maps of OHT trends
lats = np.linspace(-89.5, 89.5, 180)
lons = np.linspace(0, 359, 360)

refoht = dict() #### THIS IS now assigned to b025, but should be pi
oht_patt = dict()

areacello = xr.load_dataset(carto + 'areacello.nc')['areacello']
regr_acello = ctl.regrid_dataset(areacello, lats, lons)

for ru, col in zip(allru[3:-1], colors[3:-1]):
    print(ru)
    filo = open(carto + 'oht_{}.p'.format(ru), 'rb')

    oht700 = []
    oht2000 = []
    ohtdeep = []
    for i in range(500):
        try:
            oht_lev_i, oht700_i, oht2000_i, ohtdeep_i = pickle.load(filo)
        except:
            break
        oht700.append(oht700_i)
        oht2000.append(oht2000_i)
        ohtdeep.append(ohtdeep_i)

    filo.close()

    oht700 = np.stack(oht700)
    oht2000 = np.stack(oht2000)
    ohtdeep = np.stack(ohtdeep)

    ### Ri-divide by areacello! I want density of OHC here
    oht700 = oht700/regr_acello.values
    oht2000 = oht2000/regr_acello.values
    ohtdeep = ohtdeep/regr_acello.values

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


        if lab not in refoht:
            refoht[lab] = var[:20].mean(axis = 0) # taking b025 as reference

        oht_patt[(ru, lab, 'change')] = var[-20:].mean(axis = 0) - refoht[lab]

pickle.dump(oht_patt, open(carto + 'temp_patt_deep.p', 'wb'))



plpa = []
subt = []
hatch = []
for lev, tit in zip([700, 2000, 'deep'], ['0-700m', '700-2000m', '> 2000 m']):
    for ru, col in zip(allru[3:-1], colors[3:-1]):
        plpa.append(oht_patt[(ru, lev)])
        hatch.append(oht_patt[(ru, lev, 'pval')] < 0.05)
        subt.append(ru + ': ' + tit)

[fig] = ctl.plot_multimap_contour(plpa, lats, lons, visualization = 'Robinson', central_lat_lon = (0., -120.), filename = carto + 'temp_patt_deep.pdf', subtitles = subt, plot_anomalies = False, cmap = 'viridis', figsize = (16,12), fix_subplots_shape = (3,3), cb_label = r'Conservative temperature trend ($K yr^{-1}$)')#, add_hatching = hatch, hatch_styles = ['///', '', ''])

for ax in fig.axes[:-1]:
    ax.set_facecolor('gainsboro')
#fig.savefig(carto + 'oht_patt_deep.pdf')
fig.savefig(carto + 'temp_patt_deep.pdf')


plpa = []
subt = []
hatch = []
for lev, tit in zip([700, 2000, 'deep'], ['0-700m', '700-2000m', '> 2000 m']):
    for ru, col in zip(allru[3:-1], colors[3:-1]):
        plpa.append(oht_patt[(ru, lev, 'rel')])
        hatch.append(oht_patt[(ru, lev, 'rel_pval')] < 0.05)
        subt.append(ru + ': ' + tit)

#divnorm = colors.TwoSlopeNorm(vmin=-1., vcenter=1., vmax=3.5)
[fig] = ctl.plot_multimap_contour(plpa, lats, lons, visualization = 'Robinson', central_lat_lon = (0., -120.), filename = carto + 'temp_patt_deep_rel.pdf', subtitles = subt, plot_anomalies = False, cmap = ctl.heatmap(), cbar_range = (-1, 3), figsize = (16,12), fix_subplots_shape = (3,3), cb_label = 'OHC trend pattern', n_color_levels = 37)#, add_hatching = hatch, hatch_styles = ['///', '', ''])
for ax in fig.axes[:-1]:
    ax.set_facecolor('gainsboro')
#fig.savefig(carto + 'oht_patt_deep_rel.pdf')
fig.savefig(carto + 'temp_patt_deep_rel.pdf')

plpa = []
subt = []
hatch = []
for lev, tit in zip([700, 2000, 'deep'], ['0-700m', '700-2000m', '> 2000 m']):
    for ru, col in zip(allru[3:-1], colors[3:-1]):
        plpa.append(oht_patt[(ru, lev, 'change')])
        subt.append(ru + ': ' + tit)

[fig] = ctl.plot_multimap_contour(plpa, lats, lons, visualization = 'Robinson', central_lat_lon = (0., -120.), filename = carto + 'temp_patt_deep_change.pdf', subtitles = subt, plot_anomalies = False, cmap = 'viridis', figsize = (16,12), fix_subplots_shape = (3,3), cb_label = r'Conservative temperature change ($K$)')

for ax in fig.axes[:-1]:
    ax.set_facecolor('gainsboro')
fig.savefig(carto + 'temp_patt_deep_change.pdf')
