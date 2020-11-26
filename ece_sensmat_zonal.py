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
import pandas as pd
import glob

import tunlib as tl

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 18

################################################################################

cart_in = '/data-hobbes/fabiano/TunECS/AMIP_exps/'
cart_out = '/home/fabiano/Research/lavori/TunECS/tuning/experiments/analysis/'

mname = '{:2s}{:1s}{:1s}'
fil = cart_in + '{:4s}/post/mon/Post_{:4d}/{:4s}_{:4d}_{}.nc'

testparams = ['ENTRORG', 'RPRCON', 'DETRPEN', 'RMFDEPS', 'RVICE', 'RSNOWLIN2', 'RCLDIFF', 'RLCRIT_UPHYS']
nums = np.arange(8)
letts = 'a b c d e f g h'.split()

uff_params = dict()
uff_params['RPRCON'] = 1.34E-3
uff_params['RVICE'] = 0.137
uff_params['RLCRITSNOW'] = 4.0E-5
uff_params['RSNOWLIN2'] = 0.035
uff_params['ENTRORG'] = 1.70E-4
uff_params['DETRPEN'] = 0.75E-4
uff_params['ENTRDD'] = 3.0E-4
uff_params['RMFDEPS'] = 0.3
uff_params['RCLDIFF'] = 3.E-6
uff_params['RCLDIFFC'] = 5.0
uff_params['RLCRIT_UPHYS'] = 0.875e-5

diseqs = dict()
diseqs['m'] = -2
diseqs['n'] = -1
diseqs['p'] = 1
diseqs['q'] = 2
diseqs['l'] = -0.5
diseqs['r'] = 0.5

valchange = dict()
valchange['ENTRORG'] = np.array([1.05, 1.35, 2.05, 2.35, 1.53, 1.87])*1e-4
valchange['RPRCON'] = np.array([2.45, 1.9, 0.8, 0.25, 1.62, 1.07])*1e-3
valchange['DETRPEN'] = np.array([0.1, 0.25, 1.25, 1.75, 0.5, 1.0])*1e-4
valchange['RMFDEPS'] = np.array([0.02, 0.16, 0.44, 0.58, 0.23, 0.37])
valchange['RVICE'] = np.array([0.06, 0.098, 0.176, 0.214, 0.118, 0.157])
valchange['RSNOWLIN2'] = np.array([0.079, 0.057, 0.013, 0.001, 0.046, 0.024])
valchange['RCLDIFF'] = np.array([5, 4, 2, 1, 3.5, 2.5])*1e-6
valchange['RLCRIT_UPHYS'] = np.array([1.02, 0.95, 0.8, 0.73, 0.91, 0.84])*1e-5

allforc = ['pi', 'c4']
allchan = 'm n p q l r'.split()

def val_ok(param, change):
    iok = allchan.index(change)
    return valchange[param][iok]

forcsym = dict()
forcsym['pi'] = 'o'
forcsym['c4'] = 'x'

forccol = dict()
forccol['pi'] = 'lightseagreen'
forccol['c4'] = 'indianred'

changecol = dict()
changecol['m'] = 'darkblue'
changecol['n'] = 'steelblue'
changecol['p'] = 'darkorange'
changecol['q'] = 'indianred'
changecol['l'] = 'pink'
changecol['r'] = 'violet'

# cosa voglio fare: per ogni membro, faccio la media del toa_long e del toa_short annuale. e metto a punto una funzione che dato il cambio dei parametri, mi da l'effetto atteso sulla toa_net. Quindi forse voglio anche qui le derivate. yes, derivate normalizzate per ogni parametro.
# prima leggo, poi faccio la media zonale?, meglio in bande? SI. [-90, -65, -40, -20, 20, 40, 65, 90]
# Salvo la media nelle bande per ogni exp
# poi faccio le derivate, banda per banda, e le plotto, normalizzate. un plot per pi e uno per c4, per ogni var.

lats = [-90, -65, -40, -20, 20, 40, 65, 90]
bands = [(la1, la2) for la1, la2 in zip(lats[:-1], lats[1:])]
lacen = np.array([np.mean(laol) for laol in bands])

allvars = ['ttr', 'tsr', 'str', 'ssr', 'tcc', 'cp', 'lsp']

# tsr: This parameter is the incoming solar radiation (also known as shortwave radiation) minus the outgoing solar radiation at the top of the atmosphere. It is the amount of radiation passing through a horizontal plane. The incoming solar radiation is the amount received from the Sun. The outgoing solar radiation is the amount reflected and scattered by the Earth's atmosphere and surface.
# ttr: The thermal (also known as terrestrial or longwave) radiation emitted to space at the top of the atmosphere is commonly known as the Outgoing Longwave Radiation (OLR). The top net thermal radiation (this parameter) is equal to the negative of OLR.
# tsr + ttr = tnr. Positive -> incoming! (downward)
# ssr: This parameter is the amount of solar radiation (also known as shortwave radiation) that reaches a horizontal plane at the surface of the Earth (both direct and diffuse) minus the amount reflected by the Earth's surface (which is governed by the albedo).
# str: This parameter is the difference between downward and upward thermal radiation at the surface of the Earth. It the amount passing through a horizontal plane.
# ssr + str = snr. Positive -> downward

# energia assorbita da atm: tnr-snr

resdic = dict()
resdic_err = dict()
for varnam in allvars:
    print(varnam)
    for forc in ['pi', 'c4']:
        for nu, let, param in zip(nums, letts, testparams):
            for iic, change in enumerate(['m', 'n', 'p', 'q', 'l', 'r']):
                mok = mname.format(forc, change, let)
                print(mok)

                if forc == 'pi' and let == 'h' and change in ['p', 'n']:
                    listafil = [fil.format(mok, ye, mok, ye, varnam) for ye in range(1851, 1860)] # ten years for these, instead of 5
                else:
                    listafil = [fil.format(mok, ye, mok, ye, varnam) for ye in range(1851, 1855)] # skipping the first year

                if tl.check_file(listafil[0]):
                    var, coords, aux_info = ctl.read_ensemble_iris(listafil, netcdf4_read = True)
                    varm, varstd = ctl.zonal_seas_climatology(var, coords['dates'], 'year')
                    varstd = varstd/np.sqrt(len(listafil)-1)

                    for band in bands:
                        resdic[(forc, change, let, varnam, band)] = ctl.band_mean_from_zonal(varm, coords['lat'], band[0], band[1])
                        resdic_err[(forc, change, let, varnam, band)] = ctl.band_mean_from_zonal(varstd, coords['lat'], band[0], band[1])

                    varm, varstd = ctl.global_seas_climatology(var, coords['lat'], coords['dates'], 'year')
                    varstd = varstd/np.sqrt(len(listafil)-1)
                    resdic[(forc, change, let, varnam, 'glob')] = varm
                    resdic_err[(forc, change, let, varnam, 'glob')] = varstd

        # Loading the control
        if forc == 'pi':
            mok = 'tpa1'
        elif forc == 'c4':
            mok = 't4a1'
        else:
            mok = 'bau'
        print(mok)

        listafil = [fil.format(mok, ye, mok, ye, varnam) for ye in range(1851, 1855)] # skipping the first year
        if tl.check_file(listafil[0]):
            var, coords, aux_info = ctl.read_ensemble_iris(listafil, netcdf4_read = True)
            varm, varstd = ctl.zonal_seas_climatology(var, coords['dates'], 'year')
            varstd = varstd/np.sqrt(len(listafil)-1)

            change = 0
            let = 0
            for band in bands:
                resdic[(forc, change, let, varnam, band)] = ctl.band_mean_from_zonal(varm, coords['lat'], band[0], band[1])
                resdic_err[(forc, change, let, varnam, band)] = ctl.band_mean_from_zonal(varstd, coords['lat'], band[0], band[1])

            varm, varstd = ctl.global_seas_climatology(var, coords['lat'], coords['dates'], 'year')
            varstd = varstd/np.sqrt(len(listafil)-1)
            resdic[(forc, change, let, varnam, 'glob')] = varm
            resdic_err[(forc, change, let, varnam, 'glob')] = varstd


for var in ['toa_net', 'srf_net']:
    firl = var[0]
    for forc in ['pi', 'c4']:
        for band in bands+['glob']:
            for nu, let, param in zip(nums, letts, testparams):
                for iic, change in enumerate(['m', 'n', 'p', 'q', 'l', 'r']):
                    if (forc, change, let, firl+'tr', band) in resdic.keys():
                        resdic[(forc, change, let, var, band)] = resdic[(forc, change, let, firl+'tr', band)]+resdic[(forc, change, let, firl+'sr', band)]
                        resdic_err[(forc, change, let, var, band)] = np.mean(resdic_err[(forc, change, let, firl+'tr', band)]+resdic_err[(forc, change, let, firl+'sr', band)])

            resdic[(forc, 0, 0, var, band)] = resdic[(forc, 0, 0, firl+'tr', band)]+resdic[(forc, 0, 0, firl+'sr', band)]
            resdic_err[(forc, 0, 0, var, band)] = np.mean(resdic_err[(forc, 0, 0, firl+'tr', band)]+resdic_err[(forc, 0, 0, firl+'sr', band)])

allvars.append('toa_net')
allvars.append('srf_net')

allforc = ['pi', 'c4']

derdic = dict()
derdic_err = dict()
## Derivata con parametro normalizzato
for var in allvars:
    figs = []
    axes = []
    for band in bands + ['glob']:
        fig, ax = plt.subplots(figsize=(16,12))

        ctrl = dict()
        ctrl['pi'] = resdic[('pi', 0, 0, var, band)]
        ctrl['c4'] = resdic[('c4', 0, 0, var, band)]

        for forc, shift in zip(allforc, [-0.05, 0.05]):
            ders = []
            err_ders = []
            for nu, let, param in zip(nums, letts, testparams):
                cose = []
                errs = []
                xval = []
                if (forc, 'l', let, var, band) not in resdic:
                    for ii, change in zip([1, -1, 2], ['n', 0, 'p']):
                        if change == 0:
                            cose.append(ctrl[forc])
                        else:
                            cose.append(resdic[(forc, change, let, var, band)])
                            errs.append(resdic_err[(forc, change, let, var, band)])

                        if ii >= 0:
                            xval.append(valchange[param][ii])
                        else:
                            xval.append(uff_params[param])

                    deriv = np.gradient(np.array(cose), np.array(xval))

                    print(param, deriv)
                    derdic[(forc, param, var, band)] = deriv[1]
                    if tl.check_increasing(xval):
                        derdic[(forc, param, var, band, 'left')] = deriv[0]
                        derdic[(forc, param, var, band, 'right')] = deriv[-1]
                    elif tl.check_decreasing(xval):
                        derdic[(forc, param, var, band, 'left')] = deriv[-1]
                        derdic[(forc, param, var, band, 'right')] = deriv[0]
                    else:
                        print(xval)
                        raise ValueError('problema tenemos')

                    ders.append(uff_params[param]*deriv[1])

                    errs.insert(1, 0.)
                    errs[2] = -errs[2]
                    deriv_err = np.gradient(np.array(errs), np.array(xval))
                    err_ders.append(uff_params[param]*np.abs(deriv_err[1]))
                    derdic_err[(forc, param, var, band)] = np.abs(deriv_err[1])
                else:
                    for change in ['n', 'l', 0, 'r', 'p']:
                        if change == 0:
                            cose.append(ctrl[forc])
                            errs.append(0)
                            xval.append(uff_params[param])
                        else:
                            cose.append(resdic[(forc, change, let, var, band)])
                            errs.append(resdic_err[(forc, change, let, var, band)])
                            xval.append(val_ok(param, change))

                    deriv = np.gradient(np.array(cose), np.array(xval))
                    derdic[(forc, param, var, band)] = deriv[2]
                    print(param, deriv)
                    if tl.check_increasing(xval):
                        derdic[(forc, param, var, band, 'left')] = deriv[0]
                        derdic[(forc, param, var, band, 'right')] = deriv[-1]
                    elif tl.check_decreasing(xval):
                        derdic[(forc, param, var, band, 'left')] = deriv[-1]
                        derdic[(forc, param, var, band, 'right')] = deriv[0]
                    else:
                        print(xval)
                        raise ValueError('problema tenemos')

                    ders.append(uff_params[param]*deriv[2])

                    errs[3] = -errs[3]
                    errs[4] = -errs[4]
                    deriv_err = np.gradient(np.array(errs), np.array(xval))
                    err_ders.append(uff_params[param]*np.abs(deriv_err[2]))
                    derdic_err[(forc, param, var, band)] = np.abs(deriv_err[2])

            ax.errorbar(nums+shift, ders, yerr = err_ders, fmt = 'none', color = forccol[forc], capsize = 2, elinewidth = 1)
            ax.scatter(nums+shift, ders, color = forccol[forc], marker = forcsym[forc], s = 100, label = forc)

        ax.legend()
        ax.set_xticks(nums)
        ax.set_xticklabels(testparams, size = 'large', rotation = 60)
        ax.grid()
        ax.axhline(0., color = 'black')
        ax.set_ylabel('uff_param * derivative of '+ var + ' (W/m2)')
        axes.append(ax)

        fig.suptitle('Derivative of {} in lat band ({}, {})'.format(var, band[0], band[1]))
        figs.append(fig)

        #fig.savefig(cart_out + var+'_scattplot_{}.pdf'.format('deriv'))

    ctl.adjust_ax_scale(axes)
    ctl.plot_pdfpages(cart_out + '{}_sensmat_zonal.pdf'.format(var), figs)


### Adding simple linear deriv
linder = dict()
linder_err = dict()
for var in allvars:
    for band in bands + ['glob']:
        ctrl = dict()
        ctrl['pi'] = resdic[('pi', 0, 0, var, band)]
        ctrl['c4'] = resdic[('c4', 0, 0, var, band)]

        for forc, shift in zip(allforc, [-0.05, 0.05]):
            for nu, let, param in zip(nums, letts, testparams):
                for iic, change in zip([1, 2], ['n', 'p']):
                    if valchange[param][iic] < uff_params[param]:
                        diff = ctrl[forc] - resdic[(forc, change, let, var, band)]
                        xdi = uff_params[param] - valchange[param][iic]
                        linder[(forc, param, var, band, 'left')] = diff/xdi
                        linder_err[(forc, param, var, band, 'left')] = resdic_err[(forc, change, let, var, band)]/xdi
                    else:
                        diff = resdic[(forc, change, let, var, band)] - ctrl[forc]
                        xdi = valchange[param][iic] - uff_params[param]
                        linder[(forc, param, var, band, 'right')] = diff/xdi
                        linder_err[(forc, param, var, band, 'right')] = resdic_err[(forc, change, let, var, band)]/xdi


with open(cart_out + 'der_sensmat_zonal.p', 'wb') as filox:
    pickle.dump([resdic, resdic_err, derdic, derdic_err, linder, linder_err], filox)


## Derivata con parametro normalizzato
for var in allvars:
    figs = []
    axes = []
    for nu, let, param in zip(nums, letts, testparams):
        fig, ax = plt.subplots(figsize=(16,12))

        for forc, shift in zip(allforc, [-0.05, 0.05]):
            ders = []
            ders_left = []
            ders_right = []
            err_ders = []
            for band in bands:
                ders.append(uff_params[param]*derdic[(forc, param, var, band)])
                ders_left.append(uff_params[param]*derdic[(forc, param, var, band, 'left')])
                ders_right.append(uff_params[param]*derdic[(forc, param, var, band, 'right')])
                err_ders.append(uff_params[param]*derdic_err[(forc, param, var, band)])

            ders = np.array(ders)
            err_ders = np.array(err_ders)
            ax.fill_between(lacen, ders-err_ders, ders+err_ders, color = forccol[forc], alpha = 0.3)
            ax.plot(lacen, ders, color = forccol[forc], label = forc)
            ax.scatter(lacen, ders, color = forccol[forc], marker = forcsym[forc], s = 100)
            ax.scatter(lacen, ders_left, color = forccol[forc], marker = '<', s = 70)
            ax.scatter(lacen, ders_right, color = forccol[forc], marker = '>', s = 70)
            #ax.errorbar(lacen, ders, yerr = err_ders, fmt = 'none', color = forccol[forc], capsize = 2, elinewidth = 1)

        ax.legend()
        ax.grid()
        ax.axhline(0., color = 'black')
        ax.set_ylabel('uff_param * derivative of '+ var + ' (W/m2)')
        ax.set_xlabel('Latitude')
        axes.append(ax)

        fig.suptitle('Zonal derivative of {} wrt {}'.format(var, param))
        figs.append(fig)
        #fig.savefig(cart_out + var+'_scattplot_{}.pdf'.format('deriv'))

    ctl.adjust_ax_scale(axes)
    ctl.plot_pdfpages(cart_out + '{}_sensmat_zonal_wparam.pdf'.format(var), figs)


############# Plot toa_net diffs for each param
forcsty = dict()
forcsty['pi'] = '-'
forcsty['c4'] = '--'

for var in allvars:
    figs = []
    axes = []
    for nu, let, param in zip(nums, letts, testparams):

        for forc, shift in zip(allforc, [-0.05, 0.05]):
            fig, ax = plt.subplots(figsize=(16,12))
            ctrl = np.array([resdic[(forc, 0, 0, var, band)] for band in bands])
            print('ctrl', ctrl)

            for change in ['n', 'l', 'r', 'p']:
                if (forc, change, let, var, bands[0]) not in resdic:
                    continue
                vals = np.array([resdic[(forc, change, let, var, band)] for band in bands])
                print(vals)
                err_vals = np.array([resdic_err[(forc, change, let, var, band)] for band in bands])
                vals = vals-ctrl

                ax.fill_between(lacen, vals-err_vals, vals+err_vals, color = changecol[change], alpha = 0.3)
                ax.plot(lacen, vals, color = changecol[change], label = change, linestyle = forcsty[forc])
                ax.scatter(lacen, vals, color = changecol[change], marker = forcsym[forc], s = 100)

            ax.legend()
            ax.grid()
            ax.axhline(0., color = 'black')
            ax.set_ylabel('change of '+ var + ' (W/m2)')
            ax.set_xlabel('Latitude')
            axes.append(ax)

            fig.suptitle('{} run: change of {} wrt {}'.format(forc, var, param))
            figs.append(fig)
            #fig.savefig(cart_out + var+'_scattplot_{}.pdf'.format('deriv'))

    ctl.adjust_ax_scale(axes)
    ctl.plot_pdfpages(cart_out + '{}_changemat_zonal.pdf'.format(var), figs)


for var in allvars:
    figs = []
    axes = []
    for nu, let, param in zip(nums, letts, testparams):

        for forc, shift in zip(allforc, [-0.05, 0.05]):
            fig, ax = plt.subplots(figsize=(16,12))
            ctrl = np.array([resdic[(forc, 0, 0, var, band)] for band in bands])
            print('ctrl', ctrl)

            for iic, change in enumerate(['m', 'n', 'p', 'q', 'l', 'r']):
                if (forc, change, let, var, bands[0]) not in resdic:
                    continue
                vals = np.array([resdic[(forc, change, let, var, band)] for band in bands])
                print(vals)
                err_vals = np.array([resdic_err[(forc, change, let, var, band)] for band in bands])
                vals = vals-ctrl

                ax.fill_between(lacen, vals-err_vals, vals+err_vals, color = changecol[change], alpha = 0.3)
                ax.plot(lacen, vals, color = changecol[change], label = change, linestyle = forcsty[forc])

                cglob, czon = tl.calc_change_var(forc, param, var, valchange[param][iic], method = 'deriv')
                ax.scatter(lacen, czon, color = changecol[change], marker = '*', s = 70)
                cglob, czon = tl.calc_change_var(forc, param, var, valchange[param][iic], method = 'deriv_edge')
                ax.scatter(lacen, czon, color = changecol[change], marker = '<', s = 70)

            ax.legend()
            ax.grid()
            ax.axhline(0., color = 'black')
            ax.set_ylabel('change of '+ var + ' (W/m2)')
            ax.set_xlabel('Latitude')
            axes.append(ax)

            fig.suptitle('{} run: change of {} wrt {}'.format(forc, var, param))
            figs.append(fig)
            #fig.savefig(cart_out + var+'_scattplot_{}.pdf'.format('deriv'))

    ctl.adjust_ax_scale(axes)
    ctl.plot_pdfpages(cart_out + '{}_changemat_zonal_wcheck.pdf'.format(var), figs)
