#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os

from matplotlib import pyplot as plt
import netCDF4 as nc
import pandas as pd
from numpy import linalg as LA
from scipy import stats
from scipy import interpolate as itrp
import itertools as itt

from datetime import datetime
import pickle

import climtools_lib as ctl
import climdiags as cd

from copy import deepcopy as cp
###############################################

# ###############################################################################
# ### constants
cp = 1005.0 # specific enthalpy dry air - J kg-1 K-1
cpw = 1840.0 # specific enthalpy water vapor - J kg-1 K-1
# per la moist air sarebbe cpm = cp + q*cpw, ma conta al massimo per un 3 %
L = 2501000.0 # J kg-1
Lsub = 2835000.0 # J kg-1
g = 9.81 # m s-2
Rearth = 6371.0e3 # mean radius
KtoC = 273.15
# ###############################################################################

# Calculates radiation balances
cart_in = '/data-hobbes/fabiano/SPHINX/AMIP/radiation/'
cart_out = '/home/fabiano/Research/lavori/SPHINX_for_lisboa/radiation_balance_AMIP/'
if not os.path.exists(cart_out): os.mkdir(cart_out)

# read masks
filo = '/data-hobbes/fabiano/SPHINX/masks.nc'
fh = nc.Dataset(filo)
land_mask = fh.variables['RnfO.msk'][:]
ocean_mask = fh.variables['RnfA.msk'][:]

varnames = ['rlds', 'rlut', 'rlus', 'rsdt', 'rsut', 'rsds', 'rsus', 'hfls', 'hfss']
namefi = '{}_mon_{}_{}.nc'

years = dict()
years['a'] = np.arange(1979,2009)
years['f'] = np.arange(2038,2069)

nuvars = ['toa_balance', 'toa_balance_land', 'toa_balance_ocean', 'surf_balance', 'surf_balance_land', 'surf_balance_ocean']

ensmems = dict()
ensmems['a'] = ['lab{}'.format(i) for i in range(10)] + ['las{}'.format(i) for i in range(10)]
ensmems['f'] = ['lfb{}'.format(i) for i in range(10)] + ['lfs{}'.format(i) for i in range(10)]

# radclim = dict()
# for lett in ['a', 'f']:
#     for exp in ensmems[lett]:
#         for var in nuvars:
#             radclim[('zonal', exp, var)] = []
#             radclim[('global', exp, var)] = []
#         for year in years[lett]:
#             vardict = dict()
#             for varna in varnames:
#                 varniuu, lat, lon, dates, time_units, var_units = ctl.read3Dncfield(cart_in+namefi.format(exp,year,varna))
#                 vardict[varna] = np.mean(varniuu, axis = 0)
#                 if year == years[lett][0]:
#                     radclim[('global', exp, varna)] = [ctl.global_mean(vardict[varna], lat)]
#                 else:
#                     radclim[('global', exp, varna)].append(ctl.global_mean(vardict[varna], lat))
#
#             toa_bal = vardict['rsdt']-vardict['rsut']-vardict['rlut']
#             surf_bal = vardict['rsds']+vardict['rlds']-vardict['rsus']-vardict['rlus']-vardict['hfls']-vardict['hfss']
#             radclim[('zonal', exp, 'toa_balance')].append(ctl.zonal_mean(toa_bal))
#             radclim[('zonal', exp, 'surf_balance')].append(ctl.zonal_mean(surf_bal))
#
#             radclim[('global', exp, 'toa_balance')].append(ctl.global_mean(toa_bal, lat))
#             radclim[('global', exp, 'surf_balance')].append(ctl.global_mean(surf_bal, lat))
#
#             radclim[('global', exp, 'toa_balance_land')].append(ctl.global_mean(toa_bal, lat, mask = land_mask))
#             radclim[('global', exp, 'surf_balance_land')].append(ctl.global_mean(surf_bal, lat, mask = land_mask))
#             radclim[('global', exp, 'toa_balance_ocean')].append(ctl.global_mean(toa_bal, lat, mask = ocean_mask))
#             radclim[('global', exp, 'surf_balance_ocean')].append(ctl.global_mean(surf_bal, lat, mask = ocean_mask))
#
#         for varna in varnames+nuvars:
#             radclim[('global', exp, varna)] = np.array(radclim[('global', exp, varna)])
#
#         radclim[('zonal', exp, 'toa_balance')] = np.stack(radclim[('zonal', exp, 'toa_balance')])
#         radclim[('zonal', exp, 'surf_balance')] = np.stack(radclim[('zonal', exp, 'surf_balance')])
#
#
# factor = 2*np.pi*Rearth**2 * np.cos(np.deg2rad(lat))
#
# for lett in ['a', 'f']:
#     for exp in ensmems[lett]:
#         radclim[('zonal', exp, 'heat_flux')] = []
#         for toa_b, surf_b in zip(radclim[('zonal', exp, 'toa_balance')], radclim[('zonal', exp, 'surf_balance')]):
#             heatflu = []
#             integrand = (toa_b - surf_b)*factor
#             for ind in range(1,len(lat)+1):
#                 heatflu.append(np.trapz(integrand[:ind], x = lat[:ind]))
#             radclim[('zonal', exp, 'heat_flux')].append(np.array(heatflu))
#
#         radclim[('zonal', exp, 'heat_flux')] = np.stack(radclim[('zonal', exp, 'heat_flux')])
#
#
#
# ensmemsall = ensmems['a']+ensmems['f']
# nunams = ['base_hist', 'stoc_hist', 'base_fut', 'stoc_fut']
# expnams = ['lab', 'las', 'lfb', 'lfs']
#
# for varna in varnames+nuvars:
#     for nunam, expnam in zip(nunams, expnams):
#         radclim[('global', nunam, varna)] = np.mean([radclim[('global', exp, varna)] for exp in ensmemsall if expnam in exp], axis = 0)
#
# for varna in ['heat_flux', 'toa_balance', 'surf_balance']:
#     for nunam, expnam in zip(nunams, expnams):
#         radclim[('zonal', nunam, varna)] = np.mean([radclim[('zonal', exp, varna)] for exp in ensmemsall if expnam in exp], axis = 0)
#
# pickle.dump(radclim, open(cart_out+'radclim_yearly_AMIP.p', 'wb'))

radclim = pickle.load(open(cart_out+'radclim_yearly_AMIP.p'))
varniuu, lat, lon, dates, time_units, var_units = ctl.read3Dncfield(cart_in+namefi.format('lab0',1988,'rsut'))
del varniuu

carttemp = '/data-hobbes/fabiano/SPHINX/AMIP/tas_mon/'
# globalme = dict()
# varna = 'tas'
#
# tempnames = dict()
# tempnames['a'] = '{}_mon_1979-2008_tas.nc'
# tempnames['f'] = '{}_mon_2038-2068_tas.nc'
# for lett in ['a','f']:
#     for ens in ensmems[lett]:
#         filena = carttemp+tempnames[lett].format(ens)
#         var, lata, lona, dates, time_units, var_units = ctl.read3Dncfield(filena)
#         varye, _ = ctl.yearly_average(var, dates)
#         global_mean = ctl.global_mean(varye, lata)
#         global_mean_land = ctl.global_mean(varye, lata, mask = land_mask)
#         global_mean_ocean = ctl.global_mean(varye, lata, mask = ocean_mask)
#         globalme[(ens, varna)] = global_mean
#         globalme[(ens, varna+'_land')] = global_mean_land
#         globalme[(ens, varna+'_ocean')] = global_mean_ocean
#
#         zonal_mean = ctl.zonal_mean(varye)
#         globalme[(ens, varna+'_trop')] = ctl.band_mean_from_zonal(zonal_mean, lata, -20., 20.)
#         globalme[(ens, varna+'_mid')] = ctl.band_mean_from_zonal(zonal_mean, lata, 40., 60.)
#         globalme[(ens, varna+'_arctic')] = ctl.band_mean_from_zonal(zonal_mean, lata, 70., 90.)
#
# for varna in ['tas', 'tas_land', 'tas_ocean', 'tas_trop', 'tas_mid', 'tas_arctic']:
#     for nunam, expnam in zip(nunams, expnams):
#         globalme[(nunam, varna)] = np.mean([globalme[(exp, varna)] for exp in ensmemsall if expnam in exp], axis = 0)
#
# for key in globalme:
#     globalme[key] = globalme[key]-KtoC
#
# pickle.dump(globalme, open(cart_out+'globalme_yearly_AMIP.p', 'wb'))
globalme = pickle.load(open(cart_out+'globalme_yearly_AMIP.p', 'rb'))

# figure

# voglio:
#- un pdf con tutti i plot heat_flux zonal (ogni dieci, base, stoc e diff)
#- un pdf con plot toa_balance yearly(base, stoc e diff), surf_balance yearly (land e ocean nello stesso)(base, stoc e diff), toa_balance vs temperature (base, stoc e diff), heat_flu in qualche punto e temperatura (per detrendizzare le differenze)


figure_file = cart_out+'heat_flux_base_vs_stoc.pdf'
figures = []
axes = []
base_hist = radclim[('zonal', 'base_hist', 'heat_flux')]
print(base_hist.shape)
stoc_hist = radclim[('zonal', 'stoc_hist', 'heat_flux')]
base_fut = radclim[('zonal', 'base_fut', 'heat_flux')]
stoc_fut = radclim[('zonal', 'stoc_fut', 'heat_flux')]
fig = plt.figure()
ax = plt.subplot(1,1,1)
plt.title('Northward heat flux - HIST')
okplobase = np.mean(base_hist, axis = 0)
okplostoc = np.mean(stoc_hist, axis = 0)
ax.plot(lat, okplobase, label = 'base hist', linewidth = 2.0)
ax.plot(lat, okplostoc, label = 'stoc hist', linewidth = 2.0)
ax.plot(lat, 10*(okplostoc-okplobase), label = '(diff s-b) x 10', linewidth = 1.0, linestyle = '--')
plt.xlabel('Latitude')
plt.ylabel('Heat flux (W)')
plt.legend()
plt.grid()
figures.append(fig)
axes.append(ax)

fig = plt.figure()
ax = plt.subplot(1,1,1)
plt.title('Northward heat flux - FUT')
okplobase = np.mean(base_fut, axis = 0)
okplostoc = np.mean(stoc_fut, axis = 0)
ax.plot(lat, okplobase, label = 'base fut', linewidth = 2.0)
ax.plot(lat, okplostoc, label = 'stoc fut', linewidth = 2.0)
ax.plot(lat, 10*(okplostoc-okplobase), label = '(diff s-b) x 10', linewidth = 1.0, linestyle = '--')
plt.xlabel('Latitude')
plt.ylabel('Heat flux (W)')
plt.legend()
plt.grid()
figures.append(fig)
axes.append(ax)

ctl.adjust_ax_scale(axes, sel_axis = 'both')

ctl.plot_pdfpages(figure_file, figures)
plt.close('all')


figure_file = cart_out + 'global_diff_net_fluxes.pdf'
figures = []
axes = []

fig = plt.figure()
ax = plt.subplot(1,1,1)
plt.title('Net flux at TOA')
for lett, tim in zip(['a','f'], ['_hist', '_fut']):
    for cosone, ls in zip(['stoc', 'base'], [':', '--']):
        coso = radclim[('global', cosone+tim, 'toa_balance')]
        ax.plot(years[lett], coso, linewidth = 0.5, color = 'grey', linestyle = ls)
        rollcoso = ctl.running_mean(coso, wnd = 10)
        ax.plot(years[lett], rollcoso, label = cosone+tim, linewidth = 2.0)
plt.xlabel('Year')
plt.ylabel('Rad. forcing (W/m^2)')
plt.legend()
plt.grid()
figures.append(fig)
axes.append(ax)

# fig = plt.figure()
# ax = plt.subplot(1,1,1)
# cosone = 'diff s-b'
# coso = radclim[('global', 'stoc', 'toa_balance')]-radclim[('global', 'base', 'toa_balance')]
# plt.title('Net flux at TOA - ({})'.format(cosone))
# ax.plot(years, coso, label = 'base', linewidth = 0.5, color = 'grey')
# rollcoso = ctl.running_mean(coso, wnd = 10)
# ax.plot(years, rollcoso, label = 'base', linewidth = 2.0, color = 'green')
# plt.xlabel('Year')
# plt.ylabel('Rad. forcing (W/m^2)')
# plt.grid()
# figures.append(fig)
# #axes.append(ax)
#
# #ctl.adjust_ax_scale(axes, sel_axis = 'both')

axes = []

fig = plt.figure()
ax = plt.subplot(1,1,1)
plt.title('transient gregory plot')
for lett, tim in zip(['a','f'], ['_hist', '_fut']):
    for cosone, ls, mark in zip(['stoc', 'base'], [':', '--'], ['o','d']):
        coso = radclim[('global', cosone+tim, 'toa_balance')]
        tas = globalme[(cosone+tim,'tas')]
        print(coso.shape)
        print(tas.shape)
        ax.plot(tas, coso, linewidth = 0.5, color = 'grey', linestyle = ls)
        ax.scatter(tas, coso, color = 'darkgrey', marker = mark, s = 1)
        rollcoso = ctl.running_mean(coso, wnd = 10)
        rolltas = ctl.running_mean(tas, wnd = 10)
        ax.plot(rolltas, rollcoso, label = cosone+tim, linewidth = 2.0)
plt.xlabel('Temp (C)')
plt.ylabel('Net TOA forcing (W/m^2)')
plt.legend()
plt.grid()
figures.append(fig)

# fig = plt.figure()
# ax = plt.subplot(1,1,1)
# cosone = 'diff s-b'
# coso = radclim[('global', 'stoc', 'toa_balance')] - radclim[('global', 'base', 'toa_balance')]
# tas = globalme[('stoc','tas')] - globalme[('base','tas')]
# plt.title('gregory plot - ({})'.format(cosone))
# ax.plot(tas, coso, label = 'stoc', linewidth = 0.5, color = 'grey')
# ax.scatter(tas, coso, label = 'stoc', linewidth = 1.0, color = 'darkgrey')
# rollcoso = ctl.running_mean(coso, wnd = 10)
# rolltas = ctl.running_mean(tas, wnd = 10)
# ax.plot(rolltas, rollcoso, label = 'base', linewidth = 2.0, color = 'green')
# plt.xlabel('Temp diff (C)')
# plt.ylabel('Net TOA forcing diff (W/m^2)')
# plt.grid()
# figures.append(fig)

fig = plt.figure()
ax = plt.subplot(1,1,1)
plt.title('Upward LW radiation at TOA')
for lett, tim in zip(['a','f'], ['_hist', '_fut']):
    for cosone, ls, mark in zip(['stoc', 'base'], [':', '--'], ['o','d']):
        coso = radclim[('global', cosone+tim, 'rlut')]
        tas = globalme[(cosone+tim,'tas')]
        ax.plot(tas, coso, linewidth = 0.5, color = 'grey', linestyle = ls)
        ax.scatter(tas, coso, color = 'darkgrey', marker = mark, s = 1)
        rollcoso = ctl.running_mean(coso, wnd = 10)
        rolltas = ctl.running_mean(tas, wnd = 10)
        ax.plot(rolltas, rollcoso, label = cosone+tim, linewidth = 2.0)
plt.xlabel('Temp (C)')
plt.ylabel('Upward LW (W/m^2)')
plt.legend()
plt.grid()
figures.append(fig)

# fig = plt.figure()
# ax = plt.subplot(1,1,1)
# cosone = 'diff s-b'
# coso = radclim[('global', 'stoc', 'rlut')] - radclim[('global', 'base', 'rlut')]
# tas = globalme[('stoc','tas')] - globalme[('base','tas')]
# plt.title('Upward LW radiation at TOA - {}'.format(cosone))
# ax.plot(tas, coso, label = 'stoc', linewidth = 0.5, color = 'grey')
# ax.scatter(tas, coso, label = 'stoc', linewidth = 1.0, color = 'darkgrey')
# rollcoso = ctl.running_mean(coso, wnd = 10)
# rolltas = ctl.running_mean(tas, wnd = 10)
# ax.plot(rolltas, rollcoso, label = 'base', linewidth = 2.0, color = 'green')
# plt.xlabel('Temp diff (C)')
# plt.ylabel('Upward LW (W/m^2)')
# plt.grid()
# figures.append(fig)
# #axes.append(ax)

fig = plt.figure()
ax = plt.subplot(1,1,1)
plt.title('Upward SW radiation at TOA')
for lett, tim in zip(['a','f'], ['_hist', '_fut']):
    for cosone, ls, mark in zip(['stoc', 'base'], [':', '--'], ['o','d']):
        coso = radclim[('global', cosone+tim, 'rsut')]
        tas = globalme[(cosone+tim,'tas')]
        ax.plot(tas, coso, linewidth = 0.5, color = 'grey', linestyle = ls)
        ax.scatter(tas, coso, color = 'darkgrey', marker = mark, s = 1)
        rollcoso = ctl.running_mean(coso, wnd = 10)
        rolltas = ctl.running_mean(tas, wnd = 10)
        ax.plot(rolltas, rollcoso, label = cosone+tim, linewidth = 2.0)
plt.xlabel('Temp (C)')
plt.ylabel('Upward SW (W/m^2)')
plt.legend()
plt.grid()
figures.append(fig)

# fig = plt.figure()
# ax = plt.subplot(1,1,1)
# cosone = 'diff s-b'
# coso = radclim[('global', 'stoc', 'rsut')] - radclim[('global', 'base', 'rsut')]
# tas = globalme[('stoc','tas')] - globalme[('base','tas')]
# plt.title('Upward SW radiation at TOA - {}'.format(cosone))
# ax.plot(tas, coso, label = 'stoc', linewidth = 0.5, color = 'grey')
# ax.scatter(tas, coso, label = 'stoc', linewidth = 1.0, color = 'darkgrey')
# rollcoso = ctl.running_mean(coso, wnd = 10)
# rolltas = ctl.running_mean(tas, wnd = 10)
# ax.plot(rolltas, rollcoso, label = 'base', linewidth = 2.0, color = 'green')
# plt.xlabel('Temp diff (C)')
# plt.ylabel('Upward SW (W/m^2)')
# plt.grid()
# figures.append(fig)

axes = []

indlat = 200

for indlat in [160, 200, 240]:
    fig = plt.figure()
    ax = plt.subplot(1,1,1)
    plt.title('heat flux at {:3.0f} vs tas'.format(lat[indlat]))
    for lett, tim in zip(['a','f'], ['_hist', '_fut']):
        for cosone, ls, mark in zip(['stoc', 'base'], [':', '--'], ['o','d']):
            coso = radclim[('zonal', cosone+tim, 'heat_flux')][:, indlat]
            tas = globalme[(cosone+tim,'tas')]
            ax.plot(tas, coso, linewidth = 0.5, color = 'grey', linestyle = ls)
            ax.scatter(tas, coso, color = 'darkgrey', marker = mark, s = 1)
            rollcoso = ctl.running_mean(coso, wnd = 10)
            rolltas = ctl.running_mean(tas, wnd = 10)
            ax.plot(rolltas, rollcoso, label = cosone+tim, linewidth = 2.0)
    plt.xlabel('Temp (C)')
    plt.ylabel('Heat flux (W)')
    plt.legend()
    plt.grid()
    figures.append(fig)

    # fig = plt.figure()
    # ax = plt.subplot(1,1,1)
    # cosone = 'diff s-b'
    # coso = radclim[('zonal', 'stoc', 'heat_flux')][:, indlat] - radclim[('zonal', 'base', 'heat_flux')][:, indlat]
    # tas = globalme[('stoc','tas')] - globalme[('base','tas')]
    # plt.title('heat flux at {:3.0f} vs tas - ({})'.format(lat[indlat], cosone))
    # ax.plot(tas, coso, label = 'stoc', linewidth = 0.5, color = 'grey')
    # ax.scatter(tas, coso, label = 'stoc', linewidth = 1.0, color = 'darkgrey')
    # rollcoso = ctl.running_mean(coso, wnd = 10)
    # rolltas = ctl.running_mean(tas, wnd = 10)
    # ax.plot(rolltas, rollcoso, label = 'base', linewidth = 2.0, color = 'green')
    # plt.xlabel('Temp diff (C)')
    # plt.ylabel('Heat flux diff (W/m^2)')
    # plt.grid()
    # figures.append(fig)
#axes.append(ax)

for zollo in ['trop', 'mid', 'arctic']:
    fig = plt.figure()
    ax = plt.subplot(1,1,1)
    plt.title('mean temp {} vs tas'.format(zollo))
    for lett, tim in zip(['a','f'], ['_hist', '_fut']):
        for cosone, ls, mark in zip(['stoc', 'base'], [':', '--'], ['o','d']):
            coso = globalme[(cosone+tim, 'tas_'+zollo)]
            tas = globalme[(cosone+tim,'tas')]
            ax.plot(tas, coso, linewidth = 0.5, color = 'grey', linestyle = ls)
            ax.scatter(tas, coso, color = 'darkgrey', marker = mark, s = 1)
            rollcoso = ctl.running_mean(coso, wnd = 10)
            rolltas = ctl.running_mean(tas, wnd = 10)
            ax.plot(rolltas, rollcoso, label = cosone+tim, linewidth = 2.0)
    plt.xlabel('Temp (C)')
    plt.ylabel('Temp (C)')
    plt.legend()
    plt.grid()
    figures.append(fig)

    # fig = plt.figure()
    # ax = plt.subplot(1,1,1)
    # cosone = 'diff s-b'
    # coso = globalme[('stoc', 'tas_'+zollo)] - globalme[('base', 'tas_'+zollo)]
    # tas = globalme[('stoc','tas')] - globalme[('base','tas')]
    # plt.title('mean temp {} vs tas - {}'.format(zollo, cosone))
    # ax.plot(tas, coso, label = 'stoc', linewidth = 0.5, color = 'grey')
    # ax.scatter(tas, coso, label = 'stoc', linewidth = 1.0, color = 'darkgrey')
    # rollcoso = ctl.running_mean(coso, wnd = 10)
    # rolltas = ctl.running_mean(tas, wnd = 10)
    # ax.plot(rolltas, rollcoso, label = 'base', linewidth = 2.0, color = 'green')
    # plt.xlabel('Temp diff (C)')
    # plt.ylabel('Temp diff (C)')
    # plt.grid()
    # figures.append(fig)

#
# for stokko in ['', '_land', '_ocean']:
#     axes = []
#     fig = plt.figure()
#     ax = plt.subplot(1,1,1)
#     plt.title('forcing at surface{}'.format(stokko))
#     for cosone, ls, mark in zip(['stoc', 'base'], [':', '--'], ['o','d']):
#         coso = radclim[('global', cosone, 'surf_balance'+stokko)]
#         tas = globalme[(cosone,'tas')]
#         ax.plot(tas, coso, linewidth = 0.5, color = 'grey', linestyle = ls)
#         ax.scatter(tas, coso, color = 'darkgrey', marker = mark)
#         rollcoso = ctl.running_mean(coso, wnd = 10)
#         rolltas = ctl.running_mean(tas, wnd = 10)
#         ax.plot(rolltas, rollcoso, label = cosone, linewidth = 2.0)
#     plt.xlabel('Temp (C)')
#     plt.ylabel('Net surface forcing (W/m^2)')
#     plt.grid()
#     plt.legend()
#     figures.append(fig)
#     axes.append(ax)
#
#     fig = plt.figure()
#     ax = plt.subplot(1,1,1)
#     cosone = 'diff s-b'
#     coso = radclim[('global', 'stoc', 'surf_balance'+stokko)] - radclim[('global', 'base', 'surf_balance'+stokko)]
#     tas = globalme[('stoc','tas'+stokko)] - globalme[('base','tas'+stokko)]
#     plt.title('forcing at surface{} - ({})'.format(stokko, cosone))
#     ax.plot(tas, coso, label = 'stoc', linewidth = 0.5, color = 'grey')
#     ax.scatter(tas, coso, label = 'stoc', linewidth = 1.0, color = 'darkgrey')
#     rollcoso = ctl.running_mean(coso, wnd = 10)
#     rolltas = ctl.running_mean(tas, wnd = 10)
#     ax.plot(rolltas, rollcoso, label = 'base', linewidth = 2.0, color = 'green')
#     plt.xlabel('Temp (C)')
#     plt.ylabel('Net surface forcing (W/m^2)')
#     plt.grid()
#     figures.append(fig)
# #    axes.append(ax)
#
# for stokko in ['', '_land', '_ocean']:
#     axes = []
#     fig = plt.figure()
#     ax = plt.subplot(1,1,1)
#     plt.title('net atmos forcing {}'.format(stokko))
#     for cosone, ls, mark in zip(['stoc', 'base'], [':', '--'], ['o','d']):
#         coso = radclim[('global', cosone, 'toa_balance'+stokko)] - radclim[('global', cosone, 'surf_balance'+stokko)]
#         tas = globalme[(cosone,'tas')]
#         ax.plot(tas, coso, linewidth = 0.5, color = 'grey', linestyle = ls)
#         ax.scatter(tas, coso, color = 'darkgrey', marker = mark)
#         rollcoso = ctl.running_mean(coso, wnd = 10)
#         rolltas = ctl.running_mean(tas, wnd = 10)
#         ax.plot(rolltas, rollcoso, label = cosone, linewidth = 2.0)
#     plt.xlabel('Temp (C)')
#     plt.ylabel('Net atmos forcing (W/m^2)')
#     plt.grid()
#     plt.legend()
#     figures.append(fig)
#     axes.append(ax)
#
#     fig = plt.figure()
#     ax = plt.subplot(1,1,1)
#     cosone = 'stoc'
#     coso1 = radclim[('global', cosone, 'toa_balance'+stokko)] - radclim[('global', cosone, 'surf_balance'+stokko)]
#     cosone = 'base'
#     coso2 = radclim[('global', cosone, 'toa_balance'+stokko)] - radclim[('global', cosone, 'surf_balance'+stokko)]
#     cosone = 'diff s-b'
#     coso = coso1-coso2
#     tas = globalme[('stoc','tas'+stokko)] - globalme[('base','tas'+stokko)]
#     plt.title('net atmos forcing {} - ({})'.format(stokko, cosone))
#     ax.plot(tas, coso, label = 'stoc', linewidth = 0.5, color = 'grey')
#     ax.scatter(tas, coso, label = 'stoc', linewidth = 1.0, color = 'darkgrey')
#     rollcoso = ctl.running_mean(coso, wnd = 10)
#     rolltas = ctl.running_mean(tas, wnd = 10)
#     ax.plot(rolltas, rollcoso, label = 'base', linewidth = 2.0, color = 'green')
#     plt.xlabel('Temp (C)')
#     plt.ylabel('Net atmos forcing (W/m^2)')
#     plt.grid()
#     figures.append(fig)
# #    axes.append(ax)

ctl.plot_pdfpages(figure_file, figures)
plt.close('all')

# poi voglio figura con mappe di toa_balance e surf_balance (diff base/stoc)
# poi balance zonal anno per anno (base/stoc/diff)
