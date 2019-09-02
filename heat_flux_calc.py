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
from scipy import interpolate as itrp
import itertools as itt

from sklearn.cluster import KMeans



from datetime import datetime
import pickle

import climtools_lib as ctl
import climdiags as cd

from copy import deepcopy as cp

##############################################
### constants
# cp = 4185.5 # J kg-1 K-1 QUESTA E L'ACQUA COJOOOO
cp = 1005.0 # specific enthalpy dry air - J kg-1 K-1
cpw = 1840.0 # specific enthalpy water vapor - J kg-1 K-1
# per la moist air sarebbe cpm = cp + q*cpw, ma conta al massimo per un 3 %
L = 2501000.0 # J kg-1
Lsub = 2835000.0 # J kg-1
g = 9.81 # m s-2
Rearth = 6371.0e3 # mean radius
#####################################

#plt.ion()
shortnam = {'mshf': 'SH', 'mpef': 'PE', 'mlhf': 'LH'}
fluxnames = ['mshf', 'mpef', 'mlhf']
fluxlongnames = {'mshf':'Sensible Heat', 'mpef':'Potential Energy', 'mlhf':'Latent Heat'}
factors = [cp/g, 1., L/g]

cart_out = '/home/fabiano/Research/lavori/SPHINX_for_lisboa/heat_flux/'

figure_file_all = cart_out+'heat_fluxes_lcs0_vs_lcb0.pdf'
figures_all = []

zonal_margins = dict()
zonal_margins['tot'] = (-3.5e16, 4.5e16)
zonal_margins['mshf'] = (-3.5e16, 4.5e16)
zonal_margins['mpef'] = (-3.5e16, 4.5e16)
zonal_margins['mlhf'] = (-7e15, 9.e15)
map_margins = dict()
map_margins['mshf'] = (-1.5e10, 1.5e10)
map_margins['mpef'] = (-4.e9, 4.e9)
map_margins['mlhf'] = (-7.e8, 7.e8)
map_margins['tot'] = (-2.e10, 2.e10)

# Loading reference pressure file
pressurefile = '/data-hobbes/fabiano/SPHINX/heat_flux/1988_daily/lcs0_day_1988_ps.nc'
press0row, lat, lon, datespress, time_units, var_units = ctl.read3Dncfield(pressurefile)
press0 = dict()
press0['DJF'] = np.mean(ctl.sel_season(press0row, datespress, 'DJF', cut = False)[0], axis = 0)
press0['JJA'] = np.mean(ctl.sel_season(press0row, datespress, 'JJA', cut = False)[0], axis = 0)
press0['year'] = np.mean(press0row, axis = 0)

# Loading ERA reference
cart_era = '/data-hobbes/fabiano/OBS/ERA/ERAInterim/'
era_fi = 'prova_heatflux_1988_MM.nc'
cpc = nc.Dataset(cart_era+era_fi)
era_lat = cpc.variables['latitude'][:]
era_lon = cpc.variables['longitude'][:]
era_zonal_factor = 2*np.pi*Rearth*np.cos(np.deg2rad(era_lat))

era_fluxes_maps = dict()
era_fluxes_zonal = dict()
era_fluxes_zonal_itrp = dict()

okmon = [0,1,11]
era_fluxes_maps[('DJF','tot')] = np.mean(list(cpc.variables.values())[-1][okmon, ...], axis = 0)
era_fluxes_maps[('DJF','mpef')] = np.mean(list(cpc.variables.values())[-2][okmon, ...], axis = 0)
era_fluxes_maps[('DJF','mshf')] = np.mean(list(cpc.variables.values())[-3][okmon, ...], axis = 0)
okmon = [5,6,7]
era_fluxes_maps[('JJA','tot')] = np.mean(list(cpc.variables.values())[-1][okmon, ...], axis = 0)
era_fluxes_maps[('JJA','mpef')] = np.mean(list(cpc.variables.values())[-2][okmon, ...], axis = 0)
era_fluxes_maps[('JJA','mshf')] = np.mean(list(cpc.variables.values())[-3][okmon, ...], axis = 0)

era_fluxes_maps[('year','tot')] = np.mean(list(cpc.variables.values())[-1][:], axis = 0)
era_fluxes_maps[('year','mpef')] = np.mean(list(cpc.variables.values())[-2][:], axis = 0)
era_fluxes_maps[('year','mshf')] = np.mean(list(cpc.variables.values())[-3][:], axis = 0)

fi = 'northward_water_flux.nc'
cpc = nc.Dataset(cart_era+fi)
okmon = [0,1,11]
era_fluxes_maps[('DJF','mlhf')] = L*np.mean(list(cpc.variables.values())[-1][okmon, ...], axis = 0)
okmon = [5,6,7]
era_fluxes_maps[('JJA','mlhf')] = L*np.mean(list(cpc.variables.values())[-1][okmon, ...], axis = 0)
era_fluxes_maps[('year','mlhf')] = L*np.mean(list(cpc.variables.values())[-1][:], axis = 0)

for fu in era_fluxes_maps:
     era_fluxes_zonal[fu] = np.mean(era_fluxes_maps[fu], axis = 1)*era_zonal_factor

# ERA ref figures:
exp = 'ERA'
year = 1988
figures_era = []
for seas in press0:
    fig = plt.figure()
    plt.title('Meridional heat fluxes - {} {} {}'.format(exp, year, seas))
    for flun in fluxnames:
        plt.plot(era_lat, era_fluxes_zonal[(seas, flun)], label = shortnam[flun])
    total = np.sum([era_fluxes_zonal[(seas, flun)] for flun in fluxnames], axis = 0)
    plt.plot(era_lat, total, label = 'Total')
    plt.legend()
    plt.grid()
    plt.xlabel('Latitude')
    plt.ylabel('Integrated Net Heat Flux (W)')
    figures_era.append(fig)

    for flun in fluxnames:
        fig = ctl.plot_map_contour(era_fluxes_maps[(seas, flun)], era_lat, era_lon, title = fluxlongnames[flun]+' - {} {} {}'.format(exp, year, seas))
        figures_era.append(fig)

    total = np.sum([era_fluxes_maps[(seas, flun)] for flun in fluxnames], axis = 0)
    fig = ctl.plot_map_contour(total, era_lat, era_lon, title = 'Total meridional flux - {} {} {}'.format(exp, year, seas))
    figures_era.append(fig)

file_era = cart_out+'ERA_reference_fluxes_1988.pdf'
ctl.plot_pdfpages(file_era, figures_era)

ally = [1975, 2000, 2025, 2050, 2075, 2090, 2100]
# exp = 'lcs0'
freq = '6hrs'
# year = 1988
# Ora per gli exps
#cart_in = '/data-hobbes/fabiano/SPHINX/heat_flux/lcs0_1988_6hrs/'
cart_in = '/data-hobbes/fabiano/SPHINX/heat_flux/data_6hrs/'

allfluxes = dict()
# for exp in ['lcb0', 'lcs0']:
#     for year in [1975, 2000, 2025, 2050, 2075, 2090, 2100]:
for exp in ['las1']:
    for year in [1988]:
        pressurefile = cart_in+'{}_day_{}_ps.nc'.format(exp, year)
        press0row, latdad, londsad, datespress, time_units, var_units = ctl.read3Dncfield(pressurefile)
        press0 = dict()
        press0['DJF'] = np.mean(ctl.sel_season(press0row, datespress, 'DJF', cut = False)[0], axis = 0)
        press0['JJA'] = np.mean(ctl.sel_season(press0row, datespress, 'JJA', cut = False)[0], axis = 0)
        press0['year'] = np.mean(press0row, axis = 0)

        figure_file_exp = cart_out+'hf_{}_{}.pdf'.format(exp, year)
        figures_exp = []
        figure_file_exp_maps = cart_out+'hf_{}_{}_maps.pdf'.format(exp, year)
        figures_exp_maps = []

        fluxes_row = dict()
        fluxes_levels = dict()
        fluxes_cross = dict()
        fluxes_maps = dict()
        fluxes_zonal = dict()

        vars = dict()
        varnames = ['va', 'ta', 'zg', 'hus']
        fils = ['{}_{}_{}_{}.nc'.format(exp, freq, year, varna) for varna in varnames]

        var, level, lat, lon, dates, time_units, var_units, time_cal = ctl.read4Dncfield(cart_in+fils[0])
        vars[varnames[0]] = var

        var, level, lat, lon, dates, time_units, var_units, time_cal = ctl.read4Dncfield(cart_in+fils[1])
        vars[varnames[1]] = var
        flun = 'mshf'
        fluxes_row[flun] = factors[0]*vars['va']*vars['ta']
        del vars['ta']
        fluxes_levels[('DJF',flun)] = np.mean(ctl.sel_season(fluxes_row[flun], dates, 'DJF', cut = False)[0], axis = 0)
        fluxes_levels[('JJA',flun)] = np.mean(ctl.sel_season(fluxes_row[flun], dates, 'JJA', cut = False)[0], axis = 0)
        fluxes_levels[('year',flun)] = np.mean(fluxes_row[flun], axis = 0)
        del fluxes_row[flun]

        var, level, lat, lon, dates, time_units, var_units, time_cal = ctl.read4Dncfield(cart_in+fils[2])
        vars[varnames[2]] = var
        flun = 'mpef'
        fluxes_row[flun] = vars['va']*vars['zg']
        del vars['zg']
        fluxes_levels[('DJF',flun)] = np.mean(ctl.sel_season(fluxes_row[flun], dates, 'DJF', cut = False)[0], axis = 0)
        fluxes_levels[('JJA',flun)] = np.mean(ctl.sel_season(fluxes_row[flun], dates, 'JJA', cut = False)[0], axis = 0)
        fluxes_levels[('year',flun)] = np.mean(fluxes_row[flun], axis = 0)
        del fluxes_row[flun]

        var, level, lat, lon, dates, time_units, var_units, time_cal = ctl.read4Dncfield(cart_in+fils[3])
        vars[varnames[3]] = var
        flun = 'mlhf'
        fluxes_row[flun] = factors[2]*vars['va']*vars['hus']
        del vars['hus']

        fluxes_levels[('DJF',flun)] = np.mean(ctl.sel_season(fluxes_row[flun], dates, 'DJF', cut = False)[0], axis = 0)
        fluxes_levels[('JJA',flun)] = np.mean(ctl.sel_season(fluxes_row[flun], dates, 'JJA', cut = False)[0], axis = 0)
        fluxes_levels[('year',flun)] = np.mean(fluxes_row[flun], axis = 0)
        del fluxes_row[flun]

        # for flun in fluxes_row:
        #     fluxes_levels[('DJF',flun)] = np.mean(ctl.sel_season(fluxes_row[flun], dates, 'DJF')[0], axis = 0)
        #     fluxes_levels[('JJA',flun)] = np.mean(ctl.sel_season(fluxes_row[flun], dates, 'JJA')[0], axis = 0)
        #     fluxes_levels[('year',flun)] = np.mean(fluxes_row[flun], axis = 0)

        del fluxes_row, vars

        zonal_factor = 2*np.pi*Rearth*np.cos(np.deg2rad(lat))
        for flun in fluxes_levels:
            fluxes_maps[flun] = np.zeros(fluxes_levels[flun].shape[1:])
            fluxes_cross[flun] = np.mean(fluxes_levels[flun], axis = -1)*zonal_factor

        print('Starting vertical integration\n')
        for seas in press0:
            for ila in range(len(lat)):
                for ilo in range(len(lon)):
                    #print(lat[ila], lon[ilo])
                    p0 = press0[seas][ila, ilo]

                    # Selecting pressure levels lower than surface pressure
                    exclude_pres = level > p0
                    if np.any(exclude_pres):
                        lev0 = np.argmax(exclude_pres)
                    else:
                        lev0 = None
                    levels_ok = np.append(level[:lev0], p0)

                    if not np.all(np.diff(levels_ok) > 0):
                        print(levels_ok)
                        raise ValueError('x not increasing')

                    for flun in fluxnames:
                        coso = fluxes_levels[(seas, flun)][:, ila, ilo]
                        coso = np.append(coso[:lev0], np.interp(p0, level, coso))
                        fluxes_maps[(seas, flun)][ila, ilo] = np.trapz(coso, x = levels_ok)
                        #print(levels_ok, coso.mean(), coso[-1], fluxes_maps[(seas, flun)][ila, ilo])

        print('done\n')
        if not os.path.exists(cart_out+'Out_fluxes_nc/'): os.mkdir(cart_out+'Out_fluxes_nc/')

        for fu in fluxes_maps:
            seas, flun = fu
            namefi = '{}_{}_flux_{}_{}.nc'.format(exp, year, flun, seas)
            ctl.save2Dncfield(lat,lon,fluxes_maps[fu],flun,cart_out+'Out_fluxes_nc/'+namefi)
        # for seas in ['DJF', 'JJA', 'year']:
        #     for flun in fluxnames:
        #         namefi = '{}_{}_flux_{}_{}.nc'.format(exp, year, flun, seas)
        #         fluxes_maps[(seas, flun)] = ctl.read2Dncfield(cart_out+'Out_fluxes_nc/'+namefi)[0]

        zonal_factor = 2*np.pi*Rearth*np.cos(np.deg2rad(lat))
        for fu in fluxes_maps:
            fluxes_zonal[fu] = np.mean(fluxes_maps[fu], axis = 1)*zonal_factor

        for fu in era_fluxes_zonal:
            #print(lat, era_lat, era_fluxes_zonal[fu])
            era_fluxes_zonal_itrp[fu] = np.interp(lat, era_lat[::-1], era_fluxes_zonal[fu][::-1])
            print('{:12.3e} - {:12.3e}\n'.format(era_fluxes_zonal_itrp[fu].min(), era_fluxes_zonal_itrp[fu].max()))

        # Single exps, single season figures:
        for seas in press0:
            fig = plt.figure()
            plt.title('SPHINX Meridional heat fluxes - {} {} {}'.format(exp, year, seas))
            for flun in fluxnames:
                plt.plot(lat, fluxes_zonal[(seas, flun)], label = shortnam[flun])
            total = np.sum([fluxes_zonal[(seas, flun)] for flun in fluxnames], axis = 0)
            fluxes_zonal[(seas, 'tot')] = total
            plt.plot(lat, total, label = 'Total')
            plt.legend()
            plt.grid()
            plt.ylim(zonal_margins['tot'])
            plt.xlabel('Latitude')
            plt.ylabel('Integrated Net Heat Flux (W)')
            figures_exp.append(fig)

            fig = plt.figure()
            plt.title('ERA Meridional heat fluxes - {} {}'.format(year, seas))
            for flun in fluxnames:
                plt.plot(lat, era_fluxes_zonal_itrp[(seas, flun)], label = shortnam[flun])
            total_era = np.sum([era_fluxes_zonal_itrp[(seas, flun)] for flun in fluxnames], axis = 0)
            plt.plot(lat, total_era, label = 'Total')
            plt.legend()
            plt.grid()
            plt.ylim(zonal_margins['tot'])
            plt.xlabel('Latitude')
            plt.ylabel('Integrated Net Heat Flux (W)')
            figures_exp.append(fig)

            fig = plt.figure()
            plt.title('SPHINX - ERA diff - {} {}'.format(year, seas))
            for flun in fluxnames:
                plt.plot(lat, fluxes_zonal[(seas, flun)]-era_fluxes_zonal_itrp[(seas, flun)], label = shortnam[flun])
            plt.plot(lat, total-total_era, label = 'Total')
            plt.legend()
            plt.grid()
            plt.ylim(zonal_margins['tot'])
            plt.xlabel('Latitude')
            plt.ylabel('Integrated Net Heat Flux (W)')
            figures_exp.append(fig)

            for flun in fluxnames:
                fig = ctl.plot_map_contour(fluxes_maps[(seas, flun)], lat, lon, title = fluxlongnames[flun]+' - {} {} {}'.format(exp, year, seas), cbar_range = map_margins[flun])
                figures_exp_maps.append(fig)

            total = np.sum([fluxes_maps[(seas, flun)] for flun in fluxnames], axis = 0)
            fig = ctl.plot_map_contour(total, lat, lon, title = 'Total meridional flux - {} {} {}'.format(exp, year, seas), cbar_range = map_margins['tot'])
            figures_exp_maps.append(fig)

            pinotot = itrp.RectBivariateSpline(lat, lon, total)
            total_eragrid = np.zeros(era_fluxes_maps[(seas, 'tot')].shape)
            for lae, eral in enumerate(era_lat):
                for loe, erol in enumerate(era_lon):
                    total_eragrid[lae, loe] = pinotot(eral, erol)
            diff = total_eragrid - era_fluxes_maps[(seas, 'tot')]
            fig = ctl.plot_map_contour(diff, era_lat, era_lon, title = 'Difference vs ERA flux - {} {} {}'.format(exp, year, seas), plot_anomalies = True)
            figures_exp_maps.append(fig)

        for flun in fluxnames:
            fig = plt.figure()
            plt.title('{} fluxes - {}'.format(shortnam[flun], year))
            cset = ctl.color_set(3, bright_thres = 1)
            for seas, col in zip(press0.keys(), cset):
                plt.plot(lat, fluxes_zonal[(seas, flun)], label = seas, color = col, linewidth = 2.)
                plt.plot(lat, era_fluxes_zonal_itrp[(seas, flun)], label = 'ERA '+seas, color = col, linewidth = 0.7, linestyle = '--')
            plt.legend()
            plt.grid()
            plt.ylim(zonal_margins[flun])
            plt.xlabel('Latitude')
            plt.ylabel('Integrated Net Heat Flux (W)')
            figures_exp.append(fig)

        flun = 'tot'
        fig = plt.figure()
        plt.title('{} fluxes - {}'.format('tot', year))
        cset = ctl.color_set(3, bright_thres = 1)
        for seas, col in zip(press0.keys(), cset):
            plt.plot(lat, fluxes_zonal[(seas, flun)], label = seas, color = col, linewidth = 2.)
            plt.plot(lat, era_fluxes_zonal_itrp[(seas, flun)], label = 'ERA '+seas, color = col, linewidth = 0.7, linestyle = '--')
        plt.legend()
        plt.grid()
        plt.ylim(zonal_margins[flun])
        plt.xlabel('Latitude')
        plt.ylabel('Integrated Net Heat Flux (W)')
        figures_exp.append(fig)

        print('pitototot')
        ctl.plot_pdfpages(figure_file_exp, figures_exp)
        ctl.plot_pdfpages(figure_file_exp_maps, figures_exp_maps)

        plt.close('all')
        allfluxes[(exp, year, 'cross')] = fluxes_cross
        allfluxes[(exp, year, 'maps')] = fluxes_maps
        allfluxes[(exp, year, 'zonal')] = fluxes_zonal

pickle.dump(allfluxes, open(cart_out+'allfluxes_las1_26112018.p','wb'))
sys.exit()

# pickle.dump(allfluxes, open(cart_out+'allfluxes_23112018.p','wb'))
allfluxes = pickle.load(open(cart_out+'allfluxes_23112018.p','rb'))
# comparison between exps

for fu in era_fluxes_zonal:
    era_fluxes_zonal_itrp[fu] = np.interp(lat, era_lat[::-1], era_fluxes_zonal[fu][::-1])

figures = []
figure_file = cart_out + 'lcb0_vs_lcs0_fluxes.pdf'
seasons = ['year', 'JJA', 'DJF']
level = np.array([92500, 85000, 70000, 60000, 50000, 40000, 25000, 10000])
level = level[::-1]

# figura yearly/JJA/DJF total con tutti i base e tutti gli stoc per tutti gli anni
for flun in fluxnames:
    for year in ally:
        okflubase = allfluxes[('lcb0', year, 'zonal')]
        okflustoc = allfluxes[('lcs0', year, 'zonal')]
        fig = plt.figure()
        plt.title('{} fluxes - {}'.format(shortnam[flun], year))
        cset = ctl.color_set(3, bright_thres = 1)
        for seas, col in zip(seasons, cset):
            plt.plot(lat, okflubase[(seas, flun)], label = 'base '+seas, color = col, linewidth = 1.5)
            plt.plot(lat, okflustoc[(seas, flun)], label = 'stoc '+seas, color = col, linewidth = 1.5, linestyle = '--')
            plt.plot(lat, era_fluxes_zonal_itrp[(seas, flun)], label = 'ERA '+seas, color = col, linewidth = 0.7, linestyle = ':')
        plt.legend(fontsize = 'small')
        plt.grid()
        plt.ylim(zonal_margins[flun])
        plt.xlabel('Latitude')
        plt.ylabel('Integrated Net Heat Flux (W)')
        figures.append(fig)

flun = 'tot'
for year in ally:
    okflubase = allfluxes[('lcb0', year, 'zonal')]
    okflustoc = allfluxes[('lcs0', year, 'zonal')]
    fig = plt.figure()
    plt.title('{} fluxes - {}'.format('total', year))
    cset = ctl.color_set(3, bright_thres = 1)
    for seas, col in zip(seasons, cset):
        plt.plot(lat, okflubase[(seas, flun)], label = 'base '+seas, color = col, linewidth = 1.5)
        plt.plot(lat, okflustoc[(seas, flun)], label = 'stoc '+seas, color = col, linewidth = 1.5, linestyle = '--')
        plt.plot(lat, era_fluxes_zonal_itrp[(seas, flun)], label = 'ERA '+seas, color = col, linewidth = 0.7, linestyle = ':')
    plt.legend(fontsize = 'small')
    plt.grid()
    plt.ylim(zonal_margins[flun])
    plt.xlabel('Latitude')
    plt.ylabel('Integrated Net Heat Flux (W)')
    figures.append(fig)

# stesso per SH e PE
# poi mappe zonal di confronto con anno base (2100 vs 1975)
y1 = 1975
y2 = 2100
for exp in ['lcb0', 'lcs0']:
    okflubase1 = allfluxes[(exp, y1, 'maps')]
    okflubase2 = allfluxes[(exp, y2, 'maps')]
    print(okflubase1.keys())
    for flun in ['mshf', 'mpef']:
        for seas in seasons:
            diff = okflubase2[(seas, flun)]-okflubase1[(seas, flun)]
            fig = ctl.plot_map_contour(diff, lat, lon, title = '{} - {} diff - {} - {} vs {}'.format(shortnam[flun], exp, seas, y1, y2), plot_anomalies = True, cbar_range = map_margins[flun])
            figures.append(fig)

for y1 in [1975, 2100]:
    okflubase = allfluxes[('lcb0', y1, 'maps')]
    okflustoc = allfluxes[('lcs0', y1, 'maps')]
    for flun in ['mshf', 'mpef']:
        for seas in seasons:
            diff = okflustoc[(seas, flun)]-okflubase[(seas, flun)]
            fig = ctl.plot_map_contour(diff, lat, lon, title = '{} lcb0 vs lcs0 diff - {} - {}'.format(shortnam[flun], seas, y1), plot_anomalies = True, cbar_range = map_margins[flun])
            figures.append(fig)

# mappa di confronto fra base e stoc tra anni base e fut
# cross di confronto stesso tra anno base e fut e tra base e stoc

y1 = 1975
y2 = 2100
for exp in ['lcb0', 'lcs0']:
    okflubase1 = allfluxes[(exp, y1, 'cross')]
    okflubase2 = allfluxes[(exp, y2, 'cross')]
    for flun in ['mshf', 'mpef']:
        for seas in seasons:
            diff = okflubase2[(seas, flun)]-okflubase1[(seas, flun)]
            fig = ctl.plot_lat_crosssection(diff, lat, level, title = '{} {} - {} cross section {} vs {}'.format(shortnam[flun], seas, exp, y1, y2), plot_anomalies = True)
            figures.append(fig)

for y1 in [1975, 2100]:
    okflubase = allfluxes[('lcb0', y1, 'cross')]
    okflustoc = allfluxes[('lcs0', y1, 'cross')]
    for flun in ['mshf', 'mpef']:
        for seas in seasons:
            diff = okflustoc[(seas, flun)]-okflubase[(seas, flun)]
            fig = ctl.plot_lat_crosssection(diff, lat, level, title = '{} {} - lcb0 vs lcs0 cross section {}'.format(shortnam[flun],seas, y1), plot_anomalies = True)
            figures.append(fig)

ctl.plot_pdfpages(figure_file, figures)
