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
print('TEST run on February 1988\n\n')

cart_in = '/data-hobbes/fabiano/SPHINX/heat_flux/1988_daily/'

cart_out = '/home/fabiano/Research/lavori/SPHINX_for_lisboa/heat_flux/'
figure_file = cart_out+'heat_flux_calc_vs_ERAref.pdf'

figure_file2 = cart_out+'hfc_ERAref_recalc.pdf'
figures = []

figure_file_cross = cart_out+'hfc_crosssections_ERAref.pdf'
figures_cross = []

fluxes = dict()
fluxnames = ['mshf', 'mpef', 'mlhf']
fluxlongnames = ['Sensible Heat', 'Potential Energy', 'Latent Heat']
factors = [cp/g, 1., L/g]

vars = dict()
varnames = ['hus', 'ta', 'va', 'zg']
fils = ['lcs0_day_1988_{}.nc'.format(varna) for varna in varnames]

for varna, fi in zip(varnames, fils):
    var, level, lat, lon, dates, time_units, var_units, time_cal = ctl.read4Dncfield(cart_in+fi)
    var, okda = ctl.sel_season(var, dates, 'Feb')
    vars[varna] = var

press0, latdad, londsad, datespress, time_units, var_units = ctl.read3Dncfield(cart_in+'lcs0_day_1988_ps.nc')
press0, _ = ctl.sel_season(press0, datespress, 'Feb')
press0 = np.mean(press0, axis = 0)

mshf = factors[0]*vars['va']*vars['ta']
mpef = vars['va']*vars['zg']
mlhf = factors[2]*vars['va']*vars['hus']

mshf = np.mean(mshf, axis = 0)
mpef = np.mean(mpef, axis = 0)
mlhf = np.mean(mlhf, axis = 0)

mshfint = np.zeros(mshf.shape[1:])
mpefint = np.zeros(mpef.shape[1:])
mlhfint = np.zeros(mlhf.shape[1:])
print(mshf.shape)

# for ila in range(len(lat)):
#     for ilo in range(len(lon)):
#         print(lat[ila], lon[ilo])
#         p0 = press0[ila, ilo]
#         exclude_pres = level > p0
#         if np.any(exclude_pres):
#             lev0 = np.argmax(exclude_pres)
#         else:
#             lev0 = None
#         levels_ok = np.append(level[:lev0], p0)
#
#         if not np.all(np.diff(levels_ok) > 0):
#             print(levels_ok)
#             raise ValueError('x not increasing')
#
#         coso = mshf[:, ila, ilo]
#         coso = np.append(coso[:lev0], np.interp(p0, level, coso))
#         mshfint[ila, ilo] = np.trapz(coso, x = levels_ok)
#         print(levels_ok, coso.mean(), coso[-1], mshfint[ila, ilo])
#
#         coso = mpef[:, ila, ilo]
#         coso = np.append(coso[:lev0], np.interp(p0, level, coso))
#         mpefint[ila, ilo] = np.trapz(coso, x = levels_ok)
#
#         coso = mlhf[:, ila, ilo]
#         coso = np.append(coso[:lev0], np.interp(p0, level, coso))
#         mlhfint[ila, ilo] = np.trapz(coso, x = levels_ok)
#
# ctl.save2Dncfield(lat,lon,mshfint,'mshf',cart_out+'mshf.nc')
# ctl.save2Dncfield(lat,lon,mpefint,'mpef',cart_out+'mpef.nc')
# ctl.save2Dncfield(lat,lon,mlhfint,'mlhf',cart_out+'mlhf.nc')

zonal_factor = 2*np.pi*Rearth*np.cos(np.deg2rad(lat))
mshfcross = np.mean(mshf, axis = -1)*zonal_factor
mpefcross = np.mean(mpef, axis = -1)*zonal_factor
mlhfcross = np.mean(mlhf, axis = -1)*zonal_factor

figures_cross.append(ctl.plot_lat_crosssection(mshfcross, lat, level, title = 'SPHINX - SH cross section', cbar_range = (-3.e12, 3.e12), plot_anomalies = True))
figures_cross.append(ctl.plot_lat_crosssection(mpefcross, lat, level, title = 'SPHINX - PE cross section', cbar_range = (-1.e12, 1.e12), plot_anomalies = True))
figures_cross.append(ctl.plot_lat_crosssection(mpefcross, lat, level, title = 'SPHINX - PE cross section', cbar_range = (-1.e12, 1.e12), plot_anomalies = True, set_logscale_levels = True))
figures_cross.append(ctl.plot_lat_crosssection(mlhfcross, lat, level, title = 'SPHINX - LH cross section', cbar_range = (-1.e11, 1.e11), plot_anomalies = True))

# mshfzon = np.mean(mshfint, axis = 1)*zonal_factor
# mpefzon = np.mean(mpefint, axis = 1)*zonal_factor
# mlhfzon = np.mean(mlhfint, axis = 1)*zonal_factor
#
# fig = plt.figure()
# plt.title('Mycoso')
# plt.plot(lat, mshfzon, label = 'SH')
# plt.plot(lat, mpefzon, label = 'PE')
# plt.plot(lat, mlhfzon, label = 'LH')
# plt.plot(lat, mshfzon+mpefzon+mlhfzon, label = 'Total')
# plt.legend()
# figures.append(fig)

pinopress = itrp.RectBivariateSpline(lat, lon, press0)

cart_in = '/data-hobbes/fabiano/OBS/ERA/ERAInterim/'

#leggo gia fatto
fi = 'prova_heatflux_1988_MM.nc'
cpc = nc.Dataset(cart_in+fi)
era_lat = cpc.variables['latitude'][:]
era_lon = cpc.variables['longitude'][:]
era_zonal_factor = 2*np.pi*Rearth*np.cos(np.deg2rad(era_lat))
era_tot_Feb = list(cpc.variables.values())[-1][1, ...]
era_mpef_Feb = list(cpc.variables.values())[-2][1, ...]
era_mshf_Feb = list(cpc.variables.values())[-3][1, ...]

fi = 'northward_water_flux.nc'
cpc = nc.Dataset(cart_in+fi)
era_mlhf_Feb = L*list(cpc.variables.values())[-1][1, ...]

era_tot_Febzon = np.mean(era_tot_Feb, axis = 1)*era_zonal_factor
era_mshf_Febzon = np.mean(era_mshf_Feb, axis = 1)*era_zonal_factor
era_mpef_Febzon = np.mean(era_mpef_Feb, axis = 1)*era_zonal_factor
era_mlhf_Febzon = np.mean(era_mlhf_Feb, axis = 1)*era_zonal_factor

# fig = plt.figure()
# plt.title('ERA ref - February 1988')
# plt.plot(era_lat, era_mshf_Febzon, label = 'SH')
# plt.plot(era_lat, era_mpef_Febzon, label = 'PE')
# plt.plot(era_lat, era_mlhf_Febzon, label = 'LH')
# plt.plot(era_lat, era_mlhf_Febzon+era_mshf_Febzon+era_mpef_Febzon, label = 'Sum')
# plt.plot(era_lat, era_tot_Febzon, label = 'Total')
# plt.legend()
# figures.append(fig)

# calcolo il mio, senza correzione per la surface pressure
cpc = nc.Dataset(cart_in+'all_vtgq_1988_daily.nc')
era_levels = 100.0*cpc.variables['level'][:]
era_zonal_factor = 2*np.pi*Rearth*np.cos(np.deg2rad(era_lat))
v = list(cpc.variables.values())[-1][:]
q = list(cpc.variables.values())[-2][:]
t = list(cpc.variables.values())[-3][:]
z = list(cpc.variables.values())[-4][:]/g

mshf = factors[0]*v*t
mpef = v*z
mlhf = factors[2]*v*q
mshf = np.mean(mshf, axis = 0)
mpef = np.mean(mpef, axis = 0)
mlhf = np.mean(mlhf, axis = 0)

mshfint1 = np.trapz(mshf, x = era_levels, axis = 0).squeeze()
mpefint1 = np.trapz(mpef, x = era_levels, axis = 0).squeeze()
mlhfint1 = np.trapz(mlhf, x = era_levels, axis = 0).squeeze()

# ctl.save2Dncfield(lat,lon,mshfint,'mshf',cart_out+'era_mshf_Feb.nc')
# ctl.save2Dncfield(lat,lon,mpefint,'mpef',cart_out+'era_mpef_Feb.nc')
# ctl.save2Dncfield(lat,lon,mlhfint,'mlhf',cart_out+'era_mlhf_Feb.nc')

mshfzon1 = np.mean(mshfint1, axis = 1)*era_zonal_factor
mpefzon1 = np.mean(mpefint1, axis = 1)*era_zonal_factor
mlhfzon1 = np.mean(mlhfint1, axis = 1)*era_zonal_factor

# plt.figure()
# plt.title('ERA calc - February 1988 - daily')
# plt.plot(era_lat, mshfzon1, label = 'SH')
# plt.plot(era_lat, mpefzon1, label = 'PE')
# plt.plot(era_lat, mlhfzon1, label = 'LH')
# plt.plot(era_lat, mshfzon1+mpefzon1+mlhfzon1, label = 'Total')
# plt.legend()


# calcolo il mio, senza correzione per la surface pressure
cpc = nc.Dataset(cart_in+'all_vtgq_1988_6hrs.nc')
era_levels = 100.0*cpc.variables['level'][:]
era_zonal_factor = 2*np.pi*Rearth*np.cos(np.deg2rad(era_lat))
v = list(cpc.variables.values())[-1][:]
q = list(cpc.variables.values())[-2][:]
t = list(cpc.variables.values())[-3][:]
z = list(cpc.variables.values())[-4][:]/g

mshf00 = factors[0]*v*t
mpef00 = v*z
mlhf00 = factors[2]*v*q
mshf = np.mean(mshf00, axis = 0)
mpef = np.mean(mpef00, axis = 0)
mlhf = np.mean(mlhf00, axis = 0)

mshfint2 = np.trapz(mshf, x = era_levels, axis = 0).squeeze()
mpefint2 = np.trapz(mpef, x = era_levels, axis = 0).squeeze()
mlhfint2 = np.trapz(mlhf, x = era_levels, axis = 0).squeeze()

# ctl.save2Dncfield(lat,lon,mshfint,'mshf',cart_out+'era_mshf_Feb.nc')
# ctl.save2Dncfield(lat,lon,mpefint,'mpef',cart_out+'era_mpef_Feb.nc')
# ctl.save2Dncfield(lat,lon,mlhfint,'mlhf',cart_out+'era_mlhf_Feb.nc')

mshfzon2 = np.mean(mshfint2, axis = 1)*era_zonal_factor
mpefzon2 = np.mean(mpefint2, axis = 1)*era_zonal_factor
mlhfzon2 = np.mean(mlhfint2, axis = 1)*era_zonal_factor


mshfcross_era = np.mean(mshf, axis = -1)*era_zonal_factor
mpefcross_era = np.mean(mpef, axis = -1)*era_zonal_factor
mlhfcross_era = np.mean(mlhf, axis = -1)*era_zonal_factor

figures_cross.append(ctl.plot_lat_crosssection(mshfcross_era, era_lat, era_levels, title = 'Era - SH cross section', cbar_range = (-3.e12, 3.e12), plot_anomalies = True))
figures_cross.append(ctl.plot_lat_crosssection(mpefcross_era, era_lat, era_levels, title = 'Era - PE cross section', cbar_range = (-1.e12, 1.e12), plot_anomalies = True))
figures_cross.append(ctl.plot_lat_crosssection(mpefcross_era, era_lat, era_levels, title = 'Era - PE cross section', cbar_range = (-1.e12, 1.e12), plot_anomalies = True, set_logscale_levels = True))
figures_cross.append(ctl.plot_lat_crosssection(mlhfcross_era, era_lat, era_levels, title = 'Era - LH cross section', cbar_range = (-1.e11, 1.e11), plot_anomalies = True))

pinomshfcross = itrp.RectBivariateSpline(np.log(era_levels), era_lat[::-1], mshfcross_era[:,::-1])
pinompefcross = itrp.RectBivariateSpline(np.log(era_levels), era_lat[::-1], mpefcross_era[:,::-1])
pinomlhfcross = itrp.RectBivariateSpline(np.log(era_levels), era_lat[::-1], mlhfcross_era[:,::-1])

mshfcross_itrp = np.zeros(mshfcross.shape)
mpefcross_itrp = np.zeros(mshfcross.shape)
mlhfcross_itrp = np.zeros(mshfcross.shape)
for ila, la in enumerate(lat):
    for ile, le in enumerate(level):
        mshfcross_itrp[ile,ila] = pinomshfcross(np.log(le), la)
        mpefcross_itrp[ile,ila] = pinompefcross(np.log(le), la)
        mlhfcross_itrp[ile,ila] = pinomlhfcross(np.log(le), la)

figures_cross.append(ctl.plot_lat_crosssection(mshfcross-mshfcross_itrp, lat, level, title = 'SPHINX - Era diff - SH cross section', plot_anomalies = True))
figures_cross.append(ctl.plot_lat_crosssection(mpefcross-mpefcross_itrp, lat, level, title = 'SPHINX - Era - PE cross section', plot_anomalies = True))
figures_cross.append(ctl.plot_lat_crosssection(mpefcross-mpefcross_itrp, lat, level, title = 'SPHINX - Era - PE cross section', plot_anomalies = True, set_logscale_levels = True))
figures_cross.append(ctl.plot_lat_crosssection(mlhfcross-mlhfcross_itrp, lat, level, title = 'SPHINX - Era - LH cross section', plot_anomalies = True))

# Ora calcolo con press0
mshfint2p0 = np.zeros(mshf.shape[1:])
mpefint2p0 = np.zeros(mpef.shape[1:])
mlhfint2p0 = np.zeros(mlhf.shape[1:])

nupress = np.zeros((len(era_lat), len(era_lon)))

for ila in range(len(era_lat)):
    for ilo in range(len(era_lon)):
        print(era_lat[ila], era_lon[ilo])
        #p0 = press0[ila, ilo]
        p0 = pinopress(era_lat[ila], era_lon[ilo])
        nupress[ila,ilo] = p0
        exclude_pres = era_levels > p0
        if np.any(exclude_pres):
            lev0 = np.argmax(exclude_pres)
        else:
            lev0 = None
        levels_ok = np.append(era_levels[:lev0], p0)
        if not np.all(np.diff(levels_ok) > 0):
            print(levels_ok)
            raise ValueError('x not increasing')

        coso = mshf[:, ila, ilo]
        coso = np.append(coso[:lev0], np.interp(p0, era_levels, coso))
        mshfint2p0[ila, ilo] = np.trapz(coso, x = levels_ok)
        #print(levels_ok, coso.mean(), coso[-1], mshfint2p0[ila, ilo])

        coso = mpef[:, ila, ilo]
        coso = np.append(coso[:lev0], np.interp(p0, era_levels, coso))
        mpefint2p0[ila, ilo] = np.trapz(coso, x = levels_ok)

        coso = mlhf[:, ila, ilo]
        coso = np.append(coso[:lev0], np.interp(p0, era_levels, coso))
        mlhfint2p0[ila, ilo] = np.trapz(coso, x = levels_ok)

mshfzon2p0 = np.mean(mshfint2p0, axis = 1)*era_zonal_factor
mpefzon2p0 = np.mean(mpefint2p0, axis = 1)*era_zonal_factor
mlhfzon2p0 = np.mean(mlhfint2p0, axis = 1)*era_zonal_factor

nupress = np.stack([nupress, nupress])
mshfzon2p0_recalc = cd.quant_flux_calc(mshf00, era_lat, era_lon, era_levels, None, nupress, None, seasons = [])['zonal']['year']
mpefzon2p0_recalc = cd.quant_flux_calc(mpef00, era_lat, era_lon, era_levels, None, nupress, None, seasons = [])['zonal']['year']
mlhfzon2p0_recalc = cd.quant_flux_calc(mlhf00, era_lat, era_lon, era_levels, None, nupress, None, seasons = [])['zonal']['year']


# plt.figure()
# plt.title('ERA calc - February 1988 - 6 hrs')
# plt.plot(era_lat, mshfzon2, label = 'SH')
# plt.plot(era_lat, mpefzon2, label = 'PE')
# plt.plot(era_lat, mlhfzon2, label = 'LH')
# plt.plot(era_lat, mshfzon2+mpefzon2+mlhfzon2, label = 'Total')
# plt.legend()


# calcolo il mio, senza correzione per la surface pressure
cpc = nc.Dataset(cart_in+'all_vtgq_1988_daily_morelev.nc')
era_levels = 100.0*cpc.variables['level'][:]
era_zonal_factor = 2*np.pi*Rearth*np.cos(np.deg2rad(era_lat))
v = list(cpc.variables.values())[-1][:]
q = list(cpc.variables.values())[-2][:]
t = list(cpc.variables.values())[-3][:]
z = list(cpc.variables.values())[-4][:]/g

mshf = factors[0]*v*t
mpef = v*z
mlhf = factors[2]*v*q
mshf = np.mean(mshf, axis = 0)
mpef = np.mean(mpef, axis = 0)
mlhf = np.mean(mlhf, axis = 0)

mshfint3 = np.trapz(mshf, x = era_levels, axis = 0).squeeze()
mpefint3 = np.trapz(mpef, x = era_levels, axis = 0).squeeze()
mlhfint3 = np.trapz(mlhf, x = era_levels, axis = 0).squeeze()

# ctl.save2Dncfield(lat,lon,mshfint,'mshf',cart_out+'era_mshf_Feb.nc')
# ctl.save2Dncfield(lat,lon,mpefint,'mpef',cart_out+'era_mpef_Feb.nc')
# ctl.save2Dncfield(lat,lon,mlhfint,'mlhf',cart_out+'era_mlhf_Feb.nc')

mshfzon3 = np.mean(mshfint3, axis = 1)*era_zonal_factor
mpefzon3 = np.mean(mpefint3, axis = 1)*era_zonal_factor
mlhfzon3 = np.mean(mlhfint3, axis = 1)*era_zonal_factor

# plt.figure()
# plt.title('ERA calc - February 1988 - daily more levels')
# plt.plot(era_lat, mshfzon3, label = 'SH')
# plt.plot(era_lat, mpefzon3, label = 'PE')
# plt.plot(era_lat, mlhfzon3, label = 'LH')
# plt.plot(era_lat, mshfzon3+mpefzon3+mlhfzon3, label = 'Total')
# plt.legend()

fig = plt.figure()
plt.title('ERA calc - February 1988 - SH diffs')
#plt.plot(era_lat, mshfzon1, label = 'SH - daily')
plt.plot(era_lat, mshfzon2, label = 'SH - 6hrs')
plt.plot(era_lat, mshfzon2p0, label = 'SH - 6hrs - int_p0')
plt.plot(era_lat, mshfzon2p0_recalc, label = 'SH - hfc recalc')
#plt.plot(era_lat, mshfzon3, label = 'SH - daily more lev')
#plt.plot(lat, mshfzon, label = 'SH - sphinx daily')
plt.plot(era_lat, era_mshf_Febzon, label = 'SH era ref', linestyle = ':')
plt.legend()
figures.append(fig)


fig = plt.figure()
plt.title('ERA calc - February 1988 - PE diffs')
#plt.plot(era_lat, mpefzon1, label = 'PE - daily')
plt.plot(era_lat, mpefzon2, label = 'PE - 6hrs')
plt.plot(era_lat, mpefzon2p0, label = 'PE - 6hrs - int_p0')
plt.plot(era_lat, mpefzon2p0_recalc, label = 'PE - hfc recalc')
#plt.plot(era_lat, mpefzon3, label = 'PE - daily more lev')
#plt.plot(lat, mpefzon, label = 'PE - sphinx daily')
plt.plot(era_lat, era_mpef_Febzon, label = 'PE era ref', linestyle = ':')
plt.legend()
figures.append(fig)

fig = plt.figure()
plt.title('ERA calc - February 1988 - LH diffs')
#plt.plot(era_lat, mlhfzon1, label = 'LH - daily')
plt.plot(era_lat, mlhfzon2, label = 'LH - 6hrs')
plt.plot(era_lat, mlhfzon2p0, label = 'LH - 6hrs - int_p0')
plt.plot(era_lat, mlhfzon2p0_recalc, label = 'LH - hfc recalc')
#plt.plot(era_lat, mlhfzon3, label = 'LH - daily more lev')
#plt.plot(lat, mlhfzon, label = 'LH - sphinx daily')
plt.plot(era_lat, era_mlhf_Febzon, label = 'LH era ref', linestyle = ':')
plt.legend()
figures.append(fig)


# fig = ctl.plot_double_sidebyside(era_mshf_Feb, mshfint, [era_lat, lat], [era_lon, lon], title = 'era_ref vs sphinx', use_different_grids = True)
# figures.append(fig)
# fig = ctl.plot_double_sidebyside(era_mshf_Feb, mshfint2, era_lat, era_lon, title = 'era_ref vs era 6 hrs')
# figures.append(fig)
# fig = ctl.plot_double_sidebyside(mshfint2, mshfint2p0, era_lat, era_lon, title = 'era 6 hrs vs era 6 hrs int p0')
# figures.append(fig)

ctl.plot_pdfpages(figure_file2, figures)
figures_cross = np.array(figures_cross)
figures_cross = figures_cross[[0,4,8,1,5,9,2,6,10,3,7,11]]
ctl.plot_pdfpages(figure_file_cross, figures_cross)
