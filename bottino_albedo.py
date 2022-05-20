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

import pygrib
import cartopy.crs as ccrs

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

cart_out = cart_out + 'albedo/'
ctl.mkdir(cart_out)

# allru = ['b050', 'b100']
#
# for ru in allru:
#     filz = glob.glob(cart_out + 'alb_{}_*nc'.format(ru))
#     filz.sort()
#
#     fields_sep = []
#     fields_mar = []
#     yeas = []
#     for fi in filz:
#         pino = xr.load_dataset(fi, use_cftime = True)
#         #ctl.plot_map_contour(pino['rsuscs'])
#         fields_sep.append(pino['rsuscs'][8])
#         fields_mar.append(pino['rsuscs'][2])
#         yeas.append(str(pino.time.values[0].year))
#
#     ctl.plot_multimap_contour(fields_sep, filename = cart_out + 'greenland_albedo_sep_{}.pdf'.format(ru), cbar_range=[0,1], plot_anomalies=False, plot_margins = (-80, -10, 60, 85), subtitles = yeas, cb_label = 'Albedo (sept)', figsize = (16,9))
#
#     ctl.plot_multimap_contour(fields_mar, filename = cart_out + 'antarctic_albedo_mar_{}.pdf'.format(ru), cbar_range=[0,1], plot_anomalies=False, visualization = 'nearside', bounding_lat = -45, central_lat_lon = (-90, 0), subtitles = yeas, cb_label = 'Albedo (march)')

cart_out = cart_out + 'new_alb/'
ctl.mkdir(cart_out)

greenland_water_volume_m = 0.708*4*np.pi*(ctl.Rearth)**2*7.2
greenland_area = 2166086000000# m2
greenland_water_column_m = greenland_water_volume_m/greenland_area

antartic_water_volume = 26.5e15
print(antartic_water_volume/greenland_water_volume_m)
antartic_water_column = antartic_water_volume/14.2e12

### Loading PI ice cover, ICMGG and albedo
# surf/module/sussoil_mod.F90 contains RSNPER = 10000. : SNOW MASS per unit area FOR PERENNIAL SNOW

# Reading pre-industrial snow cover
zucu = pygrib.open(cart_out + 'ICMGGECE3INIT_ice_t255')
for grb in zucu:
    grb.expand_grid(False)
    if grb.name == 'Snow depth':
        print(grb.name)
        sd = grb.values
        latg, long = grb.latlons()

zucu.rewind()
lice = sd == 10.

# Reading initial snow cover for b100
b100_ini = pygrib.open(cart_out + '../c585-21000101/ICMGGc585INIT')

fig, ax = ctl.get_cartopy_fig_ax()
plt.scatter(long, latg, transform = ccrs.PlateCarree(), s = 1)

b100_ini.rewind()
for grb in b100_ini:
    grb.expand_grid(False)
    if grb.name == 'Snow depth':
        print(grb.name)
        sd100 = grb.values

# masking non-greenland glaciers
lice[(long < 180) & (latg > 0)] = False
lice[(long < 290) & (latg > 0)] = False
lice[(long < 300) & (latg > 0) & (latg < 73)] = False
lice[(long > 334) & (latg > 0) & (latg < 66.5)] = False

plt.scatter(long[lice], latg[lice], transform = ccrs.PlateCarree(), s = 1, color = 'orange')

sd100_mod_unif = sd100.copy()
sd100_mod_unif[lice] = 100. * sd100[lice]

sd100_mod_unif1000 = sd100.copy()
sd100_mod_unif1000[lice] = 1000.

fig, ax = ctl.get_cartopy_fig_ax()
plt.scatter(long, latg, c = sd100_mod_unif, transform = ccrs.PlateCarree(), s = 2)
plt.colorbar()

b100_ini.rewind()
allstr = []
for grb in b100_ini:
    grb.expand_grid(False)
    if grb.name == 'Snow depth':
        grb.values = sd100_mod_unif
    allstr.append(grb.tostring())

grbout = open(cart_out + 'new_b100_init_landicex100.grb','wb')
for msg in allstr:
    grbout.write(msg)
grbout.close()


b100_ini.rewind()
allstr = []
for grb in b100_ini:
    grb.expand_grid(False)
    if grb.name == 'Snow depth':
        grb.values = sd100_mod_unif1000
    allstr.append(grb.tostring())

grbout = open(cart_out + 'new_b100_init_landice1000.grb','wb')
for msg in allstr:
    grbout.write(msg)
grbout.close()

##########################################################
# Modifying bare soil albedo
bsa = pygrib.open(cart_out + 'bare_soil_albedos.grb')

allstr2 = []
for ii, grb in enumerate(bsa):
    grb.expand_grid(False)
    zuki = grb.values
    if ii < 2:
        zuki[lice] = 0.5
    else:
        zuki[lice] = 0.6
    grb.values = zuki
    allstr2.append(grb.tostring())

grbout = open(cart_out + 'new_bare_soil_albedo_landice05.grb','wb')
for msg in allstr2:
    grbout.write(msg)
grbout.close()
