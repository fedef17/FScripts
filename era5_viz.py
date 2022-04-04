
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

import matplotlib.animation as animation
from matplotlib.animation import ImageMagickFileWriter, PillowWriter

import cartopy.crs as ccrs

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 22
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 28
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.axisbelow'] = True

#############################################################################

# if os.uname()[1] == 'hobbes':
#     cart = '/home/fabiano/Research/lavori/BOTTINO/'
# elif os.uname()[1] == 'xaru':
#     cart = '/home/fedef/Research/lavori/BOTTINO/'
# elif os.uname()[1] == 'tintin':
#     cart = '/home/fabiano/work/lavori/BOTTINO/'
#
# cart_in = cart + 'seasmean/'
# cart_out = cart + 'visual_proj/'
# ctl.mkdir(cart_out)

cart = '/home/fedef/Research/lavori/era5_viz/'

# gigiuv = xr.load_dataset(cart + 'ERA5_uv850_NDJFM_2021_r25.nc')
# gigizg = xr.load_dataset(cart + 'ERA5_zg500_NDJFM_2021_r25.nc')
gigiuv = xr.load_dataset(cart + 'ERA5_uv850_NDJFM_2021.nc')
gigizg = xr.load_dataset(cart + 'ERA5_zg500_NDJFM_2021.nc')

zgok = np.array(gigizg['z'].sel(expver = 1.))[0]
uok = np.array(gigiuv['u'].sel(expver = 1.))[0]
vok = np.array(gigiuv['v'].sel(expver = 1.))[0]

# ok1 = np.where(~np.all(np.isnan(np.array(gigiuv['u'].sel(expver = 1.))), axis = (1,2)))[0]
# ok5 = np.where(~np.all(np.isnan(np.array(gigiuv['u'].sel(expver = 5.))), axis = (1,2)))[0]
#
# zg = np.concatenate([np.array(gigizg['z'].sel(expver = 1.))[ok1], np.array(gigizg['z'].sel(expver = 5.))[ok5]])
# u = np.concatenate([np.array(gigiuv['u'].sel(expver = 1.))[ok1], np.array(gigiuv['u'].sel(expver = 5.))[ok5]])
# v = np.concatenate([np.array(gigiuv['v'].sel(expver = 1.))[ok1], np.array(gigiuv['v'].sel(expver = 5.))[ok5]])
# zgok = zg[0]
# uok = u[0]
# vok = v[0]

#fig, ax = ctl.get_cartopy_fig_ax()

lat = gigizg.latitude
lon = gigizg.longitude
fig, ax = ctl.plot_map_contour(zgok, lat, lon, visualization = 'nearside', central_lat_lon = (70., 0.), return_ax = True)# add_vector_field = [uok,vok])

lame, lome = np.meshgrid(lon, lat)
ax.streamplot(lame, lome, uok, vok, transform = ccrs.PlateCarree())
# ax.quiver(lame, lome, uok, vok, transform = ccrs.PlateCarree(), regrid_shape = 10)
# ax.quiver(lame, lome, uok, vok, transform = ccrs.PlateCarree(), regrid_shape = 1)
