
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

if os.uname()[1] == 'hobbes':
    cart_in = '/data-hobbes/fabiano/ERA5_last/'
    cart_out = '/home/fabiano/Research/lavori/era5_viz/'
elif os.uname()[1] == 'xaru':
    cart_in = '/home/fedef/Research/lavori/era5_viz/'
    cart_out = '/home/fedef/Research/lavori/era5_viz/'

filregs = '/data-hobbes/fabiano/WRtool_output/ERA_ref_r25_v4/out_ERA_NDJFM_EAT_4clus_4pcs_1964-2014_dtr.p'
pino = ctl.load_wrtool(filregs)


# gigiuv = xr.load_dataset(cart_in + 'ERA5_uv850_NDJFM_2021_r25.nc')
# gigizg = xr.load_dataset(cart_in + 'ERA5_zg500_NDJFM_2021_r25.nc')
gigiuv = xr.load_dataset(cart_in + 'ERA5_uv850_NDJFM_2021.nc')
gigizg = xr.load_dataset(cart_in + 'ERA5_zg500_NDJFM_2021.nc')

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
cmappa = 'RdBu_r'

#fig, ax = ctl.get_cartopy_fig_ax()
def animate(i, ax, plot_margins = None):
    # i = iys[ii]
    proj = ctl.def_projection('nearside', central_lat_lon = (50., -20.))
    ax.clear()
    #fig, ax = ctl.plot_map_contour(zgok, lat, lon, visualization = 'nearside', central_lat_lon = (50., -20.), return_ax = True)
    map_plot = ctl.plot_mapc_on_ax(ax, zgok, lat, lon, proj, cmappa, cbar_range)
    year = anni[i]
    #print(year)
    color = cset[i]
    tam = tahiss[i] - pimean['tas']
    ax.set_title(r'{} $\to$ {:+5.1f} $^\circ$C wrt PI'.format(year, tam))
    return


lat = gigizg.latitude
lon = gigizg.longitude
fig, ax = ctl.plot_map_contour(zgok, lat, lon, visualization = 'nearside', central_lat_lon = (50., -20.), return_ax = True)# add_vector_field = [uok,vok])

lame, lome = np.meshgrid(lon, lat)
ax.streamplot(lame, lome, uok, vok, transform = ccrs.PlateCarree(), density = 3, linewidth = 0.5, integration_direction = 'backward')
# ax.quiver(lame, lome, uok, vok, transform = ccrs.PlateCarree(), regrid_shape = 10)
# ax.quiver(lame, lome, uok, vok, transform = ccrs.PlateCarree(), regrid_shape = 1)


from matplotlib.colors import LightSource
ls = LightSource(azdeg=315, altdeg=45)
zuku = ls.hillshade(zgok, vert_exag=1)

plt.figure()
plt.imshow(zuku, cmap='gray')

#fig, ax = ctl.get_cartopy_fig_ax(visualization='nearside', central_lat_lon=(50.,-20.))
#ax.imshow(zuku2, transform = ccrs.PlateCarree(), cmap = 'gray')

zgan = zgok - zgok.mean(axis = -1)[:, np.newaxis]
zuku3 = ls.hillshade(zgan, vert_exag=1)

fig, ax = ctl.get_cartopy_fig_ax(visualization='nearside', central_lat_lon=(50.,-20.))
ax.imshow(zuku3, transform = ccrs.PlateCarree(), cmap = 'gray')
plt.figure()
plt.imshow(zuku3, cmap = 'gray')

rgb = ls.shade(zgok, cm.get_cmap('viridis'), vert_exag=1, blend_mode='hsv')
plt.figure()
plt.imshow(rgb)

rgb2 = ls.shade(zgok, cm.get_cmap('viridis'), vert_exag=0.5, blend_mode='hsv')
plt.figure()
plt.imshow(rgb2)

rgb3 = ls.shade(zgok, cmap=cm.get_cmap('viridis'), vert_exag=1, blend_mode='overlay')
plt.figure()
plt.imshow(rgb3)
