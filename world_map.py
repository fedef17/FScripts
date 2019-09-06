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
import itertools as itt

from sklearn.cluster import KMeans

from datetime import datetime
import pickle

import climtools_lib as ctl

import cartopy.crs as ccrs
import cartopy.util as cartopy_util
import cartopy.feature as cfeature
from matplotlib.colors import LightSource

import seaborn as sns

##############################################
############ INPUTS #########

elev_file = '/home/fabiano/Research/temp/elev.0.25-deg.nc'

elev, lat, lon, var_units = ctl.readxDncfield(elev_file)

elev = elev.squeeze()
# map = ctl.plot_map_contour(elev, lat, lon, visualization = 'Robinson')

### Trying hillshades

proj = ccrs.Robinson(central_longitude = 0.)

# Plotting figure
fig = plt.figure(figsize = (32, 24))
ax = plt.subplot(projection = proj)

cmappa = cm.get_cmap('gist_earth')
cbar_range = None
clevels = [-10000, -7000, -5000, -4000, -3000, -2000, -1000, -500, 0, 500, 1000, 1500, 2000, 3000, 4000, 6000, 7000]
ax.set_global()
ax.coastlines(linewidth = 2)

# Shade from the northwest, with the sun 45 degrees from horizontal
ls = LightSource(azdeg=315, altdeg=45)

hillsh = ls.hillshade(elev, vert_exag = 100, dx = 10000, dy = 10000)
ax.pcolormesh(lon, lat, hillsh, cmap = 'gray', transform = ccrs.PlateCarree(), alpha = 1.0)

map = ctl.plot_mapc_on_ax(ax, elev, lat, lon, proj, cmappa, cbar_range, clevels = clevels, draw_grid = True, alphamap = 0.7)
#cb = plt.colorbar(map, orientation='horizontal')

# #### N/S Polar
proj = ctl.def_projection('nearside', (88, 0), bounding_lat = 0)
fig2 = plt.figure(figsize = (32, 24))
ax2 = plt.subplot(projection = proj)
ax2.set_global()
ax2.coastlines(linewidth = 2)

ax2.pcolormesh(lon, lat, hillsh, cmap = 'gray', transform = ccrs.PlateCarree(), alpha = 1.0)
map = ctl.plot_mapc_on_ax(ax2, elev, lat, lon, proj, cmappa, cbar_range, clevels = clevels, draw_grid = True, alphamap = 0.7)


proj = ctl.def_projection('nearside', (-88, 0), bounding_lat = 0)
fig3 = plt.figure(figsize = (32, 24))
ax3 = plt.subplot(projection = proj)
ax3.set_global()
ax3.coastlines(linewidth = 2)

ax3.pcolormesh(lon, lat, hillsh, cmap = 'gray', transform = ccrs.PlateCarree(), alpha = 1.0)
map = ctl.plot_mapc_on_ax(ax3, elev, lat, lon, proj, cmappa, cbar_range, clevels = clevels, draw_grid = True, alphamap = 0.7)

# #### Using imshow
# fig2 = plt.figure(figsize = (32, 24))
# ax2 = plt.subplot(projection = ccrs.PlateCarree())
# ax2.set_global()
# ax2.coastlines(linewidth = 2)
#
# rgb = ls.shade(elev, cmap=cmappa, blend_mode='soft', vert_exag=10, dx = 10000, dy = 10000)
#
# extent = [lon[0], lon[-1], lat[0], lat[-1]]
# ax2.imshow(rgb, origin = 'lower', extent = extent, transform = ccrs.PlateCarree())
