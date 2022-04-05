
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

from matplotlib.colors import LightSource

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
    cartreg = '/home/fabiano/Research/lavori/WeatherRegimes/ERA_ref/'
elif os.uname()[1] == 'xaru':
    cart_in = '/home/fedef/Research/lavori/era5_viz/'
    cart_out = '/home/fedef/Research/lavori/era5_viz/'
    cartreg = '/home/fedef/Research/lavori/WRtool_tests/ERA_ref/'

# filregs = '/data-hobbes/fabiano/WRtool_output/ERA_ref_r25_v4/out_ERA_NDJFM_EAT_4clus_4pcs_1964-2014_dtr.p'
# pino = ctl.load_wrtool(filregs)

fireg = 'out_ERA_NDJFM_EAT_4clus_4pcs_1979-2019.p'
gigi = ctl.load_wrtool(cartreg + fireg)
era5_all = gigi#[0]['ERA5']

cents = era5_all['centroids']
#labs = era5_all['labels']
#pcs = era5_all['pcs']


####
# fizg = cart_in + 'ERA5_zg500_NDJFM_2021_r25.nc'
# fiu = cart_in + 'ERA5_uv850_NDJFM_2021_r25.nc'
# gigizg, coords, aux_ = ctl.read_xr(fizg)
# gigiuv, coords, aux_ = ctl.read_xr(fiu)

gigiuv = xr.load_dataset(cart_in + 'ERA5_uv850_NDJFM_2021_r25.nc', use_cftime = True)
gigizg = xr.load_dataset(cart_in + 'ERA5_zg500_NDJFM_2021_r25.nc', use_cftime = True)
# gigiuv = xr.load_dataset(cart_in + 'ERA5_uv850_NDJFM_2021.nc', use_cftime = True)
# gigizg = xr.load_dataset(cart_in + 'ERA5_zg500_NDJFM_2021.nc', use_cftime = True)

# zgok = np.array(gigizg['z'].sel(expver = 1.))[0]
# uok = np.array(gigiuv['u'].sel(expver = 1.))[0]
# vok = np.array(gigiuv['v'].sel(expver = 1.))[0]

ok1 = np.where(~np.all(np.isnan(np.array(gigiuv['u'].sel(expver = 1.))), axis = (1,2)))[0]
ok5 = np.where(~np.all(np.isnan(np.array(gigiuv['u'].sel(expver = 5.))), axis = (1,2)))[0]

zg = np.concatenate([np.array(gigizg['z'].sel(expver = 1.))[ok1], np.array(gigizg['z'].sel(expver = 5.))[ok5]])[90:]
u = np.concatenate([np.array(gigiuv['u'].sel(expver = 1.))[ok1], np.array(gigiuv['u'].sel(expver = 5.))[ok5]])[90:]
v = np.concatenate([np.array(gigiuv['v'].sel(expver = 1.))[ok1], np.array(gigiuv['v'].sel(expver = 5.))[ok5]])[90:]

zg = zg/9.80655 # to meters

# zg = zg.swapaxes(1,2)
# u = u.swapaxes(1,2)
# v = v.swapaxes(1,2)
dates = gigizg['time'].values[90:]
zgan = ctl.anomalies_daily(zg, dates, era5_all['climate_mean'], era5_all['climate_mean_dates'])

lat = gigizg.lat
lon = gigizg.lon

lastwin = cd.WRtool_core(zgan, lat, lon, dates, 'EAT', numpcs = 4, numclus = 4, remove_climate_mean = False, use_reference_eofs = True, use_reference_clusters = True, ref_solver = era5_all['solver'], ref_clusters_centers = era5_all['centroids'])

labs = lastwin['labels']
pcs = lastwin['pcs']

zgok = zgan[0]
uok = u[0]
vok = v[0]

lame, lome = np.meshgrid(lon, lat)
# TEST
# lat = gigizg.latitude
# lon = gigizg.longitude
# fig, ax = ctl.plot_map_contour(zgok, lat, lon, visualization = 'nearside', central_lat_lon = (50., -20.), return_ax = True)# add_vector_field = [uok,vok])
#
# #ax.streamplot(lame, lome, uok, vok, transform = ccrs.PlateCarree(), density = 3, linewidth = 0.5, integration_direction = 'backward')
# ax.quiver(lame, lome, uok, vok, transform = ccrs.PlateCarree(), regrid_shape = 50)
# # ax.quiver(lame, lome, uok, vok, transform = ccrs.PlateCarree(), regrid_shape = 1)

#### ANIMATION
cmappa = 'RdBu_r'
cbar_range = (-300, 300)

proj = ctl.def_projection('nearside', central_lat_lon = (50., -20.))
fig = plt.figure(figsize = (16,9))
ax1 = plt.subplot(1, 2, 1, projection = proj)
ax1.set_global()
ax1.coastlines(linewidth = 1)
gl = ax1.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, linewidth = 1, color = 'gray', alpha = 0.5, linestyle = ':')

ax2 = plt.subplot(1, 2, 2, projection='3d')
colorz = ctl.color_set(4)
colorz = ['steelblue', 'orange', 'indianred', 'forestgreen']

for iii, col in enumerate(colorz):
    ax2.scatter(cents[iii][0], cents[iii][1], cents[iii][2], color = col, s = 50)

pcs_low = ctl.butter_filter(pcs, 5)

from scipy.interpolate import PchipInterpolator as spline
zuku = np.linspace(0, len(pcs), 10*len(pcs))
pcspl = spline(np.arange(len(pcs)), pcs_low)
pcs_low_hr = pcspl(zuku)

ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_zticks([])
ax2.set_xlim((-3000.,2500.))
ax2.set_ylim((-2000.,2000.))
ax2.set_zlim((-2000.,2000.))


# create light source object.
ls = LightSource(azdeg=90, altdeg=65)
# Set of all spherical angles:
alu = np.linspace(0, 2 * np.pi, 100)
alv = np.linspace(0, np.pi, 100)

# Cartesian coordinates that correspond to the spherical angles:
# (this is the equation of an ellipsoid):
surf = None
linecoso = None

ellips = []
for i, mpz in zip(range(4), [cm.Blues_r, cm.Oranges_r, cm.Purples_r, cm.Greens_r]):
    rx, ry, rz = np.std(pcs[labs == i, :3], axis = 0)
    cx, cy, cz = cents[i][:3]
    x = cx + rx * np.outer(np.cos(alu), np.sin(alv))
    y = cy + ry * np.outer(np.sin(alu), np.sin(alv))
    z = cz + rz * np.outer(np.ones_like(alu), np.cos(alv))
    ellips.append([x,y,z])

    rgb = ls.shade(z, mpz)
    surf = ax2.plot_surface(x, y, z, rstride=1, cstride=1, linewidth=0, antialiased=False, facecolors=rgb, alpha = 0.3)

ax2.view_init(30, 60)
plt.draw()


texto = fig.text(0.45, 0.05, str(dates[0])[:10], fontsize = 20)


def plot_ellipsoids(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.set_xlim((-3000.,2500.))
    ax.set_ylim((-2000.,2000.))
    ax.set_zlim((-2000.,2000.))

    for i, mpz in zip(range(4), [cm.Blues_r, cm.Oranges_r, cm.Purples_r, cm.Greens_r]):
        x,y,z = ellips[i]
        rgb = ls.shade(z, mpz)
        surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, linewidth=0, antialiased=False, facecolors=rgb, alpha = 0.3)

    ax.view_init(30, 60)
    plt.draw()
    return


def animate(i, ax1, ax2, windtype = 'quiver'):
    print(i)
    ax1.clear()
    ax2.clear()

    zgok = zgan[i]
    uok = u[i]
    vok = v[i]

    map_plot = ctl.plot_mapc_on_ax(ax1, zgok, lat, lon, proj, cmappa, cbar_range)

    if windtype == 'quiver':
        ax1.quiver(lame, lome, uok, vok, transform = ccrs.PlateCarree(), regrid_shape = 50)
    elif windtype == 'streamplot':
        ax1.streamplot(lame, lome, uok, vok, transform = ccrs.PlateCarree(), density = 3, linewidth = 0.5)

    # for iii, col in enumerate(colorz):
    #     ax2.scatter(cents[iii][0], cents[iii][1], cents[iii][2], color = col, s = 50)
    plot_ellipsoids(ax2)

    col = colorz[labs[i]]

    cx, cy, cz = pcs_low[i, :3]
    #ax2.scatter(pcs_low[i, 0], pcs_low[i, 1], pcs_low[i, 2], color = col, s = 30)
    #if surf is not None: surf.remove()
    R = 100.
    x = cx + R * np.outer(np.cos(alu), np.sin(alv))
    y = cy + R * np.outer(np.sin(alu), np.sin(alv))
    z = cz + R * np.outer(np.ones_like(alu), np.cos(alv))
    rgb = ls.shade(z, cm.Greys_r)
    surf = ax2.plot_surface(x, y, z, rstride=1, cstride=1, linewidth=0, antialiased=False, facecolors=rgb)

    # for gu in range(1, 5):
    #     ax2.scatter(pcs_low[i-gu, 0], pcs_low[i-gu, 1], pcs_low[i-gu, 2], color = colorz[labs[i-gu]], s = 5)

    #if linecoso is not None: linecoso.remove()
    linecoso = ax2.plot(pcs_low_hr[10*i-100:10*i+1, 0], pcs_low_hr[10*i-100:10*i+1, 1], pcs_low_hr[10*i-100:10*i+1, 2], color = 'gray', alpha = 0.4)
    # ax2.plot(pcs_low_hr[:10*i-100, 0], pcs_low_hr[:10*i-100, 1], pcs_low_hr[:10*i-100, 2], color = 'gray', alpha = 0.1, lw = 0.1)
    #ax2.plot(pcs_low[i-10:i+1, 0], pcs_low[i-10:i+1, 1], pcs_low[i-10:i+1, 2], color = 'gray', alpha = 0.4)

    texto.set_text(str(dates[i])[:10])

    return

windtype = 'quiver'
fps = 10
numframes = 20 # len(zgan)
line_ani = animation.FuncAnimation(fig, animate, frames = numframes, fargs = (ax1, ax2, windtype, ), interval=100, blit=False)
filename = cart_out + 'wr_anim_nearside_21-22_{}.gif'.format(windtype)
writer = PillowWriter(fps = fps)
line_ani.save(filename, writer = writer)

# windtype = 'streamplot'
# line_ani = animation.FuncAnimation(fig, animate, frames = len(zgan), fargs = (ax1, ax2, windtype, ), interval=100, blit=False)
# filename = cart_out + 'wr_anim_nearside_21-22_{}.gif'.format(windtype)
# writer = PillowWriter(fps = fps)
# line_ani.save(filename, writer = writer)

sys.exit()

def animate_hs(i, ax1, ax2, windtype = 'quiver'):
    ax1.clear()
    ax2.clear()

    zgok = zgan[i]
    uok = u[i]
    vok = v[i]

    map_plot = ctl.plot_mapc_on_ax(ax1, zgok, lat, lon, proj, cmappa, cbar_range)

    if windtype == 'quiver':
        ax1.quiver(lame, lome, uok, vok, transform = ccrs.PlateCarree(), regrid_shape = 50)
    elif windtype == 'streamplot':
        ax1.streamplot(lame, lome, uok, vok, transform = ccrs.PlateCarree(), density = 3, linewidth = 0.5)

    # year = anni[i]
    # #print(year)
    # color = cset[i]
    # tam = tahiss[i] - pimean['tas']
    # ax.set_title(r'{} $\to$ {:+5.1f} $^\circ$C wrt PI'.format(year, tam))
    return



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
