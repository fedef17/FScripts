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
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['axes.axisbelow'] = True

#############################################################################

if os.uname()[1] == 'hobbes':
    cart = '/home/fabiano/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'xaru':
    cart = '/home/fedef/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'tintin':
    cart = '/home/fabiano/work/lavori/BOTTINO/'

cart_in = cart + 'seasmean/'
cart_out = cart + 'visual_proj/'
ctl.mkdir(cart_out)

allru = ['pi', 'b025', 'b050', 'b100']
allnams = ['piControl', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']

colors = ['black', 'forestgreen', 'orange', 'violet']

allnams2 = allnams + ['ssp585', 'historical']
allru2 = allru + ['ssp585', 'hist']
colors2 = colors + ['indianred', 'steelblue']
####################################################################################################

### New colormap for temperature
col = cm.RdBu_r(np.linspace(0.,1,256))
col2 = cm.PuRd_r(np.linspace(0.,0.75,128))
col3 = cm.bone_r(np.linspace(0.25,0.75,128))
colors = np.concatenate([col3, col, col2])
from matplotlib import colors as mcolors
heatmap = mcolors.LinearSegmentedColormap.from_list('heat_strong', colors)

###
dpi = 150
save = False

glomeans, pimean, yeamean, mapmean = pickle.load(open(cart_in + 'bottino_seasmean_2D.p', 'rb'))

tahiss = np.concatenate([glomeans[('hist', 'tas')][1], glomeans[('ssp585', 'tas')][1]])
tahiss = ctl.butter_filter(tahiss, 5)

for var in ['tas', 'pr']:
    cosopi = yeamean[('pi', var)]
    pime = cosopi.mean('year')
    cosohist = yeamean[('hist', var)]
    cosossp = yeamean[('ssp585', var)]
    coso = xr.concat([cosohist, cosossp], dim = 'year')
    presme = coso.sel(year = slice(1960, 1990)).mean('year')


    if os.path.exists(cart_out + 'histssp_{}_lowpass_LR.p'.format(var)):
        cosolow, presme = pickle.load(open(cart_out + 'histssp_{}_lowpass_LR.p'.format(var), 'rb'))
    else:
        coso = ctl.regrid_dataset(coso, regrid_to_deg = 2.)
        presme = ctl.regrid_dataset(presme, regrid_to_deg = 2.)
        if var == 'tas':
            cosolow = ctl.butter_filter(coso, 5)
        elif var == 'pr':
            cosolow = ctl.butter_filter(coso, 10)
        pickle.dump([cosolow, presme], open(cart_out + 'histssp_{}_lowpass_LR.p'.format(var), 'wb'))

    anni = coso.year.values

    if var == 'tas':
        #cmappa = cm.get_cmap('RdBu_r')#'gist_ncar')
        cmappa = heatmap
        cosoanom = cosolow-presme.values[np.newaxis, ...]
        #cbar_range = ctl.get_cbar_range(cosoanom, symmetrical = True)
        cbar_range = (-8, 8)
    elif var == 'pr':
        cmappa = cm.get_cmap('BrBG')
        cosoanom = 100*(cosolow-presme.values[np.newaxis, ...])/presme.values[np.newaxis, ...]
        #cbar_range = ctl.get_cbar_range(cosoanom, symmetrical = True)
        cbar_range = (-30, 30)

    clevels = np.linspace(cbar_range[0], cbar_range[1], 21)
    cset = ctl.color_set(len(anni), bright_thres = 0., full_cb_range = True)


    def animate(ii, ax):
        i = iys[ii]
        proj = ccrs.PlateCarree()
        ax.clear()
        map_plot = ctl.plot_mapc_on_ax(ax, cosoanom[i, ...], presme.lat, presme.lon, proj, cmappa, cbar_range, draw_grid = True)
        year = anni[i]
        print(year)
        color = cset[i]
        tam = tahiss[i] - pimean['tas']
        ax.set_title(r'{} $\to$ {:+5.1f} $\circ$C wrt PI'.format(year, tam))
        return


    fig, ax = ctl.get_cartopy_fig_ax(visualization = 'standard', central_lat_lon = (0, 0), bounding_lat = None, figsize = (9, 5), coast_lw = 1)
    #tit = plt.title('1850')

    # Plotting figure
    proj = ccrs.PlateCarree()
    map_plot = ctl.plot_mapc_on_ax(ax, cosoanom[0], presme.lat, presme.lon, proj, cmappa, cbar_range)

    cax = plt.axes([0.1, 0.11, 0.8, 0.05]) #horizontal
    cb = plt.colorbar(map_plot, cax=cax, orientation='horizontal')#, labelsize=18)
    cb.ax.tick_params(labelsize=14)
    if var == 'tas':
        cb.set_label('Temperature anomaly wrt 1960-1990 (K)', fontsize=14)
    elif var == 'pr':
        cb.set_label('Relative precipitation anomaly wrt 1960-1990 (%)', fontsize=14)

    top    = 0.88  # the top of the subplots
    bottom = 0.20    # the bottom of the subplots
    left   = 0.02    # the left side
    right  = 0.98  # the right side
    hspace = 0.20   # height reserved for white space
    wspace = 0.05    # width reserved for blank space
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

    #showdate = ax.text(0.5, 0.95, '1850', fontweight = 'bold', color = cset[0], bbox=dict(facecolor='lightsteelblue', edgecolor='black', boxstyle='round,pad=1'))

    if save:
        metadata = dict(title='Temperature anomaly (EC-Earth CMIP6 - r4)', artist='F. Fabiano (ISAC - CNR)')
        #writer = ImageMagickFileWriter(fps = 10)#, frame_size = (600, 300))#, metadata = metadata,
        writer = PillowWriter(fps = 10)
        with writer.saving(fig, cart_out + "{}_anomaly_animation_flat.gif".format(var), dpi):
            for i, (year, col) in enumerate(zip(anni, cset)):
                print(year)
                animate(i, ax)
                writer.grab_frame()
                if year in [1950, 2000, 2025, 2050, 2075]:
                    for jj in range(10): writer.grab_frame()
                elif year == 2100:
                    for jj in range(50): writer.grab_frame()
    else:
        iys = [0]
        clon = 0.
        for i, year in enumerate(anni):
            iys.append(i)
            if year in [1950, 2000, 2025, 2050, 2075]:
                for jj in range(20):
                    iys.append(i)
            elif year == 2100:
                for jj in range(50):
                    iys.append(i)

        line_ani = animation.FuncAnimation(fig, animate, frames = len(iys), fargs = (ax,), interval=100, blit=False)
        filename = cart_out + "{}_anomaly_animation_flat.gif".format(var)
        writer = PillowWriter(fps = 10)
        line_ani.save(filename, writer = writer)

    ## Focus on Europe
    fig, ax = ctl.get_cartopy_fig_ax(visualization = 'nearside', central_lat_lon = (50, 20), bounding_lat = 10., figsize = (8, 6), coast_lw = 1)
    #tit = plt.title('1850')

    # Plotting figure
    proj = ccrs.PlateCarree()
    map_plot = ctl.plot_mapc_on_ax(ax, cosoanom[0], presme.lat, presme.lon, proj, cmappa, cbar_range)

    cax = plt.axes([0.1, 0.11, 0.8, 0.05]) #horizontal
    cb = plt.colorbar(map_plot, cax=cax, orientation='horizontal')#, labelsize=18)
    cb.ax.tick_params(labelsize=14)
    if var == 'tas':
        cb.set_label('Temperature anomaly wrt 1960-1990 (K)', fontsize=14)
    elif var == 'pr':
        cb.set_label('Relative precipitation anomaly wrt 1960-1990 (%)', fontsize=14)

    top    = 0.88  # the top of the subplots
    bottom = 0.20    # the bottom of the subplots
    left   = 0.02    # the left side
    right  = 0.98  # the right side
    hspace = 0.20   # height reserved for white space
    wspace = 0.05    # width reserved for blank space
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

    #showdate = ax.text(0.5, 0.95, '1850', fontweight = 'bold', color = cset[0], bbox=dict(facecolor='lightsteelblue', edgecolor='black', boxstyle='round,pad=1'))

    if save:
        metadata = dict(title='Temperature anomaly (EC-Earth CMIP6 - r4)', artist='F. Fabiano (ISAC - CNR)')
        writer = PillowWriter(fps = 10)
        #writer = ImageMagickFileWriter(fps = 10)#, frame_size = (600, 300))#, metadata = metadata,
        with writer.saving(fig, cart_out + "{}_anomaly_animation_nearside.gif".format(var), dpi):
            for i, (year, col) in enumerate(zip(anni, cset)):
                print(year)
                animate(i, ax)
                writer.grab_frame()
                if year in [1950, 2000, 2025, 2050, 2075]:
                    for jj in range(20): writer.grab_frame()
                elif year == 2100:
                    for jj in range(50): writer.grab_frame()
    else:
        iys = [0]
        for i, year in enumerate(anni):
            iys.append(i)
            if year in [1950, 2000, 2025, 2050, 2075]:
                for jj in range(20):
                    iys.append(i)
            elif year == 2100:
                for jj in range(50):
                    iys.append(i)

        line_ani = animation.FuncAnimation(fig, animate, frames = len(iys), fargs = (ax, ), interval=100, blit=False)
        filename = cart_out + "{}_anomaly_animation_nearside.gif".format(var)
        writer = PillowWriter(fps = 10)
        line_ani.save(filename, writer = writer)

    ## Rotating with focus on Northern mid-latitudes
    def animate_rotate(ii):
        i = iys[ii]
        clon = clons[ii]
        proj = ctl.def_projection('nearside', (clat, clon), bounding_lat = blat)
        fig.clear()
        ax = plt.subplot(projection = proj)
        ax.set_global()
        ax.coastlines(linewidth = 1)

        pc = ccrs.PlateCarree()
        map_plot = ctl.plot_mapc_on_ax(ax, cosoanom[i, ...], presme.lat, presme.lon, pc, cmappa, cbar_range, draw_grid = True)
        year = anni[i]
        print(year)
        color = cset[i]
        tam = tahiss[i] - pimean['tas']

        ax.set_title(r'{} $\to$ {:+5.1f} $\circ$C wrt PI'.format(year, tam))

        ## Colorbar
        cax = plt.axes([0.1, 0.11, 0.8, 0.05]) #horizontal
        cb = plt.colorbar(map_plot, cax=cax, orientation='horizontal')#, labelsize=18)
        cb.ax.tick_params(labelsize=14)
        if var == 'tas':
            cb.set_label('Temperature anomaly wrt 1960-1990 (K)', fontsize=14)
        elif var == 'pr':
            cb.set_label('Relative precipitation anomaly wrt 1960-1990 (%)', fontsize=14)

        top    = 0.88  # the top of the subplots
        bottom = 0.20    # the bottom of the subplots
        left   = 0.02    # the left side
        right  = 0.98  # the right side
        hspace = 0.20   # height reserved for white space
        wspace = 0.05    # width reserved for blank space
        plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

        return

    ##
    clat = 30
    blat = -10

    fig, ax = ctl.get_cartopy_fig_ax(visualization = 'nearside', central_lat_lon = (clat, 0), bounding_lat = blat, figsize = (8, 6), coast_lw = 1)
    #tit = plt.title('1850')

    # Plotting figure
    proj = ccrs.PlateCarree()
    map_plot = ctl.plot_mapc_on_ax(ax, cosoanom[0], presme.lat, presme.lon, proj, cmappa, cbar_range)

    if save:
        metadata = dict(title='Temperature anomaly (EC-Earth CMIP6 - r4)', artist='F. Fabiano (ISAC - CNR)')
        writer = PillowWriter(fps = 5)
        #writer = ImageMagickFileWriter(fps = 10)#, frame_size = (600, 300))#, metadata = metadata,
        clon = 0.
        with writer.saving(fig, cart_out + "{}_anomaly_animation_nearside_rotating.gif".format(var), dpi):
            for i, (year, col) in enumerate(zip(anni, cset)):
                clon = (clon-7.2)%360
                proj = ctl.def_projection('nearside', (clat, clon), bounding_lat = blat)
                fig.clear()
                ax = plt.subplot(projection = proj)
                ax.set_global()
                ax.coastlines(linewidth = 1)

                print(year)
                animate(i, ax)
                writer.grab_frame()
                if year in [1950, 2000, 2025, 2050, 2075]:
                    for jj in range(20):
                        clon = (clon-18)%360
                        proj = ctl.def_projection('nearside', (clat, clon), bounding_lat = blat)
                        fig.clear()
                        ax = plt.subplot(projection = proj)
                        ax.set_global()
                        ax.coastlines(linewidth = 1)

                        print(year)
                        animate(i, ax)
                        writer.grab_frame()
                elif year == 2100:
                    for jj in range(50):
                        clon = (clon-7.2)%360
                        proj = ctl.def_projection('nearside', (clat, clon), bounding_lat = blat)
                        fig.clear()
                        ax = plt.subplot(projection = proj)
                        ax.set_global()
                        ax.coastlines(linewidth = 1)

                        print(year)
                        animate(i, ax)
                        writer.grab_frame()
    else:
        clons = [0.]
        iys = [0]
        clon = 0.
        for i, year in enumerate(anni):
            clon = (clon-7.2)%360
            clons.append(clon)
            iys.append(i)
            if year in [2050, 2075]:
                for jj in range(50):
                    clon = (clon-7.2)%360
                    clons.append(clon)
                    iys.append(i)
            elif year == 2100:
                for jj in range(50):
                    clon = (clon-7.2)%360
                    clons.append(clon)
                    iys.append(i)

        line_ani = animation.FuncAnimation(fig, animate_rotate, frames = len(iys), interval=100, blit=False)

        filename = cart_out + "{}_anomaly_animation_nearside_rotating.gif".format(var)
        writer = PillowWriter(fps = 5)
        line_ani.save(filename, writer = writer)
