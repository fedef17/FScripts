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
col2 = cm.PuRd_r(np.linspace(0.25,0.65,100))
col3 = cm.bone_r(np.linspace(0.25,0.65,100))
colors = np.concatenate([col3, col, col2])
from matplotlib import colors as mcolors
heatmap = mcolors.LinearSegmentedColormap.from_list('heat_strong', colors)

###
dpi = 150
save = False
fps = 7

######################## LEGGO 245
filna = '/nas/archive_CMIP6/CMIP6/model-output/EC-Earth-Consortium/EC-Earth3/ssp245/atmos/Amon/r4i1p1f1/{}/{}*nc'

yeamean_245 = dict()
glomeans_245 = dict()
ru = 'ssp245'
for var in ['tas', 'pr']:
    print(var)
    fils = glob.glob(filna.format(var, var))

    kose = xr.open_mfdataset(fils, use_cftime = True)
    kose = kose.drop_vars('time_bnds')

    cosoye = kose[var].groupby("time.year").mean().compute()
    yeamean_245[(ru, var)] = cosoye

    glomeans_245[(ru, var)] = (cosoye.year.values, ctl.global_mean(kose[var]))

pickle.dump([glomeans_245, yeamean_245], open(cart_out + 'yeamean_245.p', 'wb'))

glomeans_245, yeamean_245 = pickle.load(open(cart_out + 'yeamean_245.p', 'rb'))

########################

glomeans, pimean, yeamean, mapmean = pickle.load(open(cart_in + 'bottino_seasmean_2D.p', 'rb'))

glomeans.update(glomeans_245)
yeamean.update(yeamean_245)

########

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

def animate_double(ii, ax1, ax2):
    i = iys[ii]

    year = anni[i]
    print(year)
    color = cset[i]
    tam = tahiss[i] - pimean['tas']
    tit = r'{} $\to$ {:+5.1f} $\circ$C wrt PI'.format(year, tam)
    fig.suptitle(tit)

    proj = ccrs.PlateCarree()
    ax1.clear()
    ax2.clear()
    map_plot1 = ctl.plot_mapc_on_ax(ax1, costas[i, ...], presme.lat, presme.lon, proj, cmappas[0], cbrangs[0], draw_grid = True)
    map_plot2 = ctl.plot_mapc_on_ax(ax2, cospr[i, ...], presme.lat, presme.lon, proj, cmappas[1], cbrangs[1], draw_grid = True)

    return

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

########

for ssp in ['ssp585', 'ssp245']:
    tahiss = np.concatenate([glomeans[('hist', 'tas')][1], glomeans[(ssp, 'tas')][1]])
    tahiss = ctl.butter_filter(tahiss, 5)

    cosdue = dict()

    for var in ['tas', 'pr']:
        cosopi = yeamean[('pi', var)]
        pime = cosopi.mean('year')
        cosohist = yeamean[('hist', var)]
        cosossp = yeamean[(ssp, var)]
        coso = xr.concat([cosohist, cosossp], dim = 'year')
        presme = coso.sel(year = slice(1960, 1990)).mean('year')

        # if os.path.exists(cart_out + 'histssp_{}_lowpass_LR.p'.format(var)):
        #     cosolow, presme = pickle.load(open(cart_out + 'histssp_{}_lowpass_LR.p'.format(var), 'rb'))
        # else:
        coso = ctl.regrid_dataset(coso, regrid_to_deg = 2.)
        presme = ctl.regrid_dataset(presme, regrid_to_deg = 2.)
        if var == 'tas':
            cosolow = ctl.butter_filter(coso, 5)
        elif var == 'pr':
            cosolow = ctl.butter_filter(coso, 10)
        pickle.dump([cosolow, presme], open(cart_out + 'histssp_{}_lowpass_LR_{}.p'.format(var, ssp), 'wb'))

        anni = coso.year.values

        if var == 'tas':
            #cmappa = cm.get_cmap('RdBu_r')#'gist_ncar')
            cmappa = heatmap
            cosoanom = cosolow-presme.values[np.newaxis, ...]
            #cbar_range = ctl.get_cbar_range(cosoanom, symmetrical = True)
            cbar_range = (-8, 8)
        elif var == 'pr':
            cmappa = cm.get_cmap('BrBG')
            #cosoanom = 100*(cosolow-presme.values[np.newaxis, ...])/presme.values[np.newaxis, ...]
            cosoanom = 86400.*(cosolow-presme.values[np.newaxis, ...])*365 # this should be mm per year
            cbar_range = ctl.get_cbar_range(cosoanom, symmetrical = True)
            cbar_range = (cbar_range[0]/10., cbar_range[1]/10.)
            #cbar_range = (-30, 30)

        clevels = np.linspace(cbar_range[0], cbar_range[1], 21)
        cset = ctl.color_set(len(anni), bright_thres = 0., full_cb_range = True)

        cosdue[var] = cosoanom


        fig, ax = ctl.get_cartopy_fig_ax(visualization = 'standard', central_lat_lon = (0, 0), bounding_lat = None, figsize = (16, 9), coast_lw = 1)
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
            #cb.set_label('Relative precipitation anomaly wrt 1960-1990 (%)', fontsize=14)
            cb.set_label('Precipitation anomaly wrt 1960-1990 (%)', fontsize=14)

        top    = 0.88  # the top of the subplots
        bottom = 0.20    # the bottom of the subplots
        left   = 0.02    # the left side
        right  = 0.98  # the right side
        hspace = 0.20   # height reserved for white space
        wspace = 0.05    # width reserved for blank space
        plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

        #showdate = ax.text(0.5, 0.95, '1850', fontweight = 'bold', color = cset[0], bbox=dict(facecolor='lightsteelblue', edgecolor='black', boxstyle='round,pad=1'))

        iys = [0]
        clon = 0.
        for i, year in enumerate(anni):
            iys.append(i)
            # if year in [1950, 2000, 2025, 2050, 2075]:
            #     for jj in range(20):
            #         iys.append(i)
            if year == 2100:
                for jj in range(100):
                    iys.append(i)

        line_ani = animation.FuncAnimation(fig, animate, frames = len(iys), fargs = (ax,), interval=100, blit=False)
        filename = cart_out + "{}_anim_flat_{}.gif".format(var, ssp)
        writer = PillowWriter(fps = fps)
        line_ani.save(filename, writer = writer)
        del line_ani

        ## Focus on Europe
        fig, ax = ctl.get_cartopy_fig_ax(visualization = 'nearside', central_lat_lon = (50, 20), bounding_lat = 10., figsize = (12, 9), coast_lw = 1)
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

        iys = [0]
        for i, year in enumerate(anni):
            iys.append(i)
            # if year in [1950, 2000, 2025, 2050, 2075]:
            #     for jj in range(20):
            #         iys.append(i)
            if year == 2100:
                for jj in range(100):
                    iys.append(i)

        line_ani = animation.FuncAnimation(fig, animate, frames = len(iys), fargs = (ax, ), interval=100, blit=False)
        filename = cart_out + "{}_anim_nearside_{}.gif".format(var, ssp)
        writer = PillowWriter(fps = fps)
        line_ani.save(filename, writer = writer)
        del line_ani

        ## Rotating with focus on Northern mid-latitudes
        ##
        clat = 30
        blat = -10

        fig, ax = ctl.get_cartopy_fig_ax(visualization = 'nearside', central_lat_lon = (clat, 0), bounding_lat = blat, figsize = (12, 9), coast_lw = 1)
        #tit = plt.title('1850')

        # Plotting figure
        proj = ccrs.PlateCarree()
        map_plot = ctl.plot_mapc_on_ax(ax, cosoanom[0], presme.lat, presme.lon, proj, cmappa, cbar_range)

        clons = [180.]
        iys = [0]
        clon = 0.
        for i, year in enumerate(anni):
            clon = (clon-3.6)%360
            clons.append(clon)
            iys.append(i)
            # if year in [2050, 2075]:
            #     for jj in range(100):
            #         clon = (clon-3.6)%360
            #         clons.append(clon)
            #         iys.append(i)
            if year == 2100:
                for jj in range(100):
                    clon = (clon-3.6)%360
                    clons.append(clon)
                    iys.append(i)

        line_ani = animation.FuncAnimation(fig, animate_rotate, frames = len(iys), interval=100, blit=False)

        filename = cart_out + "{}_anim_nearside_rotating_{}.gif".format(var, ssp)
        writer = PillowWriter(fps = fps)
        line_ani.save(filename, writer = writer)
        del line_ani

    #### Double anim
    ## Focus on Europe

    # cmappas = [heatmap, 'BrBG']
    # cbrangs = [(-8, 8), (-30, 30)]
    #
    # costas = cosdue['tas']
    # cospr = cosdue['pr']
    #
    # figs = ctl.plot_multimap_contour([costas[i], cospr[i]], lat = presme.lat, lon = presme.lon, visualization = 'nearside', central_lat_lon = (50, 20), cmap = cmappas, title = tit, cb_label = ['Temperature anomaly wrt 1960-1990 (K)', 'Relative precipitation anomaly wrt 1960-1990 (%)'], cbar_range = cbrangs, plot_anomalies = True, fix_subplots_shape = (1,2), figsize = (16,9), bounding_lat = 0, use_different_cbars = True, use_different_cmaps = True)
    # fig = figs[0]
    #
    # ax1 = fig.get_axes()[0]
    # ax2 = fig.get_axes()[1]
    #
    # iys = [0]
    # for i, year in enumerate(anni):
    #     iys.append(i)
    #     # if year in [1950, 2000, 2025, 2050, 2075]:
    #     #     for jj in range(20):
    #     #         iys.append(i)
    #     if year == 2100:
    #         for jj in range(100):
    #             iys.append(i)
    #
    # line_ani = animation.FuncAnimation(fig, animate_double, frames = len(iys), fargs = (ax1, ax2, ), interval=100, blit=False)
    # filename = cart_out + "taspr_anim_nearside_{}.gif".format(ssp)
    # writer = PillowWriter(fps = fps)
    # line_ani.save(filename, writer = writer)
    # del line_ani
