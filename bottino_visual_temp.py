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

glomeans, pimean, yeamean, mapmean = pickle.load(open(cart_in + 'bottino_seasmean_2D.p', 'rb'))

for var in ['tas', 'pr']:
    cosopi = yeamean[('pi', var)]
    cosohist = yeamean[('hist', var)]
    cosossp = yeamean[('ssp585', var)]
    coso = xr.concat([cosohist, cosossp], dim = 'year')

    anni = coso.year.values

    fig, ax = ctl.get_cartopy_fig_ax(visualization = 'standard', central_lat_lon = (0, 0), bounding_lat = None, figsize = (16, 12), coast_lw = 1)

    if var == 'temp':
        cb.set_label('Temp anomaly (K)', fontsize=16)
        cmappa = cm.get_cmap('gist_ncar')#('RdBu_r')
        cosoanom = coso-cosopi.mean('year')
        cbar_range = ctl.get_cbar_range(cosoanom, symmetrical = True)
        #cbar_range = (-5, 5)
    elif var == 'pr':
        cb.set_label('Pr anomaly ()', fontsize=16)
        cmappa = cm.get_cmap('BrBG')
        cosoanom = (coso-cosopi.mean('year'))/cosopi.mean('year')
        #cbar_range = ctl.get_cbar_range(cosoanom, symmetrical = True)
        cbar_range = (-50, 50)

    clevels = np.linspace(cbar_range[0], cbar_range[1], 21)
    cset = ctl.color_set(len(anni), bright_thres = 0., full_cb_range = True)

    def animate(i, ax):
        proj = ccrs.PlateCarree()
        ax.clear()
        map_plot = ctl.plot_mapc_on_ax(ax, coso.values[i, ...], coso.lat, coso.lon, proj, cmappa, cbar_range)
        year = anni[i]
        color = cset[i]
        tit.set_text('{}'.format(year))
        #showdate.set_text('{}'.format(year))#, color = color)
        #showdate.update(color = color)
        ax.relim()
        ax.autoscale_view()
        return

    # Plotting figure
    proj = ccrs.PlateCarree()
    map_plot = ctl.plot_mapc_on_ax(ax, coso.values[0], coso.lat, coso.lon, proj, cmappa, cbar_range)

    cax = plt.axes([0.1, 0.11, 0.8, 0.05]) #horizontal
    cb = plt.colorbar(map_plot, cax=cax, orientation='horizontal')#, labelsize=18)
    cb.ax.tick_params(labelsize=14)
    if var == 'temp':
        cb.set_label('Temp anomaly (K)', fontsize=16)
    elif var == 'pr':
        cb.set_label('Pr anomaly ()', fontsize=16)

    top    = 0.88  # the top of the subplots
    bottom = 0.20    # the bottom of the subplots
    left   = 0.02    # the left side
    right  = 0.98  # the right side
    hspace = 0.20   # height reserved for white space
    wspace = 0.05    # width reserved for blank space
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

    tit = plt.title('1850')
    #showdate = ax.text(0.5, 0.95, '1850', fontweight = 'bold', color = cset[0], bbox=dict(facecolor='lightsteelblue', edgecolor='black', boxstyle='round,pad=1'))

    save = True
    if save:
        metadata = dict(title='Temperature anomaly (EC-Earth CMIP6 - r4)', artist='F. Fabiano (ISAC - CNR)')
        writer = ImageMagickFileWriter(fps = 10, metadata = metadata)#, frame_size = (1200, 900))
        with writer.saving(fig, cart + "{}_anomaly_animation_flat.gif".format(var), 150):
            for i, (year, col) in enumerate(zip(anni, cset)):
                print(year)
                animate(i)
                writer.grab_frame()
    else:
        line_ani = animation.FuncAnimation(fig, animate, len(anni), interval=100, blit=False)
