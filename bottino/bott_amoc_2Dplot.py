#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
from matplotlib import pyplot as plt
from matplotlib import cm

import pickle

import climtools_lib as ctl
import xarray as xr
import glob

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 22
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['axes.axisbelow'] = True

#############################################################################

if os.uname()[1] == 'hobbes':
    cart_out = '/home/fabiano/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'xaru':
    cart_out = '/home/fedef/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'tintin':
    cart_out = '/home/fabiano/work/lavori/BOTTINO/'

cart_out = cart_out + 'amoc/'

allru = ['b990', 'b025', 'b050', 'b065', 'b080', 'b100']
allsyear = [1990, 2025, 2050, 2065, 2080, 2100]
allnams = ['stabilization-hist-1990', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2065', 'stabilization-ssp585-2080', 'stabilization-ssp585-2100']

#colors = ['teal', 'forestgreen', 'orange', 'chocolate', 'maroon', 'violet']
colors = ['lightslategray', 'forestgreen', 'orange', 'chocolate', 'maroon', 'violet']


allru_wI = ['b065', 'b65I', 'b080', 'b80I', 'b100', 'b00I']
colors_wI = ['chocolate', 'peru', 'maroon', 'firebrick', 'violet', 'plum']

amoc_2D = pickle.load(open('/home/fabiano/Research/lavori/BOTTINO/amoc/amoc_2D_1000.p', 'rb'))

######################################################################

print('figures')
# carto = '/nas/archive_CMIP6/CMIP6/model-output/EC-Earth-Consortium/EC-Earth3/piControl/ocean/Omon/r1i1p1f1/'
# moc_pi = xr.open_mfdataset(carto + 'msftyz/msftyz*nc', use_cftime = True)['msftyz']
# moc_pi = moc_pi.sel(time = slice('2330-01-01', '2380-01-01'))
# moc_pi = moc_pi.mean('time')

# amoc_2D['pi'] = moc_pi
# pickle.dump(amoc_2D, open('/home/fabiano/Research/lavori/BOTTINO/amoc/amoc_2D_1000.p', 'wb'))

cbran = (-19, 19)
cmap = ctl.heatmap()
nlevs = 20

fig, axs = plt.subplots(1, 3, sharey = True, figsize = (18,6))

map_full = ctl.plot_lat_crosssection(amoc_2D['pi'].sel(basin = 0)/1.e9, ax = axs[0], cbar_range = cbran, ylabel = 'Depth (m)', cmap = cmap, add_contour_field = amoc_2D['pi'].sel(basin = 0)/1.e9, add_contour_same_levels = True, n_lines = nlevs-1, n_color_levels = nlevs, plot_anomalies = True)

cbran = (-13, 13)
nlevs = 14
cmap = 'RdBu_r'
# ctl.plot_lat_crosssection(amoc_2D[('b100', 'ini')].sel(basin = 0)/1.e9, ax = axs[1], cbar_range = cbran, xlabel = 'Latitude', cmap = cmap)
# map_plot = ctl.plot_lat_crosssection(amoc_2D[('b100', 'fin')].sel(basin = 0)/1.e9, ax = axs[2], cbar_range = cbran, cmap = cmap)
map_diff = ctl.plot_lat_crosssection(amoc_2D[('b100', 'ini')].sel(basin = 0)/1.e9 - amoc_2D['pi'].sel(basin = 0)/1.e9, ax = axs[1], cbar_range = cbran, xlabel = 'Latitude', cmap = cmap, add_contour_field = amoc_2D[('b100', 'ini')].sel(basin = 0)/1.e9, lw_contour = 0.1, add_contour_lines_step = 3, n_color_levels = nlevs, plot_anomalies = True)
ctl.plot_lat_crosssection(amoc_2D[('b100', 'fin')].sel(basin = 0)/1.e9 - amoc_2D['pi'].sel(basin = 0)/1.e9, ax = axs[2], cbar_range = cbran, cmap = cmap, add_contour_field = amoc_2D[('b100', 'fin')].sel(basin = 0)/1.e9, lw_contour = 0.1, add_contour_lines_step = 3, n_color_levels = nlevs, plot_anomalies = True)

for ax, tit in zip(axs, ['pi', 'b100 - ini', 'b100 - fin']):
    ax.set_title(tit)

# title_obj = plt.title(title, fontsize=20, fontweight='bold')
# title_obj.set_position([.5, 1.05])


cbar1_ax = fig.add_axes([0.1, 0.11, 0.3, 0.05])  # [x, y, width, height]
cb = plt.colorbar(map_full, orientation='horizontal', fraction = 0.1, cax = cbar1_ax)
cb.ax.tick_params(labelsize=14)
cb.set_label('MOC streamfunction (Sv)', fontsize=16)

cbar2_ax = fig.add_axes([0.55, 0.11, 0.3, 0.05])  # [x, y, width, height]
cb = plt.colorbar(map_diff, orientation='horizontal', fraction = 0.1, cax = cbar2_ax)
cb.ax.tick_params(labelsize=14)
cb.set_label('MOC streamfunction change (Sv)', fontsize=16)

# cax = plt.axes([0.3, 0.11, 0.5, 0.05])
# cb = plt.colorbar(map_plot, orientation='horizontal', fraction = 0.1, cax = cax)
# cb.ax.tick_params(labelsize=14)
# cb.set_label('MOC streamfunction (Sv)', fontsize=16)

#plt.gca().invert_yaxis()

top    = 0.88  # the top of the subplots
bottom = 0.3    # the bottom of the subplots
left   = 0.1    # the left side
right  = 0.98  # the right side
hspace = 0.20   # height reserved for white space
wspace = 0.05    # width reserved for blank space
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

fig.savefig(cart_out + 'moc2d_pi_vs_b100.pdf')

####################################################### ATLANTIC
print('figures 2')


basin_ok = 1

fig, axs = plt.subplots(1, 3, sharey = True, figsize = (18,6))

cbran = (-19, 19)
nlevs = 20
cmap = ctl.heatmap()
map_full = ctl.plot_lat_crosssection(amoc_2D['pi'].sel(basin = basin_ok)/1.e9, ax = axs[0], cbar_range = cbran, ylabel = 'Depth (m)', cmap = cmap, add_contour_field = amoc_2D['pi'].sel(basin = basin_ok)/1.e9, add_contour_same_levels = True, n_lines = nlevs-1, n_color_levels = nlevs, plot_anomalies = True)

cbran = (-13, 13)
nlevs = 14
cmap = 'RdBu_r'
map_diff = ctl.plot_lat_crosssection(amoc_2D[('b100', 'ini')].sel(basin = basin_ok)/1.e9 - amoc_2D['pi'].sel(basin = basin_ok)/1.e9, ax = axs[1], cbar_range = cbran, xlabel = 'Latitude', cmap = cmap, add_contour_field = amoc_2D[('b100', 'ini')].sel(basin = basin_ok)/1.e9, lw_contour = 0.1, add_contour_lines_step = 3, n_color_levels = nlevs, plot_anomalies = True)
ctl.plot_lat_crosssection(amoc_2D[('b100', 'fin')].sel(basin = basin_ok)/1.e9 - amoc_2D['pi'].sel(basin = basin_ok)/1.e9, ax = axs[2], cbar_range = cbran, cmap = cmap, add_contour_field = amoc_2D[('b100', 'fin')].sel(basin = basin_ok)/1.e9, lw_contour = 0.1, add_contour_lines_step = 3, n_color_levels = nlevs, plot_anomalies = True)

for ax, tit in zip(axs, ['pi', 'b100 - ini', 'b100 - fin']):
    ax.set_title(tit)

cbar1_ax = fig.add_axes([0.1, 0.11, 0.3, 0.05])  # [x, y, width, height]
cb = plt.colorbar(map_full, orientation='horizontal', fraction = 0.1, cax = cbar1_ax)
cb.ax.tick_params(labelsize=14)
cb.set_label('MOC streamfunction (Sv)', fontsize=16)

cbar2_ax = fig.add_axes([0.55, 0.11, 0.3, 0.05])  # [x, y, width, height]
cb = plt.colorbar(map_diff, orientation='horizontal', fraction = 0.1, cax = cbar2_ax)
cb.ax.tick_params(labelsize=14)
cb.set_label('MOC streamfunction change (Sv)', fontsize=16)

top    = 0.88  # the top of the subplots
bottom = 0.3    # the bottom of the subplots
left   = 0.1    # the left side
right  = 0.98  # the right side
hspace = 0.20   # height reserved for white space
wspace = 0.05    # width reserved for blank space
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

fig.savefig(cart_out + 'moc2d_pi_vs_b100_atl.pdf')


####################################################### INDO-PAC
print('figures 3')

basin_ok = 2

fig, axs = plt.subplots(1, 3, sharey = True, figsize = (18,6))

cbran = (-19, 19)
nlevs = 20
cmap = ctl.heatmap()
map_full = ctl.plot_lat_crosssection(amoc_2D['pi'].sel(basin = basin_ok)/1.e9, ax = axs[0], cbar_range = cbran, ylabel = 'Depth (m)', cmap = cmap, add_contour_field = amoc_2D['pi'].sel(basin = basin_ok)/1.e9, add_contour_same_levels = True, n_lines = nlevs-1, n_color_levels = nlevs, plot_anomalies = True)

cbran = (-13, 13)
nlevs = 14
cmap = 'RdBu_r'
# ctl.plot_lat_crosssection(amoc_2D[('b100', 'ini')].sel(basin = basin_ok)/1.e9, ax = axs[1], cbar_range = cbran, xlabel = 'Latitude', cmap = cmap)
# map_plot = ctl.plot_lat_crosssection(amoc_2D[('b100', 'fin')].sel(basin = basin_ok)/1.e9, ax = axs[2], cbar_range = cbran, cmap = cmap)
map_diff = ctl.plot_lat_crosssection(amoc_2D[('b100', 'ini')].sel(basin = basin_ok)/1.e9 - amoc_2D['pi'].sel(basin = basin_ok)/1.e9, ax = axs[1], cbar_range = cbran, xlabel = 'Latitude', cmap = cmap, add_contour_field = amoc_2D[('b100', 'ini')].sel(basin = basin_ok)/1.e9, lw_contour = 0.1, add_contour_lines_step = 3, n_color_levels = nlevs, plot_anomalies = True)
ctl.plot_lat_crosssection(amoc_2D[('b100', 'fin')].sel(basin = basin_ok)/1.e9 - amoc_2D['pi'].sel(basin = basin_ok)/1.e9, ax = axs[2], cbar_range = cbran, cmap = cmap, add_contour_field = amoc_2D[('b100', 'fin')].sel(basin = basin_ok)/1.e9, lw_contour = 0.1, add_contour_lines_step = 3, n_color_levels = nlevs, plot_anomalies = True)

for ax, tit in zip(axs, ['pi', 'b100 - ini', 'b100 - fin']):
    ax.set_title(tit)

# title_obj = plt.title(title, fontsize=20, fontweight='bold')
# title_obj.set_position([.5, 1.05])

cbar1_ax = fig.add_axes([0.1, 0.11, 0.3, 0.05])  # [x, y, width, height]
cb = plt.colorbar(map_full, orientation='horizontal', fraction = 0.1, cax = cbar1_ax)
cb.ax.tick_params(labelsize=14)
cb.set_label('MOC streamfunction (Sv)', fontsize=16)

cbar2_ax = fig.add_axes([0.55, 0.11, 0.3, 0.05])  # [x, y, width, height]
cb = plt.colorbar(map_diff, orientation='horizontal', fraction = 0.1, cax = cbar2_ax)
cb.ax.tick_params(labelsize=14)
cb.set_label('MOC streamfunction change (Sv)', fontsize=16)

#plt.gca().invert_yaxis()

top    = 0.88  # the top of the subplots
bottom = 0.3    # the bottom of the subplots
left   = 0.1    # the left side
right  = 0.98  # the right side
hspace = 0.20   # height reserved for white space
wspace = 0.05    # width reserved for blank space
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

fig.savefig(cart_out + 'moc2d_pi_vs_b100_indopac.pdf')