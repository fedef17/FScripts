#!/usr/bin/python
# -*- coding: utf-8 -*-


import numpy as np
import sys
import os
from matplotlib import pyplot as plt
from matplotlib import cm

import pickle
#import netCDF4 as nc

import climtools_lib as ctl
#import climdiags as cd
#from tunlib import gregplot_on_ax

#from matplotlib.colors import LogNorm
#from datetime import datetime

#from scipy import stats
import xarray as xr
import glob
#import xclim

import multiprocessing as mp
import psutil

plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
titlefont = 22
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.axisbelow'] = True

#############################################################################

# ru = sys.argv[1]
#
# # open our log file
# logname = 'log_oceall_{}.log'.format(ru)
# logfile = open(logname,'w') #self.name, 'w', 0)
#
# # re-open stdout without buffering
# sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1)
#
# # redirect stdout and stderr to the log file opened above
# os.dup2(logfile.fileno(), sys.stdout.fileno())
# os.dup2(logfile.fileno(), sys.stderr.fileno())
#
# print('total RAM memory used 0:', psutil.virtual_memory()[3]/1.e9)

# if os.uname()[1] == 'hobbes':
#     cart_out = '/home/fabiano/Research/lavori/BOTTINO/'
# elif os.uname()[1] == 'xaru':
#     cart_out = '/home/fedef/Research/lavori/BOTTINO/'
# elif os.uname()[1] == 'tintin':
#     cart_out = '/home/fabiano/work/lavori/BOTTINO/'

cart_out = '/g100_work/IscrB_QUECLIM/BOTTINO/bottino_an/'
cart_out = cart_out + 'landice/'
ctl.mkdir(cart_out)

filna = '/g100_work/IscrB_QUECLIM/BOTTINO/{}/cmorized/cmor_*/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/stabilization-ssp585-2100/r1i1p1{}/{}/{}/*/v*/{}_*nc'

# filna = '/nas/BOTTINO/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/{}/{}i1p1f1/{}/{}/*nc'

allru = ['b100', 'b00A', 'b00I']#['pi',
#allnams = ['stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']#'piControl',
#nam = allnams[allru.index(ru)]
allmems = ['f1', 'f2', 'f3']
colors = ['violet', 'chocolate', 'steelblue']

####################################################################################################
gr_latsli = (60., 85.)
gr_lonsli = (288., 350.)

snowco = dict()

miptab = 'LImon'
var = 'snw'

fig, axs = plt.subplots(2, 1, figsize = (16,9))

for ru, mem, col in zip(allru, allmems, colors):
    print(ru)
    filz = glob.glob(filna.format(ru, mem, miptab, var, var))
    filz.sort()

    snw = xr.open_mfdataset(filz[:100], use_cftime = True)[var]

    gr_snw = snw.sel(lat = slice(*gr_latsli), lon = slice(*gr_lonsli)).groupby('time.year').mean()

    # total volume over greenland. need cell area/land-sea mask. can also approximate
    lamesh, lomesh = np.meshgrid(gr_snw.lat, gr_snw.lon)

    dtheta = np.deg2rad(gr_snw.lat[1]-gr_snw.lat[0]).values
    dphi = np.deg2rad(gr_snw.lon[1]-gr_snw.lon[0]).values
    cell_area = ctl.Rearth**2*np.cos(np.deg2rad(lamesh))*dtheta*dphi
    cell_area = np.swapaxes(cell_area,0,1)

    print(np.mean(cell_area), np.max(cell_area), np.min(cell_area))

    water_vol = gr_snw*cell_area[np.newaxis, ...]
    water_vol = water_vol.sum(['lat', 'lon'])
    snowco[(ru, 'water_vol (m3)')] = water_vol

    axs[0].plot(water_vol.year, water_vol, color = col, label = ru)

    # total freshwater flux: derivative of volume
    frw_flu = -np.diff(water_vol, axis = 0)/(1.e6*3.1e7)
    snowco[(ru, 'melted_snpw (Sv)')] = frw_flu
    axs[1].plot(water_vol.year[1:], frw_flu, color = col, label = ru)

axs[0].set_yscale('log')
axs[0].set_ylabel(r'Total water volume ($m^3$)')
axs[1].set_ylabel('Water from snow melting (Sv)')#', yearly average)')
axs[0].set_xlabel('year')
axs[0].legend()
#fig.savefig(cart_out + 'check_greenland_snw_melt.pdf')

axs[0].set_xlim(2100, 2200)
axs[1].set_xlim(2100, 2200)
fig.savefig(cart_out + 'check_greenland_snw_melt_zoom.pdf')

pickle.dump(snowco, open(cart_out + 'snowcover_{}.p'.format(ru), 'wb'))

##### Reading friver
cart_run = '/g100_scratch/userexternal/ffabiano/ece3/b00I/runtime/'
areas = xr.load_dataset(cart_run + 'areas.nc')
areavar = 'O1t0.srf'
okarea = areas[areavar]
okarea = okarea.rename({'x_3' : 'i', 'y_3' : 'j'})

miptab = 'Omon'
var = 'friver'

gr_latsli = (55., 90.)
gr_lonsli = (270., 355.)

fig, axs = plt.subplots(2, 1, figsize = (16,9))

for ru, mem, col in zip(allru, allmems, colors):
    print(ru)
    filz = glob.glob(filna.format(ru, mem, miptab, var, var))
    filz.sort()

    gigi = xr.open_mfdataset(filz[:100], use_cftime = True)[var]

    cond = (gigi.latitude > gr_latsli[0]) & (gigi.latitude < gr_latsli[1]) & (gigi.longitude > gr_lonsli[0]) & (gigi.longitude < gr_lonsli[1]) & (gigi.i > 235)
    gogo = gigi.where(cond)
    okar = okarea.where(cond)

    frwat = (gogo*okar).sum(['i', 'j'])*0.001/1.e6 # kg * m2 * s-1 to m3 * s-1 to Sv

    pino = frwat.groupby('time.year').mean()
    pimax = frwat.groupby('time.year').max()

    snowco[(ru, 'freshwater flux (Sv)')] = frwat

    axs[0].plot(pimax.year, pimax, color = col, label = ru)
    axs[1].plot(pino.year, pino, color = col, label = ru)

#axs[0].set_yscale('log')
axs[0].set_ylabel(r'Monthly max freshwater flux (Sv)')
axs[1].set_ylabel('Yearly mean freshwater flux (Sv)')
axs[1].set_xlabel('year')
axs[0].legend()
axs[0].set_xlim(2100, 2200)
axs[1].set_xlim(2100, 2200)
fig.savefig(cart_out + 'check_greenland_friver.pdf')

pickle.dump(snowco, open(cart_out + 'snowcover_{}.p'.format(ru), 'wb'))
#pickle.dump(friver, open(cart_out + 'friver_{}.p'.format(ru), 'wb'))

###
miptab = 'Amon'
var = 'tas'

fig, axs = plt.subplots(2, 1, figsize = (16,9))

fig2, ax2 = plt.subplots()

for ru, mem, col in zip(allru, allmems, colors):
    print(ru)
    filz = glob.glob(filna.format(ru, mem, miptab, var, var))
    filz.sort()

    gigi = xr.open_mfdataset(filz[:100], use_cftime = True)[var]
    gigi = gigi-273.15

    gtas_mon = ctl.global_mean(gigi[:60])
    ax2.plot(gtas_mon, color = col, label = ru)

    ygigi = gigi.groupby('time.year').mean()
    gtas = ctl.global_mean(ygigi)

    gr_tas = gigi.sel(lat = slice(*gr_latsli), lon = slice(*gr_lonsli))#.groupby('time.year').mean()

    gr_tas_max = gr_tas.groupby('time.year').max()
    gr_max = gr_tas_max.max(['lat', 'lon'])
    gr_mean = gr_tas_max.mean(['lat', 'lon'])
    gr_min = gr_tas_max.min(['lat', 'lon'])

    snowco[(ru, 'gtas')] = gtas
    snowco[(ru, 'Greenland mean temp')] = gr_mean

    axs[0].plot(ygigi.year, gtas, color = col, label = ru)

    axs[1].plot(gr_max.year, gr_mean, color = col, label = ru)
    axs[1].plot(gr_max.year, gr_max, color = col, label = ru, ls = ':')
    axs[1].plot(gr_min.year, gr_min, color = col, label = ru, ls = ':')

axs[0].set_ylabel(r'Global tas')
axs[1].set_ylabel('Mean (max/min) summer peak temp over Greenland')
axs[0].set_xlabel('year')
axs[0].legend()
axs[0].set_xlim(2100, 2200)
axs[1].set_xlim(2100, 2200)
fig.savefig(cart_out + 'check_greenland_temp.pdf')

ax2.legend()
fig2.savefig(cart_out + 'check_monthly_temp_start.pdf')

pickle.dump(snowco, open(cart_out + 'snowcover_{}.p'.format(ru), 'wb'))


###
masfi = cart_run + 'masks.nc'
cose = xr.load_dataset(masfi)
land_mask = ~cose['RnfA.msk'].values.astype('bool') # 1 over land

miptab = 'Amon'
var1 = 'rsdscs'
var2 = 'rsuscs'

fig, ax = plt.subplots(figsize = (16,9))

alb_maps = dict()

for ru, mem, col in zip(allru, allmems, colors):
    print(ru)
    filz = glob.glob(filna.format(ru, mem, miptab, var1, var1))
    filz.sort()
    gigi1 = xr.open_mfdataset(filz[:100], use_cftime = True)[var1]

    filz = glob.glob(filna.format(ru, mem, miptab, var2, var2))
    filz.sort()
    gigi2 = xr.open_mfdataset(filz[:100], use_cftime = True)[var2]

    gigi = gigi2/gigi1
    alb_maps[ru] = gigi[0]
    gigi = gigi.where(land_mask)
    gr_gigi = gigi.sel(lat = slice(*gr_latsli), lon = slice(*gr_lonsli))#.groupby('time.year').mean()

    # ygigi = gr_gigi.sel('time.month' == 9)
    ygigi = gr_gigi[gr_gigi.groupby('time.month').groups[9]]
    mean_albedo = ygigi.mean(['lat', 'lon']).groupby('time.year').mean()
    #mean_albedo = ygigi.values[:, land_mask].mean(axis = 1)

    ax.plot(mean_albedo.year, mean_albedo, color = col, label = ru)

ax.set_ylabel(r'Mean greenland albedo')
ax.legend()
ax.set_xlim(2100, 2200)

fig.savefig(cart_out + 'check_greenland_albedo.pdf')

ctl.plot_multimap_contour([alb_maps['b00A']-alb_maps['b100'], alb_maps['b00I']-alb_maps['b100']], filename = cart_out + 'check_albedo_global_month0.pdf', subtitles = ['b00A-b100', 'b00I-b100'], figsize = (16,9), plot_anomalies = True, cbar_range = (-0.2, 0.2))

#pickle.dump(snowco, open(cart_out + 'snowcover_{}.p'.format(ru), 'wb'))
