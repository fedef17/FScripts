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

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 22
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 15
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

filna = '/g100_work/IscrB_QUECLIM/BOTTINO/{}/cmorized/cmor_*/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/stabilization-ssp585-2100/r1i1p1{}/LImon/snw/gr/v*/snw_*nc'
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

fig, axs = plt.subplots(2, 1, figsize = (16,9))

for ru, mem, col in zip(allru, allmems, colors):
    filz = glob.glob(filna.format(ru, mem))
    filz.sort()

    snw = xr.open_mfdataset(filz, use_cftime = True)['snw']

    gr_snw = snw.sel(lat = slice(gr_latsli), lon = slice(gr_lonsli)).groupby('time.year').mean()

    # total volume over greenland. need cell area/land-sea mask. can also approximate
    lamesh, lomesh = np.meshgrid(gr_snw.lat, gr_snw.lon)

    dtheta = np.deg2rad(gr_snw.lat[1]-gr_snw.lat[0])
    dphi = np.deg2rad(gr_snw.lon[1]-gr_snw.lon[0])
    cell_area = ctl.Rearth**2*np.cos(np.deg2rad(lamesh))*dtheta*dphi

    print(np.mean(cell_area), np.max(cell_area), np.min(cell_area))

    water_vol = gr_snw*cell_area
    water_vol = water_vol.sum(['lat', 'lon'])
    snowco[(ru, 'water_vol (m3)')] = water_vol

    axs[0].plot(water_vol.year[:50], water_vol[:50], color = col, label = ru)

    # total freshwater flux: derivative of volume
    frw_flu = -np.diff(water_vol, axis = 0)/(1.e6*3.1e7)
    snowco[(ru, 'freshwater_flux (Sv)')] = frw_flu
    axs[1].plot(frw_flu.year[:50], frw_flu[:50], color = col, label = ru)

axs[0].set_yscale('log')
axs[0].set_ylabel(r'Total water volume ($m^3$)')

axs[1].set_ylabel('Mean meltwater flux (Sv)')
axs[0].set_xlabel('year')

axs[0].legend()

fig.savefig(cart_out + 'check_greenland_snw_melt.pdf')

pickle.dump(snowco, open(cart_out + 'snowcover_{}.p'.format(ru), 'wb'))
