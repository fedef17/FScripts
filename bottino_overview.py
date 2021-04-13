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

from matplotlib.colors import LogNorm
from datetime import datetime

from scipy import stats
import xarray as xr
import xclim as xcl
import glob

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 28
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.axisbelow'] = True

#############################################################################
import glob
import xarray as xr
import xclim

filna = '/nas/BOTTINO/{}/cmorized/cmor_*/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/stabilization-ssp585-2{}/r1i1p1f1/{}/{}/gr/v20210315/*nc'

# filist = glob.glob('/nas/BOTTINO/b100/cmorized/cmor_2*/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/stabilization-ssp585-2100/r1i1p1f1/SImon/siconc/gn/v20210315/*nc')

allru = ['b025', 'b050', 'b100']
colors = ['forestgreen', 'orange', 'violet']

#############################################################################
## SEA ICE
areacelfi = '/nas/BOTTINO/areas.nc'
acel = xr.open_dataset(areacelfi)
areaT = np.array(acel['O1t0.srf'].data)
okslat = lat > 40.

miptab = 'SImon'
var = 'siconc'

fig = plt.figure(figsize = (24,12))
ax1 = plt.subplot(1, 2, 1)
ax2 = plt.subplot(1, 2, 2)

for ru, col in zip(allru, colors):
    filist = glob.glob(filna.format(ru, ru[1:], miptab, var))
    gigi = xr.open_mfdataset(filist, use_cftime=True)

    seaice = np.array(gigi.siconc.data)
    lat = np.array(gigi.latitude.data)

    oksi = seaice[:, okslat]
    oksi[oksi < 15.] = 0.
    oksi[oksi > 15.] = 1.
    oksiarea = oksi*areaok[np.newaxis, :]
    seaicearea = np.nansum(oksiarea, axis = 1)

    date = np.array(gigi.time.data)
    okmarch = np.array([da.month == 3 for da in date])
    oksept = np.array([da.month == 9 for da in date])

    ax1.plot_date(date[okmarch], seaicearea[okmarch], linestyle='solid', marker = 'None', color = col, label = ru)
    ax2.plot_date(date[oksept], seaicearea[oksept], linestyle='solid', marker = 'None', color = col, label = ru)

ax1.set_title('March')
ax2.set_title('September')
ax2.legend()
fig.savefig(cart_out + 'bottseaice.pdf')

#pinuc = xclim.indices.sea_ice_area(gigi, acel.data)

### Mean state temperature e varianza?

### Mean state precipitazione e varianza

### Mean state wind e varianza
