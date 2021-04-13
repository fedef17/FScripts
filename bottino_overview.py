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

# filna = '/nas/BOTTINO/b{}/cmorized/cmor_*/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/stabilization-ssp585-2{}/r1i1p1f1/{}/{}/gr/v20210315/*nc'
#
# filist = glob.glob(filna.format('025', '025', 'Amon', 'pr'))
#
# gigi = xr.open_mfdataset(filist)

import glob
import xarray as xr
import xclim

filist = glob.glob('/nas/BOTTINO/b100/cmorized/cmor_2*/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/stabilization-ssp585-2100/r1i1p1f1/SImon/siconc/gn/v20210315/*nc')
gigi = xr.open_mfdataset(filist, use_cftime=True)

areacelfi = '/nas/BOTTINO/b025/cmorized/cmor_2025/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/stabilization-ssp585-2025/r1i1p1f1/fx/areacello/gn/v20210315/areacella_fx_EC-Earth3_stabilization-ssp585-2025_r1i1p1f1_gn.nc'
acel = xr.open_dataset(areacelfi)
pinuc = xclim.indices.sea_ice_area(gigi, acel.data)
