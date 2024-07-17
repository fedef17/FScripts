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
import glob
import xclim

#############################################

cart_out = '/home/fabiano/Research/lavori/BOTTINO/cmip6_data/'
#cart_out = '/home/fedef/Research/lavori/BOTTINO/forcing/'

histtag = 'CMIP_UoM-CMIP-1-2-0_gr1-GMNHSH_0000-2014'
ssptag = 'ScenarioMIP_UoM-REMIND-MAGPIE-ssp585-1-2-1_gr1-GMNHSH_2015-2500'

for gas, gnam in zip(['carbon-dioxide', 'methane', 'nitrous-oxide'], ['carbon_dioxide', 'methane', 'nitrous_oxide']):
    print(gas)
    ppms = []
    relcon = []

    # historical
    tag = histtag
    ye = 1990
    filo = glob.glob(cart_out + '*{}*{}.nc'.format(gas, tag))[0]
    co2 = xr.load_dataset(filo, decode_times = False)['mole_fraction_of_{}_in_air'.format(gnam)]

    #if picon is None:
    picon = co2.values[1850, 0]
    print('PI', picon)
    ppms.append(picon)

    cocos = co2.values[ye, 0]
    print('{} - {:6.2f}'.format(ye, cocos/picon))
    relcon.append(cocos/picon)
    ppms.append(cocos)

    # scenario
    tag = ssptag
    yeas = [2025, 2050, 2065, 2080, 2100]
    filo = glob.glob(cart_out + '*{}*{}.nc'.format(gas, tag))[0]
    co2 = xr.load_dataset(filo, use_cftime = True)['mole_fraction_of_{}_in_air'.format(gnam)]
    for ye in yeas:
        cocos = co2.sel(time = slice('{}-01-01'.format(ye), '{}-12-31'.format(ye)), sector = 0).values[0]
        print('{} - {:6.2f}'.format(ye, cocos/picon))
        relcon.append(cocos/picon)
        ppms.append(cocos)

    print((7*'{:6.0f} &').format(*ppms))
    print((6*'{:7.2f} &').format(*relcon))


aers = xr.load_dataset(cart_out + 'MACv2.0-SP_REMIND-MAGPIE-SSP5-Ref.nc')
for var in aers.data_vars:
    print(aers[var].long_name)
aers['aod_spmx'].shape
