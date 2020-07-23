#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
import glob
from matplotlib import pyplot as plt
from matplotlib import cm

import pickle
import netCDF4 as nc
import climtools_lib as ctl
import climdiags as cd

from matplotlib.colors import LogNorm
from datetime import datetime

from scipy import stats
import pandas as pd

#################################

cart_orig = '/badc/cmip6/data/CMIP6/ScenarioMIP/'
listacarts = glob.glob(cart_orig + '*/*/ssp585/r1i1*/Amon/ua/*/latest/')
cartou = '/work/scratch-nopw/fedef17/prima_outgoing/cmip6/ua_mean_strat/'

for cart in listacarts:
    mod = cart.split('/')[5]
    mem = cart.split('/')[7]
    print('----------------------------\n')
    print(mod, mem)
    print('----------------------------\n')
    cartut = cartou + '{}_{}_ssp585/'.format(mod, mem)
    ctl.mkdir(cartut)

    file_list = [co.split('/')[-1] for co in glob.glob(cart+'ua*nc')]

    for filenam in file_list:
        file_in = cart + filenam
        indpo = filenam.index('.nc')

        filepart = filenam[:indpo] + '_Smean.nc'
        file_out = cartut + filepart
        command = 'cdo sellevel,25000,20000,15000,10000,7000,5000,3000 {} {}prov.nc'.format(file_in, cartou)
        os.system(command)
        command = 'cdo vertmean {}prov.nc {}'.format(cartou, file_out)
        os.system(command)

    if len(file_list) > 1:
        command = 'cdo cat {}*Smean.nc {}ua_{}_{}_2015-2100_Smean.nc'.format(cartut, cartut, mod, mem)
        os.system(command)
    else:
        command = 'cp {}*Smean.nc {}ua_{}_{}_2015-2100_Smean.nc'.format(cartut, cartut, mod, mem)
        os.system(command)

    # regrid
    command = 'cdo remapbil,r360x180 {}ua_{}_{}_2015-2100_Smean.nc {}ua_{}_{}_2015-2100_Smean_r1.nc'.format(cartut, mod, mem, cartou, mod, mem)
    os.system(command)
