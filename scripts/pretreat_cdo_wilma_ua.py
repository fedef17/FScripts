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

for ssp in ['ssp585', 'rcp85']:
    if ssp == 'rcp85':
        cart_orig = '/archive/paolo/cmip5/CMIP5/output1/'
        listacarts = glob.glob(cart_orig + '*/*/rcp85/mon/atmos/Amon/r1i1*/ua/')
    else:
        cart_orig = '/archive/paolo/cmip6/CMIP6/model-output/'
        listacarts = glob.glob(cart_orig + '*/*/ssp585/atmos/Amon/r1i1*/ua/')

    cartou = '/home/federico/work/CMIP6/data_mon_ua/'
    cart_g_s = cartou + 'Smean/'
    ctl.mkdir(cart_g_s)

    for cart in listacarts:
        mod = cart.split('/')[7]
        if ssp == 'rcp85':
            mem = cart.split('/')[12]
        else:
            mem = cart.split('/')[11]

        print('----------------------------\n')
        print(mod, mem)
        print('----------------------------\n')

        cartstrat = cart_g_s + '{}_{}_{}/'.format(mod, mem, ssp)
        ctl.mkdir(cartstrat)

        file_list = [co.split('/')[-1] for co in glob.glob(cart+'ua*nc')]

        for filenam in file_list:
            if float(filenam.split('_')[-1][:-3].split('-')[1][:4]) <= 2100:
                pass
            else:
                continue

            file_in = cart + filenam
            indpo = filenam.index('.nc')

            filepart = filenam[:indpo] + '_Smean.nc'
            file_out = cartstrat + filepart
            command = 'cdo sellevel,25000,20000,15000,10000,7000,5000,3000 {} {}prov.nc'.format(file_in, cartstrat)
            os.system(command)
            command = 'cdo vertmean {}prov.nc {}'.format(cartstrat, file_out)
            os.system(command)

        if len(file_list) > 1:
            command = 'cdo cat {}*Smean.nc {}ua_{}_{}_2015-2100_Smean.nc'.format(cartstrat, cartstrat, mod, mem)
            os.system(command)
        else:
            command = 'cp {}*Smean.nc {}ua_{}_{}_2015-2100_Smean.nc'.format(cartstrat, cartstrat, mod, mem)
            os.system(command)

        # regrid
        command = 'cdo remapbil,r360x180 {}ua_{}_{}_2015-2100_Smean.nc {}ua_{}_{}_{}_2015-2100_Smean_r1.nc'.format(cartstrat, mod, mem, cart_g_s, mod, mem, ssp)
        os.system(command)
