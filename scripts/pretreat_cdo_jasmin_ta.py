#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
import glob
#from matplotlib import pyplot as plt
#from matplotlib import cm

#import pickle
#import netCDF4 as nc
#import climtools_lib as ctl
#import climdiags as cd

#from matplotlib.colors import LogNorm
#from datetime import datetime

#from scipy import stats
#import pandas as pd

#################################

cart_orig = '/badc/cmip6/data/CMIP6/ScenarioMIP/'
listacarts = glob.glob(cart_orig + '*/*/ssp585/r1i1*/Amon/ta/*/latest/')
cartou = '/work/scratch-nopw/fedef17/prima_outgoing/cmip6/ta_mean/'
cart_g_ut = cartou + 'UTmean/'
cart_g_lt = cartou + 'LTmean/'
cart_g_s = cartou + 'Smean/'

def mkdirif(cart):
    if not os.path.exists(cart):
        os.mkdir(cart)
    return

for cart in listacarts:
    mod = cart.split('/')[7]
    mem = cart.split('/')[9]
    print('----------------------------\n')
    print(mod, mem)
    print('----------------------------\n')

    cartut = cart_g_ut + '{}_{}_ssp585/'.format(mod, mem)
    mkdirif(cartut)
    cartlt = cart_g_lt + '{}_{}_ssp585/'.format(mod, mem)
    mkdirif(cartlt)
    cartstrat = cart_g_s + '{}_{}_ssp585/'.format(mod, mem)
    mkdirif(cartstrat)

    file_list = [co.split('/')[-1] for co in glob.glob(cart+'ta*nc')]

    for filenam in file_list:
        file_in = cart + filenam
        indpo = filenam.index('.nc')

        filepart = filenam[:indpo] + '_UTmean.nc'
        file_out = cartut + filepart
        command = 'cdo sellevel,40000,30000,25000,20000,15000 {} {}prov.nc'.format(file_in, cartou)
        os.system(command)
        command = 'cdo vertmean {}prov.nc {}'.format(cartou, file_out)
        os.system(command)

        filepart = filenam[:indpo] + '_Smean.nc'
        file_out = cartstrat + filepart
        command = 'cdo sellevel,25000,20000,15000,10000,7000,5000,3000 {} {}prov.nc'.format(file_in, cartou)
        os.system(command)
        command = 'cdo vertmean {}prov.nc {}'.format(cartou, file_out)
        os.system(command)

        filepart = filenam[:indpo] + '_LTmean.nc'
        file_out = cartlt + filepart
        command = 'cdo sellevel,100000,92500,85000,70000 {} {}prov.nc'.format(file_in, cartou)
        os.system(command)
        command = 'cdo vertmean {}prov.nc {}'.format(cartou, file_out)
        os.system(command)

    if len(file_list) > 1:
        command = 'cdo cat {}*UTmean.nc {}ta_{}_{}_2015-2100_UTmean.nc'.format(cartut, cartut, mod, mem)
        os.system(command)
        command = 'cdo cat {}*Smean.nc {}ta_{}_{}_2015-2100_Smean.nc'.format(cartstrat, cartstrat, mod, mem)
        os.system(command)
        command = 'cdo cat {}*LTmean.nc {}ta_{}_{}_2015-2100_LTmean.nc'.format(cartlt, cartlt, mod, mem)
        os.system(command)
    else:
        command = 'cp {}*UTmean.nc {}ta_{}_{}_2015-2100_UTmean.nc'.format(cartut, cartut, mod, mem)
        os.system(command)
        command = 'cp {}*Smean.nc {}ta_{}_{}_2015-2100_Smean.nc'.format(cartstrat, cartstrat, mod, mem)
        os.system(command)
        command = 'cp {}*LTmean.nc {}ta_{}_{}_2015-2100_LTmean.nc'.format(cartlt, cartlt, mod, mem)
        os.system(command)

    # regrid
    command = 'cdo remapbil,r360x180 {}ta_{}_{}_2015-2100_UTmean.nc {}ta_{}_{}_2015-2100_UTmean_r1.nc'.format(cartut, mod, mem, cart_g_ut, mod, mem)
    os.system(command)

    command = 'cdo remapbil,r360x180 {}ta_{}_{}_2015-2100_Smean.nc {}ta_{}_{}_2015-2100_Smean_r1.nc'.format(cartstrat, mod, mem, cart_g_s, mod, mem)
    os.system(command)

    command = 'cdo remapbil,r360x180 {}ta_{}_{}_2015-2100_LTmean.nc {}ta_{}_{}_2015-2100_LTmean_r1.nc'.format(cartlt, mod, mem, cart_g_lt, mod, mem)
    os.system(command)
