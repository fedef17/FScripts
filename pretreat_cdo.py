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

cart_orig = '/nas/archive_CMIP6/CMIP6_jasmin/'
listacarts = glob.glob(cart_orig + '*/*/*/*/*/*/*/*/')

for cart in listacarts:
    mod = cart.split('/')[5]
    mem = cart.split('/')[7]

    cartut = '/'.join(cart.split('/')[:-2]) + '/UTmean/'
    os.mkdir(cartut)
    cartstrat = '/'.join(cart.split('/')[:-2]) + '/Smean/'
    os.mkdir(cartstrat)
    cartlt = '/'.join(cart.split('/')[:-2]) + '/LTmean/'
    os.mkdir(cartlt)

    file_list = [co.split('/')[-1] for co in glob.glob(cart+'ta*nc')]

    for filenam in file_list:
        file_in = cart + filenam
        indpo = filenam.index('.nc')

        filepart = filenam[:indpo] + '_UTmean.nc'
        file_out = cartut + filepart
        command = 'cdo sellevel,40000,30000,25000,20000,15000 {} prov.nc'.format(file_in)
        os.system(command)
        command = 'cdo vertmean prov.nc {}'.format(file_out)
        os.system(command)

        filepart = filenam[:indpo] + '_Smean.nc'
        file_out = cartstrat + filepart
        command = 'cdo sellevel,25000,20000,15000,10000,7000,5000,3000 {} prov.nc'.format(file_in)
        os.system(command)
        command = 'cdo vertmean prov.nc {}'.format(file_out)
        os.system(command)

        filepart = filenam[:indpo] + '_LTmean.nc'
        file_out = cartlt + filepart
        command = 'cdo sellevel,100000,92500,85000,70000 {} prov.nc'.format(file_in)
        os.system(command)
        command = 'cdo vertmean prov.nc {}'.format(file_out)
        os.system(command)

    if len(file_list) > 1:
        command = 'cdo cat {}*UTmean.nc {}ta_{}_{}_2015-2100_UTmean.nc'.format(cartut, cartut, mod, mem)
        os.system(command)
        command = 'cdo cat {}*Smean.nc {}ta_{}_{}_2015-2100_Smean.nc'.format(cartstrat, cartstrat, mod, mem)
        os.system(command)
        command = 'cdo cat {}*LTmean.nc {}ta_{}_{}_2015-2100_LTmean.nc'.format(cartlt, cartlt, mod, mem)
        os.system(command)
