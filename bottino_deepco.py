#!/usr/bin/python
# -*- coding: utf-8 -*-


import numpy as np
import sys
import os
#from matplotlib import pyplot as plt
#from matplotlib import cm

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

#############################################################################

#ru = sys.argv[1]
ru = 'pi'

# open our log file
# logname = 'log_ocetrans_{}.log'.format(ru)
# logfile = open(logname,'w') #self.name, 'w', 0)

# # re-open stdout without buffering
# sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1)

# # redirect stdout and stderr to the log file opened above
# os.dup2(logfile.fileno(), sys.stdout.fileno())
# os.dup2(logfile.fileno(), sys.stderr.fileno())

print('total RAM memory used 0:', psutil.virtual_memory()[3]/1.e9)

cart_out = '/g100_work/IscrB_QUECLIM/BOTTINO/bottino_an/'
cart_out = cart_out + 'ocean3d/'
ctl.mkdir(cart_out)

####################################################################################################
process = psutil.Process(os.getpid())

lats = np.linspace(-89.5, 89.5, 180)
lons = np.linspace(0, 359, 360)

yeamean = dict()

print('total RAM memory used 0:', process.memory_info().rss/1e9)

miptab = 'Omon'
var1 = 'mlotstmax'

fia = '/g100_scratch/userexternal/ffabiano/ece3/b050/cmorized/cmor_2222/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/stabilization-ssp585-2050/r1i1p1f1/Ofx/areacello/gn/v20210315/areacello_Ofx_EC-Earth3_stabilization-ssp585-2050_r1i1p1f1_gn.nc' # areacello 
gigi_a = xr.load_dataset(fia, use_cftime = True)['areacello']

#def do_cross(fils, fils2, fils_area, fil_out):
def do_cross(fils, gigi_a, fil_out):
    print("I'm process", os.getpid())
    # Getting % usage of virtual_memory ( 3rd field)

    #cose = []
    #for fi1, fi2, fia in zip(fils, fils2, fils_area):
    #for fi1, fi2 in zip(fils, fils2):
    for fi1 in fils:
        print(fi1)
        print('total RAM memory used 1', process.memory_info().rss/1e9)

        gigi = xr.load_dataset(fi1, use_cftime = True)[var1]

        gigi_2000 = (gigi > 2000.).any('time')
        gigi_1000 = (gigi > 1000.).any('time')
        gigi_500 = (gigi > 500.).any('time')

        gigi_2000_area = (gigi_2000*gigi_a).sum(('i', 'j'))
        gigi_1000_area = (gigi_1000*gigi_a).sum(('i', 'j'))
        gigi_500_area = (gigi_500*gigi_a).sum(('i', 'j'))

        print('total RAM memory used 2', process.memory_info().rss/1e9)

        pickle.dump([gigi_500_area, gigi_1000_area, gigi_2000_area], fil_out)

        fil_out.flush()

        del gigi#, zuki100, zuki700, zuki2000, oht100, oht700, oht2000

    #coda.put(cose)
    #return cose
    return

n_proc = 10

## pi

print(ru)

print(ru)

if ru in ['b050', 'b100', 'b080']:
    datadir = '/g100_scratch/userexternal/ffabiano/ece3/{}/cmorized/'.format(ru)
    filna1 = datadir+'cmor_*/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/*/r1i1p1f1/{}/{}/g*/v*/{}*nc'.format(miptab, var1, var1)
else:
    datadir = '/g100_scratch/userexternal/ffabiano/irods_data/{}/'.format(ru)
    filna1 = datadir+'{}/{}/{}*nc'.format(miptab, var1, var1)

allfils = glob.glob(filna1)
allfils.sort()

filo = open(cart_out + 'deepco_{}.p'.format(ru), 'wb')

print(allfils[0])
do_cross(allfils, gigi_a, filo)

filo.close()
