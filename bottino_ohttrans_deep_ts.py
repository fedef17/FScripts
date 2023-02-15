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

ru = sys.argv[1]

# open our log file
logname = 'log_ocetrans_{}.log'.format(ru)
logfile = open(logname,'w') #self.name, 'w', 0)

# re-open stdout without buffering
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1)

# redirect stdout and stderr to the log file opened above
os.dup2(logfile.fileno(), sys.stdout.fileno())
os.dup2(logfile.fileno(), sys.stderr.fileno())

print('total RAM memory used 0:', psutil.virtual_memory()[3]/1.e9)

cart_out = '/g100_work/IscrB_QUECLIM/BOTTINO/bottino_an/'
cart_out = cart_out + 'ocean3d/'
ctl.mkdir(cart_out)

# if ru in ['b100', 'b050', 'b080', 'b065', 'b990']:
#     cartbase = '/g100_scratch/userexternal/ffabiano/ece3/'
# else:
#     cartbase = '/g100_work/IscrB_QUECLIM/BOTTINO/'

# filna = cartbase + '{}/cmorized/cmor_*/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/*/*/{}/{}/g*/v*/{}*nc'

####################################################################################################
process = psutil.Process(os.getpid())

basnames = ['atlmsk', 'indmsk', 'pacmsk']
subbas = xr.load_dataset(cart_out + 'subbasins.nc')

lats = np.linspace(-89.5, 89.5, 180)
lons = np.linspace(0, 359, 360)

yeamean = dict()

print('total RAM memory used 0:', process.memory_info().rss/1e9)

miptab = 'Omon'
var1 = 'bigthetao'
var2 = 'wo'
# miptab2 = 'Ofx'
# var3 = 'areacello'

fia = '/g100_scratch/userexternal/ffabiano/ece3/b050/cmorized/cmor_2222/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/stabilization-ssp585-2050/r1i1p1f1/Ofx/areacello/gn/v20210315/areacello_Ofx_EC-Earth3_stabilization-ssp585-2050_r1i1p1f1_gn.nc' # areacello is the same for all
gigi_a = xr.load_dataset(fia, use_cftime = True)['areacello']

def roundlat(ds, ndec = 10):
    ds = ds.assign_coords(lat = ds.lat.round(ndec))
    return ds

def do_cross(fils, fils2, gigi_a, fil_out):
    print("I'm process", os.getpid())
    # Getting % usage of virtual_memory ( 3rd field)

    gigi = xr.open_mfdataset(fils, use_cftime = True, chunks={'time' : 1200}, preprocess = roundlat)['bigthetao']
    gigi2 = xr.open_mfdataset(fils2, use_cftime = True, chunks={'time' : 1200}, preprocess = roundlat)['wo']

    print('Read lazy')

    pino = gigi2.interp(lev = gigi.lev)

    oht = pino*gigi ## T*w

    oht_lev = (oht*gigi_a.values[np.newaxis, np.newaxis, ...]).mean('time').sum(['i', 'j']) # integrating horizontally and averaging over time

    dtdz = gigi.differentiate('lev') ## dT/dz
    dtdz_lev = (dtdz*gigi_a.values[np.newaxis, np.newaxis, ...]).mean('time').sum(['i', 'j']) # integrating horizontally and averaging over time

    print('total RAM memory used 2', process.memory_info().rss/1e9)

    oht_lev = oht_lev.compute()
    print('Done oht_lev')
    dtdz_lev = dtdz_lev.compute()
    print('Done dtdz_lev')

    pickle.dump([oht_lev, dtdz_lev], fil_out)

    return

n_proc = 10
#for ru, nam in zip(allru, allnams):

print(ru)

if ru in ['b025', 'b050', 'b100', 'b080']:
    datadir = '/g100_scratch/userexternal/ffabiano/ece3/{}/cmorized/'.format(ru)
    filna1 = datadir+'cmor_*/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/*/r1i1p1f1/{}/{}/g*/v*/{}*nc'.format(miptab, var1, var1)
    filna2 = datadir+'cmor_*/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/*/r1i1p1f1/{}/{}/g*/v*/{}*nc'.format(miptab, var2, var2)
else:
    datadir = '/g100_scratch/userexternal/ffabiano/ece3/{}/cmorized/'.format(ru)
    filna1 = datadir+'cmor_*/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/*/r1i1p1f1/{}/{}/g*/v*/{}*nc'.format(miptab, var1, var1)
    datadir = '/g100_scratch/userexternal/ffabiano/irods_data/{}/'.format(ru)
    filna2 = datadir+'{}/{}/{}*nc'.format(miptab, var2, var2)

allfils = glob.glob(filna1)
allfils.sort()
allfils2 = glob.glob(filna2)
allfils2.sort()

filo = open(cart_out + 'ohtrans_ts_{}.p'.format(ru), 'wb')

#cose = do_cross(allfils)
print(allfils[0])
print(allfils2[0])

do_cross(allfils, allfils2, gigi_a, filo)

filo.close()

# allchu = np.array_split(allfils, n_proc)
#
# pool = mp.Pool(n_proc)
# cose = pool.map(do_cross, allchu, 1)

#cose = ctl.run_parallel(do_cross, n_proc, args = allchu)

#pickle.dump(cose, )
