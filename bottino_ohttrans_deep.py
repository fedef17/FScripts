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

if ru in ['b100', 'b050', 'b080', 'b065', 'b990']:
    cartbase = '/g100_scratch/userexternal/ffabiano/ece3/'
else:
    cartbase = '/g100_work/IscrB_QUECLIM/BOTTINO/'

filna = cartbase + '{}/cmorized/cmor_*/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/*/*/{}/{}/g*/v*/{}*nc'

####################################################################################################

basnames = ['atlmsk', 'indmsk', 'pacmsk']
subbas = xr.load_dataset(cart_out + 'subbasins.nc')

lats = np.linspace(-89.5, 89.5, 180)
lons = np.linspace(0, 359, 360)

yeamean = dict()

print('total RAM memory used 0:', psutil.virtual_memory()[3]/1.e9)

miptab = 'Omon'
var1 = 'bigthetao'
var2 = 'wo'
# miptab2 = 'Ofx'
# var3 = 'areacello'

fia = '/g100_work/IscrB_QUECLIM/BOTTINO/b050/cmorized/cmor_2222/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/stabilization-ssp585-2050/r1i1p1f1/Ofx/areacello/gn/v20210315/areacello_Ofx_EC-Earth3_stabilization-ssp585-2050_r1i1p1f1_gn.nc' # areacello is the same for all
gigi_a = xr.load_dataset(fia, use_cftime = True)['areacello']

#def do_cross(fils, fils2, fils_area, fil_out):
def do_cross(fils, fils2, gigi_a, fil_out):
    print("I'm process", os.getpid())
    # Getting % usage of virtual_memory ( 3rd field)

    #cose = []
    #for fi1, fi2, fia in zip(fils, fils2, fils_area):
    for fi1, fi2 in zip(fils, fils2):
        print(fi1)
        print('total RAM memory used 1:', psutil.virtual_memory()[3]/1.e9)

        gigi = xr.load_dataset(fi1, use_cftime = True)['bigthetao']
        gigi2 = xr.load_dataset(fi2, use_cftime = True)['wo']

        pino = gigi2.interp(lev = gigi.lev)

        oht = pino*gigi2*gigi_a.values[np.newaxis, np.newaxis, ...]
        # to have real oht must multiply by cp

        oht_lev = oht.mean('time').sum(['i', 'j'])

        oht100 = oht.sel(lev = slice(98., 99.)).mean('time').mean('lev')
        oht700 = oht.sel(lev = slice(697., 698.)).mean('time').mean('lev')
        oht2000 = oht.sel(lev = slice(1945., 1946.)).mean('time').mean('lev')

        zuki100 = ctl.regrid_dataset(oht100, lats, lons)
        zuki700 = ctl.regrid_dataset(oht700, lats, lons)
        zuki2000 = ctl.regrid_dataset(oht2000, lats, lons)

        print('total RAM memory used 2:', psutil.virtual_memory()[3]/1.e9)

        pickle.dump([oht_lev, zuki100, zuki700, zuki2000], fil_out)

        print('total RAM memory used 3:', psutil.virtual_memory()[3]/1.e9)
        fil_out.flush()

        del gigi, gigi2, oht, oht_lev, zuki100, zuki700, zuki2000, oht100, oht700, oht2000

    #coda.put(cose)
    #return cose
    return

n_proc = 10
#for ru, nam in zip(allru, allnams):

print(ru)

allfils = glob.glob(filna.format(ru, miptab, var1, var1))
allfils.sort()
allfils2 = glob.glob(filna.format(ru, miptab, var2, var2))
allfils2.sort()

filo = open(cart_out + 'ohtrans_{}.p'.format(ru), 'wb')

#cose = do_cross(allfils)
print(allfils[0])
print(allfils2[0])

do_cross(allfils[:500], allfils2[:500], gigi_a, filo)

filo.close()

# allchu = np.array_split(allfils, n_proc)
#
# pool = mp.Pool(n_proc)
# cose = pool.map(do_cross, allchu, 1)

#cose = ctl.run_parallel(do_cross, n_proc, args = allchu)

#pickle.dump(cose, )
