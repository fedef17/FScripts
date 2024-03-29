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

# plt.rcParams['xtick.labelsize'] = 15
# plt.rcParams['ytick.labelsize'] = 15
# titlefont = 22
# plt.rcParams['figure.titlesize'] = titlefont
# plt.rcParams['axes.titlesize'] = 18
# plt.rcParams['axes.labelsize'] = 15
# plt.rcParams['axes.axisbelow'] = True

#############################################################################

#ru = sys.argv[1]
ru = 'piControl'

# open our log file
logname = 'log_oceall_{}.log'.format(ru)
logfile = open(logname,'w') #self.name, 'w', 0)

# re-open stdout without buffering
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1)

# redirect stdout and stderr to the log file opened above
os.dup2(logfile.fileno(), sys.stdout.fileno())
os.dup2(logfile.fileno(), sys.stderr.fileno())

print('total RAM memory used 0:', psutil.virtual_memory()[3]/1.e9)

if os.uname()[1] == 'hobbes':
    cart_out = '/home/fabiano/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'xaru':
    cart_out = '/home/fedef/Research/lavori/BOTTINO/'
elif os.uname()[1] == 'tintin':
    cart_out = '/home/fabiano/work/lavori/BOTTINO/'

#cart_out = '/g100_work/IscrB_QUECLIM/BOTTINO/bottino_an/'
cart_out = cart_out + 'ocean3d/'
ctl.mkdir(cart_out)

filna = '/nas/BOTTINO/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/{}/{}i1p1f1/Omon/{}/{}*nc'

areafil = '/nas/BOTTINO/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/stabilization-ssp585-2025/r1i1p1f1/Ofx/areacello/areacello_Ofx_EC-Earth3_stabilization-ssp585-2025_r1i1p1f1_gn.nc'

####################################################################################################

miptab = 'Omon'
mem = 'r1'
#var = 'thetao'
var = 'bigthetao'
mvar = 'masscello'

basnames = ['atlmsk', 'indmsk', 'pacmsk']
subbas = xr.load_dataset(cart_out + 'subbasins.nc')

lats = np.linspace(-89.5, 89.5, 180)
lons = np.linspace(0, 359, 360)

yeamean = dict()

print('total RAM memory used 0:', psutil.virtual_memory()[3]/1.e9)

def do_cross(fils, fils2, fia, fil_out):#, coda):
    print("I'm process", os.getpid())
    # Getting % usage of virtual_memory ( 3rd field)
    gigi_a = xr.load_dataset(fia, use_cftime = True)['areacello']

    #cose = []
    for fi1, fi2 in zip(fils, fils2):
        print(fi1)
        print('total RAM memory used 1:', psutil.virtual_memory()[3]/1.e9)

        gigi = xr.load_dataset(fi1, use_cftime = True)['bigthetao']
        gigi2 = xr.load_dataset(fi2, use_cftime = True)['masscello']

        oht = gigi*gigi2*gigi_a.values[np.newaxis, np.newaxis, ...]
        # for the pattern, do not multiply for area!!

        # areacello is m2
        # masscello is kg/m2
        # bigthetao is C
        # to have real oht must multiply by cp

        oht_lev = oht.mean('time').sum(['i', 'j'])

        # oht700 = oht.sel(lev = slice(0., 700.)).mean('time').sum('lev')
        # oht2000 = oht.sel(lev = slice(700., 2000.)).mean('time').sum('lev')
        # oht_deep = oht.sel(lev = slice(2000., 6000.)).mean('time').sum('lev')

        # oht700 = gigi.sel(lev = slice(0., 700.)).mean('time').mean('lev')
        # oht2000 = gigi.sel(lev = slice(700., 2000.)).mean('time').mean('lev')
        # oht_deep = gigi.sel(lev = slice(2000., 6000.)).mean('time').mean('lev')

        oht700 = gigi.sel(lev = slice(697., 698.)).mean('time').mean('lev')
        oht2000 = gigi.sel(lev = slice(1945., 1946.)).mean('time').mean('lev')
        oht_deep = gigi.sel(lev = slice(4093., 4094.)).mean('time').mean('lev')

        print('nans', np.sum(np.isnan(oht_deep.values)))

        zuki700 = ctl.regrid_dataset(oht700, lats, lons)
        zuki2000 = ctl.regrid_dataset(oht2000, lats, lons)
        zuki_deep = ctl.regrid_dataset(oht_deep, lats, lons)

        print(zuki_deep)
        print(zuki_deep.max(), zuki_deep.min(), zuki_deep.mean())
        print('nans', np.sum(np.isnan(zuki_deep.values)))

        # nuvarz = dict()
        # for basnam in basnames:
        #     pinzu = np.array(subbas[basnam])
        #     pinzu[pinzu == 0] = np.nan
        #
        #     goggolo = gigi['thetao']*pinzu[np.newaxis, np.newaxis, ...]
        #     nuvarz['thetao_'+basnam[:3]] = goggolo
        # gigi = gigi.assign(nuvarz)

        print('total RAM memory used 2:', psutil.virtual_memory()[3]/1.e9)
        #
        # del nuvarz, goggolo
        #gigi = gigi.drop(['vertices_longitude', 'vertices_latitude'])

        # zuki = ctl.regrid_dataset(oht, lats, lons)
        # #gigi = gigi.drop_dims('vertices')
        #
        # #### Now the zonal cross section
        # gogcross = zuki.mean(('time', 'lon'))
        # #cose.append(gogcross)
        # pickle.dump(gogcross, fil_out)
        pickle.dump([oht_lev, zuki700, zuki2000, zuki_deep], fil_out)

        print('total RAM memory used 3:', psutil.virtual_memory()[3]/1.e9)

        fil_out.flush()

        del gigi, gigi2, oht, oht_lev, zuki700, zuki2000, zuki_deep, oht700, oht2000, oht_deep

    #coda.put(cose)
    #return cose
    return

n_proc = 10
#for ru, nam in zip(allru, allnams):

print(ru)
#nam = allnams[allru.index(ru)]
nam = ru

print(filna.format(nam, mem, var, var))
allfils = glob.glob(filna.format(nam, mem, var, var))
allfils.sort()

allfils2 = glob.glob(filna.format(nam, mem, mvar, mvar))
allfils2.sort()

print(allfils)
print(allfils2)

# allfils_a = glob.glob(filna2.format(ru, nam, mem))
# allfils_a.sort()


filo = open(cart_out + 'oht_{}.p'.format(ru), 'wb')

#cose = do_cross(allfils)

do_cross(allfils, allfils2, areafil, filo)

filo.close()

# allchu = np.array_split(allfils, n_proc)
#
# pool = mp.Pool(n_proc)
# cose = pool.map(do_cross, allchu, 1)

#cose = ctl.run_parallel(do_cross, n_proc, args = allchu)

#pickle.dump(cose, )
