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

ru = sys.argv[1]

# open our log file
logname = 'log_oceall_{}.log'.format(ru)
logfile = open(logname,'w') #self.name, 'w', 0)

# re-open stdout without buffering
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1)

# redirect stdout and stderr to the log file opened above
os.dup2(logfile.fileno(), sys.stdout.fileno())
os.dup2(logfile.fileno(), sys.stderr.fileno())

print('total RAM memory used 0:', psutil.virtual_memory()[3]/1.e9)
# if os.uname()[1] == 'hobbes':
#     cart_out = '/home/fabiano/Research/lavori/BOTTINO/'
# elif os.uname()[1] == 'xaru':
#     cart_out = '/home/fedef/Research/lavori/BOTTINO/'
# elif os.uname()[1] == 'tintin':
#     cart_out = '/home/fabiano/work/lavori/BOTTINO/'

cart_out = '/g100_work/IscrB_QUECLIM/BOTTINO/bottino_an/'
cart_out = cart_out + 'ocean3d/'
ctl.mkdir(cart_out)

filna = '/g100_work/IscrB_QUECLIM/BOTTINO/{}/cmorized/cmor_*/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/{}/{}i1p1f1/Omon/{}/gn/v20210315/{}*nc'
filna2 = '/g100_work/IscrB_QUECLIM/BOTTINO/{}/cmorized/cmor_*/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/{}/{}i1p1f1/Ofx/areacello/gn/v20210315/areacello*nc'
# filna = '/nas/BOTTINO/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/{}/{}i1p1f1/{}/{}/*nc'

allru = ['b025', 'b050', 'b100']#['pi',
allnams = ['stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']#'piControl',
#nam = allnams[allru.index(ru)]

colors = ['black', 'forestgreen', 'orange', 'violet']

####################################################################################################

miptab = 'Omon'
mem = 'r1'
#var = 'thetao'
var = 'bigthetao'
mvar = 'masscello'

# allnams2 = allnams + ['ssp585']
# allru2 = allru + ['ssp585']
# colors2 = colors + ['indianred']

# read subbasin.nc

basnames = ['atlmsk', 'indmsk', 'pacmsk']
subbas = xr.load_dataset(cart_out + 'subbasins.nc')

lats = np.linspace(-89.5, 89.5, 180)
lons = np.linspace(0, 359, 360)

yeamean = dict()

print('total RAM memory used 0:', psutil.virtual_memory()[3]/1.e9)

def do_cross(fils, fils2, fils_area, fil_out):#, coda):
    print("I'm process", os.getpid())
    # Getting % usage of virtual_memory ( 3rd field)

    #cose = []
    for fi1, fi2, fia in zip(fils, fils2, fils_area):
        print(fi1)
        print('total RAM memory used 1:', psutil.virtual_memory()[3]/1.e9)

        gigi = xr.load_dataset(fi1, use_cftime = True)['bigthetao']
        gigi2 = xr.load_dataset(fi2, use_cftime = True)['masscello']

        gigi_a = xr.load_dataset(fia, use_cftime = True)['areacello']

        oht = gigi*gigi2*gigi_a.values[np.newaxis, np.newaxis, ...]

        # areacello is m2
        # masscello is kg/m2
        # bigthetao is C
        # to have real oht must multiply by cp

        oht_lev = oht.mean('time').sum(['i', 'j'])

        oht700 = oht.sel(lev = slice(0., 700.)).mean('time').sum('lev')
        oht2000 = oht.sel(lev = slice(700., 2000.)).mean('time').sum('lev')
        oht_deep = oht.sel(lev = slice(2000., 6000.)).mean('time').sum('lev')

        zuki700 = ctl.regrid_dataset(oht700, lats, lons)
        zuki2000 = ctl.regrid_dataset(oht2000, lats, lons)
        zuki_deep = ctl.regrid_dataset(oht_deep, lats, lons)

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
nam = allnams[allru.index(ru)]

allfils = glob.glob(filna.format(ru, nam, mem, var, var))
allfils.sort()

allfils2 = glob.glob(filna.format(ru, nam, mem, mvar, mvar))
allfils2.sort()

allfils_a = glob.glob(filna2.format(ru, nam, mem))
allfils_a.sort()


filo = open(cart_out + 'oht_{}.p'.format(ru), 'wb')

#cose = do_cross(allfils)

do_cross(allfils, allfils2, allfils_a, filo)

filo.close()

# allchu = np.array_split(allfils, n_proc)
#
# pool = mp.Pool(n_proc)
# cose = pool.map(do_cross, allchu, 1)

#cose = ctl.run_parallel(do_cross, n_proc, args = allchu)

#pickle.dump(cose, )