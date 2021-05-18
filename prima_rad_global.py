#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os

from matplotlib import pyplot as plt
import matplotlib.cm as cm
import matplotlib.patheffects as PathEffects
import matplotlib.animation as animation
from matplotlib.animation import ImageMagickFileWriter

import netCDF4 as nc
import cartopy.crs as ccrs
import cartopy.util as cutil
import pandas as pd

from numpy import linalg as LA
from eofs.standard import Eof
from scipy import stats, optimize
import itertools as itt

from sklearn.cluster import KMeans

from datetime import datetime
import pickle
import iris

import glob

import climtools_lib as ctl
import climdiags as cd

import xarray as xr
import xesmf as xe
import xclim

from importlib import reload

#######################################################

cart_in = '/work/scratch-nopw/fedef17/prima_outgoing/highresSST-present/'
cart_out = '/home/users/fedef17/work/prima_rad/'
### RLUT

vars = 'clt  hfls  hfss  pr  prsn  rlds  rlus  rlut  rsds  rsdt  rsus  rsut'.split()
varstot = vars + ['net_toa', 'net_srf', 'in_toa']

mods = os.listdir(cart_in)
mods.sort()
fil = cart_in + '{}/r1*/Amon/{}/{}*nc'

fres = open(cart_out + 'prima_rad_mean.dat', 'w')
fres_std = open(cart_out + 'prima_rad_std.dat', 'w')

kosetttit = '{:15s}' + len(varstot) * '{:10s}' + '\n'
kosett = '{:15s}' + len(varstot) * '{:10.4e}' + '\n'
fres.write(kosetttit.format('model', *varstot))
fres_std.write(kosetttit.format('model', *vars))

resu = dict()
resu_std = dict()
for mod in mods:
    print(mod)
    for var in vars:
        print(var)
        fikoj = fil.format(mod, var, var)
        filok = glob.glob(fikoj)[0]

        kos, coords, auxi = ctl.read_xr(filok)
        kosme = ctl.global_mean(kos, coords['lat'])

        resu[(mod, var)] = np.mean(kosme)
        resu_std[(mod, var)] = np.std(kosme)


    resu[(mod, 'net_toa')] = resu[(mod, 'rsdt')] - resu[(mod, 'rsut')] - resu[(mod, 'rlut')]
    resu[(mod, 'net_srf')] = resu[(mod, 'rsds')] + resu[(mod, 'rlds')] - resu[(mod, 'rsus')] - resu[(mod, 'rlus')] - resu[(mod, 'hfss')] - resu[(mod, 'hfls')]

    resu[(mod, 'in_atm')] = resu[(mod, 'net_toa')] - resu[(mod, 'net_srf')]

    resulis = [resu[(mod, var)] for var in varstot]
    fres.write(kosett.format(mod, *resulis))
    resulis = [resu[(mod, var)] for var in vars]
    fres_std.write(kosett.format(mod, *resulis))

fres.close()
fres_std.close()

pickle.dump([resu, resu_std], open(cart_out + 'glomean.p', 'wb'))
