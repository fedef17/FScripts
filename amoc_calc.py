#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os

import pickle
import climtools_lib as ctl
import climdiags as cd

import xarray as xr
import glob

#############################################################################

# open log file
logname = 'log_amoc_calc.log'
logfile = open(logname,'w')
# re-open stdout without buffering
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1)
# redirect stdout and stderr to the log file opened above
os.dup2(logfile.fileno(), sys.stdout.fileno())
os.dup2(logfile.fileno(), sys.stderr.fileno())

cart_out = '/g100_work/IscrB_QUECLIM/BOTTINO/bottino_an/'
cart_out = cart_out + 'amoc/'
ctl.mkdir(cart_out)

ro = 1.025e9
####################################################################################################
allru = ['b990', 'b025', 'b050', 'b100', 'b065', 'b080', 'b00I', 'b80I', 'b65I']
miptab = 'Omon'
var = 'msftyz'

amoc_all = dict()
for ru in allru:
    print(ru)
    if ru in ['b990', 'b050', 'b100', 'b080', 'b065']:
        datadir = '/g100_scratch/userexternal/ffabiano/ece3/{}/cmorized/'.format(ru)
    else:
        datadir = '/g100_work/IscrB_QUECLIM/BOTTINO/{}/cmorized/'.format(ru)

    fis = 'f1'
    if 'I' in ru:
        fis = 'f3'

    filna = datadir+'cmor_*/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/*/r1i1p1{}/{}/{}/g*/v*/{}*nc'.format(fis, miptab, var, var)
    listafi = glob.glob(filna)
    listafi.sort()

    amoc_max = []
    amoc_wid = []
    for fi in listafi:
        print(fi)
        coso = xr.load_dataset(fi, use_cftime = True)[var]
        amax = coso.sel(basin = 1, rlat = slice(30, 50), lev = slice(500., 2000.)).mean('time').max(['rlat', 'lev']).values # basin = 1 should be Atlantic, 0 global, 2 indian/pacific
        amoc_max.append(amax)

        gino = coso.sel(basin = 1, rlat = slice(30, 50)).mean('time')
        zino = np.all(gino.values < amax/2., axis = 1)
        for i, lev in enumerate(coso.lev):
            if np.all(zino[i:]):
                awid = lev.values
                break

        amoc_wid.append(awid)
        print(amax, awid)

    amoc_max = np.stack(amoc_max)
    amoc_wid = np.stack(amoc_wid)
    amoc_all[(ru, 'max')] = amoc_max
    amoc_all[(ru, 'wid')] = amoc_wid

pickle.dump(amoc_all, open(cart_out + 'amoc_all.p', 'wb'))
