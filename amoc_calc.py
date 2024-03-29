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
allru = ['b990', 'b025', 'b050', 'b100', 'b065', 'b080']#, 'b00I', 'b80I', 'b65I']
miptab = 'Omon'
var = 'msftyz'

amoc_all = dict()

for ru in allru:
    print(ru)
    if ru in ['b025', 'b050', 'b100', 'b080']:
        datadir = '/g100_scratch/userexternal/ffabiano/ece3/{}/cmorized/'.format(ru)
        filna = datadir+'cmor_*/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/*/r1i1p1{}/{}/{}/g*/v*/{}*nc'.format(fis, miptab, var, var)
    else:
        #datadir = '/g100_work/IscrB_QUECLIM/BOTTINO/{}/cmorized/'.format(ru)
        datadir = '/g100_scratch/userexternal/ffabiano/irods_data/{}/'.format(ru)
        filna = datadir+'{}/{}/{}*nc'.format(miptab, var, var)

    fis = 'f1'
    if 'I' in ru: fis = 'f3'

    listafi = glob.glob(filna)
    listafi.sort()

    amoc_max = []
    amoc_max_lev = []
    amoc_wid = []
    aabw_max = []
    aabw_max_lev = []
    aby_max = []
    aby_max_lev = []
    smoc_max = []
    smoc_max_lev = []

    for fi in listafi:
        print(fi)
        coso = xr.load_dataset(fi, use_cftime = True)[var]
        amax = coso.sel(basin = 1, rlat = slice(30, 50), lev = slice(500., 2000.)).mean('time').max(['rlat', 'lev']).values # basin = 1 should be Atlantic, 0 global, 2 indian/pacific
        amoc_max.append(amax)

        zuki = coso.sel(basin = 1, rlat = slice(30, 50), lev = slice(500., 2000.)).mean('time').argmax(['rlat', 'lev'])
        maxlev = coso.sel(lev = slice(500., 2000.)).lev[zuki['lev']].values
        amoc_max_lev.append(maxlev)

        gino = coso.sel(basin = 1, rlat = slice(30, 50)).mean('time')
        zino = np.all(gino.values < amax/2., axis = 1)
        for i, lev in enumerate(coso.lev):
            if np.all(zino[i:]):
                awid = lev.values
                break

        amoc_wid.append(awid)
        print(amax, awid)

        ### AABW
        amax = coso.sel(basin = 0, lev = slice(500., 3000.), rlat = slice(-90, -60)).mean('time').min(['rlat', 'lev']).values
        aabw_max.append(amax)

        zuki = coso.sel(basin = 0, lev = slice(500., 3000.), rlat = slice(-90, -60)).mean('time').argmin(['rlat', 'lev'])
        maxlev = coso.sel(lev = slice(500., 3000.)).lev[zuki['lev']].values
        aabw_max_lev.append(maxlev)

        ### SMOC
        # amax = coso.sel(basin = 0, lev = slice(500., 3000.), rlat = slice(-90, -40)).mean('time').max(['rlat', 'lev']).values
        #amax = coso.sel(basin = 0, lev = slice(2000., None), rlat = slice(-90, -35)).mean('time').min(['rlat', 'lev']).values
        amax = coso.sel(basin = 0, lev = slice(3000., 4000.), rlat = slice(-50, -30)).mean(['time', 'lev']).min('rlat').values
        smoc_max.append(amax)

        zuki = coso.sel(basin = 0, lev = slice(3000., 4500.), rlat = slice(-50, -30)).mean('time').argmin(['rlat', 'lev'])
        maxlev = coso.sel(lev = slice(3000., 4500.)).lev[zuki['lev']].values
        smoc_max_lev.append(maxlev)

        ### abyssal amoc
        amax = coso.sel(basin = 0, lev = slice(2500., 4500.)).mean('time').min(['rlat', 'lev']).values
        aby_max.append(amax)

        zuki = coso.sel(basin = 0, lev = slice(2500., 4500.)).mean('time').argmin(['rlat', 'lev'])
        maxlev = coso.sel(lev = slice(2500., 4500.)).lev[zuki['lev']].values
        aby_max_lev.append(maxlev)

    amoc_all[(ru, 'amoc_max')] = np.stack(amoc_max)
    amoc_all[(ru, 'amoc_wid')] = np.stack(amoc_wid)
    amoc_all[(ru, 'amoc_maxlev')] = np.stack(amoc_max_lev)
    amoc_all[(ru, 'aabw_max')] = np.stack(aabw_max)
    amoc_all[(ru, 'aabw_maxlev')] = np.stack(aabw_max_lev)
    amoc_all[(ru, 'aby_max')] = np.stack(aby_max)
    amoc_all[(ru, 'aby_maxlev')] = np.stack(aby_max_lev)
    amoc_all[(ru, 'smoc_max')] = np.stack(smoc_max)
    amoc_all[(ru, 'smoc_maxlev')] = np.stack(smoc_max_lev)

pickle.dump(amoc_all, open(cart_out + 'amoc_all_1000.p', 'wb'))

#######
sys.exit()

amoc_2D = dict()
n_yea = 30

for ru in allru:
    print(ru)
    if ru in ['b025', 'b050', 'b100', 'b080']:
        datadir = '/g100_scratch/userexternal/ffabiano/ece3/{}/cmorized/'.format(ru)
        filna = datadir+'cmor_*/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/*/r1i1p1{}/{}/{}/g*/v*/{}*nc'.format(fis, miptab, var, var)
    else:
        #datadir = '/g100_work/IscrB_QUECLIM/BOTTINO/{}/cmorized/'.format(ru)
        datadir = '/g100_scratch/userexternal/ffabiano/irods_data/{}/'.format(ru)
        filna = datadir+'{}/{}/{}*nc'.format(miptab, var, var)

    fis = 'f1'
    if 'I' in ru: fis = 'f3'

    listafi = glob.glob(filna)
    listafi.sort()

    fi_ini = listafi[:n_yea]
    print(fi_ini[-1])
    coso_ini = xr.open_mfdataset(fi_ini, use_cftime = True)[var]
    mocstr_ini = coso_ini.mean('time')

    fi_fin = listafi[-n_yea:]
    print(fi_fin[0])
    coso_fin = xr.open_mfdataset(fi_fin, use_cftime = True)[var]
    mocstr_fin = coso_fin.mean('time')

    amoc_2D[(ru, 'ini')] = mocstr_ini.compute()
    amoc_2D[(ru, 'fin')] = mocstr_fin.compute()

pickle.dump(amoc_2D, open(cart_out + 'amoc_2D_1000.p', 'wb'))