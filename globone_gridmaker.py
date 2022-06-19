#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
from matplotlib import pyplot as plt
from matplotlib import cm

import pickle
import netCDF4 as nc

import climtools_lib as ctl
import climdiags as cd
#from tunlib import gregplot_on_ax

from matplotlib.colors import LogNorm
from datetime import datetime

from scipy import stats
import xarray as xr
import glob
import xclim

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 22
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['axes.axisbelow'] = True

#############################################################################
cart = '/home/fedef/Research/lavori/globone/grid_maker/'

all_res = ['KM078', 'KM156', 'KM312']
n_latlon = [(362, 514), (182, 258), (92, 130)]

for res, nlalo in zip(all_res, n_latlon):
    print(res, nlalo)
    cart_in = cart + 'globo_{}/'.format(res)

    nlat, nlon = nlalo

    fmask = ctl.read_globo_grid(cart_in + 'mask_globo.txt', nlon, nlat)
    latg = ctl.read_globo_grid(cart_in + 'latg_globo.txt', nlon, nlat)
    long = ctl.read_globo_grid(cart_in + 'long_globo.txt', nlon, nlat)

    fmask = fmask.squeeze()
    latg = latg.squeeze()
    long = long.squeeze()

    latg = latg*180/np.pi
    long = long*180/np.pi


    # globo: 1 ocean, 0 land
    # define oasis integer mask (0 over ocean and fractional points, 1 over land)
    fmask_oasis = fmask < 0.5
    fmask_oasis = fmask_oasis.astype('int32')

    # latvert = np.stack([latg-dlat/2, latg-dlat/2, latg+dlat/2, latg+dlat/2])
    # lonvert = np.stack([long-dlon/2, long+dlon/2, long+dlon/2, long-dlon/2])

    # plt.ion()
    plt.figure()
    plt.imshow(fmask)
    plt.colorbar()

    plt.figure()
    plt.imshow(fmask_oasis)
    plt.colorbar()

    plt.figure()
    plt.imshow(latg)
    plt.colorbar()

    plt.figure()
    plt.imshow(long)
    plt.colorbar()

    grids_ece_file = '/home/fedef/Research/git/EC-Earth/ece_runtime/oasis/AMIP/grids.nc'
    masks_ece_file = '/home/fedef/Research/git/EC-Earth/ece_runtime/oasis/AMIP/masks.nc'

    grids_ece = xr.load_dataset(grids_ece_file)
    masks_ece = xr.load_dataset(masks_ece_file)

    amip_lat = np.array(grids_ece['AMIP.lat'])
    amip_lon = np.array(grids_ece['AMIP.lon'])
    amip_msk = np.array(masks_ece['AMIP.msk'])

    dims = tuple(['x_amip', 'y_amip', 'x_glog', 'y_glog'])
    coords = dict([('x_amip', np.arange(360)), ('y_amip', np.arange(180)), ('x_glog', np.arange(nlon)), ('y_glog', np.arange(nlat))])
    #coords = dict([('x_amip', np.arange(360)), ('y_amip', np.arange(180)), ('x_glog', np.arange(nlon)), ('y_glog', np.arange(nlat)), ('iver', np.arange(4)+1)])

    data_vars = dict(zip(['AMIP.lat', 'AMIP.lon', 'glog.lon', 'glog.lat'], [(('y_amip', 'x_amip'), amip_lat), (('y_amip', 'x_amip'), amip_lon), (('y_glog', 'x_glog'), long), (('y_glog', 'x_glog'), latg)]))
    #data_vars = dict(zip(['AMIP.lat', 'AMIP.lon', 'glog.lon', 'glog.lat', 'glon.cla', 'glog.clo'], [(('y_amip', 'x_amip'), amip_lat), (('y_amip', 'x_amip'), amip_lon), (('y_glog', 'x_glog'), long), (('y_glog', 'x_glog'), latg), (('y_glog', 'x_glog', 'iver'), latvert), (('y_glog', 'x_glog', 'iver'), lonvert)]))
    grids = xr.Dataset(data_vars=data_vars, coords = coords)
    grids.to_netcdf(cart + 'grids_globo_{}.nc'.format(res))

    data_vars = dict(zip(['AMIP.msk', 'glog.msk'], [(('y_amip', 'x_amip'), amip_msk), (('y_glog', 'x_glog'), fmask_oasis)]))
    masks = xr.Dataset(data_vars=data_vars, coords = coords)
    masks.to_netcdf(cart + 'masks_globo_{}.nc'.format(res))
