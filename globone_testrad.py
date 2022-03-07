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
cart = '/home/fedef/Research/lavori/globone/test_output/'
# rsync -auv --progress -e "ssh" fabiano@tintin.bo.isac.cnr.it:/work/users/davini/globone/runtime/GLOBO_KM078_TEST/o* .

# surface fluxes
rad_flds = dict()
for var in 'chflux clwfl cqflux cswfl'.split():
    var_name, nlon, nlat, cose = ctl.read_globo_plotout(cart + var)
    rad_flds[var] = cose/86400.
    print(nlon, nlat)

lats = np.linspace(-90, 90, nlat)
lons = np.linspace(0., 360., nlon+1)[:-1]

rad_flds['net_srf'] = rad_flds['chflux']
for var in 'clwfl cqflux cswfl'.split():
    rad_flds['net_srf'] += rad_flds[var]

#ctl.plot_multimap_contour(rad_flds['net_srf'], lats, lons, plot_anomalies=False, color_percentiles=(1,99), title='SRF_NET', cmap='viridis', plot_type='pcolormesh')

ok_coso = np.mean(rad_flds['net_srf'][1:], axis = 0)
ctl.plot_map_contour(np.mean(rad_flds['net_srf'][1:], axis = 0), lats, lons, plot_anomalies=True, color_percentiles=(1,99), title='SRF_NET', cmap='viridis', plot_type='pcolormesh', filename = cart + 'map_net_srf.pdf')

print('NET SRF', ctl.global_mean(ok_coso, lats))

for var in 'osrtotc olrtotc'.split():
    var_name, nlon, nlat, cose = ctl.read_globo_plotout(cart + var)
    rad_flds[var] = cose/86400.
    print(nlon, nlat)

rad_flds['net_toa'] = rad_flds['osrtotc'] + rad_flds['olrtotc']

#ctl.plot_multimap_contour(rad_flds['net_toa'], lats, lons, plot_anomalies=False, color_percentiles=(1,99), title='TOA_NET', cmap='viridis', plot_type='pcolormesh')

ok_coso = np.mean(rad_flds['net_toa'][1:], axis = 0)
ctl.plot_map_contour(ok_coso, lats, lons, plot_anomalies=False, color_percentiles=(1,99), title='TOA_NET', cmap='viridis', plot_type='pcolormesh', filename = cart + 'map_toa_net.pdf')

print('NET TOA', ctl.global_mean(ok_coso, lats))

rad_flds['net_atm'] = rad_flds['net_toa']-rad_flds['net_srf']

ok_coso = np.mean(rad_flds['net_atm'][1:], axis = 0)
ctl.plot_map_contour(ok_coso, lats, lons, plot_anomalies=True, color_percentiles=(1,99), title='NET_ATM', plot_type='pcolormesh', filename = cart + 'map_net_atm.pdf')

print('NET ATM', ctl.global_mean(ok_coso, lats))

#####################################

for var in 'qf_accum totprec'.split():
    var_name, nlon, nlat, cose = ctl.read_globo_plotout(cart + var)
    rad_flds[var] = cose#/86400.
    print(nlon, nlat)

rad_flds['p-e'] = rad_flds['totprec'] + rad_flds['qf_accum']

ok_coso = np.mean(rad_flds['p-e'][1:], axis = 0)
ctl.plot_map_contour(ok_coso, lats, lons, plot_anomalies=False, cbar_range = (-7., 7.), title='P-E', plot_type='pcolormesh', filename = cart + 'map_p-e.pdf')


print('P-E', ctl.global_mean(ok_coso, lats))


# [1]   Terminated              ncview fml1_2013_snr_map.nc
# [2]-  Terminated              ncview fml1_2013_tnr_map.nc
# [3]+  Terminated              ncview fml1_2013_net_atm_map.nc
# (base) fabiano@hobbes:/data-hobbes/fabiano/TunECS/AMIP_worlds/fml1/post/mon/Post_2013
