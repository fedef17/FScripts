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

expnam = 'GLOBO_KM156L70'
#expnam = 'GLOBO_KM078_1YEAR'
cart = '/home/{}/Research/lavori/globone/test_output/{}/'.format(os.getlogin(), expnam)
ctl.mkdir(cart)

mask_file = '/home/{}/Research/lavori/globone/grid_maker/masks_globo_514x362.nc'.format(os.getlogin())
masks = xr.load_dataset(mask_file)
okmask = np.array(masks['glog.msk'])[:, 1:-1]

#######################################################################

#init = 45
init = 1

avfld = dict()

# surface fluxes
rad_flds = dict()
for var in 'chflux clwfl cqflux cswfl'.split():
    var_name, nlon, nlat, cose = ctl.read_globo_plotout(cart + var)
    #rad_flds[var] = cose/86400.
    rad_flds[var] = cose[init:]/86400.
    print(nlon, nlat)

lats = np.linspace(-90, 90, nlat)
lons = np.linspace(0., 360., nlon+1)[:-1]

## regrid mask to right shape
lats512 = np.linspace(-90, 90, 362)
lons512 = np.linspace(0., 360., 512+1)[:-1]
pino512 = ctl.create_xr_grid(lats512, lons512)
okpino = pino512.assign(dict([('okmask', (('lat', 'lon'), okmask))]))
numask = ctl.regrid_dataset(okpino, lats, lons)
okmask = np.array(numask['okmask']).astype(bool)

rad_flds['net_srf'] = rad_flds['chflux']
for var in 'clwfl cqflux cswfl'.split():
    rad_flds['net_srf'] += rad_flds[var]

#ctl.plot_multimap_contour(rad_flds['net_srf'], lats, lons, plot_anomalies=False, color_percentiles=(1,99), title='SRF_NET', cmap='viridis', plot_type='pcolormesh')

ok_coso = np.mean(rad_flds['net_srf'][1:], axis = 0)
avfld[(expnam, 'net_srf')] = ok_coso
ctl.plot_map_contour(np.mean(rad_flds['net_srf'][1:], axis = 0), lats, lons, plot_anomalies=True, color_percentiles=(1,99), title='SRF_NET', cmap='viridis', plot_type='pcolormesh', filename = cart + 'map_net_srf.pdf')

print('NET SRF', ctl.global_mean(ok_coso, lats))

## over land and ocean
ok_coso = avfld[(expnam, 'net_srf')]
ok_land = ctl.global_mean(ok_coso, lats, okmask)
ok_ocean = ctl.global_mean(ok_coso, lats, ~okmask)
print('NET SRF - LAND', ok_land)
print('NET SRF - OCEAN', ok_ocean)


for var in 'osrtotc olrtotc'.split():
    var_name, nlon, nlat, cose = ctl.read_globo_plotout(cart + var)
    rad_flds[var] = cose[init:]/86400.
    print(nlon, nlat)

rad_flds['net_toa'] = rad_flds['osrtotc'] + rad_flds['olrtotc']

#ctl.plot_multimap_contour(rad_flds['net_toa'], lats, lons, plot_anomalies=False, color_percentiles=(1,99), title='TOA_NET', cmap='viridis', plot_type='pcolormesh')

ok_coso = np.mean(rad_flds['net_toa'][1:], axis = 0)
avfld[(expnam, 'net_toa')] = ok_coso
ctl.plot_map_contour(ok_coso, lats, lons, plot_anomalies=False, color_percentiles=(1,99), title='TOA_NET', cmap='viridis', plot_type='pcolormesh', filename = cart + 'map_toa_net.pdf')

print('NET TOA', ctl.global_mean(ok_coso, lats))

## over land and ocean
ok_coso = avfld[(expnam, 'net_toa')]
ok_land = ctl.global_mean(ok_coso, lats, okmask)
ok_ocean = ctl.global_mean(ok_coso, lats, ~okmask)
print('NET TOA - LAND', ok_land)
print('NET TOA - OCEAN', ok_ocean)

rad_flds['net_atm'] = rad_flds['net_toa']-rad_flds['net_srf']

ok_coso = np.mean(rad_flds['net_atm'][1:], axis = 0)
avfld[(expnam, 'net_atm')] = ok_coso
ctl.plot_map_contour(ok_coso, lats, lons, plot_anomalies=True, color_percentiles=(1,99), title='NET_ATM', plot_type='pcolormesh', filename = cart + 'map_net_atm.pdf')

print('NET ATM', ctl.global_mean(ok_coso, lats))

## over land and ocean
ok_coso = avfld[(expnam, 'net_atm')]
ok_land = ctl.global_mean(ok_coso, lats, okmask)
ok_ocean = ctl.global_mean(ok_coso, lats, ~okmask)
print('NET ATM - LAND', ok_land)
print('NET ATM - OCEAN', ok_ocean)

#####################################

for var in 'qf_accum totprec'.split():
    var_name, nlon, nlat, cose = ctl.read_globo_plotout(cart + var)
    rad_flds[var] = cose[init:]#/86400.
    print(nlon, nlat)

rad_flds['p-e'] = rad_flds['totprec'] + rad_flds['qf_accum']

ok_coso = np.mean(rad_flds['p-e'][1:], axis = 0)
avfld[(expnam, 'p-e')] = ok_coso
ctl.plot_map_contour(ok_coso, lats, lons, plot_anomalies=False, cbar_range = (-7., 7.), title='P-E', plot_type='pcolormesh', filename = cart + 'map_p-e.pdf')


print('P-E', ctl.global_mean(ok_coso, lats))

for nam in rad_flds:
    if (expnam, nam) not in avfld:
        avfld[(expnam, nam)] = np.mean(rad_flds[nam][1:], axis = 0)

        print(nam, '{:f6.2}'.format(ctl.global_mean(avfld[(expnam, nam)], lats)))

        ctl.plot_map_contour(avfld[(expnam, nam)], lats, lons, plot_anomalies=False, color_percentiles=(1,99), title=nam, plot_type='pcolormesh', cmap = 'gist_heat', filename = cart + 'map_{}.pdf'.format(nam))

pickle.dump(avfld, open(cart + 'avfld.p', 'wb'))

##### stampo timeseries globali
# lowpass di 30 giorni

filt = np.where(np.any(rad_flds['cswfl'] != 0., axis = (1,2)))[0]

fig = plt.figure()
for var in 'chflux clwfl cqflux cswfl net_srf'.split():
    #pino = ctl.running_mean(ctl.global_mean(rad_flds[var][filt], lats), 30)
    pino = ctl.global_mean(rad_flds[var][filt], lats)
    plt.plot(pino, label = var)
plt.legend()
plt.grid()
plt.ylabel('W/m2')
plt.title('surface fluxes in time')
fig.savefig(cart + 'time_evmon_srf.pdf')

print(ctl.global_mean(np.mean(rad_flds['net_srf'][filt], axis = 0), lats))

fig = plt.figure()
for var in 'olrtotc osrtotc net_toa'.split():
    #pino = ctl.running_mean(ctl.global_mean(rad_flds[var], lats), 30)
    pino = ctl.global_mean(rad_flds[var][filt], lats)
    plt.plot(pino, label = var)

plt.legend()
plt.grid()
plt.ylabel('W/m2')
plt.title('TOA fluxes in time')
fig.savefig(cart + 'time_evmon_toa.pdf')


fig = plt.figure()
for var in 'qf_accum totprec p-e'.split():
    #pino = ctl.running_mean(ctl.global_mean(rad_flds[var], lats), 30)
    pino = ctl.global_mean(rad_flds[var][filt], lats)
    plt.plot(pino, label = var)

plt.legend()
plt.grid()
plt.ylabel('mm/day')
plt.title('evap/prec in time')
fig.savefig(cart + 'time_evmon_prec.pdf')


###### zonal plots

fig = plt.figure()
for var in 'chflux clwfl cqflux cswfl net_srf'.split():
    pino = ctl.zonal_mean(np.mean(rad_flds[var], axis = 0))
    plt.plot(lats, pino, label = var)

plt.legend()
plt.grid()
plt.ylabel('W/m2')
plt.title('surface fluxes - zonal mean')
fig.savefig(cart + 'zonal_srf.pdf')

ok_coso = np.mean(rad_flds['net_toa'][1:], axis = 0)

fig = plt.figure()
for var in 'olrtotc osrtotc net_toa'.split():
    pino = ctl.zonal_mean(np.mean(rad_flds[var], axis = 0))
    plt.plot(lats, pino, label = var)

plt.legend()
plt.grid()
plt.ylabel('W/m2')
plt.title('TOA fluxes - zonal mean')
fig.savefig(cart + 'zonal_toa.pdf')


fig = plt.figure()
for var in 'qf_accum totprec p-e'.split():
    pino = ctl.zonal_mean(np.mean(rad_flds[var], axis = 0))
    plt.plot(lats, pino, label = var)

plt.legend()
plt.grid()
plt.ylabel('mm/day')
plt.title('evap/prec - zonal mean')
fig.savefig(cart + 'zonal_prec.pdf')

#######

# load ece net surface shortwave for comparison
cart_ece = '/home/fedef/Research/lavori/globone/test_output/ref_ece/'
cosec = xr.load_dataset(cart_ece + 'fml1_clim_ssr.nc')

cosec_ok = ctl.regrid_dataset(cosec, lats, lons)
ssr_ece = np.array(np.mean(cosec_ok['ssr'], axis = 0))

ctl.plot_map_contour(avfld[(expnam, 'cswfl')]-ssr_ece, lats, lons, plot_anomalies=True, color_percentiles=(1,99), title='shortwave net surface: diff GLOBO - ECE', plot_type='pcolormesh', cmap = 'RdBu_r', filename = cart + 'diff_ECE_map_{}.pdf'.format('cswfl'))

# load ece net toa shortwave for comparison
cart_ece = '/home/fedef/Research/lavori/globone/test_output/ref_ece/'
cosec = xr.load_dataset(cart_ece + 'fml1_clim_tsr.nc')

cosec_ok = ctl.regrid_dataset(cosec, lats, lons)
ssr_ece = np.array(np.mean(cosec_ok['tsr'], axis = 0))

ctl.plot_map_contour(avfld[(expnam, 'osrtotc')]-ssr_ece, lats, lons, plot_anomalies=True, color_percentiles=(1,99), title='shortwave net TOA: diff GLOBO - ECE', plot_type='pcolormesh', cmap = 'RdBu_r', filename = cart + 'diff_ECE_map_{}.pdf'.format('osrtotc'))


# [1]   Terminated              ncview fml1_2013_snr_map.nc
# [2]-  Terminated              ncview fml1_2013_tnr_map.nc
# [3]+  Terminated              ncview fml1_2013_net_atm_map.nc
# (base) fabiano@hobbes:/data-hobbes/fabiano/TunECS/AMIP_worlds/fml1/post/mon/Post_2013
