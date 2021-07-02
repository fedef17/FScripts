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
from tunlib import gregplot_on_ax

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

if os.uname()[1] == 'hobbes':
    cart_out = '/home/fabiano/Research/lavori/BOTTINO/regdist/'
elif os.uname()[1] == 'xaru':
    cart_out = '/home/fedef/Research/lavori/BOTTINO/regdist/'

ctl.mkdir(cart_out)

filna = '/nas/BOTTINO/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/{}/{}i1p1f1/{}/{}/*nc'

allru = ['pi', 'b025', 'b050', 'b100']
allnams = ['piControl', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']

colors = ['black', 'forestgreen', 'orange', 'violet']

####################################################################################################
masfi = '/nas/BOTTINO/masks.nc'
cose = xr.load_dataset(masfi)
oce_mask = cose['RnfA.msk'].values # ocean mask: 1 over ocean, 0 over land

#miptab = 'day'
#allvars_2D = 'pr tasmin tasmax sfcWindmax'.split()
miptab = 'Amon'
allvars_2D = 'pr tas'.split()

seasons = ['DJFM', 'JJAS']

areas_big_names = ['Tro', 'NML', 'SML', 'NP', 'SP']
areas_big = [(-180, 180, 30, 60), (-180, 180, -30, 30), (-180, 180, -60, -30), (-180, 180, 60, 90), (-180, 180, -90, -60)]

areas_land_sect = np.concatenate([[(-20, 60, l1, l2), (60, 180, l1, l2), (-180, -20, l1, l2)] for l1, l2 in [(30, 60), (-30, 30), (-30, -60)]])
areas_ls_names = 'Euro NAsia NAme Afr SAsia CAme SAfr Aus SAme'.split()

#area_reg = ['Med', 'Eu']

# areadist = dict()
#
# for na, ru, col in zip(allnams, allru, colors):
# #for na, ru, col in zip(allnams2, allru2, colors2):
#     print(ru)
#     mem = 'r1'
#     if na == 'ssp585': mem = 'r4'
#
#     for varna in allvars_2D:
#         fils = glob.glob(filna.format(na, mem, miptab, varna))
#
#         var, coords, _ = ctl.read_xr(fils, select_var = varna)
#
#         for seas in seasons:
#             var_sea, dates_sea = ctl.sel_season(var, coords['dates'], seas, cut = True)
#
#             for area, anam in zip(areas_big, areas_big_names):
#                 var_area, lat_area, lon_area = ctl.sel_area(coords['lat'], coords['lon'], var_sea, area)
#
#                 omask_area, _, _ = ctl.sel_area(coords['lat'], coords['lon'], oce_mask, area)
#
#                 for cos, ver in zip(['oce', 'land'], [1, 0]):
#                     okpo = omask_area == ver
#
#                     var_ok = var_area[:, okpo]
#                     ## area average
#
#                     var_integ = np.mean(var_ok, axis = 1)
#
#                     areadist[(varna, seas, anam, cos, ru)] = var_integ
#
#
# pickle.dump(areadist, open(cart_out + 'bottino_mon_areadist.p', 'wb'))

areadist = pickle.load(open(cart_out + 'bottino_mon_areadist.p', 'rb'))

#### FIGURONA
## Excluding poles (no sense in land/oce over pole)

for varna in allvars_2D:
    for seas in seasons:
        # for pio in ['mon', 'sea']:
        pio = 'mon'
        fig, axs = plt.subplots(2, 3, figsize = (16,8))
        for ia, anam in enumerate(areas_big_names[:3]):
            for ic, cos in enumerate(['oce', 'land']):
                ax = axs[ic, ia]
                rugi = []
                for na, ru, col in zip(allnams, allru, colors):
                    if ru == 'pi':
                        gigi = areadist[(varna, seas, anam, cos, ru)]-np.mean(areadist[(varna, seas, anam, cos, ru)])
                    else:
                        gigi = areadist[(varna, seas, anam, cos, ru)][-800:]-np.mean(areadist[(varna, seas, anam, cos, ru)][-800:])

                    if pio == 'sea':
                        gigi = gigi.reshape(-1, 4).mean(axis = 1)

                    rugi.append(gigi)

                rucaz = np.concatenate(rugi)
                rumin, rumax = (np.min(rucaz), np.max(rucaz))

                ruvec = np.linspace(rumin-0.1*abs(rumin), rumax+0.1*abs(rumax), 100)
                for ru, gi, col in zip(allru, rugi, colors):
                    pdf = ctl.calc_pdf(gi)
                    pdfvec = pdf(ruvec)
                    pdfvec = pdfvec/np.sum(pdfvec)
                    ax.plot(ruvec, pdfvec, color = col)

                ax.grid()
                ax.set_title('{} - {}'.format(anam, cos))

        fig.suptitle('{} - {}'.format(varna, seas))
        ctl.custom_legend(fig, colors, allru, ncol = 4, add_space_below = 0.06)

        fig.savefig(cart_out + '{}_{}_{}dist.pdf'.format(varna, seas, pio))
