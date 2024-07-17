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

if os.uname()[1] == 'hobbes':
    cart_out = '/home/fabiano/Research/lavori/BOTTINO/regdist/'
elif os.uname()[1] == 'tintin':
    cart_out = '/home/fabiano/work/lavori/BOTTINO/regdist/'
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

seasons = ['DJFM', 'JJAS', 'year']

areas_big_names = ['NML', 'Tro', 'SML', 'NP', 'SP']
areas_big = [(-180, 180, 30, 60), (-180, 180, -30, 30), (-180, 180, -60, -30), (-180, 180, 60, 90), (-180, 180, -90, -60)]

areas_land_sect = np.concatenate([[(-20, 60, l1, l2), (60, 180, l1, l2), (-180, -20, l1, l2)] for l1, l2 in [(60, 80), (40, 60), (20, 40), (-20, 20), (-40, -20), (-60, -40)]])
areas_ls_names = 'NEuro Sib Can Euro CAsia NAme NAfr SAsia CNAme CAfr Indo CSAme SAfr Aus SAme Soce1 Soce2 Pata'.split()

areas_land_sect = list(areas_land_sect)

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
#         for pio, gipio in zip(['mon', 'sea'], [800, 200]):
#             for seas in seasons:
#                 if pio == 'mon' and seas == 'year': continue
#
#                 if pio == 'mon':
#                     var_sea, dates_sea = ctl.sel_season(var, coords['dates'], seas, cut = True)
#                 elif pio == 'sea':
#                     var_sea, dates_sea = ctl.seasonal_set(var, coords['dates'], seas, cut = True, seasonal_average = True)
#
#                 for area, anam in zip(areas_big+areas_land_sect, areas_big_names+areas_ls_names):
#                     var_area, lat_area, lon_area = ctl.sel_area(coords['lat'], coords['lon'], var_sea, area)
#
#                     omask_area, _, _ = ctl.sel_area(coords['lat'], coords['lon'], oce_mask, area)
#
#                     for cos, ver in zip(['oce', 'land'], [1, 0]):
#                         okpo = omask_area == ver
#
#                         var_ok = var_area[:, okpo]
#                         ## area average
#
#                         var_integ = np.mean(var_ok, axis = 1)
#
#                         areadist[(varna, pio, seas, anam, cos, ru)] = var_integ
#
#                         #### Some measure of the spatial variability. Trying this: consider the relative seasonal shift from the climatology. Relative in terms of what? osjoajoa. NO. Two possibilities:
#                         # 1 - Do this on smaller boxes, but always for area-averaged anomalies.
#                         # 2 - Do this on a map, measuring a ratio of the pdfs extremes (eg p90-p10), or a stat test of the shift of the pdf wrt pi.
#
#
# pickle.dump(areadist, open(cart_out + 'bottino_monsea_areadist.p', 'wb'))

areadist = pickle.load(open(cart_out + 'bottino_monsea_areadist.p', 'rb'))

#### FIGURONA
## Excluding poles (no sense in land/oce over pole)
for typet in ['rel', 'abs']:
    for varna in allvars_2D:
        for seas in seasons:
            for pio, gipio in zip(['mon', 'sea'], [800, 200]):
                if pio == 'mon' and seas == 'year': continue
                fig, axs = plt.subplots(2, 3, figsize = (16,8))
                for ia, anam in enumerate(areas_big_names[:3]):
                    for ic, cos in enumerate(['oce', 'land']):
                        ax = axs[ic, ia]
                        rugi = []
                        for na, ru, col in zip(allnams, allru, colors):
                            if typet == 'rel':
                                if ru == 'pi':
                                    gigi = areadist[(varna, pio, seas, anam, cos, ru)]-np.mean(areadist[(varna, pio, seas, anam, cos, ru)])
                                else:
                                    gigi = areadist[(varna, pio, seas, anam, cos, ru)][-gipio:]-np.mean(areadist[(varna, pio, seas, anam, cos, ru)][-gipio:])
                            elif typet == 'abs':
                                if ru == 'pi':
                                    gigi = areadist[(varna, pio, seas, anam, cos, ru)]
                                else:
                                    gigi = areadist[(varna, pio, seas, anam, cos, ru)][-gipio:]

                            rugi.append(gigi)

                        rucaz = np.concatenate(rugi)
                        rumin, rumax = (np.min(rucaz), np.max(rucaz))

                        ruvec = np.linspace(rumin-0.2*abs(rumin), rumax+0.2*abs(rumax), 100)
                        for ru, gi, col in zip(allru, rugi, colors):
                            pdf = ctl.calc_pdf(gi)
                            pdfvec = pdf(ruvec)
                            pdfvec = pdfvec/np.sum(pdfvec)
                            ax.plot(ruvec, pdfvec, color = col)

                        ax.grid()
                        ax.set_title('{} - {}'.format(anam, cos))

                fig.suptitle('{} - {}'.format(varna, seas))
                ctl.custom_legend(fig, colors, allru, ncol = 4, add_space_below = 0.06)

                fig.savefig(cart_out + '{}_{}_{}dist_{}.pdf'.format(varna, seas, pio, typet))

    #### Now with smaller regions
    for varna in allvars_2D:
        for seas in seasons:
            for pio, gipio in zip(['mon', 'sea'], [800, 200]):
                if pio == 'mon' and seas == 'year': continue

                fig, axs = plt.subplots(3, 5, figsize = (24,12))
                for ia, anam in enumerate(areas_ls_names[:-3]):
                    ax = axs.T.flatten()[ia]
                    rugi = []
                    for na, ru, col in zip(allnams, allru, colors):
                        if typet == 'rel':
                            if ru == 'pi':
                                gigi = areadist[(varna, pio, seas, anam, cos, ru)]-np.mean(areadist[(varna, pio, seas, anam, cos, ru)])
                            else:
                                gigi = areadist[(varna, pio, seas, anam, cos, ru)][-gipio:]-np.mean(areadist[(varna, pio, seas, anam, cos, ru)][-gipio:])
                        elif typet == 'abs':
                            if ru == 'pi':
                                gigi = areadist[(varna, pio, seas, anam, cos, ru)]
                            else:
                                gigi = areadist[(varna, pio, seas, anam, cos, ru)][-gipio:]

                        rugi.append(gigi)

                    rucaz = np.concatenate(rugi)
                    rumin, rumax = (np.min(rucaz), np.max(rucaz))

                    ruvec = np.linspace(rumin-0.2*abs(rumin), rumax+0.2*abs(rumax), 100)
                    for ru, gi, col in zip(allru, rugi, colors):
                        pdf = ctl.calc_pdf(gi)
                        pdfvec = pdf(ruvec)
                        pdfvec = pdfvec/np.sum(pdfvec)
                        ax.plot(ruvec, pdfvec, color = col)

                    ax.grid()
                    ax.set_title('{} - {}'.format(anam, cos))

                fig.suptitle('{} - {}'.format(varna, seas))
                ctl.custom_legend(fig, colors, allru, ncol = 4, add_space_below = 0.06)

                fig.savefig(cart_out + 'smalreg_{}_{}_{}dist_{}.pdf'.format(varna, seas, pio, typet))
