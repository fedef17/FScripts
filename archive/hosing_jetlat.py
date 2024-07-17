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

from matplotlib.colors import LogNorm
from datetime import datetime

from scipy import stats
import xarray as xr
import glob
import xclim

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 28
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.axisbelow'] = True

#############################################################################

# cart = '/home/federico/work/Tipes/tipes_hosing/'
#
# fils = ['uaday_cntrl_1850-1999_850hPa_remap25.nc', 'uaday_ho03_1850-1899_850hPa_remap25.nc']#, 'uaday_c3r5_1900-1999_850hPa_remap25.nc']
#
allru = ['pi', 'ho03', 'ho03b']#, 'c3r5']
colors = ['steelblue', 'indianred', 'orange']#, 'forestgreen']

# cart_out = '/home/fabiano/work/lavori/hosing/'
#
# cart = '/home/bellomo/work/ec_hosing/ho03/day/r2/'
# filone = 'ua_day_ho03b_1850-1899.nc'
#
# resdict = dict()
# jli, jspeed, dates = cd.jli_from_files(cart + filone, orogfile = '/home/federico/work/Tipes/geopot_vegcover_25.nc')
# ru = 'ho03b'
# resdict[(ru, 'jli')] = jli
# resdict[(ru, 'jspeed')] = jspeed
# resdict[(ru, 'dates')] = dates
#
# with open(cart_out + 'res_jli_hosing2.p', 'wb') as filox:
#     pickle.dump(resdict, filox)
#
# sys.exit()

#

# #############################################################################
# resdict = dict()
#
# for filone, ru, col in zip(fils, allru, colors):
#     jli, jspeed, dates = cd.jli_from_files(cart + filone, orogfile = '/home/federico/work/Tipes/geopot_vegcover_25.nc')
#
#     resdict[(ru, 'jli')] = jli
#     resdict[(ru, 'jspeed')] = jspeed
#     resdict[(ru, 'dates')] = dates
#
# with open(cart + 'res_jli_hosing.p', 'wb') as filox:
#     pickle.dump(resdict, filox)

cart_out = '/home/fedef/Research/lavori/tipes/'
with open(cart_out + 'res_jli_hosing.p', 'rb') as filox:
    resdict = pickle.load(filox)

with open(cart_out + 'res_jli_hosing2.p', 'rb') as filox:
    resdict2 = pickle.load(filox)

resdict.update(resdict2)

fig = plt.figure(figsize = (24,12))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

latsel = np.arange(20, 71, 0.5)
kos = np.concatenate([resdict[(ru, 'jspeed')] for ru in allru])
vmin, vmax = (np.min(kos), np.max(kos))
print(vmin, vmax)
vsel = np.linspace(vmin, vmax, 100)

def dopdf(var, xi, bnd_width):
    pdf = ctl.calc_pdf(var, bnd_width = bnd_width)
    pdfok = pdf(xi)
    pdfok /= np.sum(pdfok)
    return pdfok

for ru, col in zip(allru, colors):
    jli = resdict[(ru, 'jli')]
    jspeed = resdict[(ru, 'jspeed')]

    if ru == 'ho03':
        jli, _ = ctl.sel_time_range(jli, resdict[(ru, 'dates')], ctl.range_years(1860, 1880), ignore_HHMM = False)
        jspeed, _ = ctl.sel_time_range(jspeed, resdict[(ru, 'dates')], ctl.range_years(1860, 1880), ignore_HHMM = False)
    elif ru == 'pi':
        jliserie = ctl.bootstrap(jli, resdict[(ru, 'dates')], None, apply_func = dopdf, func_args = [latsel, 0.22], n_choice = 20, n_bootstrap = 1000)
        jspedserie = ctl.bootstrap(jspeed, resdict[(ru, 'dates')], None, apply_func = dopdf, func_args = [vsel, None], n_choice = 20, n_bootstrap = 1000)

    pdf = ctl.calc_pdf(jli, bnd_width = 0.22)
    pdfok = pdf(latsel)
    pdfok /= np.sum(pdfok)
    if ru == 'pi':
        jlimin = np.percentile(jliserie, 10, axis = 0)
        jlimax = np.percentile(jliserie, 90, axis = 0)
        ax1.fill_between(latsel, jlimin, jlimax, color = col, alpha = 0.3)

    ax1.plot(latsel, pdfok, label = ru, color = col, linewidth = 3)

    pdf = ctl.calc_pdf(jspeed)
    pdfok = pdf(vsel)
    pdfok /= np.sum(pdfok)
    if ru == 'pi':
        jlimin = np.percentile(jspedserie, 10, axis = 0)
        jlimax = np.percentile(jspedserie, 90, axis = 0)
        ax2.fill_between(vsel, jlimin, jlimax, color = col, alpha = 0.3)

    ax2.plot(vsel, pdfok, label = ru, color = col, linewidth = 3)

ax1.grid()
ax2.grid()
ax1.set_xlabel('Latitude')
ax2.set_xlabel('u wind (m/s)')
ax1.set_title('Jet latitude index')
ax2.set_title('Jet speed')
ax1.legend()
ax2.legend()
fig.savefig(cart_out + 'jlinspeed_hosing2.pdf')
