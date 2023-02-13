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
# import tunlib as tl

from matplotlib.colors import LogNorm
from datetime import datetime

from scipy import stats
import xarray as xr
import glob
import xclim

import psutil

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 24
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.axisbelow'] = True

#############################################################################
# cart_out = '/home/fabiano/Research/lavori/BOTTINO/seasmean/'
# ctl.mkdir(cart_out)
#

tip = sys.argv[1]

logname = 'log_overview_seaice.log'
logfile = open(logname,'w') #self.name, 'w', 0)

# re-open stdout without buffering
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1)

# redirect stdout and stderr to the log file opened above
os.dup2(logfile.fileno(), sys.stdout.fileno())
os.dup2(logfile.fileno(), sys.stderr.fileno())


if os.uname()[1] == 'hobbes':
    cart_out = '/home/{}/Research/lavori/BOTTINO/seasmean/'.format(os.getlogin())
    cart_run = '/home/{}/Research/git/ece_runtime/run/'.format(os.getlogin())
elif os.uname()[1] == 'xaru':
    cart_out = '/home/{}/Research/lavori/BOTTINO/seasmean/'.format(os.getlogin())
elif os.uname()[1] == 'tintin':
    cart_out = '/home/{}/work/lavori/BOTTINO/seasmean/'.format(os.getlogin())
    cart_run = '/home/{}/Research/git/ece_runtime/run/'.format(os.getlogin())
else:
    user = 'ffabiano'
    cart_out = '/g100_work/IscrB_QUECLIM/BOTTINO/bottino_an/seasmean/'
    cart_run = '/g100_scratch/userexternal/{}/ece3/b00I/runtime/'.format(user)

allruadd2 = ['b065', 'b080']
allnadd2 = ['stabilization-ssp585-2065', 'stabilization-ssp585-2080']
colorsadd2 = ['chocolate', 'maroon']

allruall = ['pi', 'hist', 'ssp585', 'b990', 'b025', 'b050', 'b065', 'b080', 'b100']
okye = [(2015, 2516), (1850, 2015), (2015, 2101), (1990, 2490), (2025, 2525), (2050, 2550), (2065, 2565), (2080, 2580), (2100, 2600)]
#colall = ['black', 'steelblue', 'indianred', 'teal', 'forestgreen', 'orange', 'chocolate', 'maroon', 'violet']
colall = ['black', 'royalblue', 'crimson', 'lightslategray', 'forestgreen', 'orange', 'chocolate', 'maroon', 'violet']

##################################################################

allru = ['pi', 'b025', 'b050', 'b100']
allnams = ['piControl', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']

colors = ['black', 'forestgreen', 'orange', 'violet']

process = psutil.Process(os.getpid())
#############################################################################
## SEA ICE
#areacelfi = '/nas/BOTTINO/areas.nc'
areacelfi = cart_run + 'areas.nc'
acel = xr.open_dataset(areacelfi)
areaT = np.array(acel['O1t0.srf'].data)

#miptab = 'SImon_r1'
miptab = 'SImon'
varnam = 'siconc'

resdict = pickle.load(open(cart_out + 'seaicearea.p', 'rb'))
#resdict = dict()

fig, axs = plt.subplots(2,2,figsize = (16,9))
#
allnams2 = ['ssp585', 'historical'] + allnams
allru2 = ['ssp585', 'hist'] + allru
colors2 = ['indianred', 'steelblue'] + colors

allnams3 = allnams2 + ['stabilization-hist-1990']
allru3 = allru2 + ['b990']
colors3 = colors2 + ['teal']

allruok = ['b990', 'b025', 'b050', 'b065', 'b080', 'b100']
colok = ['lightslategray', 'forestgreen', 'orange', 'chocolate', 'maroon', 'violet']

na = ''

#for ru, col in zip(allruadd2, colorsadd2):
if tip == 'calc':
    for ru, col in zip(allruok, colok):
        print(ru)
        mem = 'r1'
        if ru in ['ssp585', 'hist']: mem = 'r4'

        if os.uname()[1] == 'hobbes':
            filna = '/nas/BOTTINO/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/{}/{}i1p1f1/{}/{}/*nc'
            filist = glob.glob(filna.format(na, mem, miptab, varnam))
        else:
            datadir = '/g100_scratch/userexternal/{}/ece3/{}/cmorized/'.format(user, ru)
            filna = datadir+'cmor_*/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/*/r1i1p1f1/{}/{}/g*/v*/{}*nc'
            filist = glob.glob(filna.format(miptab, varnam, varnam))

        gigi = xr.open_mfdataset(filist, use_cftime=True)

        try:
            lat = np.array(gigi.lat.data)
        except:
            print('lat name is latitude')
            lat = np.array(gigi.latitude.data)

        seaice = np.array(gigi.siconc.data)

        print('total RAM memory used for {}: '.format(ru), process.memory_info().rss/1e9)


        #okslat = lat > 40.
        for ii, emi, okslat in zip([0,1], ['N', 'S'], [lat > 40, lat < -40]):
            areaok = areaT[okslat]
            oksi = seaice[:, okslat]
            oksi[oksi < 15.] = 0.
            oksi[oksi > 15.] = 1.
            oksiarea = oksi*areaok[np.newaxis, :]
            seaicearea = np.nansum(oksiarea, axis = 1)

            dates = np.array(gigi.time.data)
            okmarch = np.array([da.month == 3 for da in dates])
            oksept = np.array([da.month == 9 for da in dates])

            resdict[(ru, varnam, 'glomean', 'mar', emi)] = seaicearea[okmarch]
            resdict[(ru, varnam, 'glomean', 'sep', emi)] = seaicearea[oksept]

        gigi.close()
        del seaice, oksi, oksiarea, seaicearea

        print('total RAM memory used after {}: '.format(ru), process.memory_info().rss/1e9)

    pickle.dump(resdict, open(cart_out + 'seaicearea_1000.p', 'wb'))

if tip == 'plot':
    resdict = pickle.load(open(cart_out + 'seaicearea.p', 'rb'))
    resdict_1000 = pickle.load(open(cart_out + 'seaicearea_1000.p', 'rb'))
    resdict.update(resdict_1000)
    print('plot')

    for ru, col, (ye1, ye2) in zip(allruall, colall, okye):
        print(ru)
        for ii, emi in enumerate(['N', 'S']):
            sia_march = resdict[(ru, varnam, 'glomean', 'mar', emi)]
            sia_sept = resdict[(ru, varnam, 'glomean', 'sep', emi)]

            if emi == 'N':
                sia_summer = sia_sept
                sia_winter = sia_march
            else:
                sia_summer = sia_march
                sia_winter = sia_sept

            yeaok = np.arange(ye1, ye1 + len(sia_march))

            siwi10 = ctl.running_mean(sia_winter, 10)
            sisu10 = ctl.running_mean(sia_summer, 10)

            axs[ii,0].plot(yeaok, siwi10, linestyle='solid', marker = 'None', color = col, label = ru, linewidth = 2)
            axs[ii,1].plot(yeaok, sisu10, linestyle='solid', marker = 'None', color = col, label = ru, linewidth = 2)
            axs[ii,0].plot(yeaok, sia_winter, linestyle='solid', color = col, alpha = 0.3, linewidth = 0.5)
            axs[ii,1].plot(yeaok, sia_summer, linestyle='solid', color = col, alpha = 0.3, linewidth = 0.5)

    axs[0,0].set_title('winter')
    axs[0,1].set_title('summer')
    axs[0,0].set_ylabel(r'Sea ice extent (m$^2$)')
    axs[1,0].set_ylabel(r'Sea ice extent (m$^2$)')
    axs[0,0].text(0.05, 0.7, 'Arctic', horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 24)
    axs[0,0].text(0.05, 0.3, 'Antarctic', horizontalalignment='center', verticalalignment='center', rotation='vertical',transform=fig.transFigure, fontsize = 24)
    #axs[1,1].legend()
    ctl.custom_legend(fig, colall, allruall, ncol = 5, add_space_below = 0.06)
    fig.savefig(cart_out + 'bottseaice_1000.pdf')

sys.exit()

#
# miptab = 'SImon_r1'
# var_map_200 = ['siconc', 'sithick']
#
# mapmean = dict()
# for na, ru, col in zip(allnams, allru, colors):
#     mem = 'r1'
#     filist = np.concatenate([glob.glob(filna.format(na, mem, miptab, varnam)) for varnam in var_map_200])
#     gigi = xr.open_mfdataset(filist, use_cftime=True)
#
#     if ru != 'pi':
#         gigi = gigi.sel(time = gigi['time.year'] >= gigi['time.year'].data[-1]-200)
#
#     for var in var_map_200:
#         print(var)
#         gigi_sclim = ctl.seasonal_climatology(gigi[var], allseasons = ['Mar', 'Sep'])
#         mapmean[(ru, var)] = gigi_sclim
#
# pickle.dump(mapmean, open(cart_out + 'seaicemap.p', 'wb'))

var_map_200 = ['siconc', 'sithick']
mapmean = pickle.load(open(cart_out + 'seaicemap.p', 'rb'))

allcopls = ['seamean', 'seastd', 'seap10', 'seap90']
###### Plots 2D
figs_map = []
for var in var_map_200:
    for copl in allcopls:
        #mappe = [mapmean[('pi', var, copl)]] + [mapmean[(ru, var, copl)]-mapmean[('pi', var, copl)] for ru in allru[1:]]
        #mappe = [mapmean[('pi', var)][copl]] + [mapmean[(ru, var)][copl]-mapmean[('pi', var)][copl] for ru in allru[1:]]
        mappe = [mapmean[(ru, var)][copl] for ru in allru]
        # for ma in mappe:
        #     ma[ma == 0.] = np.nan

        #mappeseas = [ma.sel(time = ma['time.season'] == seasok). for seasok in ['DJF', 'MAM', 'JJA', 'SON'] for ma in mappe]
        mappeseas = [ma.sel(season = seasok) for seasok in ['Mar', 'Sep'] for ma in mappe]
        mapcont = [None if ru == 'pi' else mapmean[('pi', var)][copl].sel(season = seasok) for seasok in ['Mar', 'Sep'] for ru in allru]
        subtitles = ['{} - {}'.format(ru, seasok) for seasok in ['Mar', 'Sep'] for ru in allru]

        cba = None
        cma = 'viridis'
        if var == 'siconc':
            cba = (0,1)
            cma = 'Blues_r'
        #for pol in ['Npolar', 'Spolar']:
        for clat in [90, -90]:
            print(var, copl, clat)
            fig = ctl.plot_multimap_contour(mappeseas, figsize = (24,12), plot_anomalies = False, subtitles = subtitles, title = var+' - '+copl, add_contour_field = mapcont, add_contour_plot_anomalies = False, visualization = 'nearside', cmap = cma, cbar_range = cba, central_lat_lon = (clat, 0), bounding_lat = 40 * np.sign(clat), fix_subplots_shape = (2,4))

            figs_map.append(fig)

figs_map = np.concatenate(figs_map)
fignames = [var+'_'+copl+'_'+pole for var in var_map_200 for copl in allcopls for pole in ['N', 'S']]
ctl.plot_pdfpages(cart_out + 'bottino_seaicemap.pdf', figs_map, True, fignames)
