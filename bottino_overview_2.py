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

import psutil

plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
titlefont = 22
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.axisbelow'] = True
plt.rcParams['legend.fontsize'] = 18

#############################################################################

tip = sys.argv[1]

if tip == 'calc':
    ru = sys.argv[2]

logname = 'log_overview_{}.log'.format(ru)
logfile = open(logname,'w') #self.name, 'w', 0)

# re-open stdout without buffering
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1)

# redirect stdout and stderr to the log file opened above
os.dup2(logfile.fileno(), sys.stdout.fileno())
os.dup2(logfile.fileno(), sys.stderr.fileno())

if tip == 'calc':
    print('Running calculation for {}'.format(ru))

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

ctl.mkdir(cart_out)

#filna = '/nas/BOTTINO/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/{}/{}i1p1f1/{}/{}/*nc'

allru = ['pi', 'b025', 'b050', 'b100']
allnams = ['piControl', 'stabilization-ssp585-2025', 'stabilization-ssp585-2050', 'stabilization-ssp585-2100']

colors = ['black', 'forestgreen', 'orange', 'violet']

####################################################################################################

#
masfi = cart_run + 'masks.nc'
cose = xr.load_dataset(masfi)
oce_mask = cose['RnfA.msk'].values.astype('bool') # 1 over ocean

miptab = 'Amon'
#allvars_2D = 'tas pr clt rlut rsdt rsut'.split()
#allvars_2D = 'rsds rsus rlds rlus hfss hfls'.split()
allvars_2D = 'clt pr rlut rsdt rsut tas rsds rsus rlds rlus hfss hfls'.split()
add_uas = False

allvars_3D = []#'ta ua'.split()

#var_glob_mean = 'tas pr clt rlut rsut net_toa'.split()  # plot global timeseries, including ssp585
#var_glob_mean = 'tas pr clt rlut rsut net_toa rsds rsus rlds rlus hfss hfls net_srf'.split()
var_glob_mean = 'tas pr clt rlut rsut net_toa rsds rsus rlds rlus hfss hfls net_srf'.split()

var_glob_ylabel = ['Temp. anomaly (K)', r'Prec. anomaly (kg m$^{-2}$ s$^{-1}$)', 'Cloud cover', 'Outgoing Longwave Radiation (W m$^{-2}$)', 'Outgoing Shortwave Radiation (W m$^{-2}$)', 'Net incoming TOA flux (W m$^{-2}$)']
var_glob_ylabel += 'rsds rsus rlds rlus hfss hfls net_srf'.split()


allruok = ['b990', 'b025', 'b050', 'b065', 'b080', 'b100']
colok = ['lightslategray', 'forestgreen', 'orange', 'chocolate', 'maroon', 'violet']

#colall = ['black', 'royalblue', 'crimson', 'lightslategray', 'forestgreen', 'orange', 'chocolate', 'maroon', 'violet']

figs_glob = []
axs_glob = []
pimean = dict()
glomeans = dict()
yeamean = dict()
mapmean = dict()

# Read already computed experiments
gigi, pimean = pickle.load(open(cart_out + 'bottino_glomeans.p', 'rb'))
for ke in gigi:
    if 'pi' in gigi:
        glomeans[ke] = gigi[ke]

del gigi

#glomeans, pimean, yeamean, mapmean = pickle.load(open(cart_out + 'bottino_seasmean_2D.p', 'rb')) # too heavy!

#for ru in allruok:
if tip == 'calc':
    print(ru)
    mem = 'r1'
    if ru in ['ssp585', 'hist']: mem = 'r4'

    datadir = '/g100_scratch/userexternal/{}/ece3/{}/cmorized/'.format(user, ru)
    # else:
    #     datadir = '/g100_work/IscrB_QUECLIM/BOTTINO/{}/cmorized/'.format(ru)

    filna = datadir+'cmor_*/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/*/r1i1p1f1/{}/{}/g*/v*/{}*nc'#.format(fis, miptab, var, var)

    fils = np.concatenate([glob.glob(filna.format(miptab, var, var)) for var in allvars_2D])

    if len(fils) == 0:
        print('no files for {}'.format(na))
        continue

    # for var in allvars_2D:
    #     print(var, len(glob.glob(filna.format(na, mem, miptab, var))))
    #     print(filna.format(na, mem, miptab, var))
    # continue

    kose = xr.open_mfdataset(fils, use_cftime = True, chunks={'time' : 1200})
    kose = kose.drop_vars('time_bnds')

    kose = kose.assign(net_toa = kose.rsdt - kose.rlut - kose.rsut) # net downward energy flux at TOA
    print('assigned net_toa')
    kose = kose.assign(net_srf = kose.rsds + kose.rlds - kose.rsus - kose.rlus - kose.hfss - kose.hfls) # net downward energy flux at srf
    print('assigned net_srf')

    # # Separate for uas
    # if add_uas:
    #     fils = glob.glob(filna.format(na, mem, miptab, 'uas'))
    #     if len(fils) > 0:
    #         kosettt = xr.open_mfdataset(fils, use_cftime = True)
    #         kosettt = kosettt.drop_vars('time_bnds')
    #         kosettt = kosettt.drop_vars('height')
    #         kose = kose.assign(uas = kosettt.uas)

    for var in var_glob_mean:
        print(var)
        if var not in kose:
            if ru == 'pi':
                pimean[var] = 0.
            continue

        cosoye = kose[var].groupby("time.year").mean().compute()

        print('total RAM memory used after {}:'.format(var), psutil.virtual_memory()[3]/1.e9)

        if var in ['tas', 'pr', 'rlut', 'rsut', 'net_toa', 'net_srf']:
            yeamean[(ru, var)] = cosoye

        if var == 'net_toa':
            for vauu in 'tas pr clt rlut rsut net_toa'.split():
                del kose[vauu]

        # coso = cosoye.mean('lon')
        # glomean = np.average(coso, weights = abs(np.cos(np.deg2rad(coso.lat))), axis = -1)
        glomean = ctl.global_mean(cosoye)
        glomean_oce = ctl.global_mean(cosoye, mask = oce_mask)
        glomean_land = ctl.global_mean(cosoye, mask = ~oce_mask)

        # if ru == 'pi':
        #     years = cosoye.year.data-2256+2015
        #     pimean[var] = np.mean(glomean)
        # else:
        #     years = cosoye.year.data

        glomeans[(ru, var)] = (years, glomean)

        glomeans[(ru, var, 'oce')] = (years, glomean_oce)
        glomeans[(ru, var, 'land')] = (years, glomean_land)


    if ru == 'ssp585':
        continue

    # if ru != 'pi':
    #     kose = kose.sel(time = kose['time.year'] >= kose['time.year'].data[-1]-200)

    # for var in var_map_200:
    #     print(var)
    #     kose_sclim = ctl.seasonal_climatology(kose[var])
    #     mapmean[(ru, var)] = kose_sclim

    pickle.dump(yeamean, open(cart_out + 'bottino_seasmean_2D_{}_1000.p'.format(ru), 'wb'))

    pickle.dump(glomeans, open(cart_out + 'bottino_glomeans_{}_1000.p'.format(ru), 'wb'))

    del yeamean
    del glomeans
    yeamean = dict()
    glomeans = dict()

#pickle.dump([glomeans, pimean], open(cart_out + 'bottino_glomeans_1000.p', 'wb'))
#pickle.dump([glomeans, pimean], open(cart_out + 'bottino_glomeans.p', 'wb'))

if tip == 'plot':
    glomeans = dict()
    for ru in allruok:
        gigi = pickle.load(open(cart_out + 'bottino_glomeans_{}_1000.p'.format(ru), 'rb'))
        glomeans.update(gigi)

    pickle.dump([glomeans, pimean], open(cart_out + 'bottino_glomeans_1000.p', 'wb'))

    for var in ['tas', 'pr', 'net_toa', 'net_srf']:
        yeamean_var = dict()
        for ru in allruok:
            gigi = pickle.load(open(cart_out + 'bottino_seasmean_2D_{}_1000.p'.format(ru), 'rb'))
            yeamean_var[(ru, var)] = gigi[(ru, var)]

        pickle.dump(yeamean_var, open(cart_out + 'bottino_seasmean_2D_{}_1000.p'.format(var), 'wb'))
        del yeamean_var


# pickle.dump([glomeans, pimean, yeamean, mapmean], open(cart_out + 'bottino_seasmean_2D_srf.p', 'wb'))

#glomeans, pimean, yeamean, mapmean = pickle.load(open(cart_out + 'bottino_seasmean_2D.p', 'rb'))

### Add b990
# ru = 'b990'
# mem = 'r1'
#
# fil = '/nas/BOTTINO/indices/b990/tas_stabilization-hist-1990_199001-218912.nc'
#
# kose = xr.open_dataset(fil, use_cftime = True)
# kose = kose.drop_vars('time_bnds')
#
# var = 'tas'
#
# cosoye = kose[var].groupby("time.year").mean().compute()
# yeamean[(ru, var)] = cosoye
#
# coso = cosoye.mean('lon')
# glomean = np.average(coso, weights = abs(np.cos(np.deg2rad(coso.lat))), axis = -1)
# years = coso.year.data
#
# glomeans[(ru, var)] = (years, glomean)

#
# pickle.dump([glomeans, pimean, yeamean, mapmean], open(cart_out + 'bottino_seasmean_2D.p', 'wb'))

# # 3D vars
# for na, ru, col in zip(allnams2, allru2, colors2):
#     mem = 'r1'
#     if na == 'ssp585': mem = 'r4'
#
#     fils = np.concatenate([glob.glob(filna.format(na, mem, miptab, var)) for var in allvars_3D])
#
#     kose = xr.open_mfdataset(fils, use_cftime = True)
#     kose = kose.drop_vars('time_bnds')
#
#     for var in allvars_3D:
#         print(var)
#         cosoye = kose[var].groupby("time.year").mean().compute()
#         yeamean[(ru, var)] = cosoye
#
#     if ru == 'ssp585':
#         continue
#
#     if ru != 'pi':
#         kose = kose.sel(time = kose['time.year'] >= kose['time.year'].data[-1]-200)
#
#     for var in allvars_3D:
#         print(var)
#         kose_sclim = ctl.seasonal_climatology(kose[var])
#         mapmean[(ru, var)] = kose_sclim
#
#     # kose_smean = kose.groupby("time.season").mean()
#     # kose_sstd = kose.groupby("time.season").std()
#     # kose_p90 = kose.groupby("time.season").percentile(90)
#     # kose_p10 = kose.groupby("time.season").percentile(10)
#     #
#     # for var in var_map_200:
#     #     mapmean[(ru, var, 'mean')] = kose_smean
#     #     mapmean[(ru, var, 'std')] = kose_sstd
#     #     mapmean[(ru, var, 'p90')] = kose_p90
#     #     mapmean[(ru, var, 'p10')] = kose_p10
#
# pickle.dump([glomeans, pimean, yeamean, mapmean], open(cart_out + 'bottino_seasmean.p', 'wb'))

# glomeans, pimean, yeamean, mapmean = pickle.load(open(cart_out + 'bottino_seasmean.p', 'rb'))

if tip == 'plot':
    figs_glob = []
    axs_glob = []
    for var in var_glob_mean:
        fig, ax = plt.subplots(figsize = (16,9))
        axs_glob.append(ax)
        figs_glob.append(fig)
        ax.set_title(var)

    fig_greg, ax_greg = plt.subplots(figsize = (16,9))

    #for na, ru, col in zip(allnams, allru, colors):
    #for na, ru, col in zip(allnams3, allru3, colors3):
    for ru, col in zip(allruall, colall):
        print(ru)

        for var, fig, ax, vylab in zip(var_glob_mean, figs_glob, axs_glob, var_glob_ylabel):
            print(var)
            if (ru, var) not in glomeans.keys():
                print('NOT found')
                continue

            # cosoye = yeamean[(ru, var)]
            years, glomean = glomeans[(ru, var)]

            ax.plot(years, glomean-pimean[var], label = ru, color = col)

            if ru == 'pi':
                ax.set_ylabel(vylab)
                ax.set_xlabel('Year')

        # gregory
        try:
            #ax_greg.plot(glomeans[(ru, 'tas')][1]-pimean['tas'], glomeans[(ru, 'net_toa')][1], label = ru, color = col)
            ctl.gregplot_on_ax(ax_greg, glomeans[(ru, 'tas')][1]-pimean['tas'], glomeans[(ru, 'net_toa')][1], color = col, label = ru, calc_ERF = False, calc_ECS = False, nyea = 10, point_dim = 20)
        except Exception as exc:
            print(ru, exc)
            pass

    ax_greg.legend()
    ax_greg.grid()

    for ax in axs_glob:
        ax.legend()
        ax.grid()

    ctl.plot_pdfpages(cart_out + 'bottino_glomeans_1000.pdf', figs_glob, True, )

    ax_greg.set_xlabel('Global mean tas (K)')
    ax_greg.set_ylabel('Global net incoming TOA flux (W/m2)')
    fig_greg.savefig(cart_out + 'bottino_gregory_1000.pdf')

    ax_greg.set_xlim((-0.5, 0.5))
    ax_greg.set_ylim((-0.5, 0.5))
    fig_greg.savefig(cart_out + 'bottino_gregory_pizoom.pdf')


sys.exit()

# var_map_200 = 'clt pr tas rlut uas'.split()  # plot last 200 mean map, stddev, low/high var wrt pi
#
# allcopls = ['seamean', 'seastd', 'seap10', 'seap90']
# for ru in allru:
#     zup = mapmean[(ru, 'tas')]
#     mapmean[(ru, 'tas_patt')] = zup.copy(deep = True) # it's back!
#     zupme = ctl.global_mean(zup['seamean'])
#     for copl in allcopls:
#         mapmean[(ru, 'tas_patt')][copl] = zup[copl]-zupme[:, np.newaxis, np.newaxis]
#
#     # zup = mapmean[(ru, 'pr')]
#     # mapmean[(ru, 'pr_perc')] = zup
#     # zupme = mapmean[('pi', 'pr')]['seamean']
#     # for copl in allcopls:
#     #     mapmean[(ru, 'pr_perc')][copl] = (zup[copl]-zupme)/zupme
#
# print('CHECK! -> ', mapmean[(ru, 'tas')] is mapmean[(ru, 'tas_patt')], mapmean[(ru, 'tas')][copl] is mapmean[(ru, 'tas_patt')][copl])
#
# var_map_200 += ['tas_patt', 'pr_perc']
# ###### Plots 2D
# figs_map = []
# for var in var_map_200:
#     for copl in allcopls:
#         #mappe = [mapmean[('pi', var, copl)]] + [mapmean[(ru, var, copl)]-mapmean[('pi', var, copl)] for ru in allru[1:]]
#         if var == 'pr_perc':
#             zupme = mapmean[('pi', 'pr')]['seamean']
#             mappe = [mapmean[('pi', 'pr')][copl]] + [(mapmean[(ru, 'pr')][copl]-zupme)/zupme for ru in allru[1:]]
#             mapcont = [None if ru == 'pi' else mapmean[('pi', 'pr')][copl].sel(season = seasok) for seasok in ['DJF', 'MAM', 'JJA', 'SON'] for ru in allru]
#         else:
#             mappe = [mapmean[('pi', var)][copl]] + [mapmean[(ru, var)][copl]-mapmean[('pi', var)][copl] for ru in allru[1:]]
#             mapcont = [None if ru == 'pi' else mapmean[('pi', var)][copl].sel(season = seasok) for seasok in ['DJF', 'MAM', 'JJA', 'SON'] for ru in allru]
#
#         plot_anoms = [False if ru == 'pi' else True for se in range(4) for ru in allru]
#         cmaps = ['viridis' if ru == 'pi' else 'RdBu_r' for se in range(4) for ru in allru]
#
#         cbr0 = ctl.get_cbar_range(mappe[0].values, False, (2,98))
#         cbr1 = ctl.get_cbar_range([ma.values for ma in mappe[1:]], True, (2,98))
#         cbar_range = [cbr0 if ru == 'pi' else cbr1 for se in range(4) for ru in allru]
#
#         #mappeseas = [ma.sel(time = ma['time.season'] == seasok). for seasok in ['DJF', 'MAM', 'JJA', 'SON'] for ma in mappe]
#         mappeseas = [ma.sel(season = seasok) for seasok in ['DJF', 'MAM', 'JJA', 'SON'] for ma in mappe]
#
#         subtitles = ['{} - {}'.format(ru, seasok) for seasok in ['DJF', 'MAM', 'JJA', 'SON'] for ru in allru]
#
#         fig = ctl.plot_multimap_contour(mappeseas, figsize = (20,12), cmap = cmaps, cbar_range = cbar_range, use_different_cbars = True, use_different_cmaps = True, subtitles = subtitles, title = var+' - '+copl, add_contour_field = mapcont, add_contour_plot_anomalies = False)
#         figs_map.append(fig)
#
# figs_map = np.concatenate(figs_map)
# fignames = [var+'_'+copl for var in var_map_200 for copl in allcopls]
# ctl.plot_pdfpages(cart_out + 'bottino_mapmeans.pdf', figs_map, True, fignames)
#
# figs_map = []
# for var in allvars_3D:
#     for copl in allcopls:
#         fig, axs = plt.subplots(4, 4, figsize = (20,12))
#         #mappe = [mapmean[('pi', var, copl)]] + [mapmean[(ru, var, copl)]-mapmean[('pi', var, copl)] for ru in allru[1:]]
#         mappe = [mapmean[('pi', var)][copl]] + [mapmean[(ru, var)][copl]-mapmean[('pi', var)][copl] for ru in allru[1:]]
#
#         #mappeseas = [ma.sel(time = ma['time.season'] == seasok). for seasok in ['DJF', 'MAM', 'JJA', 'SON'] for ma in mappe]
#         mappeseas = [ma.sel(season = seasok) for seasok in ['DJF', 'MAM', 'JJA', 'SON'] for ma in mappe]
#         subtitles = ['{} - {}'.format(ru, seasok) for seasok in ['DJF', 'MAM', 'JJA', 'SON'] for ru in allru]
#
#         for ma, ax, subt in zip(mappeseas, axs.flatten(), subtitles):
#             guplo = ma.mean('lon').plot.contourf(x = 'lat', y = 'plev', ax = ax, levels = 11, ylim = (1.e5, 1.e3), yscale = 'log')#vmax = )
#             try:
#                 guplo.colorbar.set_label('')
#             except Exception as exc:
#                 print(exc)
#                 pass
#
#             if 'pi' not in subt:
#                 seasok = subt.split('-')[1].strip()
#                 guplo2 = mapmean[('pi', var)][copl].sel(season = seasok).mean('lon').plot.contour(x = 'lat', y = 'plev', ax = ax, levels = 11, ylim = (1.e5, 1.e3), yscale = 'log', colors="k", add_colorbar=False)#vmax = )
#             #guplo.set_titles(template='{value}', maxchar = 13, fontsize = 12)
#             ax.set_title(subt)
#
#         for i in range(4):
#             for j in range(4):
#                 if j > 0:
#                     ax.set_ylabel('')
#
#         fig.suptitle(var+' - '+copl)
#         plt.tight_layout()
#
#         figs_map.append(fig)
#
# fignames = [var+'_'+copl for var in var_map_200 for copl in allcopls]
# ctl.plot_pdfpages(cart_out + 'bottino_crossmeans.pdf', figs_map, True, fignames)
