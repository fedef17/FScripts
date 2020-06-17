#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
from matplotlib import pyplot as plt
from matplotlib import cm

import pickle

import climtools_lib as ctl
import climdiags as cd

from matplotlib.colors import LogNorm
from datetime import datetime

from scipy import stats

#######################################
cart_out = '/home/fabiano/Research/lavori/prima_regimes_KS/'
if not os.path.exists(cart_out): os.mkdir(cart_out)

# apro il file originale
cart = '/home/fabiano/Research/lavori/WeatherRegimes/prima_coup_v7/'
filogen = cart + 'out_prima_coup_v7_DJF_EAT_4clus_4pcs_1957-2014_refEOF.p'

# apro il file con le pcs nuove
with open(cart_out + 'residual_pc_clustering_as_dict.pkl', 'rb') as fil:
    nures = pickle.load(fil)
with open(cart_out + 'raw_pc_clustering_as_dict.pkl', 'rb') as fil:
    nures_nofil = pickle.load(fil)
with open(cart_out + 'ERA20C_r_clusters.pkl', 'rb') as fil:
    refres = pickle.load(fil)
with open(cart_out + 'ERA20C_u_clusters.pkl', 'rb') as fil:
    refres_nofil = pickle.load(fil)

refres['time'] = refres_nofile['time']

refcen = gigiref[4]['centroids'][:,:4]/9.81 # these should be the counterparts of the original regimes in the residual phase space
refpcs, refdat = ctl.sel_time_range(refres['pcs'], refres['time'], ctl.range_years(1954, 2010))
refpcs = refpcs/9.81
reflab, _ = ctl.sel_time_range(refres[4]['states'], refres['time'], ctl.range_years(1954, 2010))

perm = np.array([0, 2, 3, 1])
refcen, reflab = ctl.change_clus_order(refcen, reflab, perm)

##############################################################################

restot = pickle.load(open(filogen, 'rb'))
results_old = restot['models']
results_ref = restot['reference']

results_ref['dates'] = refdat
results_ref['pcs'] = refpcs
results_ref['centroids'] = refcen
results_ref['labels'] = reflab

results = dict()
for ke in results_old:
    if ke in nures:
        check = np.all(results_old[ke]['pcs'] == nures_nofil[ke]['pcs'])
        print('OFUND: ', ke, check)
        if check:
            results[ke] = results_old[ke]
            results[ke]['pcs'] = nures[ke]['pcs']

            perm = ctl.match_pc_sets(refcen, nures[ke][4]['centroids'])
            centroids, labels = ctl.change_clus_order(nures[ke][4]['centroids'], nures[ke][4]['states'], perm)

            results[ke]['labels'] = labels
            results[ke]['centroids'] = centroids
    else:
        print('{} NOT FOUND!!!! Skipping.....'.format(ke))

results['ERA_0'] = results_ref

n_choice = 30
n_bootstrap = 200

#filo = open(cart_out + 'res_bootstrap_v7_500_relent2_nosig.p', 'wb')
filo = open(cart_out + 'res_bootstrap_v7_KJ.p', 'wb')

model_names = ['AWI-CM-1-1-LR', 'AWI-CM-1-1-HR', 'CMCC-CM2-HR4', 'CMCC-CM2-VHR4', 'CNRM-CM6-1', 'CNRM-CM6-1-HR', 'EC-Earth3P', 'EC-Earth3P-HR', 'ECMWF-IFS-LR', 'ECMWF-IFS-MR', 'ECMWF-IFS-HR', 'MPI-ESM1-2-HR', 'MPI-ESM1-2-XR', 'HadGEM3-GC31-LL', 'HadGEM3-GC31-MM', 'HadGEM3-GC31-HM', 'HadGEM3-GC31-HH']

vers = ['LR', 'HR', 'LR', 'HR', 'LR', 'HR', 'LR', 'HR', 'LR', 'MR', 'HR', 'LR', 'HR', 'LR', 'MR', 'MR', 'HR'] + ['OBS']
model_coups = ['AWI-CM-1-1', 'CMCC-CM2', 'CNRM-CM6-1', 'EC-Earth3P', 'ECMWF-IFS', 'MPI-ESM1-2', 'HadGEM3-GC31']

model_names_all = model_names + ['ERA']

colors = ctl.color_set(len(model_coups)*2, sns_palette = 'Paired')
mr_cos = np.where(np.array(vers) == 'MR')[0]
for gi in mr_cos:
    colors.insert(gi, np.mean(colors[gi-1:gi+1], axis = 0))

colors_wERA = colors + ['black']

regnam = ['NAO +', 'Sc. BL', 'AR', 'NAO -']

################################################################################

# results, results_ref = pickle.load(open(filogen, 'rb'))
# restot = pickle.load(open(filogen, 'rb'))
# results = restot['models']
# results_ref = restot['reference']


all_mods = np.array([ke.split('_')[0] for ke in results.keys()])
all_mems = np.array([ke.split('_')[1] for ke in results.keys()])

ref_cen = results['ERA_0']['centroids']
labels = results['ERA_0']['labels']
pcs = results['ERA_0']['pcs']
#dates = results['ERA_0']['dates']

grid_i = []
for i in range(3):
    (x0, x1) = (np.percentile(pcs[:, i], 1), np.percentile(pcs[:, i], 99))
    xss = np.linspace(x0, x1, 50)
    grid_i.append(xss)

xi_grid, yi_grid, zi_grid = np.meshgrid(*grid_i)

zi_ref = []
for reg in range(4):
    okclu = labels == reg
    okpc = pcs[okclu, :]
    kufu = ctl.calc_pdf(okpc[:, :3].T)
    zi = kufu(np.vstack([xi_grid.flatten(), yi_grid.flatten(), zi_grid.flatten()]))
    zi = zi/np.max(zi)

    zi_ref.append(zi)

# filt_labels = ctl.regime_filter_long(labels, dates, days_thres = 5)
# ref_cen_filt = []
# zi_ref_filt = []
# for reg in range(4):
#     okclu = filt_labels == reg
#     okpc = pcs[okclu, :]
#     kufu = ctl.calc_pdf(okpc[:, :3].T)
#     zi = kufu(np.vstack([xi_grid.flatten(), yi_grid.flatten(), zi_grid.flatten()]))
#     zi = zi/np.max(zi)
#
#     zi_ref_filt.append(zi)
#
#     centr = np.mean(okpc, axis = 0)
#     ref_cen_filt.append(centr)


##
#allkeysss = ['significance', 'varopt', 'autocorr', 'freq', 'dist_cen', 'resid_times_av', 'resid_times_90', 'trans_matrix', 'centroids', 'relative_entropy', 'RMS', 'patcor', 'filt_relative_entropy', 'filt_RMS', 'filt_patcor']
allkeysss = ['significance', 'varopt', 'autocorr', 'freq', 'dist_cen', 'resid_times_av', 'resid_times_90', 'trans_matrix', 'centroids', 'relative_entropy', 'patcor']#, 'filt_dist_cen', 'filt_relative_entropy', 'filt_patcor']
#allkeysss = ['varopt', 'autocorr', 'freq', 'dist_cen', 'resid_times_av', 'resid_times_90', 'trans_matrix', 'centroids', 'relative_entropy', 'patcor', 'filt_dist_cen', 'filt_relative_entropy', 'filt_patcor']

for mod in model_names_all:
    print(mod)
    whos_mod = all_mods == mod
    ok_mems = np.sort(all_mems[whos_mod])

    bootstraps_all = dict()
    for mem in ok_mems:
        print(mod, mem)
        modmem = mod + '_' + mem
        if modmem not in results.keys():
            print('NO data for {}, continuing...'.format(modmem))
            continue
        results[modmem]['varopt'] = ctl.calc_varopt_molt(results[modmem]['pcs'], results[modmem]['centroids'], results[modmem]['labels'])

        results[modmem]['autocorr_wlag'] = ctl.calc_autocorr_wlag(results[modmem]['pcs'], results[modmem]['dates'], out_lag1 = True)

        pcs_seas_set, dates_seas_set = ctl.seasonal_set(results[modmem]['pcs'], results[modmem]['dates'], 'DJF')
        labels_seas_set, dates_seas_set = ctl.seasonal_set(results[modmem]['labels'], results[modmem]['dates'], 'DJF')

        n_seas = len(dates_seas_set)
        years_set = np.array([dat[0].year for dat in dates_seas_set])
        bootstraps = dict()

        for nam in allkeysss:
            bootstraps[nam] = []

        t0 = datetime.now()
        for i in range(n_bootstrap):
            #if i % 10 == 0:
            print(i)
            ok_yea = np.sort(np.random.choice(list(range(n_seas)), n_choice))
            pcs = np.concatenate(pcs_seas_set[ok_yea])
            labels = np.concatenate(labels_seas_set[ok_yea])
            dates = np.concatenate(dates_seas_set[ok_yea])

            centroids = []
            for iclu in range(4):
                okla = labels == iclu
                centroids.append(np.mean(pcs[okla], axis = 0))
            centroids = np.stack(centroids)

            sig = ctl.clusters_sig(pcs, centroids, labels, dates, nrsamp = 200)
            # if sig < 10:
            #     sig2 = ctl.clusters_sig(pcs, centroids, labels, dates, nrsamp = 1000)
            #     print('RECHECK ', sig, sig2, sig3)

            varopt = ctl.calc_varopt_molt(pcs, centroids, labels)
            autocorr = ctl.calc_autocorr_wlag(pcs, dates, out_lag1 = True)
            bootstraps['significance'].append(sig)
            bootstraps['varopt'].append(varopt)
            bootstraps['autocorr'].append(autocorr)

            bootstraps['freq'].append(ctl.calc_clus_freq(labels, 4))

            centdist = np.array([ctl.distance(centroids[iclu], ref_cen[iclu]) for iclu in range(4)])
            bootstraps['dist_cen'].append(centdist)
            bootstraps['centroids'].append(centroids)

            resid_times = ctl.calc_regime_residtimes(labels, dates = dates)[0]
            av_res = np.array([np.mean(resid_times[reg]) for reg in range(4)])
            av_res_90 = np.array([np.percentile(resid_times[reg], 90) for reg in range(4)])
            bootstraps['resid_times_av'].append(av_res)
            bootstraps['resid_times_90'].append(av_res_90)

            bootstraps['trans_matrix'].append(ctl.calc_regime_transmatrix(1, labels, dates))

            # relative entropy, RMS, patcor
            relent_all = []
            for reg in range(4):
                okclu = labels == reg
                okpc = pcs[okclu, :]
                for comp in range(4):
                    okpc[:, comp] = okpc[:, comp] - centroids[reg, comp] + ref_cen[reg, comp]
                kufu = ctl.calc_pdf(okpc[:, :3].T)
                zi = kufu(np.vstack([xi_grid.flatten(), yi_grid.flatten(), zi_grid.flatten()]))
                zi = zi/np.max(zi)

                relent = stats.entropy(zi, zi_ref[reg])
                relent_all.append(relent)

            bootstraps['relative_entropy'].append(relent_all)

            #bootstraps['RMS'].append([ctl.distance(ce, refce) for ce, refce in zip(centroids, ref_cen)])
            bootstraps['patcor'].append([ctl.Rcorr(ce, refce) for ce, refce in zip(centroids, ref_cen)])

            # redo the same for filtered regimes
            # filt_labels = ctl.regime_filter_long(labels, dates, days_thres = 5)
            #
            # relent_all = []
            # filt_centroids = []
            # for reg in range(4):
            #     okclu = filt_labels == reg
            #     okpc = pcs[okclu, :]
            #     kufu = ctl.calc_pdf(okpc[:, :3].T)
            #     zi = kufu(np.vstack([xi_grid.flatten(), yi_grid.flatten(), zi_grid.flatten()]))
            #     zi = zi/np.max(zi)
            #
            #     relent = stats.entropy(zi, zi_ref_filt[reg])
            #     relent_all.append(relent)
            #
            #     centr = np.mean(okpc, axis = 0)
            #     filt_centroids.append(centr)
            #
            # bootstraps['filt_relative_entropy'].append(relent_all)
            #
            # centdist = np.array([ctl.distance(filt_centroids[iclu], ref_cen_filt[iclu]) for iclu in range(4)])
            # bootstraps['filt_dist_cen'].append(centdist)
            #
            # #bootstraps['filt_RMS'].append([ctl.E_rms(ce, refce) for ce, refce in zip(filt_centroids, ref_cen_filt)])
            # bootstraps['filt_patcor'].append([ctl.Rcorr(ce, refce) for ce, refce in zip(filt_centroids, ref_cen_filt)])

        for nam in allkeysss:
            bootstraps[nam] = np.array(bootstraps[nam])

        bootstraps_all[mem] = bootstraps

    # Calculate ensemble mean, std of:
        # - bootstrap mean, std, p10, p90 for each quantity;

    allmems = list(bootstraps_all.keys())
    bootstraps_all['boot_mean'] = dict() # list of all boot means
    bootstraps_all['boot_std'] = dict()
    bootstraps_all['boot_p10'] = dict()
    bootstraps_all['boot_p90'] = dict()

    bootstraps_all['ens_mean'] = dict() # mean of all boot means
    bootstraps_all['ens_std'] = dict()

    for ke in allkeysss:
        bootstraps_all['boot_mean'][ke] = np.array([np.mean(bootstraps_all[mem][ke], axis = 0) for mem in allmems])
        bootstraps_all['ens_mean'][ke] = np.mean(bootstraps_all['boot_mean'][ke], axis = 0)
        bootstraps_all['ens_std'][ke] = np.std(bootstraps_all['boot_mean'][ke], axis = 0)

        bootstraps_all['boot_mean'][ke] = np.mean(bootstraps_all['boot_mean'][ke], axis = 0)
        bootstraps_all['boot_std'][ke] = np.std(np.concatenate([bootstraps_all[mem][ke] for mem in allmems]), axis = 0)

        if bootstraps_all[allmems[0]][ke].ndim == 1:
            bootstraps_all['boot_p10'][ke] = np.percentile(np.concatenate([bootstraps_all[mem][ke] for mem in allmems]), 10)
            bootstraps_all['boot_p90'][ke] = np.percentile(np.concatenate([bootstraps_all[mem][ke] for mem in allmems]), 90)
        elif bootstraps_all[allmems[0]][ke].ndim == 2:
            bootstraps_all['boot_p10'][ke] = np.array([np.percentile(np.concatenate([bootstraps_all[mem][ke][:, reg] for mem in allmems]), 10) for reg in range(4)])
            bootstraps_all['boot_p90'][ke] = np.array([np.percentile(np.concatenate([bootstraps_all[mem][ke][:, reg] for mem in allmems]), 90) for reg in range(4)])
        else:
            bootstraps_all['boot_p10'][ke] = np.zeros((4,4))
            bootstraps_all['boot_p90'][ke] = np.zeros((4,4))
            for i in range(4):
                for j in range(4):
                    bootstraps_all['boot_p10'][ke][i, j] = np.percentile(np.concatenate([bootstraps_all[mem][ke][:, i, j] for mem in allmems]), 10)
                    bootstraps_all['boot_p90'][ke][i, j] = np.percentile(np.concatenate([bootstraps_all[mem][ke][:, i, j] for mem in allmems]), 90)

    t1 = datetime.now()
    print('Performed in {:10.1f} min\n'.format((t1-t0).total_seconds()/60.))
    pickle.dump(bootstraps_all, filo)

filo.close()
