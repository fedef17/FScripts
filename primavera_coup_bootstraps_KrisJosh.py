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
with open(cart_out + 'residual_pc_clustering_as_dict_v2.pkl', 'rb') as fil:
    nures = pickle.load(fil)
with open(cart_out + 'raw_pc_clustering_as_dict_v2.pkl', 'rb') as fil:
    nures_nofil = pickle.load(fil)

for numclus in [3,4,5]:
    refres = nures['reference']
    refres_nofil = nures_nofil['reference']
    refcen = refres[numclus]['centroids'][:,:4]#/9.81 # these should be the counterparts of the original regimes in the residual phase space
    refpcs = refres['pcs']
    reflab = refres[numclus]['states']
    refdat = refres['time']

    ##############################################################################

    results_old, results_ref = ctl.load_wrtool(filogen)

    results_ref['dates'] = refdat
    results_ref['pcs'] = refpcs[:, :4]
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

                perm = ctl.match_pc_sets(refcen, nures[ke][numclus]['centroids'])
                centroids, labels = ctl.change_clus_order(nures[ke][numclus]['centroids'], nures[ke][numclus]['states'], perm)

                results[ke]['labels'] = labels
                results[ke]['centroids'] = centroids
            else:
                results[ke] = results_old[ke]
                results[ke]['pcs'] = nures[ke]['pcs']
                if 'HadGEM3-GC31' in ke:
                    results[ke]['dates'] = ctl.adjust_360day_dates(nures[ke]['time'])
                elif ke in ['CNRM-CM6-1-HR_r1i1p1f2', 'EC-Earth3P-HR_r1i1p1f1']:
                    results[ke]['dates'] = nures[ke]['time']
                else:
                    print(ke)
                    print(len(nures[ke]['time']), len(results_old[ke]['dates']))
                    raise ValueError('no rule')

                perm = ctl.match_pc_sets(refcen, nures[ke][numclus]['centroids'])
                centroids, labels = ctl.change_clus_order(nures[ke][numclus]['centroids'], nures[ke][numclus]['states'], perm)

                results[ke]['labels'] = labels
                results[ke]['centroids'] = centroids
        else:
            print('{} NOT FOUND!!!! Skipping.....'.format(ke))

    for ke in results:
        results[ke]['var_ratio'] = ctl.calc_varopt_molt(results[ke]['pcs'], results[ke]['centroids'], results[ke]['labels'])

        results[ke]['freq_clus'] = ctl.calc_clus_freq(results[ke]['labels'], numclus)

        rms, patcor = ctl.calc_RMS_and_patcor(refcen, results[ke]['centroids'])
        results[ke]['RMS'] = rms
        results[ke]['patcor'] = patcor

        results[ke]['resid_times'] = ctl.calc_regime_residtimes(results[ke]['labels'], dates = results[ke]['dates'])[0]

    results_ref['var_ratio'] = ctl.calc_varopt_molt(results_ref['pcs'], results_ref['centroids'], results_ref['labels'])
    results_ref['freq_clus'] = ctl.calc_clus_freq(results_ref['labels'], numclus)
    results_ref['resid_times'] = ctl.calc_regime_residtimes(results_ref['labels'], dates = results_ref['dates'])[0]

    pickle.dump([results, results_ref], open(cart_out + 'out_prima_coup_v7_DJF_EAT_4clus_4pcs_1957-2014_refEOF_FILTEREDKJ_k{}.p'.format(numclus), 'wb'))
    #sys.exit()

    results['ERA_0'] = results_ref


    n_choice = 30
    n_bootstrap = 100

    filo = open(cart_out + 'res_bootstrap_v7_KJ_ref_k{}.p'.format(numclus), 'wb')

    model_names = ['AWI-CM-1-1-LR', 'AWI-CM-1-1-HR', 'CMCC-CM2-HR4', 'CMCC-CM2-VHR4', 'CNRM-CM6-1', 'CNRM-CM6-1-HR', 'EC-Earth3P', 'EC-Earth3P-HR', 'ECMWF-IFS-LR', 'ECMWF-IFS-MR', 'ECMWF-IFS-HR', 'MPI-ESM1-2-HR', 'MPI-ESM1-2-XR', 'HadGEM3-GC31-LL', 'HadGEM3-GC31-MM', 'HadGEM3-GC31-HM', 'HadGEM3-GC31-HH']

    vers = ['LR', 'HR', 'LR', 'HR', 'LR', 'HR', 'LR', 'HR', 'LR', 'MR', 'HR', 'LR', 'HR', 'LR', 'MR', 'MR', 'HR'] + ['OBS']
    model_coups = ['AWI-CM-1-1', 'CMCC-CM2', 'CNRM-CM6-1', 'EC-Earth3P', 'ECMWF-IFS', 'MPI-ESM1-2', 'HadGEM3-GC31']

    model_names_all = model_names + ['ERA']

    colors = ctl.color_set(len(model_coups)*2, sns_palette = 'Paired')
    mr_cos = np.where(np.array(vers) == 'MR')[0]
    for gi in mr_cos:
        colors.insert(gi, np.mean(colors[gi-1:gi+1], axis = 0))

    colors_wERA = colors + ['black']

    #regnam = ['NAO +', 'Sc. BL', 'AR', 'NAO -']
    regnam = ['clus_{}'.format(i) for i in range(numclus)]

    ################################################################################
    all_mods = np.array([ke.split('_')[0] for ke in results.keys()] + ['ERA'])
    all_mems = np.array([ke.split('_')[1] for ke in results.keys()] + ['0'])

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
    for reg in range(numclus):
        okclu = labels == reg
        okpc = pcs[okclu, :]
        kufu = ctl.calc_pdf(okpc[:, :3].T)
        zi = kufu(np.vstack([xi_grid.flatten(), yi_grid.flatten(), zi_grid.flatten()]))
        zi = zi/np.max(zi)

        zi_ref.append(zi)

    allkeysss = ['varopt', 'autocorr', 'freq', 'dist_cen', 'resid_times_av', 'resid_times_90', 'centroids', 'trans_matrix', 'relative_entropy', 'patcor']

    #for mod in model_names_all:
    for mod in ['ERA']:
        print(mod)
        if mod not in all_mods:
            print('Skipping....')
            continue
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
                for iclu in range(numclus):
                    okla = labels == iclu
                    centroids.append(np.mean(pcs[okla], axis = 0))
                centroids = np.stack(centroids)

                #sig = ctl.clusters_sig(pcs, centroids, labels, dates, nrsamp = 200)
                # if sig < 10:
                #     sig2 = ctl.clusters_sig(pcs, centroids, labels, dates, nrsamp = 1000)
                #     print('RECHECK ', sig, sig2, sig3)

                varopt = ctl.calc_varopt_molt(pcs, centroids, labels)
                autocorr = ctl.calc_autocorr_wlag(pcs, dates, out_lag1 = True)
                #bootstraps['significance'].append(sig)
                bootstraps['varopt'].append(varopt)
                bootstraps['autocorr'].append(autocorr)

                bootstraps['freq'].append(ctl.calc_clus_freq(labels, numclus))

                centdist = np.array([ctl.distance(centroids[iclu], ref_cen[iclu]) for iclu in range(numclus)])
                bootstraps['dist_cen'].append(centdist)
                bootstraps['centroids'].append(centroids)

                resid_times = ctl.calc_regime_residtimes(labels, dates = dates)[0]
                av_res = np.array([np.mean(resid_times[reg]) for reg in range(numclus)])
                av_res_90 = np.array([np.percentile(resid_times[reg], 90) for reg in range(numclus)])
                bootstraps['resid_times_av'].append(av_res)
                bootstraps['resid_times_90'].append(av_res_90)

                bootstraps['trans_matrix'].append(ctl.calc_regime_transmatrix(1, labels, dates))

                # relative entropy, RMS, patcor
                relent_all = []
                for reg in range(numclus):
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
                # for reg in range(numclus):
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
                # centdist = np.array([ctl.distance(filt_centroids[iclu], ref_cen_filt[iclu]) for iclu in range(numclus)])
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
                bootstraps_all['boot_p10'][ke] = np.array([np.percentile(np.concatenate([bootstraps_all[mem][ke][:, reg] for mem in allmems]), 10) for reg in range(numclus)])
                bootstraps_all['boot_p90'][ke] = np.array([np.percentile(np.concatenate([bootstraps_all[mem][ke][:, reg] for mem in allmems]), 90) for reg in range(numclus)])
            else:
                print(ke)
                if ke == 'centroids':
                    J2 = 4
                elif ke == 'trans_matrix':
                    J2 = numclus
                else:
                    raise ValueError('not specified')
                bootstraps_all['boot_p10'][ke] = np.zeros((numclus, J2))
                bootstraps_all['boot_p90'][ke] = np.zeros((numclus, J2))
                for i in range(numclus):
                    for j in range(J2):
                        bootstraps_all['boot_p10'][ke][i, j] = np.percentile(np.concatenate([bootstraps_all[mem][ke][:, i, j] for mem in allmems]), 10)
                        bootstraps_all['boot_p90'][ke][i, j] = np.percentile(np.concatenate([bootstraps_all[mem][ke][:, i, j] for mem in allmems]), 90)

        t1 = datetime.now()
        print('Performed in {:10.1f} min\n'.format((t1-t0).total_seconds()/60.))
        pickle.dump(bootstraps_all, filo)

    filo.close()
