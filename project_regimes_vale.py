#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os

from matplotlib import pyplot as plt
import pickle

import climtools_lib as ctl
import climdiags as cd

from scipy import stats

####################################

#cart = '/home/fabiano/Research/lavori/WeatherRegimes/ERA5/'
cart = '/home/fedef/Research/lavori/valembo_era5/'
filmask = cart + 'mask_1D_{}_{}.p'

ofidpc = open(cart + 'rel_delta_pcs.txt', 'w')

reg_events = dict()

figs_scatter = []
figs_clouds = []

with open(cart+'regimes_ref.p', 'rb') as figi:
    regimes_ref = pickle.load(figi)

ttests = dict()
deltadist_all = dict()

for area in ['EAT', 'PNA', 'NML']:
    for season in ['DJF', 'JJA']:
        koze_ref = regimes_ref[(season, area)]

        # conto i giorni
        if season == 'DJF':
            nye = 33 # from 1979 to 2012
            lensea = 90
            add = 9+31
            skip = 59
        elif season == 'JJA':
            nye = 34
            lensea = 92
            add = 0
            skip = 0

        for tip in ['pos', 'neg']:
            print(area, season, tip)
            with open(filmask.format(season.lower(), tip), 'rb') as filok:
                mask = pickle.load(filok)
                print(len(mask)) # ci sono anche i leap!

            nmiss = np.nansum(mask[:skip*4])
            print('missing {} events'.format(nmiss))
            mask = mask[skip*4:] # skip first 2 months (JF), these are 6-hourly data
            print('ok events: {}'.format(np.nansum(mask)))

            print(len(koze_ref['labels']), len(mask))

            koze = dict()

            for ke in ['lat_area', 'lon_area', 'model_eofs']:
                koze[ke] = koze_ref[ke]

            # salvo le serie temporali filtrate e non (_all)
            for ke in ['labels', 'dist_centroid', 'pcs', 'dates']:
                koze[ke] = koze_ref[ke][mask == 1]

            # ricostruisco il pattern degli stati selezionati
            patts = ctl.reconstruct_from_pcs(koze['pcs'], koze['model_eofs'])
            cluspatt_ok = []
            for reg in range(4):
                okpo = koze['labels'] == reg
                cluspatt_ok.append(np.mean(patts[okpo, ...], axis = 0))

            koze['cluspattern_area'] = np.stack(cluspatt_ok)
            koze['freq_clus'] = ctl.calc_clus_freq(koze['labels'], 4)
            koze['centroids'] = ctl.calc_effective_centroids(koze['pcs'], koze['labels'], 4)

            if season == 'DJF':
                cbar_range = (-140., 140.)
            else:
                cbar_range = (-70., 70.)

            if area == 'EAT':
                clatlo = (65, -20)
            elif area == 'PNA':
                clatlo = (65, -140)
            else:
                clatlo = (89, 0)

            cd.plot_regimes(koze['lat_area'], koze['lon_area'], koze['cluspattern_area'], cart + 'regpatt_{}_{}_{}.pdf'.format(season, area, tip), clatlo = clatlo, cbar_range = cbar_range, draw_contour_lines = False, cmappa = 'RdBu_r', n_color_levels = 21)
            print(koze['lon_area'])

            reg_events[(season, tip, area)] = koze

        #### scatter plot piano EOF0/1
        fig = plt.figure()

        coso_pos = reg_events[(season, 'pos', area)]
        coso_neg = reg_events[(season, 'neg', area)]

        # check change in pcs
        pcpos = coso_pos['pcs'].mean(axis = 0)
        pcneg = coso_neg['pcs'].mean(axis = 0)
        pcref = koze_ref['pcs'].mean(axis = 0)

        testyp = []
        testyn = []
        for pcp, pcn, pcr in zip(coso_pos['pcs'].T, coso_neg['pcs'].T, koze_ref['pcs'].T):
            testyp.append(stats.ttest_ind(pcp, pcr, equal_var = False)[1])
            testyn.append(stats.ttest_ind(pcn, pcr, equal_var = False)[1])
        ttests[(area, season, 'pos')] = np.stack(testyp)
        ttests[(area, season, 'neg')] = np.stack(testyn)

        deltapos = (pcpos-pcref)/koze_ref['pcs'].std(axis = 0)
        deltaneg = (pcneg-pcref)/koze_ref['pcs'].std(axis = 0)
        deltadist_all[(area, season, 'pos')] = deltapos
        deltadist_all[(area, season, 'neg')] = deltaneg

        ofidpc.write('------- {} - {} ----------\n'.format(area, season))
        ofidpc.write(('pcs deltapos: '+'{:10.2e}'*len(deltapos) + '\n').format(*deltapos))
        ofidpc.write(('pcs deltaneg: '+'{:10.2e}'*len(deltaneg) + '\n').format(*deltaneg))
        ofidpc.write('---------------------------\n')

        colors = ctl.color_set(4)
        for ii, col in enumerate(colors):
            # eventi positivi
            okpo = coso_pos['labels'] == ii
            plt.scatter(coso_pos['pcs'][okpo, 0], coso_pos['pcs'][okpo, 1], color = col, s = 10, marker = 'o')

            # eventi negativi
            okpo = coso_neg['labels'] == ii
            plt.scatter(coso_neg['pcs'][okpo, 0], coso_neg['pcs'][okpo, 1], color = col, s = 10, marker = 'x')

            # centroidi globali
            plt.scatter(koze_ref['centroids'][ii][0], koze_ref['centroids'][ii][1], color = col, s = 100, label = str(ii))

        plt.legend()
        plt.title(area + ' - ' + season)
        #fig.savefig(cart + 'scatter_{}_{}.pdf'.format(season, area))
        figs_scatter.append(fig)

        # #### regime clouds
        # clouds = dict()
        # clouds['reference'] = koze_ref
        # clouds['pos'] = coso_pos
        # clouds['neg'] = coso_neg
        #
        # xlims = (np.percentile(clouds['reference']['pcs'], 1), np.percentile(clouds['reference']['pcs'], 99))
        #
        # fig = ctl.plot_multimodel_regime_pdfs(clouds, model_names = ['reference', 'pos', 'neg'], eof_proj = [(0,1), (2,3)], reference = 'reference', check_for_eofs = False, colors = ['black', 'forestgreen', 'indianred'], eof_axis_lim = xlims)
        # fig.suptitle(area + ' - ' + season)
        # figs_clouds.append(fig)

ctl.plot_pdfpages(cart + 'scatter_regimes.pdf', figs_scatter)
ctl.plot_pdfpages(cart + 'clouds_regimes.pdf', figs_clouds)

pickle.dump(reg_events, open(cart + 'regimes_masked.p', 'wb'))

pickle.dump([ttests, deltadist_all], open(cart + 'delta_pcs.p', 'wb'))

ofidpc.close()

ofidpc = open(cart + 'rel_delta_pcs_latex.txt', 'w')

from matplotlib import colors as mcolors

colo = '#d73027 #f46d43 #fdae61 #fee090 #ffffff #e0f3f8 #abd9e9 #74add1 #4575b4'
colo = colo.split()
colo = colo[::-1]
cmappa = mcolors.ListedColormap(colo)
cmappa.set_over('#800026') #662506
cmappa.set_under('#023858') #542788

fig, axs = plt.subplots(2, 2, figsize = (12, 12))

for ise, season in enumerate(['DJF', 'JJA']):
    all_pos = []
    all_neg = []
    sig_pos = []
    sig_neg = []
    for area in ['EAT', 'PNA', 'NML']:
        #ofidpc.write('------- {} - {} ----------\n'.format(area, season))
        deltapos = deltadist_all[(area, season, 'pos')]
        all_pos.append(deltapos[:4])
        ofidpc.write(('{:4s}, {:4s}, pos: &'.format(area,season) + '{:10.2e} & '*(len(deltapos)-1) + '{:10.2e}' + '\\ \n').format(*deltapos))
        tist = ['*' if ti < 0.05 else '' for ti in ttests[(area, season, 'pos')]]
        sig_pos.append(ttests[(area, season, 'pos')][:4])
        ofidpc.write(('{:15s} & '.format('')+'{:10s} &'*(len(deltapos)-1) +'{:10s}' + '\\ \n').format(*tist))

        deltaneg = deltadist_all[(area, season, 'neg')]
        all_neg.append(deltaneg[:4])
        ofidpc.write(('{:4s}, {:4s}, neg: &'.format(area,season) + '{:10.2e} & '*(len(deltaneg)-1) + '{:10.2e}' + '\\ \n').format(*deltaneg))
        tist = ['*' if ti < 0.05 else '' for ti in ttests[(area, season, 'neg')]]
        sig_neg.append(ttests[(area, season, 'neg')][:4])
        ofidpc.write(('{:15s} & '.format('') + '{:10s} &'*(len(deltaneg)-1) +'{:10s}' + '\\ \n').format(*tist))
        ofidpc.write('\hline \n')
        #ofidpc.write('---------------------------\n')

    ax = axs[ise, 0]
    all_pos = np.stack(all_pos)
    sig_pos = np.stack(sig_pos)
    ext = [0, 4, 0, 3]
    gigifig = ax.imshow(all_pos, vmin = -0.1125, vmax = 0.1125, cmap = cmappa, origin = 'lower',  extent = ext, aspect = 1)

    for ix in range(4):
        for iy in range(3):
            if sig_pos[iy, ix] < 0.01:
                ax.scatter(ix+0.5, iy+0.5, s=60, edgecolors = 'black', facecolors='white')
            elif sig_pos[iy, ix] < 0.05:
                ax.scatter(ix+0.5, iy+0.5, s=20, edgecolors = 'black', facecolors='white')

    # ora il neg
    ax = axs[ise, 1]
    all_neg = np.stack(all_neg)
    sig_neg = np.stack(sig_neg)
    ext = [0, 4, 0, 3]
    gigifig = ax.imshow(all_neg, vmin = -0.1125, vmax = 0.1125, cmap = cmappa, origin = 'lower',  extent = ext, aspect = 1)

    for ix in range(4):
        for iy in range(3):
            if sig_neg[iy, ix] < 0.01:
                ax.scatter(ix+0.5, iy+0.5, s=60, edgecolors = 'black', facecolors='white')
            elif sig_neg[iy, ix] < 0.05:
                ax.scatter(ix+0.5, iy+0.5, s=20, edgecolors = 'black', facecolors='white')

    axs[ise, 0].set_yticks(0.5+np.arange(3), minor = False)
    axs[ise, 0].set_yticklabels(['EAT', 'PAC', 'NH'], va='center')


axs[1, 0].set_xticks(0.5+np.arange(4), minor = False)
axs[1, 0].set_xticklabels(['pc{}'.format(nu) for nu in range(1,5)], ha='center')
axs[1, 1].set_xticks(0.5+np.arange(4), minor = False)
axs[1, 1].set_xticklabels(['pc{}'.format(nu) for nu in range(1,5)], ha='center')

axs[0, 0].set_xticks([])
axs[0, 1].set_xticks([])
axs[0, 1].set_yticks([])
axs[1, 1].set_yticks([])

# axs[0, 0].set_title('POS')
# axs[0, 1].set_title('NEG')

for ax in axs.flatten():
    for co in np.arange(1, 4):
        ax.axvline(co, color = 'white', linewidth = 0.1)
    for co in np.arange(1, 3):
        ax.axhline(co, color = 'white', linewidth = 0.1)

showdate = fig.text(0.3, 0.9, 'POS', rotation = 'horizontal', color = 'black', fontsize = 20)
showdate = fig.text(0.75, 0.9, 'NEG', rotation = 'horizontal', color = 'black', fontsize = 20)

showdate = fig.text(0.02, 0.35, 'JJA', rotation = 'vertical', color = 'black', fontsize = 20)

showdate = fig.text(0.02, 0.35, 'JJA', rotation = 'vertical', color = 'black', fontsize = 20)#, bbox=dict(facecolor='lightsteelblue', edgecolor='black', boxstyle='round,pad=1'))
showdate = fig.text(0.02, 0.7, 'DJF', rotation = 'vertical', color = 'black', fontsize = 20)

#cax = fig.add_subplot(gs[6, :])
cax = plt.axes([0.1, 0.1, 0.8, 0.05])
cb = plt.colorbar(gigifig, cax=cax, orientation='horizontal', extend = 'both')
cb.ax.tick_params(labelsize=18)
cb.set_label('Fractional change in PC mean', fontsize=20)
plt.subplots_adjust(left=0.1, bottom=0.2, right=0.98, top=0.86, wspace=0.05, hspace=0.1)

fig.savefig(cart + 'pc_shift.pdf')

ofidpc.close()



# data, datacoords, aux_info = ctl.read_xr('', regrid_to_deg = 2.5, extract_level_hPa = 500.)
#
# var_sea, dates_sea = ctl.sel_season(data, datacoords['dates'], 'DJF', cut = False)
#
# lat = resu['ERA5']['lat']
# lon = resu['ERA5']['lon']
#
# res2 = cd.WRtool_core(var_sea, lat, lon, dates_sea, 'EAT', ref_solver = resu['ERA5']['solver'], ref_patterns_area = resu['ERA5']['cluspattern_area'], run_significance_calc = False, use_reference_eofs = True, use_reference_clusters = True, ref_clusters_centers = resu['ERA5']['centroids'], climate_mean = resu['ERA5']['climate_mean'], dates_climate_mean = resu['ERA5']['climate_mean_dates'])
