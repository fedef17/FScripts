#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os

from matplotlib import pyplot as plt
import pickle

import climtools_lib as ctl
import climdiags as cd

####################################

#cart = '/home/fabiano/Research/lavori/WeatherRegimes/ERA5/'
cart = '/home/fedef/Research/lavori/valembo_era5/'
fil = 'out_ERA5_{}_{}_4clus_55perc_allyrs.p'
filmask = cart + 'mask_1D_{}_{}.p'

reg_events = dict()
regime_ref = dict()

figs_scatter = []
figs_clouds = []

for area in ['EAT', 'PNA', 'NML']:
    for season in ['DJF', 'JJA']:
        resu, resu_ref = ctl.load_wrtool(cart+fil.format(season, area)) # carico i regimi del wrtool
        gigi = resu['ERA5']



        for ke in ['labels', 'dist_centroid', 'pcs', 'dates']:
            cose = []
            gigi[ke] = gigi[ke][:nye*lensea+add]

            for co in gigi[ke]:
                cose.append([co]*4) # moltiplico per 4 per avere lo stesso numero dei dati 6hr (qui stiamo approssimando ma credo cambi veramente poco)

            gigi[ke] = np.concatenate(cose, axis = 0)

        koze_ref = dict()
        # salvo un po' di cose generali dei clusters
        for ke in ['cluspattern', 'cluspattern_area', 'lat', 'lat_area', 'lon', 'lon_area', 'centroids', 'model_eofs', 'model_eofs_varfrac', 'freq_clus', 'var_ratio', 'labels', 'dist_centroid', 'pcs', 'dates']:
            koze_ref[ke] = gigi[ke]

        regime_ref[(season, area)] = koze_ref

        for tip in ['pos', 'neg']:
            print(area, season, tip)
            with open(filmask.format(season.lower(), tip), 'rb') as filok:
                mask = pickle.load(filok)
                print(len(mask)) # ci sono anche i leap!

            nmiss = np.nansum(mask[:skip*4])
            print('missing {} events'.format(nmiss))
            mask = mask[skip*4:] # skip first 2 months (JF), these are 6-hourly data
            print('ok events: {}'.format(np.nansum(mask)))

            print(len(gigi['labels']), len(mask))

            koze = dict()
            for ke in ['lat_area', 'lon_area', 'model_eofs']:
                koze[ke] = koze_ref[ke]

            # salvo le serie temporali filtrate e non (_all)
            for ke in ['labels', 'dist_centroid', 'pcs', 'dates']:
                koze[ke] = gigi[ke][mask == 1]

            # ricostruisco il pattern degli stati selezionati
            patts = ctl.reconstruct_from_pcs(koze['pcs'], koze['model_eofs'])
            cluspatt_ok = []
            for reg in range(4):
                okpo = koze['labels'] == reg
                cluspatt_ok.append(np.mean(patts[okpo, ...], axis = 0))

            koze['cluspattern_area'] = np.stack(cluspatt_ok)
            koze['freq_clus'] = ctl.calc_clus_freq(koze['labels'], 4)
            koze['centroids'] = ctl.calc_effective_centroids(koze['pcs'], koze['labels'], 4)

            cd.plot_regimes(koze['lat_area'], koze['lon_area'], koze['cluspattern_area'], cart + 'regpatt_{}_{}_{}.pdf'.format(season, area, tip), plot_type = 'pcolormesh', clatlo = (89,0))
            print(koze['lon_area'])

            reg_events[(season, tip, area)] = koze

        #### scatter plot piano EOF0/1
        fig = plt.figure()

        coso_pos = reg_events[(season, 'pos', area)]
        coso_neg = reg_events[(season, 'neg', area)]

        colors = ctl.color_set(4)
        for ii, col in enumerate(colors):
            # eventi positivi
            okpo = coso_pos['labels'] == ii
            plt.scatter(coso_pos['pcs'][okpo, 0], coso_pos['pcs'][okpo, 1], color = col, s = 10, marker = 'o')

            # eventi negativi
            okpo = coso_neg['labels'] == ii
            plt.scatter(coso_neg['pcs'][okpo, 0], coso_neg['pcs'][okpo, 1], color = col, s = 10, marker = 'x')

            # centroidi globali
            plt.scatter(coso_pos['centroids'][ii][0], coso_pos['centroids'][ii][1], color = col, s = 100, label = str(ii))

        plt.legend()
        plt.title(area + ' - ' + season)
        #fig.savefig(cart + 'scatter_{}_{}.pdf'.format(season, area))
        figs_scatter.append(fig)

        #### regime clouds
        clouds = dict()
        clouds['reference'] = resu['ERA5']
        clouds['pos'] = reg_events[(season, 'pos', area)]
        clouds['neg'] = reg_events[(season, 'neg', area)]

        xlims = (np.percentile(clouds['reference']['pcs'], 1), np.percentile(clouds['reference']['pcs'], 99))

        fig = ctl.plot_multimodel_regime_pdfs(clouds, model_names = ['reference', 'pos', 'neg'], eof_proj = [(0,1), (2,3)], reference = 'reference', check_for_eofs = False, colors = ['black', 'forestgreen', 'indianred'], eof_axis_lim = xlims)
        fig.suptitle(area + ' - ' + season)
        figs_clouds.append(fig)

ctl.plot_pdfpages(cart + 'scatter_regimes.pdf', figs_scatter)
ctl.plot_pdfpages(cart + 'clouds_regimes.pdf', figs_clouds)

pickle.dump(reg_events, open(cart + 'regimes_masked.p', 'wb'))
pickle.dump(regime_ref, open(cart + 'regimes_ref.p', 'wb'))

# data, datacoords, aux_info = ctl.read_xr('', regrid_to_deg = 2.5, extract_level_hPa = 500.)
#
# var_sea, dates_sea = ctl.sel_season(data, datacoords['dates'], 'DJF', cut = False)
#
# lat = resu['ERA5']['lat']
# lon = resu['ERA5']['lon']
#
# res2 = cd.WRtool_core(var_sea, lat, lon, dates_sea, 'EAT', ref_solver = resu['ERA5']['solver'], ref_patterns_area = resu['ERA5']['cluspattern_area'], run_significance_calc = False, use_reference_eofs = True, use_reference_clusters = True, ref_clusters_centers = resu['ERA5']['centroids'], climate_mean = resu['ERA5']['climate_mean'], dates_climate_mean = resu['ERA5']['climate_mean_dates'])
