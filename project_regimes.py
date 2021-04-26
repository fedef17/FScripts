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

for area in ['EAT', 'PNA', 'NML']:
    for season in ['DJF', 'JJA']:
        resu, resu_ref = ctl.load_wrtool(cart+fil.format(season, area))
        gigi = resu['ERA5']

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

        for ke in ['labels', 'dist_centroid', 'pcs', 'dates']:
            cose = []
            gigi[ke] = gigi[ke][:nye*lensea+add]

            for co in gigi[ke]:
                cose.append([co]*4)

            gigi[ke] = np.concatenate(cose, axis = 0)

        for tip in ['pos', 'neg']:
            with open(filmask.format(season.lower(), tip), 'rb') as filok:
                mask = pickle.load(filok) # ci sono anche i leap!
                print(len(mask))

            nmiss = np.nansum(mask[:skip*4])
            print('missing {} events'.format(nmiss))
            mask = mask[skip*4:] # skip first 2 months (JF), these are 6-hourly data

            print(len(gigi['labels']), len(mask))

            koze = dict()
            for ke in ['cluspattern', 'cluspattern_area', 'lat', 'lat_area', 'lon', 'lon_area', 'centroids', 'model_eofs', 'model_eofs_varfrac', 'freq_clus', 'var_ratio']:
                koze[ke] = gigi[ke]

            for ke in ['labels', 'dist_centroid', 'pcs', 'dates']:
                koze[ke] = gigi[ke][mask == 1]

            reg_events[(season, tip, area)] = koze

pickle.dump(reg_events, open(cart + 'regimes_masked.p', 'wb'))

# data, datacoords, aux_info = ctl.read_xr('', regrid_to_deg = 2.5, extract_level_hPa = 500.)
#
# var_sea, dates_sea = ctl.sel_season(data, datacoords['dates'], 'DJF', cut = False)
#
# lat = resu['ERA5']['lat']
# lon = resu['ERA5']['lon']
#
# res2 = cd.WRtool_core(var_sea, lat, lon, dates_sea, 'EAT', ref_solver = resu['ERA5']['solver'], ref_patterns_area = resu['ERA5']['cluspattern_area'], run_significance_calc = False, use_reference_eofs = True, use_reference_clusters = True, ref_clusters_centers = resu['ERA5']['centroids'], climate_mean = resu['ERA5']['climate_mean'], dates_climate_mean = resu['ERA5']['climate_mean_dates'])
