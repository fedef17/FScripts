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

cart = '/home/fabiano/Research/lavori/WeatherRegimes/ERA5/'
fil = 'out_ERA5_DJFM_EAT_4clus_4pcs_allyrs.p'

resu, resu_ref = ctl.load_wrtool(cart+fil)
gigi = resu['ERA5']

nye = 33 # from 1979 to 2012
for ke in ['labels', 'dist_centroid', 'pcs', 'dates']:
    gigi[ke] = gigi[ke][:nye*90]

filmask = '/home/fabiano/Downloads/mask_1D_djf_neg.p'
with open(filmask, 'rb') as filok:
    mask = pickle.load(filok)

# ci sono anche i leap!


# data, datacoords, aux_info = ctl.read_xr('', regrid_to_deg = 2.5, extract_level_hPa = 500.)
#
# var_sea, dates_sea = ctl.sel_season(data, datacoords['dates'], 'DJF', cut = False)
#
# lat = resu['ERA5']['lat']
# lon = resu['ERA5']['lon']
#
# res2 = cd.WRtool_core(var_sea, lat, lon, dates_sea, 'EAT', ref_solver = resu['ERA5']['solver'], ref_patterns_area = resu['ERA5']['cluspattern_area'], run_significance_calc = False, use_reference_eofs = True, use_reference_clusters = True, ref_clusters_centers = resu['ERA5']['centroids'], climate_mean = resu['ERA5']['climate_mean'], dates_climate_mean = resu['ERA5']['climate_mean_dates'])
