#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
from matplotlib import pyplot as plt
from matplotlib import cm

import pickle
import numpy as np
from scipy import stats
import climtools_lib as ctl
import climdiags as cd
import glob

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 18
#############################################################################

cart_out = '/home/fedef/Research/lavori/CMIP6/Clusters_of_regime_hist/'
ctl.mkdir(cart_out)

for area in ['EAT', 'PNA']:
    filox = '/home/fedef/Research/lavori/CMIP6/cmip6_hist/out_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refEOF_dtr.p'.format(area)
    coso, coso_ref = ctl.load_wrtool(filox)

    lat = coso['BCC-CSM2-MR_r1i1p1f1']['lat']
    lon = coso['BCC-CSM2-MR_r1i1p1f1']['lon']

    for reg in range(4):
        cluspatterns = [coso[mo]['cluspattern'][reg] for mo in coso.keys()]

        centroids, labels, patts, repres, distances = cd.EnsClus_light(cluspatterns, lat, numclus = 4, numpcs = 4, flag_perc = False, perc = None)

        filename = cart_out + 'clu2_{}_{}.pdf'.format(area, reg)
        ctl.plot_multimap_contour(patts, lat, lon, filename, plot_margins = area, title = '{} - reg {}'.format(area, reg))
