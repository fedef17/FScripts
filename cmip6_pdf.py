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
import pandas as pd

from scipy.ndimage import gaussian_filter as gfilt
#############################################################################

area = 'EAT'
cart_out = '/home/fedefab/Scrivania/Research/Post-doc/lavori/CMIP6/Results_v2/EAT_NDJFM/'

erafi = '/home/fedefab/Scrivania/Research/Post-doc/lavori/ecmwf_seasonal/ERA_ref/out_ERA_DJF_EAT_4clus_4pcs_1957-2018.p'
results_ref = pickle.load(open(erafi, 'rb'))

#last20 = True

pdfssp = pickle.load(open(cart_out + 'pdfs_refCLUS_last20_{}.p'.format(area), 'rb'))

for filt in [False, True]:
    for last20 in [False, True]:
        if last20:
            ke = 'all_last20'
        else:
            ke = 'all'

        xss = np.linspace(-2500., 2500., 201)
        xi_grid, yi_grid = np.meshgrid(xss, xss)

        zi = pdfssp[('hist', ke)]
        if filt:
            zi = gfilt(zi, 10)

        zi = zi/np.sum(zi)
        levzi = np.linspace(np.percentile(zi, 10), np.percentile(zi, 90), 9)

        figs = []

        levels = np.linspace(-1e-5, 1e-5, 10)
        if filt: levels = np.linspace(-8e-6, 8e-6, 10)

        for ssp in ['ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']:
            zi2 = pdfssp[(ssp, ke)]
            if filt:
                zi2 = gfilt(zi2, 10)

            zi2 = zi2/np.sum(zi2)

            zidif = zi2 - zi

            fig = plt.figure(figsize = (16,12))
            ax = fig.add_subplot(111)

            cont = ax.contourf(xi_grid, yi_grid, zidif.reshape(xi_grid.shape), levels, cmap =
             'RdBu_r', extend = 'both' )#, linewidths = lw)
            plt.colorbar(cont)
            ax.axhline(0., color = 'grey', alpha = 0.6)
            ax.axvline(0., color = 'grey', alpha = 0.6)

            #cont = ax.contour(xi_grid, yi_grid, zi.reshape(xi_grid.shape), levzi, color = 'black')#, linewidths = lw)

            for reg, col in zip(range(4), ['blue', 'orange', 'green', 'red']):
                centr = results_ref['centroids'][reg]
                zireg = pdfssp[('hist', reg)]
                ax.scatter(centr[0], -1*centr[1], color = col, marker = 'x', s = 100)
                cmappa = ctl.custom_alphagradient_cmap(col)
                cont = ax.contour(xi_grid, yi_grid, zireg.reshape(xi_grid.shape), [0.2, 0.5], cmap = cmappa)

            fig.suptitle(ssp)
            figs.append(fig)
            fig.savefig(cart_out + 'PDF_diff_{}_{}.pdf'.format(ssp, area))

        nam = 'PDF_diff_allssps'
        if last20: nam += '_last20'
        if filt: nam += '_wfilt'
        nam += '_{}.pdf'.format(area)

        ctl.plot_pdfpages(cart_out + nam, figs)

#ax.scatter(cent[0], cent[1], color = colsim[0], s = 10, marker = 'x')

#ctl.custom_legend(fig, colsim, allsims, ncol = 3)
