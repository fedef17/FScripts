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

cart_in = '/home/fabiano/Research/lavori/CMIP6/'
cart_out_orig = cart_in + 'Results_v2_rebase/'

cart_in = '/data-hobbes/fabiano/WR_CMIP6/'
file_hist_refEOF = cart_in + 'out_NEW_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refEOF_dtr.p'
file_hist = cart_in + 'out_NEW_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refCLUS_dtr_light.p'
gen_file_ssp = cart_in + 'out_NEW_cmip6_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr_histrebase.p'

xss3 = np.linspace(-3000., 3000., 301)
xi_grid3, yi_grid3 = np.meshgrid(xss3, xss3)

ref_pdf = dict()
for area in ['EAT', 'PNA']:
    results_hist, results_ref = ctl.load_wrtool(file_hist_refEOF.format(area))

    dat1 = pd.Timestamp('09-01-1995').to_pydatetime()
    dat2 = pd.Timestamp('04-01-2014').to_pydatetime()
    okpc = results_ref['pcs'][:, :2]
    okpcok, dats = ctl.sel_time_range(okpc, results_ref['dates'], (dat1, dat2))

    kufu = ctl.calc_pdf(okpcok.T)

    zi = kufu(np.vstack([xi_grid3.flatten(), yi_grid3.flatten()]))
    zi = zi/np.max(zi)
    ref_pdf[(area, 'last20')] = zi

    kufu = ctl.calc_pdf(okpc.T)

    zi = kufu(np.vstack([xi_grid3.flatten(), yi_grid3.flatten()]))
    zi = zi/np.max(zi)
    ref_pdf[(area, 'tot50')] = zi

pickle.dump(ref_pdf, open(cart_out_orig + 'refpdf_ERA.p', 'wb'))

area = 'EAT'

for area in ['EAT', 'PNA']:
    cart_out = cart_out_orig + '{}_NDJFM/'.format(area)
    results_hist, results_ref = ctl.load_wrtool(file_hist_refEOF.format(area))

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

            for ssp in ['ssp126', 'ssp245', 'ssp370', 'ssp585']:
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
                    ax.scatter(centr[0], centr[1], color = col, marker = 'x', s = 100)
                    cmappa = ctl.custom_alphagradient_cmap(col)
                    cont = ax.contour(xi_grid, yi_grid, zireg.reshape(xi_grid.shape), [0.2, 0.5], cmap = cmappa)

                fig.suptitle(ssp)
                figs.append(fig)
                fig.savefig(cart_out + 'PDF_diff_{}_{}.pdf'.format(ssp, area))

            nam = 'PDF_diff_allssps'
            if last20:
                nam += '_last20'
            else:
                nam += '_tot50'
            if filt: nam += '_wfilt'
            nam += '_{}.pdf'.format(area)

            ctl.plot_pdfpages(cart_out + nam, figs)
            plt.close('all')

    ## Now for single models!
    modpdf = pickle.load(open(cart_out + 'modpdfs_refCLUS_{}.p'.format(area), 'rb'))
    modall = dict()
    for cos in ['hist', 'ssp126', 'ssp245', 'ssp370', 'ssp585']:
        modall[cos] = np.unique([ke[1] for ke in modpdf if ke[0] == cos])

    for filt in [False, True]:
        for last20 in [False, True]:
            if last20:
                ke = 'last20'
            else:
                ke = 'tot50'

            figall = plt.figure(figsize = (28, 12))

            xss3 = np.linspace(-3000., 3000., 301)
            xi_grid3, yi_grid3 = np.meshgrid(xss3, xss3)

            histzi = []
            for mod in modall['hist']:
                zi = modpdf[('hist', mod, ke)]
                zi = zi/np.sum(zi)
                histzi.append(zi)

            zi = np.mean(histzi, axis = 0)
            zistd = np.std(histzi, axis = 0)

            refzi = ref_pdf[(area, ke)]
            refzi = refzi/np.sum(refzi)
            zidif = zi-refzi

            if filt:
                zi = gfilt(zi, 10)
                zistd = gfilt(zistd, 10)

            levzi = np.linspace(np.percentile(zi, 10), np.percentile(zi, 90), 9)

            figs = []

            levels = np.linspace(-5e-6, 5e-6, 10)
            if filt: levels = np.linspace(-4e-6, 4e-6, 10)

            icos = 1
            ax2 = figall.add_subplot(2, 3, icos)

            if filt:
                zidif = gfilt(zidif, 10)

            contuno = ax2.contourf(xi_grid3, yi_grid3, zidif.reshape(xi_grid3.shape), levels, cmap =
             'RdBu_r', extend = 'both' )#, linewidths = lw)

            ax2.axhline(0., color = 'grey', alpha = 0.6)
            ax2.axvline(0., color = 'grey', alpha = 0.6)

            for reg, col in zip(range(4), ['blue', 'orange', 'green', 'red']):
                centr = results_ref['centroids'][reg]
                zireg = pdfssp[('hist', reg)]
                ax2.scatter(centr[0], centr[1], color = col, marker = 'x', s = 100)
                cmappa = ctl.custom_alphagradient_cmap(col)
                cont = ax2.contour(xi_grid, yi_grid, zireg.reshape(xi_grid.shape), [0.2, 0.5], cmap = cmappa)

            ax2.set_title('hist vs ref')

            for ssp in ['ssp126', 'ssp245', 'ssp370', 'ssp585']:
                icos += 1
                ax2 = figall.add_subplot(2, 3, icos)

                histzi = []
                for mod in modall[ssp]:
                    zi2 = modpdf[(ssp, mod, ke)]
                    zi2 = zi2/np.sum(zi2)
                    histzi.append(zi2)

                zi2 = np.mean(histzi, axis = 0)
                zi2std = np.std(histzi, axis = 0)

                if filt:
                    zi2 = gfilt(zi2, 10)
                    zi2std = gfilt(zi2std, 10)

                zidif = zi2 - zi
                #zidifstd = np.max([zistd, zi2std], axis = 0)

                add_hatching = np.abs(zidif) >= 0.5*zistd

                fig = plt.figure(figsize = (16,12))
                ax = fig.add_subplot(111)

                cont = ax.contourf(xi_grid3, yi_grid3, zidif.reshape(xi_grid3.shape), levels, cmap = 'RdBu_r', extend = 'both' )#, linewidths = lw)
                hatch = ax.contourf(xi_grid3, yi_grid3, add_hatching.reshape(xi_grid3.shape), levels = [0.2, 0.8], hatches = ['', '', '...'], colors = 'none')
                plt.colorbar(cont)
                ax.axhline(0., color = 'grey', alpha = 0.6)
                ax.axvline(0., color = 'grey', alpha = 0.6)

                cont = ax2.contourf(xi_grid3, yi_grid3, zidif.reshape(xi_grid3.shape), levels, cmap = 'RdBu_r', extend = 'both')
                ax2.axhline(0., color = 'grey', alpha = 0.6)
                ax2.axvline(0., color = 'grey', alpha = 0.6)
                ax2.set_title(ssp)

                #cont = ax.contour(xi_grid, yi_grid, zi.reshape(xi_grid.shape), levzi, color = 'black')#, linewidths = lw)

                for reg, col in zip(range(4), ['blue', 'orange', 'green', 'red']):
                    centr = results_ref['centroids'][reg]
                    zireg = pdfssp[('hist', reg)]
                    ax.scatter(centr[0], centr[1], color = col, marker = 'x', s = 100)
                    ax2.scatter(centr[0], centr[1], color = col, marker = 'x', s = 100)
                    cmappa = ctl.custom_alphagradient_cmap(col)
                    cont = ax.contour(xi_grid, yi_grid, zireg.reshape(xi_grid.shape), [0.2, 0.5], cmap = cmappa)
                    cont = ax2.contour(xi_grid, yi_grid, zireg.reshape(xi_grid.shape), [0.2, 0.5], cmap = cmappa)

                fig.suptitle(ssp)
                figs.append(fig)
                fig.savefig(cart_out + 'PDF_diff_models_{}_{}.pdf'.format(ssp, area))

            nam = 'PDF_diff_models'
            if last20:
                cos = '_last20'
            else:
                cos = '_tot50'
            nam += cos
            if filt: nam += '_wfilt'
            nam += '_{}.pdf'.format(area)

            ctl.plot_pdfpages(cart_out + nam, figs)

            plt.figure(figall.number)
            cax = plt.axes([0.9, 0.15, 0.03, 0.7]) #horizontal
            cb = plt.colorbar(contuno, cax=cax, orientation='vertical', format='%.1e')#, labelsize=18)
            cb.ax.tick_params(labelsize=14)
            # cb.set_label(cb_label, fontsize=16)
            plt.subplots_adjust(left=0.04, bottom=0.05, right=0.86, top=0.9, wspace=0.05, hspace=0.2)
            # cax = plt.axes([0.1, 0.11, 0.8, 0.05]) #horizontal
            # cb = plt.colorbar(contuno, cax=cax, orientation='horizontal', format='%.1e')#, labelsize=18)
            # cb.ax.tick_params(labelsize=14)
            # # cb.set_label(cb_label, fontsize=16)
            # plt.subplots_adjust(left=0.02, bottom=0.2, right=0.98, top=0.88, wspace=0.05, hspace=0.2)
            figall.savefig(cart_out + 'Allssp_PDFdiff_models_{}{}.pdf'.format(area, cos))

            plt.close('all')

    ## Now for single models element-wise
    modpdf = pickle.load(open(cart_out + 'modpdfs_refCLUS_{}.p'.format(area), 'rb'))
    modall = dict()
    for cos in ['hist', 'ssp126', 'ssp245', 'ssp370', 'ssp585']:
        modall[cos] = np.unique([ke[1] for ke in modpdf if ke[0] == cos])

    for filt in [False, True]:
        for last20 in [False, True]:
            if last20:
                ke = 'last20'
            else:
                ke = 'tot50'

            figall = plt.figure(figsize = (28, 12))

            xss3 = np.linspace(-3000., 3000., 301)
            xi_grid3, yi_grid3 = np.meshgrid(xss3, xss3)

            figs = []

            levels = np.linspace(-5e-6, 5e-6, 10)
            if filt: levels = np.linspace(-4e-6, 4e-6, 10)

            icos = 1
            ax2 = figall.add_subplot(2, 3, icos)

            refzi = ref_pdf[(area, ke)]
            refzi = refzi/np.sum(refzi)

            alldiffs = []
            for mod in modall['hist']:
                zi2 = modpdf[('hist', mod, ke)]
                zi2 = zi2/np.sum(zi2)

                zidif = zi2-refzi
                alldiffs.append(zidif)

            zidif = np.mean(alldiffs, axis = 0)
            if filt:
                zidif = gfilt(zidif, 10)

            contuno = ax2.contourf(xi_grid3, yi_grid3, zidif.reshape(xi_grid3.shape), levels, cmap =
             'RdBu_r', extend = 'both' )#, linewidths = lw)

            ax2.axhline(0., color = 'grey', alpha = 0.6)
            ax2.axvline(0., color = 'grey', alpha = 0.6)

            for reg, col in zip(range(4), ['blue', 'orange', 'green', 'red']):
                centr = results_ref['centroids'][reg]
                zireg = pdfssp[('hist', reg)]
                ax2.scatter(centr[0], centr[1], color = col, marker = 'x', s = 100)
                cmappa = ctl.custom_alphagradient_cmap(col)
                cont = ax2.contour(xi_grid, yi_grid, zireg.reshape(xi_grid.shape), [0.2, 0.5], cmap = cmappa)

            ax2.set_title('hist vs ref')

            for ssp in ['ssp126', 'ssp245', 'ssp370', 'ssp585']:
                icos += 1
                ax2 = figall.add_subplot(2, 3, icos)

                histzi = []
                signzi = []
                modok = []
                for mod in modall[ssp]:
                    if mod not in modall['hist']: continue

                    zi2 = modpdf[(ssp, mod, ke)]
                    zi2 = zi2/np.sum(zi2)

                    zi = modpdf[('hist', mod, ke)]
                    zi = zi/np.sum(zi)

                    histzi.append(zi2 - zi)
                    signzi.append(np.sign(zi2 - zi))
                    modok.append(mod)

                zidif = np.mean(histzi, axis = 0)
                signzi_all = np.sum(signzi == np.sign(zidif), axis = 0)

                add_hatching = signzi_all > 0.8*len(modall[ssp])
                print(area, filt, last20)
                print(ssp, len(modall[ssp]), np.sum(add_hatching))

                if filt:
                    zidif = gfilt(zidif, 10)

                    if ssp == 'ssp585':
                        modfigs = []
                        for ziii, mod in zip(histzi, modok):
                            fig = plt.figure(figsize = (16,12))
                            ax = fig.add_subplot(111)

                            ziii = gfilt(ziii, 10)

                            cont = ax.contourf(xi_grid3, yi_grid3, ziii.reshape(xi_grid3.shape), levels, cmap =
                             'RdBu_r', extend = 'both' )#, linewidths = lw)
                            hatch = ax.contour(xi_grid3, yi_grid3, zidif.reshape(xi_grid3.shape), levels, colors = 'k')

                            plt.colorbar(cont)
                            ax.axhline(0., color = 'grey', alpha = 0.6)
                            ax.axvline(0., color = 'grey', alpha = 0.6)

                            plt.title(mod)
                            modfigs.append(fig)

                        nam = 'MODPDF_ssp585_diff'
                        if last20:
                            nam += '_last20'
                        else:
                            nam += '_tot50'
                        if filt: nam += '_wfilt'
                        nam += '_{}.pdf'.format(area)
                        ctl.plot_pdfpages(cart_out + nam, modfigs)
                        plt.close('all')

                fig = plt.figure(figsize = (16,12))
                ax = fig.add_subplot(111)

                cont = ax.contourf(xi_grid3, yi_grid3, zidif.reshape(xi_grid3.shape), levels, cmap =
                 'RdBu_r', extend = 'both' )#, linewidths = lw)
                hatch = ax.contourf(xi_grid3, yi_grid3, add_hatching.reshape(xi_grid3.shape), levels = [0.2, 0.8], hatches = ['', '', '...'], colors = 'none')

                plt.colorbar(cont)
                ax.axhline(0., color = 'grey', alpha = 0.6)
                ax.axvline(0., color = 'grey', alpha = 0.6)

                cont = ax2.contourf(xi_grid3, yi_grid3, zidif.reshape(xi_grid3.shape), levels, cmap = 'RdBu_r', extend = 'both')
                ax2.axhline(0., color = 'grey', alpha = 0.6)
                ax2.axvline(0., color = 'grey', alpha = 0.6)
                ax2.set_title(ssp)

                #cont = ax.contour(xi_grid, yi_grid, zi.reshape(xi_grid.shape), levzi, color = 'black')#, linewidths = lw)

                for reg, col in zip(range(4), ['blue', 'orange', 'green', 'red']):
                    centr = results_ref['centroids'][reg]
                    zireg = pdfssp[('hist', reg)]
                    cmappa = ctl.custom_alphagradient_cmap(col)
                    ax.scatter(centr[0], centr[1], color = col, marker = 'x', s = 100)
                    cont = ax.contour(xi_grid, yi_grid, zireg.reshape(xi_grid.shape), [0.2, 0.5], cmap = cmappa)
                    ax2.scatter(centr[0], centr[1], color = col, marker = 'x', s = 100)
                    cont = ax2.contour(xi_grid, yi_grid, zireg.reshape(xi_grid.shape), [0.2, 0.5], cmap = cmappa)

                fig.suptitle(ssp)
                figs.append(fig)
                fig.savefig(cart_out + 'PDF_diff_singlemodels_{}_{}.pdf'.format(ssp, area))

            nam = 'PDF_diff_singlemodels'
            if last20:
                cos = '_last20'
            else:
                cos = '_tot50'
            nam += cos
            if filt: nam += '_wfilt'
            nam += '_{}.pdf'.format(area)

            ctl.plot_pdfpages(cart_out + nam, figs)

            plt.figure(figall.number)
            cax = plt.axes([0.9, 0.15, 0.03, 0.7]) #horizontal
            cb = plt.colorbar(contuno, cax=cax, orientation='vertical', format='%.1e')#, labelsize=18)
            cb.ax.tick_params(labelsize=14)
            # cb.set_label(cb_label, fontsize=16)
            plt.subplots_adjust(left=0.04, bottom=0.05, right=0.86, top=0.9, wspace=0.05, hspace=0.2)
            figall.savefig(cart_out + 'Allssp_PDFdiff_singlemodels_{}{}.pdf'.format(area, cos))

            plt.close('all')

#ax.scatter(cent[0], cent[1], color = colsim[0], s = 10, marker = 'x')

#ctl.custom_legend(fig, colsim, allsims, ncol = 3)
