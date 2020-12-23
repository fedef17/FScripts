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

plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 28
plt.rcParams['axes.labelsize'] = 22

#############################################################################

yr10 = 10 # length of running mean
#dtrtyp = 'light'
dtrtyp = 'histrebase'
cart_out_orig = '/home/fabiano/Research/lavori/CMIP6/Results_v6_eceens/'
ctl.mkdir(cart_out_orig)

cart_in = '/data-hobbes/fabiano/WR_CMIP6/'
#file_hist = cart_in + 'out_NEW_cmip6_hist_NDJFM_{}_4clus_4pcs_1964-2014_refCLUS_dtr_light.p'
file_hist = cart_in + 'out_eceens_hist_NDJFM_{}_4clus_4pcs_1964-2014_refCLUS_dtr.p'
file_hist_refEOF = cart_in + 'out_eceens_hist_NDJFM_{}_4clus_4pcs_1964-2014_refEOF_dtr.p'
gen_file_ssp = cart_in + 'out_eceens_{}_NDJFM_{}_4clus_4pcs_2015-2100_refCLUS_dtr_reb.p'

area = 'EAT'
ssp = 'ssp585'

reshist, resref = ctl.load_wrtool(file_hist.format(area))
reshist_re, _ = ctl.load_wrtool(file_hist.format(area))
resssp, _ = ctl.load_wrtool(gen_file_ssp.format(ssp, area))
