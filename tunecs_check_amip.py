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
import glob

import tunlib as tl

from scipy.optimize import Bounds, minimize, least_squares
import itertools as itt
from multiprocessing import Process, Queue, Pool
import random
#from multiprocessing import set_start_method
#from multiprocessing import get_context
#set_start_method("spawn")


plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
titlefont = 24
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 18

################################################################################

cart_in = '/data-hobbes/fabiano/TunECS/AMIP_exps/'
cart_out = '/home/fabiano/Research/lavori/TunECS/tuning/experiments/analysis/'
cart_out_2 = '/home/fabiano/Research/lavori/TunECS/tuning/coupled/'

testparams = ['ENTRORG', 'RPRCON', 'DETRPEN', 'RMFDEPS', 'RVICE', 'RSNOWLIN2', 'RCLDIFF', 'RLCRIT_UPHYS']
nums = np.arange(8)
letts = 'a b c d e f g h'.split()

uff_params = dict()
uff_params['RPRCON'] = 1.34E-3
uff_params['RVICE'] = 0.137
uff_params['RLCRITSNOW'] = 4.0E-5
uff_params['RSNOWLIN2'] = 0.035
uff_params['ENTRORG'] = 1.70E-4
uff_params['DETRPEN'] = 0.75E-4
uff_params['ENTRDD'] = 3.0E-4
uff_params['RMFDEPS'] = 0.3
uff_params['RCLDIFF'] = 3.E-6
uff_params['RCLDIFFC'] = 5.0
uff_params['RLCRIT_UPHYS'] = 0.875e-5

bounds = (np.array([1e-4, 9e-4, 0.4e-4, 0.1, 0.08, 0.015, 1.e-6, 0.7e-5]), np.array([2.4e-4, 2.e-3, 1.3e-4, 0.5, 0.2, 0.055, 5.e-6, 1.1e-5]))

##########################################

[parsets, pich, c4ch] = pickle.load(open(cart_out_2 + 'altparams.p', 'rb'))
with open(cart_out + 'der_sensmat_global.p', 'rb') as filox:
    [resdic_mean, resdic_err, derdic, derdic_err] = pickle.load(filox)

ctrl_pi = resdic_mean[('pi', 0, 0, 'toa_net')]
ctrl_pi_err = resdic_err[('pi', 0, 0, 'toa_net')]
ctrl_c4 = resdic_mean[('c4', 0, 0, 'toa_net')]
ctrl_c4_err = resdic_err[('c4', 0, 0, 'toa_net')]

parsets = parsets[1:]
pich = pich[1:]
c4ch = c4ch[1:]

carttab = '/home/fabiano/Research/lavori/TunECS/tuning/experiments/table/'

fig1, ax1 = plt.subplots(figsize=(16,12))
fig2, ax2 = plt.subplots(figsize=(16,12))
fig3, ax3 = plt.subplots(figsize=(16,12))
fig4, ax4 = plt.subplots(figsize=(16,12))
for parset, chanpi, chanc4, num in zip(parsets, pich, c4ch, range(1, 10)):
    expnam = 'pia{}'.format(num)
    pi_score = tl.read_PI(carttab, expnam)
    toa_pi, err_toa_pi = tl.read_toa_net(carttab, expnam)

    expnam = 'c4a{}'.format(num)
    toa_c4, err_toa_c4 = tl.read_toa_net(carttab, expnam)

    labp = None
    laba = None
    if num == 1:
        labp = 'predicted'
        laba = 'actual'

    ax1.scatter(num, chanpi, color = 'steelblue', marker = 'o', s = 100, label = labp)
    ax1.errorbar(num, toa_pi - ctrl_pi, yerr = err_toa_pi, fmt = 'none', color = 'indianred', capsize = 2, elinewidth = 1)
    ax1.scatter(num, toa_pi - ctrl_pi, color = 'indianred', marker = 'x', s = 100, label = laba)

    ax2.scatter(num, chanc4, color = 'steelblue', marker = 'o', s = 100, label = labp)
    ax2.errorbar(num, toa_c4 - ctrl_c4, yerr = err_toa_c4, fmt = 'none', color = 'indianred', capsize = 2, elinewidth = 1)
    ax2.scatter(num, toa_c4 - ctrl_c4, color = 'indianred', marker = 'x', s = 100, label = laba)

    ax3.scatter(num, chanc4-chanpi, color = 'steelblue', marker = 'o', s = 100, label = labp)
    ax3.errorbar(num, (toa_c4 - ctrl_c4) - (toa_pi - ctrl_pi), yerr = (err_toa_pi+err_toa_c4)/2., fmt = 'none', color = 'indianred', capsize = 2, elinewidth = 1)
    ax3.scatter(num, (toa_c4 - ctrl_c4) - (toa_pi - ctrl_pi), color = 'indianred', marker = 'x', s = 100, label = laba)

    ax4.scatter(num, pi_score, color = 'steelblue', marker = 'o', s = 100)

pi_score = tl.read_PI(carttab, 'tpa1')
ax4.axhline(pi_score)

for ax in [ax1, ax2, ax3]:
    ax.grid()
    ax.axhline(0.)
    ax.legend()
ax1.set_ylabel('Change in pi toa_net (W/m2)')
ax2.set_ylabel('Change in 4xCO2 toa_net (W/m2)')
ax3.set_ylabel('Change in sensitivity at toa_net (W/m2)')

fig1.savefig(cart_out_2 + 'check_change_pi.pdf')
fig2.savefig(cart_out_2 + 'check_change_c4.pdf')
fig3.savefig(cart_out_2 + 'check_change_sens.pdf')
fig4.savefig(cart_out_2 + 'check_pi_score.pdf')
