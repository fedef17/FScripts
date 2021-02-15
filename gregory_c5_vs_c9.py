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

cart_in = '/home/fabiano/Research/ecearth/TunECS/'
cart_out = '/home/fabiano/Research/lavori/TunECS/results/'

exps = ['pic0', 'pic5', 'pic6', 'pic8', 'pic9', 'c4c5', 'c4c9']
colors = ctl.color_set(5)
colors += [colors[1], colors[-1]]
markers = 5*['o'] + 2*['D']

fig, ax = plt.subplots(figsize=(16,12))
for exp, col, mar in zip(exps, colors, markers):
    gregfi = cart_in + '{}/ecmean/gregory_{}.txt'.format(exp, exp)
    anni, toa_net, srf_net, tas = tl.read_gregory(gregfi)

    toa5 = []
    tas5 = []
    for i in range(0, len(anni), 5):
        if len(anni) - i < 2: continue
        toa5.append(np.mean(toa_net[i:i+5]))
        tas5.append(np.mean(tas[i:i+5]))

    ax.scatter(tas5, toa5, color = col, marker = mar, label = exp)
    ax.plot(tas5, toa5, color = col, linewidth = 0.5)

ax.grid()
ax.legend()
fig.savefig(cart_out + 'gregory_5y.pdf')


fig, ax = plt.subplots(figsize=(16,12))
for exp, col, mar in zip(exps[:5], colors[:5], markers[:5]):
    gregfi = cart_in + '{}/ecmean/gregory_{}.txt'.format(exp, exp)
    anni, toa_net, srf_net, tas = tl.read_gregory(gregfi)

    toa5 = []
    tas5 = []
    for i in range(0, len(anni), 5):
        if len(anni) - i < 2: continue
        toa5.append(np.mean(toa_net[i:i+5]))
        tas5.append(np.mean(tas[i:i+5]))

    ax.scatter(tas5[1:-1], toa5[1:-1], color = col, marker = mar, label = exp)
    ax.scatter(tas5[0], toa5[0], color = col, marker = '>')
    ax.scatter(tas5[-1], toa5[-1], color = col, marker = '<')
    ax.plot(tas5, toa5, color = col, linewidth = 0.5)

ax.grid()
ax.legend()
fig.savefig(cart_out + 'gregory_5y_pi.pdf')

fig, ax = plt.subplots(figsize=(16,12))
for exp, col, mar in zip(exps[5:], colors[5:], markers[5:]):
    gregfi = cart_in + '{}/ecmean/gregory_{}.txt'.format(exp, exp)
    anni, toa_net, srf_net, tas = tl.read_gregory(gregfi)

    toa5 = []
    tas5 = []
    for i in range(0, len(anni), 5):
        if len(anni) - i < 2: continue
        toa5.append(np.mean(toa_net[i:i+5]))
        tas5.append(np.mean(tas[i:i+5]))

    #ax.scatter(tas5, toa5, color = col, marker = mar, label = exp)
    ax.scatter(tas5[1:-1], toa5[1:-1], color = col, marker = mar, label = exp)
    ax.scatter(tas5[0], toa5[0], color = col, marker = '>')
    ax.scatter(tas5[-1], toa5[-1], color = col, marker = '<')
    ax.plot(tas5, toa5, color = col, linewidth = 0.5)

gregc4co = '/home/fabiano/Research/lavori/TunECS/gregory_c4co.txt'
anni, toa_net, srf_net, tas = tl.read_gregory(gregc4co)
col = 'black'

toa5 = []
tas5 = []
for i in range(0, len(anni), 5):
    if len(anni) - i < 2: continue
    toa5.append(np.mean(toa_net[i:i+5]))
    tas5.append(np.mean(tas[i:i+5]))

#ax.scatter(tas5, toa5, color = col, marker = mar, label = exp)
ax.scatter(tas5[1:-1], toa5[1:-1], color = col, marker = mar, label = 'c4co (cmip6 r4)')
ax.scatter(tas5[0], toa5[0], color = col, marker = '>')
ax.scatter(tas5[-1], toa5[-1], color = col, marker = '<')
ax.plot(tas5, toa5, color = col, linewidth = 0.5)

gregb100 = '/home/fabiano/Research/lavori/TunECS/gregory_b100.txt'
anni, toa_net, srf_net, tas = tl.read_gregory(gregb100)

col = 'violet'

toa5 = []
tas5 = []
for i in range(0, len(anni), 5):
    if len(anni) - i < 2: continue
    toa5.append(np.mean(toa_net[i:i+5]))
    tas5.append(np.mean(tas[i:i+5]))

ax.scatter(tas5[1:-1], toa5[1:-1], color = col, marker = mar, label = 'b100 (bottino)')
ax.scatter(tas5[0], toa5[0], color = col, marker = '>')
ax.scatter(tas5[-1], toa5[-1], color = col, marker = '<')
ax.plot(tas5, toa5, color = col, linewidth = 0.5)

ax.grid()
ax.legend()
fig.savefig(cart_out + 'gregory_5y_c4.pdf')
