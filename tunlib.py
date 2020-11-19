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

def read_gregory(filnam):
    with open(filnam, 'r') as filoz:
        filoz.readline()
        linee = filoz.readlines()
        anni = np.array([int(lin.split(':')[0].split('(')[1][:4]) for lin in linee])

        cose = np.stack([lin.rstrip().split(':')[1].split() for lin in linee])

        toa_net = np.array([float(co) for co in cose[:, 0]])
        srf_net = np.array([float(co) for co in cose[:, 1]])
        tas = np.array([float(co) for co in cose[:, 2]])

        gigi = np.argsort(anni)
        anni = anni[gigi]
        toa_net = toa_net[gigi]
        srf_net = srf_net[gigi]
        tas = tas[gigi]

    return anni, toa_net, srf_net, tas

def check_file(filnam):
    if os.path.exists(filnam):
        return True
    else:
        return False
