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
import xarray as xr
import glob
import xclim

plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
titlefont = 22
plt.rcParams['figure.titlesize'] = titlefont
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.axisbelow'] = True
plt.rcParams['legend.fontsize'] = 18

#############################################################################

# simple 2-layer ebm (geoffroy2013, held2010, gregory2000)
# constants
oce_mass = 1.381107e+21 # global and vertical sum of masscello*areacello
#ml_mass = 3.7352024e+19 # first 100 m globally
#bulk_mass = 1.3436264e+21 # below 100 m globally
ml_mass = 5.595727e+19 # first 150 m
bulk_mass = 7.689922e+20 # btw 150 and 2500 m
oce_area = 3.6481962e+14 # m2
cp0 = 3989.245 # J/kg/K

Cs = cp0*ml_mass/oce_area
Cd = cp0*bulk_mass/oce_area
# parameters
lamb = 1. # feedback param, W/m2/K
gamma0 = 1. # efficiency of heat transfer to deep ocean
epsilon = 0.6 # efficacy of heat removal from surface
forc = 7.4 # W/m2 (set to 4xCO2?)

dt = 86400*10
dt_year = 86400*365
n_year = 1000
n_step = int(n_year * dt_year/dt)

gammafu = np.append(np.linspace(1., 0.1, n_step//2), 0.1*np.ones(n_step//2))
# initial conditions
Ts_in = 0.
Td_in = 0.

def step(Ts, Td, gamma = gamma0):
    dTs = dt*(-lamb*Ts -epsilon*gamma*(Ts-Td) + forc)/Cs
    dTd = dt*(gamma*(Ts-Td))/Cd
    return dTs, dTd

plt.ion()

fig = plt.figure(figsize = (16,9))

# with constant gamma = 1
Ts_all = []
Td_all = []
years = []

dTs, dTd = step(Ts_in, Td_in)
Ts = Ts_in + dTs
Td = Td_in + dTd
for i in range(n_step):
    dTs, dTd = step(Ts, Td)
    Ts += dTs
    Td += dTd
    Ts_all.append(Ts)
    Td_all.append(Td)
    years.append(i*dt/dt_year)

plt.plot(years, Ts_all, color = 'blue', lw = 2, label = r'$\gamma = 1$')
plt.plot(years, Td_all, color = 'orange', lw = 2)

# with constant gamma = 0.1
gamma2=0.1
Ts_all = []
Td_all = []
years = []

dTs, dTd = step(Ts_in, Td_in)
Ts = Ts_in + dTs
Td = Td_in + dTd
for i in range(n_step):
    dTs, dTd = step(Ts, Td, gamma = gamma2)
    Ts += dTs
    Td += dTd
    Ts_all.append(Ts)
    Td_all.append(Td)
    years.append(i*dt/dt_year)

plt.plot(years, Ts_all, color = 'blue', lw = 2, ls = '--', label = r'$\gamma = 0.1$')
plt.plot(years, Td_all, color = 'orange', lw = 2, ls = '--')

### with varying gamma
Ts2_all = []
Td2_all = []
years = []

dTs, dTd = step(Ts_in, Td_in, gammafu[0])
Ts = Ts_in + dTs
Td = Td_in + dTd
#gamma = np.linspace(1., 0., 100000)
for i in range(n_step):
    dTs, dTd = step(Ts, Td, gammafu[i])
    Ts += dTs
    Td += dTd
    Ts2_all.append(Ts)
    Td2_all.append(Td)
    years.append(i*dt/dt_year)

plt.plot(years, Ts2_all, color = 'blue', lw = 2, ls = ':', label = r'$\gamma = 1 \to 0.1$')
plt.plot(years, Td2_all, color = 'orange', lw = 2, ls = ':')

plt.legend()

plt.xlabel('Years')
plt.ylabel(r'$\Delta T\,$(K)')
