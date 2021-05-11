#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
import glob

import climtools_lib as ctl
import climdiags as cd

################################################################################

bcart = '/nas/BOTTINO/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/'

okvars = ['zg', 'ua']
oklevs = dict()
oklevs['zg'] = 50000.
oklevs['ua'] = 85000.

taglev = dict()
taglev['zg'] = '_500hPa'
taglev['ua'] = '_850hPa'

allexps = os.listdir(bcart)

for exp in allexps:
    expcart = bcart + exp + '/'
    allmems = os.listdir(expcart)
    print(exp)

    for mem in allmems:
        memcart = expcart + mem + '/'
        print(mem)

        daycart = memcart + 'day/'
        r25cart = memcart + 'day_r25/'
        ctl.mkdir(r25cart)

        for var in okvars:
            print(var)
            cart_in = daycart + var + '/'
            cart_out = r25cart + var + '/'
            cd.preprocess_cdo(cart_in, cart_out, sel_levels = oklevs[var], taglev = taglev[var])
