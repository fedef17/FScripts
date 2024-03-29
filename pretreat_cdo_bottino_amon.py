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

# okvars = ['zg', 'ua']
# oklevs = dict()
# oklevs['zg'] = 50000.
# oklevs['ua'] = 85000.
#
# taglev = dict()
# taglev['zg'] = '_500hPa'
# taglev['ua'] = '_850hPa'


filerr = open('/nas/BOTTINO/log_Amon_err.txt', 'w')

allexps = os.listdir(bcart)

#for exp in ['piControl']:
for exp in allexps:
    expcart = bcart + exp + '/'
    allmems = os.listdir(expcart)
    print(exp)

    for mem in allmems:
        memcart = expcart + mem + '/'
        print(mem)

        # daycart = memcart + 'Amon/'
        # r25cart = memcart + 'Amon_r25/'
        daycart = memcart + 'SImon/'
        r25cart = memcart + 'SImon_r1/'
        ctl.mkdir(r25cart)

        #for var in ['clt', 'rsut']:
        for var in os.listdir(daycart):
            print(var)
            cart_in = daycart + var + '/'
            cart_out = r25cart + var + '/'
            if os.path.exists(cart_out):
                if len(os.listdir(cart_in)) == len(os.listdir(cart_out)):
                    print('Already processed\n')
                    continue

            try:
                #cd.preprocess_cdo(cart_in, cart_out, skip_existing = True)
                cd.preprocess_cdo(cart_in, cart_out, skip_existing = True, grid = 'r360x180', gridtag = '1')
            except Exception as exc:
                filerr.write('Error for {} {} {}\n'.format(exp, mem, var))
                filerr.write(str(exc))
                filerr.write('\n')

filerr.close()
