#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os

cart = '/data-hobbes/fabiano/Medscope/seasonal_forecasts_1d5/'

lista = [fi for fi in os.listdir(cart) if fi[-3:] == '.nc']

for fil in lista:
    fil = cart+fil
    nunam = fil[:-3]+'.grb'
    print(fil, nunam)
    command = 'mv '+fil+' '+nunam
    print(command)
    os.system(command)
    command = 'cdo -f nc -R copy '+nunam+' '+fil
    print(command)
    os.system(command)
