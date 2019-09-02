#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os

from matplotlib import pyplot as plt

#import climtools_lib as ctl

from scipy import io

#############################################################

def trova_spip(ifile, hasha = '#', read_past = False):
    """
    Trova il '#' nei file .dat
    """
    gigi = 'a'
    while gigi != hasha :
        linea = ifile.readline()
        gigi = linea[0]
    else:
        if read_past:
            return linea[1:]
        else:
            return

def read_input_atm_man(filename):
    """
    Reads input atmosphere in manuel standard.
    :param filename:
    :return:
    """
    infile = open(filename,'r')
    trova_spip(infile,hasha='$')
    n_alt = int(infile.readline())
    trova_spip(infile,hasha='$')
    prof = []
    while len(prof) < n_alt:
        line = infile.readline()
        prof += list(map(float, line.split()))
    alts = np.array(prof)

    trova_spip(infile,hasha='$')
    prof = []
    while len(prof) < n_alt:
        line = infile.readline()
        prof += list(map(float, line.split()))
    pres = np.array(prof)

    trova_spip(infile,hasha='$')
    prof = []
    while len(prof) < n_alt:
        line = infile.readline()
        prof += list(map(float, line.split()))
    temp = np.array(prof)

    return alts, temp, pres


def read_input_vmr_man(filename):
    """
    Reads input atmosphere in manuel standard.
    :param filename:
    :return:
    """
    infile = open(filename,'r')
    trova_spip(infile,hasha='$')
    n_alt = int(infile.readline())
    trova_spip(infile,hasha='$')
    prof = []
    while len(prof) < n_alt:
        line = infile.readline()
        prof += list(map(float, line.split()))
    alts = np.array(prof)

    trova_spip(infile,hasha='$')
    n_mols = int(infile.readline().strip())
    mol_vmrs = dict()

    for i in range(n_mols):
        molnam = trova_spip(infile, hasha='$', read_past = True).strip()
        prof = []
        while len(prof) < n_alt:
            line = infile.readline()
            prof += list(map(float, line.split()))
        mol_vmrs[molnam] = np.array(prof)

    return alts, mol_vmrs

###################################################################

cartsav = '/home/fedefab/Scrivania/Research/Post-doc/CO2_cooling/new_param/sent2/sav/'
cartatm = '/home/fedefab/Scrivania/Research/Post-doc/CO2_cooling/new_param/sent2/atm/'

filsav = 'data_cira_mle_co2_1.sav'
filvmr = 'vmr_cira_mle_co2_1.prf'
filatm = 'pt_cira_mle.prf'

alts, temp, pres = read_input_atm_man(cartatm+filatm)
alts, mol_vmr = read_input_vmr_man(cartatm+filvmr)

coso = io.readsav(cartsav+filsav)['data']
nomi = coso.dtype.names
#('WNUMLOW', 'WNUMHIGH','CM', 'CSURF', 'LSPACE', 'HR_KDAY', 'SRCFCTN', 'RSURF', 'LOG_PRESS', 'TEMPERATURE', 'PRESSURE', 'ALTITUDE', 'CO2_VMR', 'O_VMR')

cm = coso.CM[0]
hr = coso.HR_KDAY[0]

co2 = coso.CO2_VMR[0]
alts = coso.ALTITUDE[0]
temp = coso.TEMPERATURE[0]
pres = coso.PRESSURE[0]

frlo = coso.WNUMLOW[0]
frhi = coso.WNUMHIGH[0]

# quindi abbiamo 66 alts, 26 freq intervals
# ora bisogna capire come usare ste cose
# ...
