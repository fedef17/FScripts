#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
import matplotlib.pyplot as plt

#################################

nomefile = '/data-hobbes/fabiano/temp/ghg_rcp85.txt'

with open(nomefile, 'r') as infile: # open the file
    infile.readline() # skip the first line
    fieldnames = infile.readline().rstrip().split() # read column names
    lines = np.array([list(map(float, lin.rstrip().split())) for lin in infile.readlines()]) # read all the file and transform it into a matrix of floats

quant = dict()
for col, nam in zip(lines.T, fieldnames): # create a dictionary with fieldnames and arrays
    quant[nam] = col

fig = plt.figure()
plt.plot(quant['YEAR'], quant['CO2']) # plot
fig.savefig('my_figure.pdf')
