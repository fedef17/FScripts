import os
import numpy as np
from scipy import signal
import netCDF4
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from pylab import *

def lowpass(x,fs = 1,fc = 1/10):
    fs = 1     # Sampling frequency of data
    fc = 1/10  # Cut-off frequency of the filter
    w = fc / (fs / 2) # Normalize the frequency
    bfil, afil = signal.butter(5, w, 'low')
    return signal.filtfilt(bfil, afil, x)

def stand(x):
    xs=(x-np.mean(x,axis=0))/np.std(x,axis=0)
    return xs

#################### CMIP6 - NAN ###############

EXPS='ssp585'
VARt='ta'
VARu='ua'
tmpDIR='/data-hobbes/fabiano/cmip6_UTW/'

INDIR_AA='/nas/archive_CMIP6/utw_aa_pvs/LTmean/'
INDIR_TA='/nas/archive_CMIP6/utw_aa_pvs/UTmean/'
INDIR_PVt='/nas/archive_CMIP6/utw_aa_pvs/Smean/'
INDIR_PVu='/nas/archive_CMIP6/utw_aa_pvs/Smean_u/'

MODELS=['ACCESS-CM2', 'BCC-CSM2-MR', 'CESM2-WACCM', 'CNRM-CM6-1-HR', 'CNRM-CM6-1', 'CNRM-ESM2-1', 'EC-Earth3', 'FGOALS-f3-L', 'FGOALS-g3', 'GISS-E2-1-G', 'HadGEM3-GC31-LL', 'HadGEM3-GC31-MM', 'INM-CM4-8', 'INM-CM5-0', 'IPSL-CM6A-LR', 'KACE-1-0-G', 'MIROC6',  'MPI-ESM1-2-HR', 'MPI-ESM1-2-LR', 'MRI-ESM2-0', 'NorESM2-LM', 'NorESM2-MM', 'UKESM1-0-LL']

ENSEMBLES=['r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f2', 'r1i1p1f2', 'r1i1p1f2', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f2', 'r1i1p1f3', 'r1i1p1f3', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1',  'r1i1p1f1', 'r1i1p1f2']

for mod, ens in zip(MODELS, ENSEMBLES):
    os.system('cdo -s -selyear,2015/2100 -fldmean -sellonlatbox,-180,180,60,90 -yearmean -selmonth,11,12,1,2,3 {}/ta_{}_{}_ssp585_2015-2100_LTmean_r1.nc {}AA_{}.nc'.format(INDIR_AA, mod, ens, tmpDIR, mod)
    os.system('cdo -s -selyear,2015/2100 -fldmean -sellonlatbox,-180,180,-30,30 -yearmean -selmonth,11,12,1,2,3 {}/ta_{}_{}_ssp585_2015-2100_UTmean_r1.nc {}TA_{}.nc'.format(INDIR_TA, mod, ens, tmpDIR, mod)
    os.system('cdo -s -selyear,2015/2100 -fldmean -sellonlatbox,-180,180,70,90 -yearmean -selmonth,11,12,1,2,3 {}/ta_{}_{}_ssp585_2015-2100_Smean_r1.nc {}PVt_{}.nc'.format(INDIR_PVt, mod, ens, tmpDIR, mod)
    os.system('cdo -s -selyear,2015/2100 -fldmean -sellonlatbox,-180,180,70,90  -yearmean -selmonth,11,12,1,2,3 {}/ta_{}_{}_ssp585_2015-2100_Smean_r1.nc {}PVu_{}.nc'.format(INDIR_PVu, mod, ens, tmpDIR, mod)

Nens=len(MODELS)
Yrs='2015/2100'

AA = {}
dAA = {}
TA = {}
dTA = {}
PVt = {}
dPVt = {}
PVu = {}
dPVu = {}
time={}

for mod in range(Nens):
    T = Dataset(tmpDIR+'/AA_'+MODELS[mod]+'.nc')
    A = np.squeeze(np.array(T.variables['ta'][:,0,0]))
    T.close()
    time[MODELS[mod]]=[datetime.datetime(2015,6,15)+relativedelta(years=1*val) for val in range(0,len(A))]
    AA[MODELS[mod]] = A
    os.system('cdo -s -trend '+tmpDIR+'/AA_'+MODELS[mod]+'.nc a.nc b.nc')
    T = Dataset('b.nc')
    A = np.squeeze(np.array(T.variables['ta'][0,0,0]))
    T.close()
    dAA[MODELS[mod]] = A
    os.system('rm a.nc b.nc')

    T = Dataset(tmpDIR+'/TA_'+MODELS[mod]+'.nc')
    A = np.squeeze(np.array(T.variables['ta'][:,0,0]))
    T.close()
    TA[MODELS[mod]] = A
    os.system('cdo -s -trend '+tmpDIR+'/TA_'+MODELS[mod]+'.nc a.nc b.nc')
    T = Dataset('b.nc')
    A = np.squeeze(np.array(T.variables['ta'][0,0,0]))
    T.close()
    dTA[MODELS[mod]] = A
    os.system('rm a.nc b.nc')

    T = Dataset(tmpDIR+'/PVt_'+MODELS[mod]+'.nc')
    A = np.squeeze(np.array(T.variables['ta'][:,0,0]))
    T.close()
    PVt[MODELS[mod]] = A
    os.system('cdo -s -trend '+tmpDIR+'/PVt_'+MODELS[mod]+'.nc a.nc b.nc')
    T = Dataset('b.nc')
    A = np.squeeze(np.array(T.variables['ta'][0,0,0]))
    T.close()
    dPVt[MODELS[mod]] = A
    os.system('rm a.nc b.nc')

    T = Dataset(tmpDIR+'/PVu_'+MODELS[mod]+'.nc')
    A = np.squeeze(np.array(T.variables['ua'][:,0,0]))
    T.close()
    PVu[MODELS[mod]] = A
    os.system('cdo -s -trend '+tmpDIR+'/PVu_'+MODELS[mod]+'.nc a.nc b.nc')
    T = Dataset('b.nc')
    A = np.squeeze(np.array(T.variables['ua'][0,0,0]))
    T.close()
    dPVu[MODELS[mod]] = A
    os.system('rm a.nc b.nc')

######################### FIGURE ######################

colors = ['slategrey','royalblue','navy','mediumblue','slateblue', 'darkslateblue', 'mediumslateblue', 'rebeccapurple', 'blueviolet', 'indigo', 'darkorchid', 'mediumorchid', 'thistle','plum', 'violet', 'purple','darkmagenta','fuchsia','orchid','mediumvioletred','deeppink', 'hotpink','palevioletred','crimson','pink']
#colors=['forestgreen','limegreen','darkgreen','green','g', 'lime','seagreen','mediumseagreen','springgreen','acquamarine','turquoise','lightseagreen','mediumturquoise', 'paleturquoise','darkslategray','darkcyan','cyan','acqua','powerblue','lightblue','deepskyblue', 'skyblue','lightskyblue','steelblue','dogerblue','slategrey', 'royalblue','navy','blue']

fig, axs = plt.subplots(2,2,figsize=(16, 10))
plt.subplots_adjust(hspace=0.25, wspace=0.25)
ax = fig.add_subplot(axs[0, 0])
for mod in range(Nens):
    ax.plot(time[MODELS[mod]],AA[MODELS[mod]], color=colors[mod], marker='', linestyle='-', linewidth=1.5, markersize=1, label=MODELS[mod])
#    ax.set_ylim(13.5, 21.5)
#    ax.yaxis.set_major_locator(plt.MaxNLocator(9))
#    ax.yaxis.set_minor_locator(plt.MaxNLocator(17))
    ax.set_ylabel('T (K)',fontsize=14)
    ax.set_title('Arctic: 60N-90N',fontsize=16)

ax = fig.add_subplot(axs[0, 1])
for mod in range(Nens):
    ax.plot(time[MODELS[mod]],TA[MODELS[mod]], color=colors[mod], marker='', linestyle='-', linewidth=1.5, markersize=1, label=MODELS[mod])
#    ax.set_ylim(-20, 5)
#    ax.yaxis.set_major_locator(plt.MaxNLocator(6))
#    ax.yaxis.set_minor_locator(plt.MaxNLocator(25))
    ax.set_title('Tropics: 20S-20N',fontsize=16)
    ax.set_ylabel('T (K)',fontsize=14)
    ax.legend(bbox_to_anchor=(1.1,1.0),ncol=1,loc='upper left',frameon=False,fontsize=14)

ax = fig.add_subplot(axs[1, 0])
for mod in range(Nens):
    ax.plot(time[MODELS[mod]],PVt[MODELS[mod]], colors[mod], marker='', linestyle='-', linewidth=1.5, markersize=1, label=MODELS[mod])
#    ax.set_ylim(8, 18)
#    ax.yaxis.set_major_locator(plt.MaxNLocator(12))
#    ax.yaxis.set_minor_locator(plt.MaxNLocator(23))
    ax.set_ylabel('T (K)',fontsize=14)
    ax.set_title('Polar Vortex (ta): 70N-90N',fontsize=16)

ax = fig.add_subplot(axs[1, 1])
for mod in range(Nens):
    ax.plot(time[MODELS[mod]],PVu[MODELS[mod]], colors[mod], marker='', linestyle='-', linewidth=1.5, markersize=1, label=MODELS[mod])
#    ax.set_ylim(12, 20)
#    ax.yaxis.set_major_locator(plt.MaxNLocator(6))
#    ax.yaxis.set_minor_locator(plt.MaxNLocator(25))
    ax.set_ylabel('U (m s-1)',fontsize=14)
    ax.set_title('Polar Vortex (ua): 70N-90N',fontsize=16)

for ax in axs.flat:
    ax.set_xlim([datetime.date(2016, 6, 15), datetime.date(2100, 6, 15)])
    ax.xaxis.set_major_locator(YearLocator(20, month=6, day=1))
    ax.xaxis.set_minor_locator(YearLocator(5, month=6, day=1))
    years_fmt = mdates.DateFormatter('%Y')
    ax.xaxis.set_major_formatter(years_fmt)
    ax.tick_params(labelsize=12)

plt.subplots_adjust(right=0.75)

plt.suptitle('Indices',fontsize=24)
plt.savefig('timeseries_Indices_'+EXPS+'.png',dpi=200)

######################### FIGURE ######################

fig, axs = plt.subplots(2,2,figsize=(16, 10))
plt.subplots_adjust(hspace=0.25, wspace=0.25)

ax = fig.add_subplot(axs[0, 0])
for mod in range(Nens):
    ax.plot(time[MODELS[mod]],lowpass(AA[MODELS[mod]]), color=colors[mod], marker='', linestyle='-', linewidth=1.5, markersize=1, label=MODELS[mod])
#    ax.set_ylim(13.5, 21.5)
#    ax.yaxis.set_major_locator(plt.MaxNLocator(9))
#    ax.yaxis.set_minor_locator(plt.MaxNLocator(17))
    ax.set_ylabel('T (K)',fontsize=14)
    ax.set_title('Arctic: 60N-90N',fontsize=16)

ax = fig.add_subplot(axs[0, 1])
for mod in range(Nens):
    ax.plot(time[MODELS[mod]],lowpass(TA[MODELS[mod]]), colors[mod], marker='', linestyle='-', linewidth=1.5, markersize=1, label=MODELS[mod])
#    ax.set_ylim(-20, 5)
#    ax.yaxis.set_major_locator(plt.MaxNLocator(6))
#    ax.yaxis.set_minor_locator(plt.MaxNLocator(25))
    ax.set_ylabel('T (K)',fontsize=14)
    ax.set_title('Tropics: 20S-20N',fontsize=16)
    ax.legend(bbox_to_anchor=(1.1,1.0),ncol=1,loc='upper left',frameon=False,fontsize=14)

ax = fig.add_subplot(axs[1, 0])
for mod in range(Nens):
    ax.plot(time[MODELS[mod]],lowpass(PVt[MODELS[mod]]), colors[mod], marker='', linestyle='-', linewidth=1.5, markersize=1, label=MODELS[mod])
#    ax.set_ylim(8, 18)
#    ax.yaxis.set_major_locator(plt.MaxNLocator(12))
#    ax.yaxis.set_minor_locator(plt.MaxNLocator(23))
    ax.set_ylabel('T (K)',fontsize=14)
    ax.set_title('Polar Vortex (ta): 70N-90N',fontsize=16)

ax = fig.add_subplot(axs[1, 1])
for mod in range(Nens):
    ax.plot(time[MODELS[mod]],lowpass(PVu[MODELS[mod]]), colors[mod], marker='', linestyle='-', linewidth=1.5, markersize=1, label=MODELS[mod])
#    ax.set_ylim(12, 20)
#    ax.yaxis.set_major_locator(plt.MaxNLocator(6))
#    ax.yaxis.set_minor_locator(plt.MaxNLocator(25))
    ax.set_ylabel('U (m s-1)',fontsize=14)
    ax.set_title('Polar Vortex (ua): 70N-90N',fontsize=16)

for ax in axs.flat:
    ax.set_xlim([datetime.date(2016, 6, 15), datetime.date(2100, 6, 15)])
    ax.xaxis.set_major_locator(YearLocator(20, month=6, day=1))
    ax.xaxis.set_minor_locator(YearLocator(5, month=6, day=1))
    years_fmt = mdates.DateFormatter('%Y')
    ax.xaxis.set_major_formatter(years_fmt)
    ax.tick_params(labelsize=12)

plt.subplots_adjust(right=0.75)

plt.suptitle('Indices - 10yrs lowpass',fontsize=24)
plt.savefig('timeseries_Indices_lowpass_'+EXPS+'.png',dpi=200)

###################### Writing a table ################

TName = [MODELS[mod] for mod in range(Nens)]
AAf = [dAA[MODELS[mod]] for mod in range(Nens)]
TAf = [dTA[MODELS[mod]] for mod in range(Nens)]
PVtf = [dPVt[MODELS[mod]] for mod in range(Nens)]
PVuf = [dPVu[MODELS[mod]] for mod in range(Nens)]

#table1 = '\n'.join('\t'.join(map(str,row)) for row in zip(TName,AAf,TAf,PVtf,PVvf))
table1 = '\n'.join('\t'.join(map(str,row)) for row in zip(TName,AAf,TAf,PVtf,PVuf))
print(table1)
