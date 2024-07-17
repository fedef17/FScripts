
# %%
#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
import glob
import climtools_lib as ctl
import xarray as xr

import pickle

cart_out = '/home/fabiano/Research/lavori/TunECS/results/feedbacks/'
ctl.mkdir(cart_out)

cart_in = '/data-hobbes/fabiano/TunECS/coupled/'
filin_pi = cart_in + 'pic{0}/cmorized/cmor_*/CMIP6/CMIP/EC-Earth-Consortium/EC-Earth3/piControl/r1i1p{0}f1/Amon/{1}/gr/v*/{1}*nc'
filin_4c = cart_in + 'c4c{0}/cmorized/cmor_*/CMIP6/CMIP/EC-Earth-Consortium/EC-Earth3/abrupt-4xCO2/r1i1p{0}f1/Amon/{1}/gr/v*/{1}*nc'


### IMPORTANT!!!

# ERROR IN depth of first layer (assuming constant depth, instead of considering varying surface pressure)

## TO FIX BEFORE publication

# %% [markdown]
# ### TODO
# - Per tutte le variabili:
#     - leggi pi e calcola pi mean (con seasonal cycle?)
#     - leggi c4 e fai anomaly con pi (con seasonal cycle)
# - Leggi kernels
# - Meglio leggere una var alla volta e moltiplicare per il kernel
# 
# Da Zelinka 2020: For each month of the 150-year experiment, spatially-resolved kernels are multiplied by the relevant climate field anomalies. These are then vertically integrated up to a time-varying tropopause and then annually averaged to produce a 150-year time series of TOA radiative flux anomalies due to each field. These are then regressed on T to yield the individual radiative feedback components.

# %%
#allvars = 'clt hfss pr rlds rlut rsds rsus rsut ta hfls hus prsn rlus rlutcs rsdscs rsuscs rsutcs tas'
allvars = 'clt rlut rsut ta hus rlutcs rsutcs tas rsus rsds'.split()

pimean = dict()
for exp in [5, 9]:
    for varnam in allvars:
        print(exp, varnam)
        filist = glob.glob(filin_pi.format(exp, varnam))
        filist.sort()

        var = xr.open_mfdataset(filist)

        var_mean = var.groupby('time.month').mean('time')

        # ctl.regrid_dataset(var_mean, kernel.lat, kernel.lon)

        pimean[(exp, varnam)] = var_mean[varnam]

for exp in [5, 9]:
    pimean[(exp, 'alb')] = pimean[(exp, 'rsus')]/pimean[(exp, 'rsds')]


# %% [markdown]
# ## Reading kernels
# 
# from Huang readme:
# 
# kernels calculated based on 5 yrs ERA-Interim reanalysis and RRTMG
# 
# the units of the atmospgeric temperature and water vapor kernels are w/m2/k/100hpa
# dR(Ta)=sum(Kt(i)*dTa(i)*dP(i)/100), i denotes each layer, dP(i) is the thickness of each layer
# dR(q) =sum(Kq(i)*dlnq(i)*Ta(i)^2*Rv/Lv*dP(i)/100)
# 
# in the file dp.nc
# player denotes the middle pressure of each layer we perturbed
# plevel denotes the boundary pressure of each layer
# dp denotes the thickness (in hPa) of each layer
# 
# 4. calculate the feedbacks
# - $dr_q =  \sum{Kq*dqn*dps/100.}$
# - $dr_{Ta} =  \sum{Kt*dTa*dps/100.}$
# - $dr_{Ts} =  Kts*dTs$
# - $dr_{alb} = Kalb*dalb$

# %%
cart_k = '/data-hobbes/fabiano/radiative_kernels/Huang/toa/'

finam = 'RRTMG_{}_toa_{}_highR.nc'

vnams = ['t', 'ts', 'wv_lw', 'wv_sw', 'alb']
tips = ['clr', 'cld']

allkers = dict()

for tip in tips:
    for vna in vnams:
        ker = xr.load_dataset(cart_k + finam.format(vna, tip))

        allkers[(tip, vna)] = ker

vlevs = xr.load_dataset(cart_k + 'dp.nc')

kernel = allkers[('cld', 't')].lwkernel
# ## Ora leggo i 4xCO2

# %%
vars_fb = ['tas', 'ta', 'hus', 'alb']
feedbacks = dict()

calc_fb = False

if calc_fb:
    for exp in [5, 9]:
        for varnam in vars_fb:
            print(exp, varnam)

            if varnam != 'alb':
                filist = glob.glob(filin_4c.format(exp, varnam))
                filist.sort()

                var = xr.open_mfdataset(filist)

                var = var[varnam]
            else:
                filist_1 = glob.glob(filin_4c.format(exp, 'rsus'))
                filist_1.sort()

                var_rsus = xr.open_mfdataset(filist_1)['rsus']

                filist_2 = glob.glob(filin_4c.format(exp, 'rsds'))
                filist_2.sort()

                var_rsds = xr.open_mfdataset(filist_2)['rsds']

                var = var_rsus/var_rsds

            ## regridding to kernel dimension
            var = ctl.regrid_dataset(var, kernel.lat, kernel.lon)
            pivar = ctl.regrid_dataset(pimean[(exp, varnam)], kernel.lat, kernel.lon)

            ## need to compute before apply_ufunc
            piok = pivar.compute()
            var = var.compute()

            # Removing inf and nan from alb
            if varnam == 'alb': 
                piok = piok.where(piok > 0., 0.)
                var = var.where(var > 0., 0.)

            ## store ta for use with hus
            if varnam == 'ta':
                ta_abs = var.interp(plev = 100*vlevs.player) 

            ## Computing log of hus
            if varnam == 'hus':
                var = np.log(var)
                piok = np.log(piok)

            anoms = xr.apply_ufunc(lambda x, mean: x - mean, var.groupby('time.month'), piok)

            # match varnam:
            #     case 'ta':
            #         pass
            #     case 'tas':
            #         pass    

            if varnam in ['ta', 'hus']:
                anoms_ok = anoms.interp(plev = 100*vlevs.player) 
            else:
                anoms_ok = anoms
            
            # Calculating the feedbacks
            for tip in ['cld', 'clr']:
                if varnam == 'ta':
                    kernel = allkers[(tip, 't')].lwkernel

                    # Here I was removing the vertical mean warming, but better removing the tas anomaly
                    anoms_mea = anoms_ok.mean('player')
                    anoms_lr = anoms_ok - anoms_mea
                    
                    # Uniform warming (equal to surface) is part of the Planck feedback
                    anoms_lr = anoms_ok - tas_anom
                    anoms_unif = anoms_ok - anoms_lr

                    dRt_unif = (xr.apply_ufunc(lambda x, ker: x*ker, anoms_unif.groupby('time.month'), kernel) * vlevs.dp / 100.).sum('player').groupby('time.year').mean('time')

                    dRt_lr = (xr.apply_ufunc(lambda x, ker: x*ker, anoms_lr.groupby('time.month'), kernel) * vlevs.dp / 100.).sum('player').groupby('time.year').mean('time')

                    dRt_unif_glob = ctl.global_mean(dRt_unif)
                    dRt_lr_glob = ctl.global_mean(dRt_lr)

                    feedbacks[(exp, tip, 'planck-atmo')] = dRt_unif_glob
                    feedbacks[(exp, tip, 'lapse-rate')] = dRt_lr_glob

                    sys.exit()

                    del dRt_unif, dRt_lr

                elif varnam == 'tas':
                    kernel = allkers[(tip, 'ts')].lwkernel

                    gtas = ctl.global_mean(anoms_ok).groupby('time.year').mean('time')
                    tas_anom = anoms_ok

                    dRt = xr.apply_ufunc(lambda x, ker: x*ker, anoms_ok.groupby('time.month'), kernel).groupby('time.year').mean('time')
                    dRt_glob = ctl.global_mean(dRt)

                    feedbacks[(exp, tip, 'planck-surf')] = dRt_glob
                    
                    del dRt

                elif varnam == 'hus':
                    # dR(q) =sum(Kq(i)*dlnq(i)*Ta(i)^2*Rv/Lv*dP(i)/100)
                    Rv = 487.5 # gas constant of water vapor
                    Lv = 2.5e+06 # latent heat of water vapor

                    kernel_lw = allkers[(tip, 'wv_lw')].lwkernel
                    kernel_sw = allkers[(tip, 'wv_sw')].swkernel
                    kernel = kernel_lw + kernel_sw

                    coso = anoms_ok * ta_abs**2
                    coso = coso.compute()

                    #dRt = (xr.apply_ufunc(lambda x, ta, ker: x*ta**2*ker, anoms_ok.groupby('time.month'), ta_abs.groupby('time.month'), kernel) * Rv/Lv * vlevs.dp / 100.).sum('player').groupby('time.year').mean('time')
                    dRt = (xr.apply_ufunc(lambda x, ker: x*ker, coso.groupby('time.month'), kernel) * Rv/Lv * vlevs.dp / 100.).sum('player').groupby('time.year').mean('time')

                    dRt_glob = ctl.global_mean(dRt)

                    feedbacks[(exp, tip, 'water-vapor')] = dRt_glob
                    
                    del dRt
                
                elif varnam == 'alb':
                    kernel = allkers[(tip, 'alb')].swkernel

                    dRt = xr.apply_ufunc(lambda x, ker: x*ker, anoms_ok.groupby('time.month'), kernel).groupby('time.year').mean('time')
                    dRt_glob = ctl.global_mean(dRt)

                    feedbacks[(exp, tip, 'albedo')] = 100*dRt_glob
                    
                    del dRt

        feedbacks[(exp, 'gtas')] = gtas

    pickle.dump(feedbacks, open(cart_out + 'feedbacks_tunecs.p', 'wb'))

else:
    feedbacks = pickle.load(open(cart_out + 'feedbacks_tunecs.p', 'rb'))

# for exp in [5, 9]:
#     for tip in ['clr', 'cld']:
#         feedbacks[(exp, tip, 'albedo')] = 100*feedbacks[(exp, tip, 'albedo')]

# %%
from scipy import stats

fbnams = ['planck-surf', 'planck-atmo', 'lapse-rate', 'water-vapor', 'albedo']

fb_coef = dict()

for tip in ['clr', 'cld']:
    print('\n -------------- {} ---------------- \n'.format(tip))
    for fbn in fbnams:
        for exp in [5, 9]:
            gtas = feedbacks[(exp, 'gtas')]
            coso = feedbacks[(exp, tip, fbn)]

            res = stats.linregress(gtas, coso)
            fb_coef[(exp, tip, fbn)] = res

            print('{} feedback, exp {} = {:6.2f} +/- {:6.2f} W/m^2/K'.format(fbn, exp, res.slope, res.stderr))


# %%
## Computing cloud feedback

for exp in [5, 9]:
    gtas = feedbacks[(exp, 'gtas')]

    varnam = 'rlut'
    filist = glob.glob(filin_4c.format(5, varnam))
    filist.sort()
    rlut = xr.open_mfdataset(filist)['rlut']

    varnam = 'rsut'
    filist = glob.glob(filin_4c.format(5, varnam))
    filist.sort()
    rsut = xr.open_mfdataset(filist)['rsut']

    varnam = 'rsutcs'
    filist = glob.glob(filin_4c.format(5, varnam))
    filist.sort()
    rsutcs = xr.open_mfdataset(filist)['rsutcs']

    varnam = 'rlutcs'
    filist = glob.glob(filin_4c.format(5, varnam))
    filist.sort()
    rlutcs = xr.open_mfdataset(filist)['rlutcs']

    ###

    N = - rlut - rsut
    N0 = - rsutcs - rlutcs

    crf = (N0 - N)
    crf = crf.groupby('time.year').mean('time')

    N = N.groupby('time.year').mean('time')
    N0 = N0.groupby('time.year').mean('time')

    crf_glob = ctl.global_mean(crf).compute()
    N_glob = ctl.global_mean(N).compute()
    N0_glob = ctl.global_mean(N0).compute()

    res_N = stats.linregress(gtas, N_glob)
    res_N0 = stats.linregress(gtas, N0_glob)
    res_crf = stats.linregress(gtas, crf_glob)

    for nam, res in zip(['N', 'N0', 'Crf'], [res_N, res_N0, res_crf]):
        print(r'{} - slope: ${:6.2f} \pm {:6.2f}$, intercept: ${:6.2f} \pm {:6.2f}$'.format(nam, res.slope, res.stderr, res.intercept, res.intercept_stderr))

    F0 = res_N0.intercept + pimean[(exp, 'rlutcs')] + pimean[(exp, 'rsutcs')]
    F = res_N.intercept + pimean[(exp, 'rlut')] + pimean[(exp, 'rsut')]
    F0.compute()
    F.compute()

    F_glob = ctl.global_mean(F).mean('month')
    F0_glob = ctl.global_mean(F0).mean('month')
    F_glob = F_glob.compute()
    F0_glob = F0_glob.compute()

    print(F0_glob-F_glob)

    # fb_cloud = res_crf.slope + np.nansum([fb_coef[(exp, 'clr', fbn)].slope - fb_coef[(exp, 'cld', fbn)].slope for fbn in fbnams]) + (F0_glob - F_glob)/gtas[-5:].mean() ## as in Soden

    fb_cloud = -res_crf.slope + np.nansum([fb_coef[(exp, 'clr', fbn)].slope - fb_coef[(exp, 'cld', fbn)].slope for fbn in fbnams])

    fb_cloud_err = np.sqrt(res_crf.stderr**2 + np.nansum([fb_coef[(exp, 'cld', fbn)].stderr**2 for fbn in fbnams]))

    fb_coef[(exp, 'cloud')] = fb_cloud
    fb_coef[(exp, 'cloud_err')] = fb_cloud_err

    print('cloud feedback, exp {} = {:6.2f} +/- {:6.2f} W/m^2/K'.format(exp, fb_cloud, fb_cloud_err))
# %%

for exp in [5, 9]:
    fb_coef[(exp, 'tot')] = np.sum([fb_coef[(exp, 'cld', fbn)].slope for fbn in fbnams]) + fb_coef[(exp, 'cloud')]

    fb_coef[(exp, 'tot_err')] = np.sqrt(np.sum([fb_coef[(exp, 'cld', fbn)].stderr**2 for fbn in fbnams]) + fb_coef[(exp, 'cloud_err')]**2)

    print('TOT feedback, exp {} = {:6.2f} +/- {:6.2f} W/m^2/K'.format(exp, fb_coef[(exp, 'tot')], fb_coef[(exp, 'tot_err')]))


# %%
pickle.dump(feedbacks, open(cart_out + 'feedbacks_tunecs.p', 'wb'))

pickle.dump(fb_coef, open(cart_out + 'fb_coef_tunecs.p', 'wb'))
# %%

### Appunti per improvements
# - usa ts invece che tas
# - separa sw e lw cloud feedback (devi separare anche water vapor)
# - guarda la distribuzione spaziale dei feedbacks
# - in generale, rivedi il cloud feedback che non mi convince tantissimo
