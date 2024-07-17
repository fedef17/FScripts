# %%
import xarray as xr
from cdo import Cdo
import numpy as np
from climtools import climtools_lib as ctl
import glob
import numpy as np
import pygrib
from scipy.interpolate import interp1d

# %% [markdown]
# - Ciao Fede, dovrei avere tutti i trend lineari che servono per fare un primo test con perturbazione in MOLOCH
# - /home/zappa/EM-moloch/data/nc
# - Perturberei queste quattro: ln_surface_pressure; skin_temperature; soil_temperature_level_1; temperature ... più l'aggiustamento sulla specific humidity assumendo relative humidity costante
# - (anche se la relative humidity non sembra costante, da trend lineare!)

# %%
# cart_init = '/home/fabiano/moloch_ini/ERA5/GRIB2/'
# cart_out = '/home/fabiano/moloch_ini/perturbed_v1/'

def fact(T1, T2):
    """
    Scaling factor for specific humidity.
    """
    alpha1 = 17.67*(T1-273.15)/(T1-29.65)
    alpha2 = 17.67*(T2-273.15)/(T2-29.65)
    return np.exp(alpha2-alpha1)

cart_init = '/data-hobbes/fabiano/alluvione/init_ERA5_GRIB2/GRIB2/'
cart_out = '/data-hobbes/fabiano/alluvione/init_ERA5_GRIB2/perturbed_v1_ok/'
ctl.mkdir(cart_out)

allfi_grib = glob.glob(cart_init + '*.grib2')
allfi_grib.sort()

allfi_nc = glob.glob(cart_init + '../orig_nc/*.nc')
allfi_nc.sort()

cart = '/data-hobbes/fabiano/alluvione/trend_v1/'
trends = xr.open_mfdataset(cart + '*amj*linTrend.nc')

### Interpolate T trend (consider surf P of first file)
inicon = xr.open_dataset(allfi_nc[0])

aka = inicon['hyam'].values
bika = inicon['hybm'].values

ps_clim = xr.load_dataset(cart + 'surface_pressure_amj_clim.nc')

pres = aka[24:, np.newaxis, np.newaxis] + bika[24:, np.newaxis, np.newaxis] * ps_clim['aps'].values #np.exp(inicon['lnsp'].values)
pres = pres.squeeze()

T_ok = np.empty_like(pres)
# Interpolate T to actual values
for lo in np.arange(len(trends.lon)):
    print(lo)
    for la in np.arange(len(trends.lat)):
        pio = interp1d(np.log(trends.plev), trends.T.values[:, la, lo], fill_value='extrapolate')
        T_ok[:, la, lo] = pio(np.log(pres[:, la, lo]))

trnew = xr.Dataset()
trnew['delta_T'] = xr.DataArray(T_ok, dims=['lev', 'lat', 'lon'], coords={'lat': trends.lat, 'lon': trends.lon, 'lev': np.arange(25, 138, dtype=float)}, name='delta_T')

trnew['delta_Q'] = xr.DataArray(fact(old_T, new_T))

for var in ['hyam', 'hybm', 'hyai', 'hybi']:
    trnew[var] = inicon[var]

trnew.lev.attrs = inicon.lev.attrs
trnew.to_netcdf(cart + 'trend_ERA5_levels.nc')


#for fil_nc, fil_grib in zip(allfi_nc, allfi_grib):
for fil_grib in allfi_grib:
    print(fil_grib)

    ### Uncomment here to consider surf P of each file differently
    # inicon = xr.open_dataset(fil_nc)

    # aka = inicon['hyam'].values
    # bika = inicon['hybm'].values

    # pres = aka[24:, np.newaxis, np.newaxis] + bika[24:, np.newaxis, np.newaxis] * np.exp(inicon['lnsp'].values)
    # pres = pres.squeeze()

    # # Interpolate T to actual values
    # for lo in np.arange(len(trends.lon)):
    #     print(lo)
    #     for la in np.arange(len(trends.lat)):
    #         pio = interp1d(np.log(trends.plev), trends.T.values[:, lo, la], fill_value='extrapolate')
    #         T_ok[:, lo, la] = pio(np.log(pres[:, lo, la]))

    # Create new inicon
    orig = pygrib.open(fil_grib)
    orig.rewind()

    allstr = []
    for grb in orig:
        #print(grb.name, grb.level)
        grb.expand_grid(False)

        if grb.name == 'Temperature':
            ### Sum trend
            old_T = grb.values
            new_T = old_T - T_ok[grb.level-25, ...].reshape(241*421)
            grb.values = new_T
        elif grb.name == 'Specific humidity':
            ### Adjust to new temp (read above)
            grb.values = grb.values * fact(old_T, new_T)
        elif grb.name == 'Soil temperature':
            if grb.level == 0:
                grb.values = grb.values + trends['STL1'].values.reshape(241*421)
        elif grb.name == 'Logarithm of surface pressure':
            grb.values = grb.values + trends['lnsp'].values.reshape(241*421)
        elif grb.name == 'Skin temperature':
            grb.values = grb.values + trends['SKT'].values.reshape(241*421)

        allstr.append(grb.tostring())

    grbout = open(cart_out + fil_grib.split('/')[-1], 'wb')
    for msg in allstr:
        grbout.write(msg)
    grbout.close()
