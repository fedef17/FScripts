# %%
import xarray as xr
from cdo import Cdo
import numpy as np

# %%
from climtools import climtools_lib as ctl

# %%
cdo = Cdo(debug=True)

inicon = xr.open_dataset('/home/fabiano/Research/lavori/alluvione/input_data/pstocchi/IFS_Era5_2023050202.nc')

# %%
aka_era5 = inicon['hyam'].values
bika_era5 = inicon['hybm'].values

# %%
cart = '/home/fabiano/Research/lavori/alluvione/input_data/'

trends = xr.open_mfdataset(cart + '*amj*.nc')

# %%
# Akabika ERA5
aka = aka_era5
bika = bika_era5

# %%
from matplotlib import pyplot as plt
#%matplotlib inline

# %%
pres = aka[24:, np.newaxis, np.newaxis] + bika[24:, np.newaxis, np.newaxis] * trends['aps'].values

# %%
pres.shape

# %%
from scipy.interpolate import interp1d

# %%
T_ok = np.empty_like(pres)
R_ok = np.empty_like(pres)

# %%
#interp1d?
# Traceback (most recent call last):
#   File "/home/fabiano/Research/git/globothon/retrieve-ic/build_moloch_IC.py", line 52, in <module>
#     print(np.log(pres[:, lo, la])[0], np.log(pres[:, lo, la])[-1])
#                  ~~~~^^^^^^^^^^^
# IndexError: index 201 is out of bounds for axis 1 with size 201

# %%
for lo in np.arange(len(pres.lon)):
    print(lo)
    for la in np.arange(len(pres.lat)):
        print(np.log(pres[:, lo, la])[0], np.log(pres[:, lo, la])[-1])
        print(np.log(trends.plev))
        pio = interp1d(np.log(trends.plev), trends.T.values[:, lo, la], fill_value='extrapolate')
        T_ok[:, lo, la] = pio(np.log(pres[:, lo, la]))

        pio = interp1d(np.log(trends.plev), trends.R.values[:, lo, la], fill_value='extrapolate')
        R_ok[:, lo, la] = pio(np.log(pres[:, lo, la]))
    

# %%
trnew = xr.Dataset()
trnew['T_trend'] = xr.DataArray(T_ok, dims=['lev', 'lat', 'lon'], coords={'lat': trends.lat, 'lon': trends.lon, 'lev': np.arange(25, 138, dtype=float)}, name='T_trend')
trnew['R_trend'] = xr.DataArray(R_ok, dims=['lev', 'lat', 'lon'], coords={'lat': trends.lat, 'lon': trends.lon, 'lev': np.arange(25, 138, dtype=float)}, name='R_trend')

for var in ['hyam', 'hybm', 'hyai', 'hybi']:
    trnew[var] = inicon[var]

trnew.lev.attrs = inicon.lev.attrs


# %%
trnew.to_netcdf(cart + 'trend_ERA5_levels.nc')
