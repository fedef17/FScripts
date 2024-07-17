# %%
from climtools import climtools_lib as ctl
import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
#%matplotlib inline
import pickle

# %%
Index_MonthRM, freqs, power_spec, stdT_MonthRM, std100_MonthRM, stdRM_MonthRM = pickle.load(open('/home/montanarini/output/Task-5/Variables/Task-5_Dataset-hist.p', 'rb'))

# %%
plt.plot(Index_MonthRM)

# %%
len(Index_MonthRM)/12

# %%
cart_out = '/home/fabiano/Research/abstractsepapers/tesi/tesi_alemont/'

# %%
# # read thetao
thetao = xr.open_mfdataset('/nas/BOTTINO/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/historical/r4i1p1f1/Omon_r25/thetao/*nc')
thetaok = thetao.sel(lat=slice(-5,5), lon=slice(120,290))
thet = thetaok.compute()
thet.to_netcdf(cart_out + 'thet_hist.nc')

# %%
thet = xr.open_dataset(cart_out + 'thet_hist.nc')

# %%
thet_ok = thet.mean('lat')
thet_clim = thet_ok.groupby('time.month').mean('time')

thet_anom = thet_ok.groupby('time.month') - thet_clim

# %%
thet_anom

# %%
oknino = Index_MonthRM > 0.4

# %%
thet_clim

# %%
thet_cut = thet_anom.sel(time = Index_MonthRM.time)

# %%
thet_cut

# %%
nino_comp = xr.where(Index_MonthRM > 0.4, thet_cut, np.nan).mean('time')

# %%
fig = plt.figure()
thet_clim['thetao'].mean('month').plot.contourf(x = 'lon', y ='lev', levels = 17, ylim = (600., 0.))
plt.savefig(cart_out + 'thetao_clim_hist.pdf')

# %%
fig = plt.figure()
nino_comp['thetao'].plot.contourf(x = 'lon', y ='lev', levels = 17, ylim = (600., 0.))
plt.savefig(cart_out + 'nino_comp_hist.pdf')

# %%



