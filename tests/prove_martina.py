# coding: utf-8
get_ipython().run_line_magic('clear', '')
cart = '/data-hobbes/fabiano/CMIP6/cmip6_rad/amip/'
fil = cart + '{}/{}/{}/{}/{}*nc'
fil_ok = fil.format('EC-Earth3', 'r1i1p1f1', 'Amon', 'rlut', 'rlut')
fil_ok
import glob
glob.glob(fil_ok)
fil1 = glob.glob(fil_ok)
fil1
import xarray as xr
import matplotlib.pyplot as plt
import climtools_lib as ctl
gigi = xr.load_dataset(fil1)
fil1
fil1 = fil1[0]
gigi = xr.load_dataset(fil1)
get_ipython().run_line_magic('clear', '')
gigi
gigi['rlut']
rlut = gigi['rlut']
rlut.sel(lat = slice(20, 40))
rlut.mean('time')
rlut.mean('time').plot()
plt.close('all')
plt.ion()
rlut.mean('time').plot()
ctl.plot_map_contour(rlut.mean('time'))
rlut_t = rlut.mean('time')
rlut_t.mean('lat', 'lon')
rlut_t.mean(['lat', 'lon'])
ctl.global_mean(rlut_t)
rlut_yr = rlut.groupby('time.year').mean('time')
rlut_yr
rlut_yr.std('year')
ctl.plot_map_contour(rlut_yr.std('year'))
ctl.plot_map_contour(rlut_yr.std('year'), cbar_range=(0,10), central_lat_lon=(0, 180), visualization='Robinson')
