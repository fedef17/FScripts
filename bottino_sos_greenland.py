# coding: utf-8
get_ipython().run_line_magic('run', 'imports.py')
cart_in = '/nas/BOTTINO/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/stabilization-ssp585-2100/r1i1p1f1/Omon/sos/'
allfi = glob.glob(cart_in + 'sos_*012.nc')
allfi.sort()
sos = xr.open_mfdataset(allfi)
sos = xr.open_mfdataset(allfi, use_cftime = True)
sos
sos = sos.drop(['vertices_longitude', 'vertices_latitude'])
zuki = ctl.regrid_dataset(sos, regrid_to_deg=1.)
plt.ion()
zuki
pinzu = zuki.sel(lon = slice(-80, 0), lat = slice(45,90))
pinzu
pinzu = zuki.sel(lon = slice(280, 360), lat = slice(45,90))
pinzu
ctl.plot_map_contour(zuki[0], plot_margins=(-80, 0, 45, 90))
ctl.plot_map_contour(zuki['sos'][0], plot_margins=(-80, 0, 45, 90))
ctl.plot_map_contour(zuki['sos'][0], plot_margins=(-80, 0, 45, 90), cbar_range=(20,40))
ctl.plot_map_contour(zuki['sos'][10], plot_margins=(-80, 0, 45, 90), cbar_range=(20,40))
ctl.plot_map_contour(zuki['sos'][10]-zuki['sos'][0], plot_margins=(-80, 0, 45, 90), cbar_range=(20,40))
ctl.plot_map_contour(zuki['sos'][10]-zuki['sos'][0], plot_margins=(-80, 0, 45, 90))
ctl.plot_map_contour(zuki['sos'][5]-zuki['sos'][0], plot_margins=(-80, 0, 45, 90))
