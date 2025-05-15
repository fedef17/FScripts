#!/usr/bin/env python
import cdsapi
c = cdsapi.Client()

import os

"""
some examples of variables names:

u_component_of_wind
v_component_of_wind
geopotential

"""

update_db = False

for var in ['total_column_cloud_liquid_water', 'total_column_cloud_ice_water']:
	if not update_db:
		for year in range(1940,2023):
			year = str(year)
			c.retrieve('reanalysis-era5-single-levels-monthly-means',
    			{
			        'format': 'netcdf',
			        'product_type': 'monthly_averaged_reanalysis',
				    'variable': var,
			    	'year': year,
				    'month':[
				        '01','02','03',
			    	    '04','05','06',
		    	    	'07','08','09',
			    	    '10','11','12'
				    ],
				    'time':'00:00',
				 },
				 'ERA5_Amon_{}_{}.nc'.format(var, year))
	else:
		print('Updating db')
		year = 2023
		year = str(year)

		if os.path.exists('ERA5_{}_{}.nc'.format(var, year)):
			print('Renaming old file')
			os.rename('ERA5_Amon_{}_{}.nc'.format(var, year), 'ERA5_Amon_{}_{}_OLD.nc'.format(var, year))

		c.retrieve('reanalysis-era5-single-levels-monthly-means',
   			{
		        'format': 'netcdf',
		        'product_type': 'monthly_averaged_reanalysis',
			    'variable': var,
		    	'year': year,
			    'month':[
			        '01','02','03',
		    	    '04','05','06',
	    	    	'07','08','09',
		    	    '10','11','12'
			    ],
			    'time':'00:00',
			 },
			 'ERA5_Amon_{}_{}.nc'.format(var, year))

