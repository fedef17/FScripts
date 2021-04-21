### appunti per post-processing di FORUM

grib_copy ICMGGfrm6+185001 ICMGG_[surface].grb

cdo setctomiss,-10000 ICMGG_surface.grb outprova.grb   # questo mette a NaN i -10000

cdo -f nc setgridtype,regularnn outprova.grb outprova2.nc # questo converte a regular gaussian ma con nearest neighbour

cdo selparam,52.126 outprova2.nc outprova2_52.nc  # seleziona il parametro
