### appunti per post-processing di FORUM

grib_copy ICMGGfrm6+185001 ICMGG_[surface].grb

cdo setctomiss,-10000 ICMGG_surface.grb outprova.grb   # questo mette a NaN i -10000

cdo -f nc setgridtype,regularnn outprova.grb outprova2.nc # questo converte a regular gaussian ma con nearest neighbour

cdo selparam,52.126 outprova2.nc outprova2_52.nc  # seleziona il parametro


# Nuova procedura
grib_copy ICMGGfor3+200701 ICMGG_[hybrid].grb

cdo selparam,51.126 ICMGG_undef.grb ICMGG_var51.grb # e 52, 53
cdo selparam,164.128 ICMGG_undef.grb ICMGG_clouds164.grb

cdo setctomiss,-10000 ICMGG_var51.grb ICMGG_var51_miss.grb

### CLEAR SKY
cdo lec,0.1 ICMGG_clouds164.grb ICMGG_clearsky.grb
cdo ifthen ICMGG_clearsky.grb ICMGG_var51_miss.grb ICMGG_var51_clearsky.grb

cdo timmean ICMGG_var51_clearsky.grb ICMGG_var51_clearsky_monme.grb #### ATTENZIONE! This only takes gridpoints which are always non missing..

cdo -f nc setgridtype,regularnn ICMGG_var51_clearsky_monme.grb ICMGG_var51_clearsky_monme.nc
#cdo zonmean ICMGG_var51_clearsky_monme.nc ICMGG_var51_clearsky_monme_zonme.nc

cdo mergetime ICMGG_var51_clearsky_monme*2007*nc ICMGG_var51_clearsky_2007.nc
