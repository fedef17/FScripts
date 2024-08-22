#!/bin/bash

set -ue

#. load_modules.sh

ECEDIR=/ec/res4/hpcperm/ccff/ecearth3/revisions/r9409-cmip6-fwf-knmi
. ${ECEDIR}/../trunk/load_modules.sh

#ECEDIR=/ec/res4/hpcperm/ccff/ecearth3/revisions/r7870
#ECEDIR=/ec/res4/hpcperm/ccff/ecearth3/revisions/r8083-waterhosing
cd $ECEDIR/sources

# run ec-conf
echo 'Linking build ec-conf'
#python2.7 ./util/ec-conf/ec-conf --platform ecmwf-tems-intel-openmpi config-build.xml
./util/ec-conf/ec-conf3 --platform ecmwf-hpc2020-intel-openmpi config-build.xml

# compile oasis
echo '============== OASIS ==================='
cd $ECEDIR/sources/oasis3-mct/util/make_dir
make realclean BUILD_ARCH=ecconf -f TopMakefileOasis3
make BUILD_ARCH=ecconf -f TopMakefileOasis3

# compile xios
echo ' ===================XIOS ==================='
cd $ECEDIR/sources/xios-2.5
./make_xios --arch ecconf --use_oasis oasis3_mct --netcdf_lib netcdf4_par --job 8

# compile ifs
echo ' ===================IFS ==================='
cd $ECEDIR/sources/ifs-36r4
make clean BUILD_ARCH=ecconf
make BUILD_ARCH=ecconf -j 4 lib
make BUILD_ARCH=ecconf master

# compile table 126
cd ${ECEDIR}/sources/util/grib_table_126
./define_table_126.sh

# compile nemo
echo ' ===================NEMO ==================='
cd $ECEDIR/sources/nemo-3.6/CONFIG
./makenemo -n ORCA1L75_LIM3 -m ecconf clean
./makenemo -n ORCA1L75_LIM3 -m ecconf -j 4
./makenemo -n ORCA1L75_LIM3_standalone -m ecconf -j 4

# compile Elpin (needed by NEMO)
cd ${ECEDIR}/sources/util/ELPiN
make clean
make

# compile runoff mapper and amip-forcing
echo ' ===================RUNOFF and AMIP ==================='
cd $ECEDIR/sources/runoff-mapper/src
make clean
make

cd $ECEDIR/sources/amip-forcing/src
make clean
make

echo 'COMPILATION COMPLETE!!'

########
echo 'Linking runtime ec-conf'

cd $ECEDIR/runtime/classic/
#python2.7 $ECEDIR/sources/util/ec-conf/ec-conf --platform ecmwf-tems-intel-openmpi config-run.xml
$ECEDIR/sources/util/ec-conf/ec-conf3 --platform ecmwf-hpc2020-intel-openmpi config-run.xml

echo 'READY TO RUN!'
