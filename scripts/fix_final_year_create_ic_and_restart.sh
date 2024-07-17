#!/bin/bash

# this create is made to conclude a CMIP6 simulation:
# it aims at several things
# 1. it creates the initial condition from the last save_ic produced
#    this is done rebuilding the nemo files and rearranging various files in new folder

expname=pi9r
year=1870

modeltype="EC-EARTH-AOGCM" 

DIR=$SCRATCH/ece3/$expname
SAVEICDIR=$DIR/save_ic/${year}0101
ICDIR=/ec/res4/hpcperm/ccff/ICs/${modeltype}_${expname}_${year}
mkdir -p $ICDIR
cd $ICDIR

echo "NEMO"
# nemo rebuild
NEMOREBUILD=$HPCPERM/ecearth3/revisions/r7870/sources/nemo-3.6/TOOLS/REBUILD_NEMO/rebuild_nemo
for kind in ice oce ; do
	cp -v ${SAVEICDIR}/nemo/${expname}_????????_restart_${kind}_????.nc $ICDIR
	nfiles=$(ls $ICDIR/${expname}_????????_restart_${kind}_????.nc | wc -l )
	tstep=$(basename $ICDIR/${expname}_????????_restart_${kind}_0000.nc | cut -f2 -d"_")
	echo $nfiles $tstep
	
 	$NEMOREBUILD -t 2 $ICDIR/${expname}_${tstep}_restart_${kind} $nfiles
	ln -sf ${expname}_${tstep}_restart_${kind}.nc restart_${kind}.nc
	rm -f $ICDIR/${expname}_????????_restart_${kind}_????.nc
done


echo "IFS"
cp -v ${SAVEICDIR}/ifs/ICM{SH,GG}${expname}INI* $ICDIR

echo "OASIS"
cp -v ${SAVEICDIR}/oasis/rst{a,o}s.nc $ICDIR


if [[ $modeltype == "EC-EARTH-Veg" ]] ; then
    echo "LPJG"
    LPJGDIR=$DIR/restart/Restart_${year}0101/LPJG/lpjg_state_${year}
    cp -rv $LPJGDIR $ICDIR
    ln -sf lpjg_state_${year} lpjg_state_init
    OASISDIR=$DIR/restart/Restart_${year}0101/OASIS
    cp -v $OASISDIR/lpjgv.nc $ICDIR
    cp -v $OASISDIR/vegin.nc $ICDIR
fi

#tarfile=$SCRATCH/PRIMAVERA/CMIP6-INIT/${modeltype}_${expname}_${year}.tar.gz
#tar cfvz $tarfile -C $SCRATCH/PRIMAVERA/CMIP6-INIT ${modeltype}_${expname}_${year}
#ecp $tarfile  ec:/ccpd/ICs/T255L91
#echmod 644 ec:/ccpd/ICs/T255L91/${modeltype}_${expname}_${year}.tar.gz
#rm -f $tarfile
