#!/bin/bash

# script to remove archived data from scratch

set -eu

expname=$1
firstyear=$2
year2=$3
dryrun=${4:-1}

DIR=$SCRATCH

# First year of simulation (the first year output needs not be removed)
#firstyear=2065

OUTDIR=${SCRATCH}/ece3/${expname}/output/

mincmorsize=13

# detect automatically the first year to be removed
year1=$( ls $OUTDIR | head -2 | tail -1 | cut -f2 -d"_" )
yearlast=$( ls $OUTDIR | tail -1 | cut -f2 -d"_" )

if [[ $year2 -ge $(( yearlast -1 )) ]] ; then
  year2=$(( yearlast-2 ))
  echo Leaving there the last two years, setting year2 = $year2
fi

# origin directory: DIR+ORIGBASE+..  ------  dest. directory: ec:/ccff/+ECBASE+..
ORIGBASE=ece3
infocmordir=${PERM}/ecearth3/infodir/cmorized/${expname}

#use ecp command with -u (update) flag
cmorcheck=0
cmorcheck2=0
nccheck=0
nccheck2=0

# loop on the year, disabled by submitter
for year in $(seq $year1 $year2) ; do
    ORIGTREE=${OUTDIR}/Output_${year}

    cmorcheck=0
    nccheck=0
    nccheck2=0

    # Check cmor prepare is there
    PREPARE=${infocmordir}/PrePARE_${expname}_${year}.txt
    if [ -f $PREPARE ]; then
       echo 'OK! CMOR successfully completed!'
       cmorcheck=1
    else
       echo CMOR has not been completed for this year!
    fi


    cmordim=$(du -sh $SCRATCH/ece3/${expname}/cmorized/$(ls $SCRATCH/ece3/${expname}/cmorized/ | head -n 1) | cut -f 1 | cut -d 'G' -f 1)

    if (( $cmordim >= $mincmorsize )); then
	    echo 'OK! CMOR folder larger than '${mincmorsize}'G'
	    cmorcheck2=1
    else
	    echo CMOR folder too small!!
    fi

    pino=$(grep 'skipped' ${infocmordir}/EC-Earth_nctcck_${expname}.txt | cut -d ':' -f 2)
    pino2=$(grep 'overlap' ${infocmordir}/EC-Earth_nctcck_${expname}.txt | cut -d ':' -f 2)
    pino3=$(grep 'broken' ${infocmordir}/EC-Earth_nctcck_${expname}.txt | cut -d ':' -f 2)

    if [[ pino -eq 0 && pino2 -eq 0 && pino3 -eq 0 ]]; then
	    nccheck=1
	    echo 'nccheck ok!'
    else
	    echo 'nccheck not ok....'
    fi

#    if [[ ( "$g" -eq 1 && "$c" = "123" ) || ( "$g" -eq 2 && "$c" = "456" ) ]]
#        then echo "g = $g; c = $c; true"
#    else echo "g = $g; c = $c; false"
#    fi
    if [[ nccheck -eq 1 && cmorcheck2 -eq 1 && cmorcheck -eq 1 && $year -ne $firstyear ]] ; then
       # OK, cancella output
       if [ $dryrun -eq 1 ]; then
           if [[ -d ${ORIGTREE} ]] ; then
             echo DRYRUN Removing output directory from scratch: ${ORIGTREE}
           else
             echo Output $year of exp $expname has already been removed!
           fi
       else
           echo Removing output directory from scratch: ${ORIGTREE}
           if [[ -d ${ORIGTREE} ]] ; then
             rm -rf ${ORIGTREE}
           else
             echo Output $year of exp $expname has already been removed!
           fi
       fi
    elif [[ archivecheck -eq 1 && dataontape -eq 1 && cmorcheck -eq 1 && $year -eq $firstyear ]] ; then
       echo Everything ok, but this is the first year, leaving it there.
    else
       # continue, leave everything there
       echo NOOOO '->' $expname $year
       echo NOOOO '->' $expname $year cmor: $cmorcheck archive: $archivecheck dataontape: $dataontape >> problems_${expname}.txt
    fi
    echo
done
