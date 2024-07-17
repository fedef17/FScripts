#!/bin/bash

# script to remove archived data from scratch

set -e

expname=$1
firstyear=$2
year2=$3
dryrun=${4:-1}

DIR=${SCRATCH}/

# First year of simulation (the first year output needs not be removed)
#firstyear=1850

INFODIR=$PERM/ecearth3/infodir/archive/${expname}

OUTDIR=${SCRATCH}/ece3/${expname}/output/

# detect automatically the first year to be removed
year1=$( ls $OUTDIR | head -2 | tail -1 | cut -f2 -d"_" )
yearlast=$( ls $OUTDIR | tail -1 | cut -f2 -d"_" )

#if [[ $year2 -ge $(( yearlast -1 )) ]] ; then
#  year2=$(( yearlast-2 ))
#  echo Leaving there the last two years, setting year2 = $year2
#fi

# origin directory: DIR+ORIGBASE+..  ------  dest. directory: ec:/ccff/+ECBASE+..
ORIGBASE=ece3
ECBASE=ece3/hosing

mincmorsize=21
nemonum=0
#use ecp command with -u (update) flag
cmorcheck=0
cmorcheck2=0
cmorcheck3=0
dataontape=0
archivecheck=0

# loop on the year, disabled by submitter
for year in $(seq $year1 $year2) ; do
    cmorcheck=0
    cmorcheck2=0
    cmorcheck3=0
    dataontape=0
    archivecheck=0

    echo "Checking archived data for exp $expname year $year..."
    TREEDIR=${ECBASE}/${expname}/output/Output_${year}
    ORIGTREE=${ORIGBASE}/${expname}/output/Output_${year}

    if ! [ -d ${DIR}${ORIGTREE} ] ; then
      echo Output $year of exp $expname has already been removed!
      continue
    fi

    # output: archive file by file for IFS and NEMO, a tar file for LPJG
    LOGFILE1=$INFODIR/archive_output_${expname}_${year}.txt
    if [ -f $LOGFILE1 ] ; then
       echo OK! archive logfile found: $LOGFILE1
       archivecheck=1
    else
       echo archive logfile NOT found!!: $LOGFILE1
    fi

    #check data on tape
    if [ $year -eq $firstyear ]; then ifsnum=26; else ifsnum=24; fi

    checkifs=$(els -l ec:/ccff/${TREEDIR}/IFS/ | wc -l)
    checknemo=$(els -l ec:/ccff/${TREEDIR}/NEMO/ | wc -l)
    if [ $checkifs -eq $ifsnum ] && [ $checknemo -eq $nemonum ]; then
       echo 'OK! Data are on tape'
       dataontape=1
    else
       echo Found only $checkifs ifs files and $checknemo nemo files on tape!!
    fi

    # Check cmor prepare is there
    PREPARE=${PERM}/ecearth3/infodir/cmorized/${expname}/PrePARE_${expname}_${year}.txt
    if [ -f $PREPARE ]; then
       echo 'OK! CMOR successfully completed!'
       cmorcheck=1
    else
       echo CMOR has not been completed for this year!
    fi

    cmordim=$(du -sh $SCRATCH/ece3/${expname}/cmorized/cmor_${year} | cut -f 1 | cut -d 'G' -f 1)

    if (( $cmordim >= $mincmorsize )); then
	    echo 'OK! CMOR folder larger than '${mincmorsize}'G'
	    cmorcheck2=1
    else
	    echo CMOR folder too small!!
    fi

    # Check cmor prepare is there also for the next year
    yearnext=$(( year+1 ))
    PREPARE=${PERM}/ecearth3/infodir/cmorized/${expname}/PrePARE_${expname}_${yearnext}.txt
    if [ -f $PREPARE ]; then
       echo 'OK! CMOR successfully completed also for next year!'
       cmorcheck3=1
    else
       echo CMOR has not been completed for next year!
    fi
#    if [[ ( "$g" -eq 1 && "$c" = "123" ) || ( "$g" -eq 2 && "$c" = "456" ) ]]
#        then echo "g = $g; c = $c; true"
#    else echo "g = $g; c = $c; false"
#    fi
#    if [[ archivecheck -eq 1 && dataontape -eq 1 && cmorcheck -eq 1 && cmorcheck2 -eq 1 && cmorcheck3 -eq 1 && $year -ne $firstyear ]] ; then
    if [[ archivecheck -eq 1 && dataontape -eq 1 && cmorcheck2 -eq 1 ]] ; then
       # OK, cancella output
       if [ $dryrun -eq 1 ]; then
           if [[ -d ${DIR}${ORIGTREE} ]] ; then
             echo DRYRUN Removing output directory from scratch: ${DIR}${ORIGTREE}
           else
             echo Output $year of exp $expname has already been removed!
           fi
       else
           echo Removing output directory from scratch: ${DIR}${ORIGTREE}
           if [[ -d ${DIR}${ORIGTREE} ]] ; then
             rm -rf ${DIR}${ORIGTREE}
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
