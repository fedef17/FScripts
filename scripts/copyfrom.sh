#!/bin/bash

#expname=b025
#year1=2025
#year2=2099

expname=${expname}
year1=${year1}
year2=${year2}

echo $expname, $year1, $year2

username=ffabiano
ipaddress=login.galileo.cineca.it

TARGETDIR=/gpfs/scratch/userexternal/pdavini0/ece3/${expname}/output/
OUTDIR=/scratch/ms/it/ccff/cinbkp/ece3/${expname}/output/

#rsync -av --append-verify --progress -e "ssh" ffabiano@login.galileo.cineca.it:/gpfs/scratch/userexternal/pdavini0/ece3/b025/output/Output_20* ece3/b025/output/

copycommand="rsync -av --append-verify --progress"

# loop until succeed for three times
function smartcopy {
        ok_year=$1
        MAX_RETRIES=20; i=0; rcheck=255
        while ( [[ $rcheck -ne 0 ]] && [[ $i -lt $MAX_RETRIES ]] ) ; do
                i=$(($i+1))
		echo "$copycommand -e "ssh -T -o Compression=no -x" ${username}@${ipaddress}:${TARGETDIR}Output_${ok_year} ${OUTDIR}"
                $copycommand -e "ssh -T -o Compression=no -x" ${username}@${ipaddress}:${TARGETDIR}Output_${ok_year} ${OUTDIR}
                rcheck=$?
        done
}

for ye in $(seq $year1 $year2); do
  echo Copying $expname year $ye ...
  smartcopy $ye
done

exit 0
