#!/bin/sh

set -e

## Script to transfer bottino data to DRES_D_QUECLI via irods

echo $PATH
# Before running
module load icommands
module list

export PATH=$PATH:/cineca/prod/opt/tools/icommands/4.2.10/binary/

echo $PATH
ils

dryrun=0

cart_orig_1=/g100_scratch/userexternal/ffabiano/ece3/
cart_orig_2=/g100_work/IscrB_QUECLIM/BOTTINO/

#for mem in b990 b025 b050 b100 b80I b65I b00I b080 b065; do
for mem in b050 b100 b025 b80I b65I b00I b080 b065; do
    echo '------------------------------------------'
    echo '----------------> Trasferring '$mem
    echo '------------------------------------------'

    check_file='check_irods_trasferred_'${mem}'.dat'
    if [ -f ${check_file} ]; then
	lyear=$(tail -n 1 ${check_file})
    else
	lyear=0
    fi

    ysho=$(echo ${mem} | cut -d'b' -f 2)
    if [ ${mem} == 'b050' ] || [ ${mem} == 'b100' ]; then
	cartbase=${cart_orig_1}/${mem}/cmorized/
    else
	cartbase=${cart_orig_2}/${mem}/cmorized/
    fi

    if [ ${mem} == 'b990' ]; then
	mem_id='r1i1p1f1'
	run_id=stabilization-hist-1${ysho}
    elif [[ ${mem} == *'I' ]]; then
	mem_id='r1i1p1f3'
	run_id=stabilization-ssp585-2${ysho}
    else
	mem_id='r1i1p1f1'
	run_id=stabilization-ssp585-2${ysho}
    fi
    cartadd=CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/${run_id}/${mem_id}/

    cart_dest=QUECLIM/${cartadd}/

    for cos in $(ls -1 ${cartbase} | cat); do
	yea=${cos:(-4)}
	echo ${yea}
	cart_orig=${cartbase}/cmor_${yea}/${cartadd}/
	if (( $yea <= $lyear )); then
	    ils ${cart_dest}/Amon/tas/tas_Amon_EC-Earth3_${run_id}_${mem_id}_gr_${yea}01-${yea}12.nc
	    echo 'Year '${yea}' already transferred'
	    continue
	else
	    # Transfer to irods
	    for miptab in $(ls -1 ${cart_orig} | cat); do
		for var in $(ls -1 ${cart_orig}/${miptab} | cat); do
		    if (( $lyear == 0 )); then imkdir -p ${cart_dest}/${miptab}/${var}; fi
		    if (( $dryrun == 1 )); then
			ls ${cart_orig}/${miptab}/${var}/*/*/${var}*nc
			echo ${cart_dest}/${miptab}/${var}/
		    else
#			iput -P -b -T --retries 3 -X checkpoint-file-${mem} ${cart_orig}/${miptab}/${var}/*/*/${var}*nc ${cart_dest}/${miptab}/${var}/
			iput -P -b -f ${cart_orig}/${miptab}/${var}/*/*/${var}*nc ${cart_dest}/${miptab}/${var}/
		    fi
		done
	    done
	    lyear=$yea
	    echo $yea >> $check_file
	fi
    done
done

