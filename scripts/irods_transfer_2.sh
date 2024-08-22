#!/bin/sh

#set -e

## Script to transfer bottino data to DRES_D_QUECLI via irods

#echo $PATH
# Before running
module purge
module load icommands
module list

#echo 'making aliases for icommands'
#ICOMMANDS_IMAGE=/cineca/prod/opt/tools/icommands/4.2.10/binary/icommands.sif

#alias iinit=`singularity exec -B/g100,/g100_scratch,/g100_work,/scratch_local,/tmp $ICOMMANDS_IMAGE iinit`
#alias ils=`singularity exec -B/g100,/g100_scratch,/g100_work,/scratch_local,/tmp $ICOMMANDS_IMAGE ils`
#alias iput=`singularity exec -B/g100,/g100_scratch,/g100_work,/scratch_local,/tmp $ICOMMANDS_IMAGE iput`
#alias imkdir=`singularity exec -B/g100,/g100_scratch,/g100_work,/scratch_local,/tmp $ICOMMANDS_IMAGE iput`

#echo 'trying an ils'
#ils

#export PATH=$PATH:/cineca/prod/opt/tools/icommands/4.2.10/binary/
#echo $PATH

function irods_transfer {
	mem=$1
	dryrun=0

	cart_orig_1=/g100_scratch/userexternal/ffabiano/ece3/
	cart_orig_2=/g100_work/IscrB_QUECLIM/BOTTINO/

	echo '------------------------------------------'
	echo '----------------> Trasferring '$mem
	echo '------------------------------------------'

	check_file='check_irods_trasferred_'${mem}'.dat'
	check_file_2='check_irods_mipvar_ok_'${mem}'.dat'

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
	    if [ ${mem} == 'b00I' ]; then
		    run_id=stabilization-ssp585-2100
	    elif [ ${mem} == 'b65I' ]; then
		    run_id=stabilization-ssp585-2065
	    elif [ ${mem} == 'b80I' ]; then
		    run_id=stabilization-ssp585-2080
	    fi
	else
	    mem_id='r1i1p1f1'
	    run_id=stabilization-ssp585-2${ysho}
	fi
	echo $run_id $mem_id
	cartadd=CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/${run_id}/${mem_id}/
	
	cart_dest=QUECLIM/${cartadd}/

	cos=$(ls -1 ${cartbase} | head -n 1)
	cos2=$(ls -1 ${cartbase} | tail -n 1)
	lyear=${cos2:(-4)}
	echo $lyear >> ${check_file}
 
	yea=${cos:(-4)}
	echo ${yea}
	cart_orig=${cartbase}/cmor_${yea}/${cartadd}/

	if [ -f $check_file_2 ]; then
		lastdone=$(tail -n 1 $check_file_2)
		lastdone=($lastdone)
		lasttab=${lastdone[0]}
		lastvar=${lastdone[1]}
		echo $lasttab
		echo $lastvar
		foundinit=0
		incomplete=0
	else
		foundinit=1
		incomplete=0
	fi
     
	# Transfer to irods
	for miptab in $(ls -1 ${cart_orig} | cat); do
	    if [[ $lasttab != $miptab ]] && (( $foundinit == 0 )); then
		echo $miptab 'already transferred'
		ils ${cart_dest}/${miptab}
		continue
	    else
		echo 'transferring ' $miptab
	    fi

	    for var in $(ls -1 ${cart_orig}/${miptab} | cat); do
		if (( $foundinit == 0 )); then
		    if [[ $lastvar != $var ]]; then
		        echo $var 'already transferred'
   		        ils ${cart_dest}/${miptab}/${var} | wc -l
		        continue
	            else
		        foundinit=1
			incomplete=1
		        echo $var 'already transferred'
   		        ils ${cart_dest}/${miptab}/${var} | wc -l
		        continue
	            fi
		fi

		echo 'transferring ' $var
	        imkdir -p ${cart_dest}/${miptab}/${var}
	        if (( $dryrun == 1 )); then
	            ls ${cart_orig}/${miptab}/${var}/*/*/${var}*nc
	            echo ${cart_dest}/${miptab}/${var}/
		elif (( $incomplete == 1 )); then
		    strin=$(ils ${cart_dest}/${miptab}/${var}/ | tail -1)
		    gigi=$( echo $strin | cut -d '_' -f 7 | cut -d '.' -f 1 )
		    laye=${gigi:(0):(4)}
		    echo $laye

	            echo 'transfer incomplete'
		    for ye in $(seq $laye $lyear); do
			iput -P -b -f ${cartbase}/cmor_${ye}/${cartadd}/${miptab}/${var}/*/*/${var}*nc ${cart_dest}/${miptab}/${var}/
		    done
		    echo $miptab' '$var >> $check_file_2
		    incomplete=0
		else
	            echo 'transfer'
		    iput -P -b -f -T --retries 3 -X checkpoint-file-${mem} ${cartbase}/cmor_*/${cartadd}/${miptab}/${var}/*/*/${var}*nc ${cart_dest}/${miptab}/${var}/
		    echo $miptab' '$var >> $check_file_2
        	fi
	    done
	done
}

function check_irods_transfer {
	mem=$1

	cart_orig_1=/g100_scratch/userexternal/ffabiano/ece3/
	cart_orig_2=/g100_work/IscrB_QUECLIM/BOTTINO/

	echo '------------------------------------------'
	echo '----------------> Checking '$mem
	echo '------------------------------------------'

	check_file='check_irods_trasferred_'${mem}'.dat'
	check_file_2='check_irods_mipvar_ok_'${mem}'.dat'

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
	    if [ ${mem} == 'b00I' ]; then
		    run_id=stabilization-ssp585-2100
	    elif [ ${mem} == 'b65I' ]; then
		    run_id=stabilization-ssp585-2065
	    elif [ ${mem} == 'b80I' ]; then
		    run_id=stabilization-ssp585-2080
	    fi
	else
	    mem_id='r1i1p1f1'
	    run_id=stabilization-ssp585-2${ysho}
	fi
	cartadd=CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/${run_id}/${mem_id}/
	
	cart_dest=QUECLIM/${cartadd}/

	cos=$(ls -1 ${cartbase} | head -n 1)
	cos2=$(ls -1 ${cartbase} | tail -n 1)
	lyear=${cos2:(-4)}
	fyear=${cos:(-4)}
 
	cart_orig=${cartbase}/cmor_${lyear}/${cartadd}/

	# Transfer to irods
	for miptab in $(ls -1 ${cart_orig} | cat); do
            echo 'checking ' $miptab

	    for var in $(ls -1 ${cart_orig}/${miptab} | cat); do
                echo 'checking ' $var

		ils ${cart_dest}/${miptab}/${var}
		if [ $? -eq 0 ]; then
		    echo $var ' folder is there'
		    nucos=$(ils ${cart_dest}/${miptab}/${var} | wc -l)
		    if (( $nucos > 1 )); then
		        strin=$(ils ${cart_dest}/${miptab}/${var}/ | tail -1)
		        gigi=$( echo $strin | cut -d '_' -f 7 | cut -d '.' -f 1 )
		        laye=${gigi:(0):(4)}
   		        if (( $laye < $lyear )); then
	                    echo 'transfer incomplete at year '$laye
		            for ye in $(seq $laye $lyear); do
			        iput -P -b -f ${cartbase}/cmor_${ye}/${cartadd}/${miptab}/${var}/*/*/${var}*nc ${cart_dest}/${miptab}/${var}/
		            done
			elif (( $nucos < $lyear - $fyear + 1 )); then
			    echo 'AAAAAAAAAAAAAAAAAAAAA missing years'
			    echo $miptab $var $nucos >> missing_${mem}.dat
		            iput -P -b -f -T --retries 3 -X checkpoint-file-${mem} ${cartbase}/cmor_*/${cartadd}/${miptab}/${var}/*/*/${var}*nc ${cart_dest}/${miptab}/${var}/
		            #for ye in $(seq $fyear $lyear); do
				#ils ${cart_dest}/${miptab}/${var}/${var}_${miptab}_EC-Earth3_${run_id}_${mem_id}_gr_${ye}01-${ye}12.nc
			    #done
			else
			    echo $var ' ok!' $nucos $fyear $lyear >> check_${mem}.dat
	    	        fi
		    else
	                echo $var ' not transferred'
		        iput -P -b -f -T --retries 3 -X checkpoint-file-${mem} ${cartbase}/cmor_*/${cartadd}/${miptab}/${var}/*/*/${var}*nc ${cart_dest}/${miptab}/${var}/
	   	    fi
		else
		    echo $var ' folder is NOT there'
	            echo $var ' not transferred'
	            imkdir -p ${cart_dest}/${miptab}/${var}
		    iput -P -b -f -T --retries 3 -X checkpoint-file-${mem} ${cartbase}/cmor_*/${cartadd}/${miptab}/${var}/*/*/${var}*nc ${cart_dest}/${miptab}/${var}/
		fi
		    
	    done
	done
}

function irods_transfer_year {
	mem=$1
	yea=$2
	dryrun=0

	cart_orig_1=/g100_scratch/userexternal/ffabiano/ece3/
	cart_orig_2=/g100_work/IscrB_QUECLIM/BOTTINO/

	echo '------------------------------------------'
	echo '----------------> Trasferring '$mem
	echo '------------------------------------------'

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
	    if [ ${mem} == 'b00I' ]; then
		    run_id=stabilization-ssp585-2100
	    elif [ ${mem} == 'b65I' ]; then
		    run_id=stabilization-ssp585-2065
	    elif [ ${mem} == 'b80I' ]; then
		    run_id=stabilization-ssp585-2080
	    fi
	else
		mem_id='r1i1p1f1'
		run_id=stabilization-ssp585-2${ysho}
	fi
	cartadd=CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/${run_id}/${mem_id}/

	cart_dest=QUECLIM/${cartadd}/

	echo '-----> Transferring '${yea}
	cart_orig=${cartbase}/cmor_${yea}/${cartadd}/
	# Transfer to irods
	for miptab in $(ls -1 ${cart_orig} | cat); do
		for var in $(ls -1 ${cart_orig}/${miptab} | cat); do
		    imkdir -p ${cart_dest}/${miptab}/${var}
		    if (( $dryrun == 1 )); then
			ls ${cart_orig}/${miptab}/${var}/*/*/${var}*nc
			echo ${cart_dest}/${miptab}/${var}/
		    else
#			iput -P -b -T --retries 3 -X checkpoint-file-${mem} ${cart_orig}/${miptab}/${var}/*/*/${var}*nc ${cart_dest}/${miptab}/${var}/
			iput -P -b -f ${cart_orig}/${miptab}/${var}/*/*/${var}*nc ${cart_dest}/${miptab}/${var}/
		    fi
		done
	done
}

function irods_restore_var {
	mem=$1
	miptab=$2
	var=$3

	cart_rest=/g100_scratch/userexternal/ffabiano/irods_data/
	mkdir -p $cart_rest

	echo '------------------------------------------'
	echo '----------------> Restoring '$mem $miptab $var
	echo '------------------------------------------'

	ysho=$(echo ${mem} | cut -d'b' -f 2)

	if [ ${mem} == 'b990' ]; then
	    mem_id='r1i1p1f1'
	    run_id=stabilization-hist-1${ysho}
	elif [[ ${mem} == *'I' ]]; then
	    mem_id='r1i1p1f3'
	    if [ ${mem} == 'b00I' ]; then
		    run_id=stabilization-ssp585-2100
	    elif [ ${mem} == 'b65I' ]; then
		    run_id=stabilization-ssp585-2065
	    elif [ ${mem} == 'b80I' ]; then
		    run_id=stabilization-ssp585-2080
	    fi
	else
	    mem_id='r1i1p1f1'
	    run_id=stabilization-ssp585-2${ysho}
	fi
	cartadd=CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/${run_id}/${mem_id}/
	
	cart_irods=QUECLIM/${cartadd}/

	# Restore from irods
	mkdir -p ${cart_rest}/${cartadd}/${miptab}
	iget -P -r -f ${cart_irods}/${miptab}/${var} ${cart_rest}/${cartadd}/${miptab}/

	if [ $? -eq 0 ]; then echo 'Done!'; else echo 'Problema!'; fi	
}


function irods_restore_var_years {
	mem=$1
	miptab=$2
	var=$3
	year1=$4
	year2=$5

	cart_rest=/g100_scratch/userexternal/ffabiano/irods_data/
	mkdir -p $cart_rest

	echo '------------------------------------------'
	echo '----------------> Restoring '$mem $miptab $var
	echo '------------------------------------------'

	ysho=$(echo ${mem} | cut -d'b' -f 2)

	if [ ${mem} == 'b990' ]; then
	    mem_id='r1i1p1f1'
	    run_id=stabilization-hist-1${ysho}
	elif [[ ${mem} == *'I' ]]; then
	    mem_id='r1i1p1f3'
	    if [ ${mem} == 'b00I' ]; then
		    run_id=stabilization-ssp585-2100
	    elif [ ${mem} == 'b65I' ]; then
		    run_id=stabilization-ssp585-2065
	    elif [ ${mem} == 'b80I' ]; then
		    run_id=stabilization-ssp585-2080
	    fi
	else
	    mem_id='r1i1p1f1'
	    run_id=stabilization-ssp585-2${ysho}
	fi
	cartadd=CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/${run_id}/${mem_id}/
	
	cart_irods=QUECLIM/${cartadd}/

	# Restore from irods
	mkdir -p ${cart_rest}/${cartadd}/${miptab}/${var}

	echo 'NOT IMPLEMENTED'

# MORE COMPLICATED: filename changes with miptab	
#	for ye in $(seq $year1 $year2); do
#	    iget -P -b -f ${cart_irods}/${miptab}/${var}/${var}_${miptab}_EC-Earth3_${run_id}_${mem_id}_gr_${ye}01-${ye}12.nc ${cart_rest}/${cartadd}/${miptab}/${var}/
#	    if [ $? -ne 0 ]; then
#		iget -P -b -f ${cart_irods}/${miptab}/${var}/${var}_${miptab}_EC-Earth3_${run_id}_${mem_id}_gn_${ye}01-${ye}12.nc ${cart_rest}/${cartadd}/${miptab}/${var}/
#	    fi
#	done
}
