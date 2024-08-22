#!/bin/bash

#SBATCH --job-name archive
#SBATCH --qos=nf
#SBATCH --time=23:57:00
#SBATCH -n 1
#SBATCH --account=spitdav2

# script to archive on ECFS data from ecearth simulations
set -uxe

expname=$1
mem=$2

cd /etc/ecmwf/nfs/dh1_perm_b/ccff/post/archive/

username=fabiano
ipaddress=tintin.bo.isac.cnr.it

INCMORDIR=/work_big/users/fabiano/irods_move/${expname}/${mem}/
OUTCMORDIR=${SCRATCH}/tmp_transfer/
TAPEDIR="ec:/ccff/ece3/bottino_archive/${expname}/"
emkdir -p ${TAPEDIR}
emkdir -p ${TAPEDIR}/${mem}

OUT_DIR=$OUTCMORDIR/${expname}/${mem}/
mkdir -p $OUT_DIR
echo $OUT_DIR

bname="b${expname: -3}"

filecheck="lista_ok_tabvars_${bname}_${mem}.txt"
touch $filecheck

#INFODIR=$PERM/ecearth3/infodir/archive/${bname}

#copycommand="rsync -av --append-verify --progress"
copycommand="rsync -av"

# loop until succeed for three times
function smartcopy_cmor {
		miptab=$1
		var=$2

        MAX_RETRIES=20; i=0; rcheck=255
        while ( [[ $rcheck -ne 0 ]] && [[ $i -lt $MAX_RETRIES ]] ) ; do
                i=$(($i+1))

				echo "$copycommand ${username}@${ipaddress}:${INCMORDIR}/${miptab}/${var} ${OUT_DIR}"
				$copycommand ${username}@${ipaddress}:${INCMORDIR}/${miptab}/${var} ${OUT_DIR}
                rcheck=$?
        done
}

for tabvar in $(cat lista_tabvars.txt); do
    echo $tabvar
    IFS='/' read -ra parts <<< $tabvar
    miptab=${parts[0]}
    var=${parts[1]}

	emkdir -p ${TAPEDIR}/${mem}/${miptab}/

	if grep -Fxq "$tabvar" "$filecheck" # Check if not done already
    then
        echo "Already transferred: $miptab $var"
    else
        echo "Transferring $miptab $var.."

		#narch=$(els ec:/ccff/ece3/${bname}/cmorized/cmor_${year} | wc -l)
		#if [ $? -ne 0 ] || [ "$narch" -eq 0 ]; then
		# not archived

		# CMOR part: archive as a tar file divided in different chunks (limit of 32 GB on ECFS)
		echo "Copying $expname: $miptab $var ..."
		smartcopy_cmor $miptab $var

		echo "Packing NetCDF to a single tar ..."
		# tar the file
		tar cfv ${OUT_DIR}/${var}_${miptab}.tar ${OUT_DIR}/${var}

		echo "Splitting to 32GB chunks ..."
		split -b 32G ${OUT_DIR}/${var}_${miptab}.tar "${OUT_DIR}/${var}_${miptab}.part."
		if [ "$?" -eq 0 ] ; then
			rm ${OUT_DIR}/${var}_${miptab}.tar
		fi

		echo "Archiving partfiles to ECFS ..."
		partfiles=$(ls ${OUT_DIR}/${var}_${miptab}.part.*)
		for partfile in $partfiles ; do
			echo $partfile
			ecp -u  $partfile ${TAPEDIR}/${mem}/${miptab}/
			if [ "$?" -eq 0 ] ; then
				rm $partfile
				echo "Cmorized packed file $partfile has been archived!"
			fi
		done

		echo "Removing all data from scratch..."
		rm -r ${OUT_DIR}/${var}
		echo "${tabvar}" >> $filecheck
	fi

done

