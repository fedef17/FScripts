#!/bin/bash

#SBATCH --job-name archive
#SBATCH --qos=nf
#SBATCH --time=23:57:00
#SBATCH -n 1
#SBATCH --account=spitdav2

# script to archive on ECFS data from ecearth simulations
set -uxe

expname=$1
year1=$2
year2=$3

username=fabiano
ipaddress=tintin.bo.isac.cnr.it

INCMORDIR=/work_big/users/fabiano/irods_move/${expname}/
OUTCMORDIR=${SCRATCH}/tmp_transfer/

bname="b${expname: -3}"

INFODIR=$PERM/ecearth3/infodir/archive/${bname}

#copycommand="rsync -av --append-verify --progress"
copycommand="rsync -avR"

# loop until succeed for three times
function smartcopy_cmor {
        ok_year=$1

		CMOR_YR_DIR=$OUTCMORDIR/cmor_${ok_year}/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/${expname}/
		mkdir -p $CMOR_YR_DIR

        MAX_RETRIES=20; i=0; rcheck=255
        while ( [[ $rcheck -ne 0 ]] && [[ $i -lt $MAX_RETRIES ]] ) ; do
                i=$(($i+1))

				echo "$copycommand ${username}@${ipaddress}:${INCMORDIR}./r1i1p1f1/*/*/*_${ok_year}01*nc ${CMOR_YR_DIR}"
                $copycommand ${username}@${ipaddress}:${INCMORDIR}./r1i1p1f1/*/*/*_${ok_year}01*nc ${CMOR_YR_DIR}
                rcheck=$?
        done

		## change dir structure

		base_dir=${CMOR_YR_DIR}/r1i1p1f1/
		# Loop through all miptab and var directories
		for miptab in "$base_dir"/*/; do
			miptab=${miptab%/}  # Remove trailing slash
			for var in "$miptab"/*/; do
				var=${var%/}  # Remove trailing slash
				for file in "$var"/*; do
					if [[ -f "$file" ]]; then  # Ensure it's a file
						# Extract the 8th field from the filename
						filename=$(basename "$file")
						grid=$(echo "$filename" | cut -d'_' -f8)
						version=$(ncdump -h $file | grep creation | grep -oP '(?<=creation_date = ")\d{4}-\d{2}-\d{2}' | sed 's/-//g')

						# Create the destination directory if it doesn't exist
						destination_dir="$var/$grid/v$version"
						mkdir -p "$destination_dir"

						# Move the file to the new directory
						mv "$file" "$destination_dir/"
					fi
				done
			done
		done

}

# loop on the year, disabled by submitter
for year in $(seq $year1 $year2) ; do

	# Check if year is already archived, if not, archive
	narch=$(els ec:/ccff/ece3/${bname}/cmorized/cmor_${year} | wc -l)
	if [ $? -ne 0 ] || [ "$narch" -eq 0 ]; then
    	# not archived

		# CMOR part: archive as a tar file divided in different chunks (limit of 32 GB on ECFS)
		echo Copying $expname year $year ...
			smartcopy_cmor $year

		echo "Archiving data for year $year..."

		LOGFILE2=$INFODIR/archive_cmorized_${expname}_${year}.txt
		[ -f $LOGFILE2 ] && { head $LOGFILE2; exit 1; }

		echo "Copying cmorized folder for year $year..."
		
		TREEDIR=ece3/${bname}/cmorized/cmor_${year}
		PACKDIR=${OUTCMORDIR}Pack_${year}
		mkdir -p $PACKDIR
		emkdir -p ec:/ccff/${TREEDIR}

		compresslist=CMIP6
		echo "Packing NetCDF to a single tar ..."
		rm -f ${PACKDIR}/${bname}_cmorized_${year}.tar
		# set striping to optimize I/O
		lfs setstripe -c 20 ${PACKDIR}
		# tar the file
		tar cfv ${PACKDIR}/${bname}_cmorized_${year}.tar -C ${OUTCMORDIR}cmor_${year} ${compresslist}

		echo "Splitting to 32GB chunks ..."
		split -b 32G ${PACKDIR}/${bname}_cmorized_${year}.tar "${PACKDIR}/${bname}_cmorized_${year}.part."
		if [ "$?" -eq 0 ] ; then
			rm ${PACKDIR}/${bname}_cmorized_${year}.tar
		fi

		echo "Archiving partfiles to ECFS ..."
		partfiles=$(ls ${PACKDIR}/${bname}_cmorized_${year}.part.*)
		for partfile in $partfiles ; do
			echo $partfile
			ecp -u  $partfile ec:/ccff/$TREEDIR/
			if [ "$?" -eq 0 ] ; then
				rm $partfile
				echo "Cmorized packed file $partfile for $year has been archived!" >> $LOGFILE2
			fi
		done
		rmdir ${PACKDIR}

		echo "Removing all data from scratch..."
		rm -rf ${OUTCMORDIR}cmor_${year}

	else
    	echo "$year already archived"
	fi

done

