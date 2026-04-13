#!/bin/bash
#SBATCH --job-name=copyssh   # Job name
#SBATCH --output=copyssh_%j.out # Output file
#SBATCH --error=copyssh_%j.err  # Error file
#SBATCH --time=48:00:00                # Wall time limit
#SBATCH --nodes=1                      # Number of nodes
#SBATCH --ntasks-per-node=1            # Number of tasks
#SBATCH --mem=8G                       # Memory per node


#expname=b025
#year1=2025
#year2=2099

expname=$1
# year1=${year1}
# year2=${year2}

echo $expname
#, $year1, $year2

username=fabiano
ipaddress=tintin.bo.isac.cnr.it

# TARGETDIR=/gpfs/scratch/userexternal/pdavini0/ece3/${expname}/output/
# OUTDIR=/scratch/ms/it/ccff/cinbkp/ece3/${expname}/output/
TARGETDIR=/ec/res4/scratch/ecme3038/ece4/${expname}/output/oifs
OUTDIR=/work/users/clima/fabiano/ece4/sim_marianna/${expname}/

#rsync -av --append-verify --progress -e "ssh" ffabiano@login.galileo.cineca.it:/gpfs/scratch/userexternal/pdavini0/ece3/b025/output/Output_20* ece3/b025/output/

copycommand="rsync -av --append-verify --progress"

# loop until succeed for three times
# function smartcopy {
#         ok_year=$1
#         MAX_RETRIES=20; i=0; rcheck=255
#         while ( [[ $rcheck -ne 0 ]] && [[ $i -lt $MAX_RETRIES ]] ) ; do
#                 i=$(($i+1))
# 		echo "$copycommand -e "ssh -T -o Compression=no -x" ${username}@${ipaddress}:${TARGETDIR}Output_${ok_year} ${OUTDIR}"
#                 $copycommand -e "ssh -T -o Compression=no -x" ${username}@${ipaddress}:${TARGETDIR}Output_${ok_year} ${OUTDIR}
#                 rcheck=$?
#         done
# }

function smartcopy {
        MAX_RETRIES=1000; i=0; rcheck=255
        while ( [[ $rcheck -ne 0 ]] && [[ $i -lt $MAX_RETRIES ]] ) ; do
                i=$(($i+1))
		echo "$copycommand -e "ssh -T -o Compression=no -x" ${TARGETDIR} ${username}@${ipaddress}:${OUTDIR}"
                $copycommand -e "ssh -T -o Compression=no -x" ${TARGETDIR} ${username}@${ipaddress}:${OUTDIR}
                rcheck=$?
        done
}


smartcopy
echo 'All done!'

exit 0
