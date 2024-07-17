#!/bin/bash

#---REQUIRED ARGUMENTS----#

usage()
{
   echo "Usage:"
   echo "      ./launch_copy.sh expname year1 year2"
   echo
}


expname=$1
year1=$2
year2=$3

if [ $# -lt 3 ]; then
   usage
   exit 1
fi

CALOG=/scratch/ms/it/ccff/cinbkp/logrsy/

MACHINE_OPT="-l EC_billing_account=spitfabi -q ns -l EC_memory_per_task=10GB -l EC_hyperthreads=1"

BASE_OPT="expname=$expname,year1=$year1,year2=$year2"
#BASE_OPT="$expname,$year1,$year2"

JOB_ATM='qsub ${MACHINE_OPT} -l walltime=47:59:00 -v ${BASE_OPT} -N cpf-${expname}-${year1}-${year2}
        -o ${CALOG}cpf_${expname}_${year1}_${year2}.out -e ${CALOG}cpf_${expname}_${year1}_${year2}.err
        ./copyfrom.sh'

JOBID=$(eval ${JOB_ATM})
set -- $JOBID
JOBID=$4
echo $JOBID
