#!/bin/bash

exp=b050
ye1=2284
ye2=2549

omonvars='so thetao vo'
amonvars='evspsbl'

eccart=ec:/ccff/ece3/${exp}/cmorized/
outgen=/scratch/ms/it/ccff/outgoing/${exp}/

#####################
ysho=$(echo ${exp} | cut -d'b' -f 2)

cd /scratch/ms/it/ccff/tmp/

#for exp in allexps
for ye in $(seq $ye1 $ye2); do
  echo $ye
  ecp ${eccart}cmor_${ye}/${exp}_cmorized_${ye}.part.aa /scratch/ms/it/ccff/tmp/

  tmpcart=/scratch/ms/it/ccff/tmp/${exp}_cmorized_${ye}
  mkdir -p ${tmpcart}

  tar xf ${exp}_cmorized_${ye}.part.aa -C ${tmpcart}

  for var in ${omonvars}; do
    miptab=Omon
    grid=gn
    outcart=${outgen}/${miptab}/${var}/
    mkdir -p ${outcart}

    mv ${tmpcart}/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/stabilization-ssp585-2${ysho}/r1i1p1f1/${miptab}/${var}/${grid}/v20210315/*nc ${outcart}
  done

  for var in ${amonvars}; do
    miptab=Amon
    grid=gr
    outcart=${outgen}/${miptab}/${var}/
    mkdir -p ${outcart}

    mv ${tmpcart}/CMIP6/LongRunMIP/EC-Earth-Consortium/EC-Earth3/stabilization-ssp585-2${ysho}/r1i1p1f1/${miptab}/${var}/${grid}/v20210315/*nc ${outcart}
  done

  rm -rf ${tmpcart}
  rm -f ${exp}_cmorized_${ye}.part.aa
done
