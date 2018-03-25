#!/bin/bash

Nproc=8

SSmin=26
SSmax=135

step=$(( ( $SSmax - $SSmin ) / $Nproc + 1 ))

for (( i = 0 ; i < $Nproc ; i++ ))

  do
    SSL=$(( $SSmin + i * step ))
    SSU=$(( $SSmin + ( i + 1 ) * step ))
    SSUf=$(( $SSU - 1 ))
    
    if [ $SSUf -gt $SSmax ]
      then
      SSUf=$SSmax
    fi

    NewDir=strain_$SSL-$SSUf

    mkdir $NewDir
    cp strain $NewDir

    cd $NewDir
    echo "./strain $SSL $SSUf &"
    cd ..
  done
