#!/bin/bash

nProc=8

SSmin=26
SSmax=135

step=$(( ( $SSmax - $SSmin ) / $nProc + 1 ))

for (( iProc = 1 ; iProc < $nProc+1 ; iProc++ ))
do

  SSL=$(( $SSmin + ( $iProc - 1 )* $step ))
  SSU=$(( $SSL + $step - 1 ))
    
  if [ $SSU -gt $SSmax ]
    then
    SSU=$SSmax
  fi

  NewDir=strain_$SSL-$SSU

  echo "Processor $iProc"
  echo "mkdir $NewDir"
  echo "cp strain $NewDir"
  echo "cd $NewDir"
  echo "./strain $SSL $SSU &"
  echo "cd .."
  echo

#  mkdir $NewDir
#  cp strain $NewDir
#  cd $NewDir
#  ./strain $SSL $SSU &
#  cd ..
  
done
