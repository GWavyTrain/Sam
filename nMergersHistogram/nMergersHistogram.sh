#!/bin/bash

nProc=8

SSmin=26
SSmax=135

step=$(( ( $SSmax - $SSmin ) / $nProc + 1 ))

for (( iProc = 1; iProc < $nProc+1; iProc++ ))
do
  SSL=$(( $SSmin + ( $iProc - 1 ) * $step ))
  SSU=$(( $SSL + $step - 1 ))
    
  if [ $SSU -gt $SSmax ]
    then
    SSU=$SSmax
  fi

  echo "Processor $iProc"
  echo "./nMergersHistogram $SSL $SSU &"
  echo "$SSL-$SSU: $! >> PIDs.txt"
  echo

#  ./nMergersHistogram $SSL $SSU &
#  echo "$SSL-$SSU: $! >> PIDs.txt"

done
