#!/bin/bash

mkdir -p confs
mkdir -p data
mkdir -p logs

for INDEX in `seq 6 1 13`

do
  for aux in `seq 0 1 0`

  do
    #aux1=`echo "0.50+0.02*$aux" | bc -l`
    #BETA=`printf %.6f $aux1`
    aux2=`echo "130+8*$aux" | bc -l`
    LEN=`printf %d $aux2`
  
    #sed "s/BETA/${BETA}/" simparam  | sed "s/RANDOM/${RANDOM}/" | sed "s/INDEX/${INDEX}/" > input_${BETA}_${INDEX}.in
    sed "s/LEN/${LEN}/" simparam  | sed "s/RANDOM/${RANDOM}/" | sed "s/INDEX/${INDEX}/" > input_${LEN}_${INDEX}.in
    #sed "s/BETA/${BETA}/" simparam  | sed "s/LEN/${LEN}/" | sed "s/RANDOM/${RANDOM}/" | sed "s/INDEX/${INDEX}/" > input_${BETA}_${INDEX}.in
  done
done
