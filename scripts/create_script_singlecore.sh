#!/bin/bash

for INDEX in 0
do
  for aux in `seq 0 1 31`
  do
    aux1=`echo "0.6+0.04*$aux" | bc -l`
    BETA=`printf %.6f $aux1`
  
    scriptname="script_cp1_${BETA}_${INDEX}.sh"
  
    echo "#!/bin/bash" > $scriptname
    echo " " >> $scriptname
    echo "cd $PWD" >> $scriptname
    
    echo "./cp1 < input_${BETA}_${INDEX}.in &" >> $scriptname
    
    echo "wait" >> $scriptname
    chmod +x $scriptname
  done
done  
  
