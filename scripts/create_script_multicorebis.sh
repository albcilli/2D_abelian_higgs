#!/bin/bash

scriptname="script_cp9.sh"
 
echo "#!/bin/bash" > $scriptname
echo " " >> $scriptname
echo "cd $PWD" >> $scriptname

#for INDEX in 1 2
#do
#  for aux in `seq 0 1 15`
#  do
#    #aux1=`echo "0.85+0.01*$aux" | bc -l`
#    #BETA=`printf %.6f $aux1`
#    aux2=`echo "5+5*$aux" | bc -l`
#    LEN=`printf %d $aux2`
#
#    #echo "./cp6 < input_${BETA}_${INDEX}a.in &" >> $scriptname
#    echo "./cp6 < input_${LEN}_${INDEX}.in &" >> $scriptname
# done
#done

for file in input_*.in
do
   echo "./cp9cool < $file &" >> $scriptname
done

#for INDEX in 0
#do
#  for aux in `seq 0 1 15`
#  do
#    #aux1=`echo "0.81+0.01*$aux" | bc -l`
#    #BETA=`printf %.6f $aux1`
#    aux1=`echo "5+5*$aux" | bc -l`
#    LEN=`printf %d $aux1`
#    echo "./cp6 < input_${LEN}_${INDEX}a.in &" >> $scriptname
#    #echo "./cp6 < input_${BETA}_${INDEX}b.in &" >> $scriptname
# done
#done

#for INDEX in 0 1 2 3
#do
#  for aux in `seq 0 1 3`
#  do
#    aux1=`echo "0.77+0.01*$aux" | bc -l`
#    BETA=`printf %.6f $aux1`
#    #aux1=`echo "5+5*$aux" | bc -l`
#    #LEN=`printf %d $aux1`
#    #echo "./cp6 < input_${LEN}_${INDEX}b.in &" >> $scriptname
#    echo "./cp6 < input_${BETA}_${INDEX}c.in &" >> $scriptname
# done
#done
  
echo "wait" >> $scriptname
chmod +x $scriptname
  
