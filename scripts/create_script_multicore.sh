#!/bin/bash

scriptname="script_cp9.sh"
 
echo "#!/bin/bash" > $scriptname
echo " " >> $scriptname
echo "cd $PWD" >> $scriptname

for INDEX in 0 1
do
  for aux in `seq 0 1 15`
  do
    #aux1=`echo "0.50+0.02*$aux" | bc -l`
    #BETA=`printf %.6f $aux1`
    aux2=`echo "5+5*$aux" | bc -l`
    LEN=`printf %d $aux2`

    #echo "./cp9new < input_${BETA}_${INDEX}a.in &" >> $scriptname
    echo "./cp9cool < input_${LEN}_${INDEX}.in &" >> $scriptname
 done
done

#for INDEX in 5 6 7 8
#do
#  for aux in `seq 0 1 3`
#  do
#    #aux1=`echo "0.80+0.03*$aux" | bc -l`
#    #BETA=`printf %.6f $aux1`
#    aux1=`echo "180+15*$aux" | bc -l`
#    LEN=`printf %d $aux1`
#    echo "./cp9new < input_${LEN}_${INDEX}.in &" >> $scriptname
# done
#done

#for INDEX in 0
#do
#  for aux in `seq 0 1 10`
#  do
#    #aux1=`echo "0.80+0.03*$aux" | bc -l`
#    #BETA=`printf %.6f $aux1`
#    aux1=`echo "4+4*$aux" | bc -l`
#    LEN=`printf %d $aux1`
#    echo "./cp9new < input_${LEN}_${INDEX}f.in &" >> $scriptname
# done
#done
  

echo "wait" >> $scriptname
chmod +x $scriptname
  
