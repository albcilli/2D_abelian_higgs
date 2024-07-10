#!/bin/bash


for file in ./data/fulldata/dati_*_0.900_-1.000.dat
do
	#tmpind=`echo "1"`
        nome0=$file

        tmp=${nome0##./data/fulldata/dati_}
        tmp1=${tmp%%.dat}

        lungh0=`wc -l < $nome0`

        index=`echo "$tmp1" | cut -d"_" -f1`
        L=`echo "$tmp1" | cut -d"_" -f2`
        beta=`echo "$tmp1" | cut -d"_" -f3`
        lambda=`echo "$tmp1" | cut -d"_" -f4`

	nomebis=`echo "./data/dati_${L}_${beta%.3}_${lambda%.3}.dat"`
	
	if (( $lungh0 <= 40000 ));
	then
	#da controllarhe
	a=1	

	else
        	sed -n 40000,${lungh0}p $nome0 >> $nomebis
		#awk -v awk_num=$lungh0 'FNR>=40000 && FNR<=awk_num' $nome >> $nomebis
	fi
	
	lunghres="$((40000-$lungh0))"
	#lunghtmp=`echo "-1"`	
	for fileaux in ./data/fullloaddata/dati_${index}*_${L}_${beta%.3}_${lambda%.3}.dat
        do
                #echo "$fileaux"
                nomeaux=$fileaux

		tmpaux=${nomeaux##./data/fullloaddata/dati_}
                tmp1aux=${tmpaux%%.dat}


                indexaux=`echo "$tmp1aux" | cut -d"_" -f1`
                indexaux2=${indexaux%?}
		indexaux3=`echo "${indexaux: -1}"`
		#echo "$index, $indexaux, $indexaux2, $indexaux3"
		if [ "$indexaux3" == "*" ]; then
                        break
                fi


                if [ "$index" != "$indexaux2" ]; then
                        continue
                fi

                #echo "$indexaux, $indexaux2"

                lunghaux=`wc -l < $nomeaux`
                #lunghres="$((40000-$lungh0))"
                if (( $lunghres <= 0 ));
                then
                        lunghres=`echo "1"`
                fi


                if (( $(($lunghaux-$lunghres)) > 0));
                then
			#lunghtmp="$(($lunghres-$lunghaux))"
			#lunghres="$(($lunghtmp))"
			#continue
			sed -n ${lunghres},${lunghaux}p $nomeaux >> $nomebis
			lunghres=`echo "1"`
		else
			lunghtmp="$(($lunghres-$lunghaux))"
                        lunghres="$(($lunghtmp))"
		
		fi

        done
done


