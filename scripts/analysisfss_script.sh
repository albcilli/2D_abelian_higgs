#!/bin/bash

for file in ./data/dati_*_0.900_-1.000.dat

do 
	#nome=$1
	nome=$file
	nomeinf=`echo ./data/dati_130_0.900_-1.000.dat`
	Linf=`echo 130`

	#echo "$nomeinf"
 
	tmp=${nome##./data/dati_}
	#tmp=${nome##dati_}
	tmp1=${tmp%%.dat}

	rand=`echo ${RANDOM}`
	#echo "$rand"
	#echo "6124"

	lungh=`wc -l < $nome`
	words=`wc -w < $nome`
	#echo "$lungh	$words"
	lunghinf=`wc -l < $nomeinf`
	wordsinf=`wc -w < $nomeinf`

	#index=`echo "$tmp1" | cut -d"_" -f1`
	L=`echo "$tmp1" | cut -d"_" -f1`
	beta=`echo "$tmp1" | cut -d"_" -f2`
	lambda=`echo "$tmp1" | cut -d"_" -f3`

	#block size for blocking and bootstrap, termalization steps and number of samples for bootstrap

	#if data comes from merged files use term = 0, else choose it
	block=`echo 3000`
	term=`echo 0`
	sample=`echo 1000`


	tot="$nome
	$nomeinf
	$rand
	$lungh  $words	$lunghinf	$wordsinf
	$block $term    $sample
	$L      $Linf	$beta   $lambda"
	
	#echo "$tot"

	./fsscpncool < <(echo "$tot")
	wait
done


