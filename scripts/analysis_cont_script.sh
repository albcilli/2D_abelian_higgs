#!/bin/bash


for file in ./data/dati_*_*_0.000.dat

do 
	#nome=$1
	nome=$file

	#echo "$nome"
 
	tmp=${nome##./data/dati_}
	#tmp=${nome##dati_}
	tmp1=${tmp%%.dat}

	rand=`echo ${RANDOM}`
	#echo "$rand"
	#echo "6124"

	lungh=`wc -l < $nome`
	words=`wc -w < $nome`
	#echo "$lungh	$words"

	#index=`echo "$tmp1" | cut -d"_" -f1`
	L=`echo "$tmp1" | cut -d"_" -f1`
	beta=`echo "$tmp1" | cut -d"_" -f2`
	lambda=`echo "$tmp1" | cut -d"_" -f3`

	#block size for blocking and bootstrap, termalization steps and number of samples for bootstrap

	#if data comes from merged files use term = 0, else choose it
	block=`echo 5000`
	term=`echo 0`
	sample=`echo 1000`


	#echo "$block	$term	$sample"

	#echo "$L $beta	$lambda"

	tot="$nome
	$rand
	$lungh  $words
	$block $term    $sample
	$L      $beta   $lambda"
	
	#echo "$tot"

	./ancontcp9 < <(echo "$tot")
	wait
done

