scriptname=$1

aux=${scriptname%%.sh}
id_string=${aux##script_}
 
bsub -q local -o out_${id_string} -e err_${id_string} -J run_${id_string} ${PWD}/${scriptname}

