scriptname=$1

aux=${scriptname%%.sh}
id_string=${aux##script_}
 
bsub -q longparallel -o out_${id_string} -e err_${id_string} -J run_${id_string} -n32 ${PWD}/${scriptname}

