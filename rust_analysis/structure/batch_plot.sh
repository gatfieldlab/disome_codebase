#!/bin/bash

declare -a FUNCTIONS=("secondary_struct_simple") # "secondary_struct" "topology")
declare -a PROFILES=( "rp" "tr_all" "tr_rp" "tr_di" "di" )
for function in "${FUNCTIONS[@]}"; do
  for profile in "${PROFILES[@]}"; do
  	out_folder="${profile}_${function}"
  	out_file="${profile}_${function}_norm.pdf"
  	meta_file=$( ls ${out_folder}/RUST*meta.txt )
  	rust_file=$( echo ${meta_file%_*} )
    ./plotstructure ${rust_file} ${out_file} ${function} && echo "Done with ${out_folder}" &
    pid_last=$!
  done
done

wait ${pid_last}
exit 0