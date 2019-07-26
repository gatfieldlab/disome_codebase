#!/bin/bash

declare -a FUNCTIONS=("codon" "6mer" "9mer" "aminoacid" "dipeptide" "tripeptide")
declare -a PROFILES=( "rp" "tr" "di" )
for function in "${FUNCTIONS[@]}"; do
  for profile in "${PROFILES[@]}"; do
  	out_folder="${profile}_${function}"
  	out_file="${profile}_${function}_norm.pdf"
  	meta_file=$( ls ${out_folder}/RUST*meta.txt )
  	rust_file=$( echo ${meta_file%_*} )
    ./plotcodon ${rust_file} ${out_file} ${function} && echo "Done with ${out_folder}" &
    pid_last=$!
  done
done

wait ${pid_last}
exit 0