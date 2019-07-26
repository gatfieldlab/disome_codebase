#!/bin/bash

declare -a FUNCTIONS=("conservation") 
# "codon_usage" "codon_usage2p" "codon_stai" "codon_stai_2p"
#                      "codon_tai" "codon_tai_2p")
declare -a PROFILES=( "rp" "tr_all" "tr_rp" "tr_di" "di" ) #"rp_long" "di_long" 
#	                    "tr_long_all" "tr_long_rp" "tr_long_di" )
for function in "${FUNCTIONS[@]}"; do
  for profile in "${PROFILES[@]}"; do
  	out_folder="${profile}_${function}"
  	out_file="${profile}_${function}_norm.pdf"
  	meta_file=$( ls ${out_folder}/RUST*meta.txt )
  	rust_file=$( echo ${meta_file%_*} )
    ./plotconservation ${rust_file} ${out_file} ${function} && echo "Done with ${out_folder}" &
    pid_last=$!
  done
done

wait ${pid_last}
exit 0
