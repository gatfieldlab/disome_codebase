#!/bin/bash

declare -a FUNCTIONS=( "charge3p_bin" "charge5p_bin" )

# Following are not used
# "charge1a" "charge2p" "charge2p_mean" "charge3p_mean" 
# "charge4p_mean" "charge5p_mean" ) # "charge4p_aa")

declare -a PROFILES=( 
    "di_long" "rp_long" "tr_long_all"
    "tr_long_rp" "tr_long_di" )

#"rp" "tr_all" "tr_rp" "tr_di" "di" "rp_long" "di_long" 
#                      "tr_long_all" "tr_long_rp" "tr_long_di" )

for function in "${FUNCTIONS[@]}"; do
  for profile in "${PROFILES[@]}"; do
  	out_folder="${profile}_${function}"
  	out_file="${profile}_${function}_norm.pdf"
  	meta_file=$( ls ${out_folder}/RUST*meta.txt )
  	rust_file=$( echo ${meta_file%_*} )
    ./plotcharge ${rust_file} ${out_file} ${function} && echo "Done with ${out_folder}" &
    pid_last=$!
  done
done

wait ${pid_last}
exit 0