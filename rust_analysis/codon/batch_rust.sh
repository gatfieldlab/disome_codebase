#!/bin/bash

declare -a FUNCTIONS=("codon" "6mer" "9mer" "aminoacid" "dipeptide" "tripeptide")
declare -A PROFILES=( ["rp"]="RP_M30-30_W30-20_T10259_Rgt_1_mean_RUST_11-Jun-18_12.47_profile.txt"
	                  ["tr"]="TR_M30-30_W30-20_T10259_Rgt_1_mean_RUST_23-Apr-18_14.40_profile.txt"
	                  ["di"]="Di_M30-30_W30-20_T10259_Rgt_1_mean_RUST_23-Apr-18_14.23_profile.txt" )
for function in "${FUNCTIONS[@]}"; do
  for profile in "${!PROFILES[@]}"; do
  	out_folder="${profile}_${function}"
    rust.py --quantum -o ${out_folder} ${function} profile ${PROFILES[${profile}]} && echo "Done with ${out_folder}" &
    pid_last=$!
  done
done

wait ${pid_last}
exit 0