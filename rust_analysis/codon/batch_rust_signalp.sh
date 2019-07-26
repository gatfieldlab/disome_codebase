#!/bin/bash

declare -a FUNCTIONS=( "aminoacid" "codon" )
declare -A PROFILES=( 
    ["rp"]="RP_M7-85_W25-15_T10259_Rgt_1_q1_RUST_01-Apr-19_14.11_profile.txt"
    ["tr"]="TR_M7-85_W25-15_T10259_Rgt_1_q1_RUST_01-Apr-19_14.31_profile.txt"
    ["di"]="Di_M7-85_W25-15_T10259_Rgt_1_q1_RUST_01-Apr-19_14.08_profile.txt" )
declare -a SIGNALP=( "include" "exclude" )
declare -a PIDS=( )

for function in "${FUNCTIONS[@]}"; do
  for profile in "${!PROFILES[@]}"; do
    for signalp_action in "${SIGNALP[@]}"; do
  	out_folder="${profile}_${function}_signalp_${signalp_action}"
	rust.py --quantum -o ${out_folder} --${signalp_action} signalp_trids.txt ${function} profile ${PROFILES[${profile}]} && echo "Done with ${out_folder}" &
	PIDS+=($!)
    done
  done
done

for pid in ${PIDS[*]}; do
    wait $pid
done

exit 0