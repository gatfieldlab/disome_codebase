#!/bin/bash

declare -a FUNCTIONS=( 
    "charge3p_bin" "charge5p_bin"
)

# The follwing functions are left out 
# charge1a charge2p charge2p_mean charge3p_mean 
# charge5p_mean charge4p_aa charge4p_mean

declare -A PROFILES=(
#    ["rp"]="RP_M30-30_W30-20_T10259_Rgt_1_mean_RUST_11-Jun-18_12.47_profile.txt"
#    ["rp_long"]="RP_M15-15_W50-20_T10259_Rgt_1_mean_RUST_12-Jun-18_10.23_profile.txt"
#    ["tr_all"]="TR_M30-30_W30-20_T10259_Rgt_1_mean_RUST_23-Apr-18_14.40_profile.txt"
#    ["tr_di"]="TR_M30-30_W30-20_T10259_Rgt_1_mean_RUST_26-Apr-18_16.16_profile.txt"
#    ["tr_rp"]="TR_M30-30_W30-20_T10259_Rgt_1_mean_RUST_11-Jun-18_12.36_profile.txt"
#    ["tr_long_all"]="TR_M15-15_W50-20_T10259_Rgt_1_mean_RUST_12-Jun-18_10.53_profile.txt"
#    ["tr_long_rp"]="TR_M15-15_W50-20_T10259_Rgt_1_mean_RUST_21-Jun-18_12.02_profile.txt"
#    ["tr_long_di"]="TR_M15-15_W50-20_T10259_Rgt_1_mean_RUST_21-Jun-18_12.06_profile.txt"
#    ["di"]="Di_M30-30_W30-20_T10259_Rgt_1_mean_RUST_23-Apr-18_14.23_profile.txt"
#    ["di_long"]="Di_M15-15_W50-20_T10259_Rgt_1_mean_RUST_12-Jun-18_09.29_profile.txt"
    ["di_long"]="Di_M15-15_W50-30_T10259_Rgt_1_mean_RUST_12-Oct-18_19.27_profile.txt"
    ["rp_long"]="RP_M15-15_W50-30_T10259_Rgt_1_mean_RUST_12-Oct-18_19.32_profile.txt"
    ["tr_long_all"]="TR_M15-15_W50-30_T10259_Rgt_1_mean_RUST_12-Oct-18_19.39_profile.txt"
    ["tr_long_di"]="TR_M15-15_W50-30_T10259_Rgt_1_mean_RUST_12-Oct-18_19.43_profile.txt"
    ["tr_long_rp"]="TR_M15-15_W50-30_T10259_Rgt_1_mean_RUST_12-Oct-18_19.48_profile.txt"
)
for function in "${FUNCTIONS[@]}"; do
  for profile in "${!PROFILES[@]}"; do
  	out_folder="${profile}_${function}"
    rust.py --quantum -e signalp_trids.txt -o ${out_folder} ${function} profile ${PROFILES[${profile}]} && echo "Done with ${out_folder}" &
    pid_last=$!
  done
  wait ${pid_last}
done


exit 0