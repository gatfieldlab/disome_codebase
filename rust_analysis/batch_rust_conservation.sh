#!/bin/bash

declare -a FUNCTIONS=( "conservation" ) # secondary_struct
declare -A OPTIONS=( ["rp"]=""
	             ["tr_all"]=""
	             ["di"]=""
		     ["tr_di"]="--offset Di:default"
		     ["tr_rp"]="--offset Mono:default" )
declare -A TYPES=( ["rp"]="Mono"
                   ["tr_all"]="Tr"
		   ["tr_di"]="Tr"
		   ["tr_rp"]="Tr"
		   ["di"]="Di" )

for function in "${FUNCTIONS[@]}"; do
    for profile in "${!OPTIONS[@]}"; do
	out_folder="conservation/${profile}_${function}"
	rust.py --quantum -o ${out_folder} ${OPTIONS[${profile}]} -w 30,20 -m 30,30 --save-profile ${function} db --treatment ZT00,ZT02,ZT12 --rep 1,2 --type ${TYPES[${profile}]} --cds ../cds_preparation/disome_prepared_cds.tsv "../config/disome_specs.yaml" ../extra_data/conscores_pos12.fa ../analysis/density_files_db_autogen.tsv && echo "Done with ${out_folder}" &
	pid_last=$!
    done
done

wait ${pid_last}
exit 0
