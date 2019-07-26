#!/bin/bash

LIMIT=32
export SIZE_DIR="insert_lengths"
export MAPPING_DIR="mapping_data/sam"
export MAP_TYPE="$1"
export SIZE_EXT="_insert_sizes.tsv"

mkdir -p ${SIZE_DIR}

cat $2 | xargs -n 1 -P $LIMIT -I INFILE bash -c \
'echo "Processing INFILE ..."; '\
'samtools view ${MAPPING_DIR}/INFILE${MAP_TYPE} |'\
'read_size.py > ${SIZE_DIR}/INFILE${SIZE_EXT} '\
'&& echo "Finished processing of INFILE"'

echo " Script completed."
exit 0
