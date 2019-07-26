#!/bin/bash
# The project name is intentionally NOT set by argv
PROJECT="disome"

LIMIT=32
export DENSITY_DIR="density_data_merged"
export BAM_DIR="mapping_merged/sam"
export BAM_EXT=".mouse_cDNA.bam"
export EXPRESSED="cds_preparation_merged/${PROJECT}_filtered_transcripts.tsv"
export PROTCODING="config/Mmusculus.GRCm38.91_prot_coding_gene_IDs.txt"

mkdir -p ${DENSITY_DIR}

cat $1 | xargs -n 1 -P $LIMIT -I INFILE bash -c \
'echo "Processing INFILE ..."; '\
'make_density.py  --expressed ${EXPRESSED} '\
'--list-prot ${PROTCODING} --output - ${BAM_DIR}/INFILE${BAM_EXT} '\
'| qsort | bgzip -c > ${DENSITY_DIR}/INFILE_5prime_count.bgz '\
'&& index_density.py ${DENSITY_DIR}/INFILE_5prime_count.bgz '\
'&& echo "Finished processing of INFILE"'

echo " Script completed."
exit 0
