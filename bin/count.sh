#!/bin/bash
# The project name is intentionally NOT set by argv
PROJECT="disome"

LIMIT=32
export COUNT_DIR="count_pri_data"
export DENSITY_DIR="density_data"
export PROJECT_SPECS="config/${PROJECT}_specs.yaml"
export GENES="cds_preparation/${PROJECT}_filtered_genes.tsv"
export CDS_FILE="cds_preparation/${PROJECT}_prepared_cds.tsv"
export BGZ="_5prime_count.bgz"
export COUNT_EXT="_counts.txt"
export READTYPE="$1"
export FLAG="unique"
[ ! -z "$3" ] && FLAG="$3"

mkdir -p ${COUNT_DIR}

cat $2 | xargs -n 1 -P $LIMIT -I INFILE bash -c \
'echo "Processing INFILE ..."; '\
'count_per_gene.py -c $CDS_FILE -s ${PROJECT_SPECS},${READTYPE} '\
'-g $GENES ${DENSITY_DIR}/INFILE${BGZ} -f ${FLAG} '\
'> ${COUNT_DIR}/INFILE${COUNT_EXT} 2> ${COUNT_DIR}/INFILE${COUNT_EXT}.log '\
'&& echo "Finished processing of INFILE"'

echo " Script completed."
exit 0
