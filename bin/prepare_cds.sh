#!/bin/bash

# The project name is intentionally NOT set by argv
PROJECT="disome"

LIMIT=4
export PROCESS=8
export STRING_DIR="stringtie_data"
export CDS_DIR="cds_preparation"
export MAP_DIR="mapping_data/STAR_genome"
export GTF_FILE="config/Mmusculus.GRCm38.91.gtf"
export OUTFILE="${CDS_DIR}/${PROJECT}_filtered_genes_transcripts.tsv"
PROT_CODING="config/Mmusculus.GRCm38.91_prot_coding_genes.tsv"
FASTA="config/Mmusculus.GRCm38.91.cdna.ensembl.fa"
TR_FILE="${CDS_DIR}/${PROJECT}_filtered_transcripts.tsv"
GENE_FILE="${CDS_DIR}/${PROJECT}_filtered_genes.tsv"
PROT_FILE="${CDS_DIR}/${PROJECT}_filtered_Mmusculus.GRCm38.91_prot_coding_genes.txt"
PREPARED_CDS="${CDS_DIR}/${PROJECT}_prepared_cds.tsv"
PREPARED_LOG="${CDS_DIR}/${PROJECT}_prepared_cds_info.txt"
ANNOT_FILE="${CDS_DIR}/${PROJECT}_gene_annotation.tsv"
P_GTF="${CDS_DIR}/${PROJECT}_filtered.GRCm38.91.gtf"
IGV_GTF="${CDS_DIR}/${PROJECT}_igv.gtf"
IGV_FASTA="${CDS_DIR}/${PROJECT}_igv.fa"
LOCI="${CDS_DIR}/${PROJECT}_loci_IDs.txt"

mkdir -p ${STRING_DIR}
mkdir -p ${CDS_DIR}
touch ${OUTFILE}

echo "Transcript analysis with StringTie"
cat $1 | xargs -n 1 -P $LIMIT -I INFILE bash -c \
'echo "... processing INFILE"; '\
'stringtie ${MAP_DIR}/INFILE/Aligned.sortedByCoord.out.bam '\
'-o ${STRING_DIR}/INFILE/out.gtf -p ${PROCESS} -G ${GTF_FILE} '\
'-A ${STRING_DIR}/gene_abund.tab -C ${STRING_DIR}/cov_refs.gtf -B -e '\
'&& filter_tdata.py ${STRING_DIR}/INFILE/t_data.ctab >> ${OUTFILE} '\
'&& echo "... finished processing of INFILE"'

echo "Sorting and removing redundancy ..."
sort ${OUTFILE} | uniq --repeated > _temp_sort_file && mv _temp_sort_file ${OUTFILE}

echo "Splitting Gene and Transcript IDs ..."
cut -f2 ${OUTFILE} | sort > ${TR_FILE}
cut -f1 ${OUTFILE} | sort | uniq > ${GENE_FILE}

echo "Filtering protein coding genes ..."
filter_by_ids.py --not-skip-firstline ${TR_FILE} ${PROT_CODING} > ${PROT_FILE}

echo "Preparing CDS models for filtered protein coding transcipts"
prepare_cds_models.py --list ${PROT_FILE} ${GTF_FILE} >${PREPARED_CDS} 2>${PREPARED_LOG}

echo "Preparing project GTF"
filter_gtf_by_ids.py --gene-id --tr-id --filter-file ${OUTFILE} ${GTF_FILE} > ${P_GTF}

echo "Preparing project annotations"
filter_summarize_gtf.py --filter-file ${OUTFILE} --gene-id --tr-id \
  -o gene_id,gene_name,gene_biotype,transcript_id,transcript_name,transcript_biotype \
  ${P_GTF} | sort | uniq > ${ANNOT_FILE}

echo "Preparing GTF and FASTA for IGV"
make_gtf_for_igv.py ${ANNOT_FILE} ${PREPARED_CDS} > ${IGV_GTF}
awk 'BEGIN {OFS="|"} {print $1,$3}' ${OUTFILE} > ${LOCI}
filter_fasta_by_ids.py --filter-file ${LOCI} ${FASTA} > ${IGV_FASTA}
samtools faidx ${IGV_FASTA}

echo "Script completed."
exit 0
