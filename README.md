# disome_codebase

This repository contains code used in the publication "Transcriptome-wide sites of collided ribosomes reveal principles of translational pausing" from the Gatfield lab at the UNIL.

Code is collected under three main categories:

1) Initial preprocessing and mapping of raw reads. For this purpose many bash and some utility Python scripts were used which can be found in 'pipeline' folder.
2) Annotation, filtering of transcript models, counting and RUST analyses were performed predominantly by Python scripts which can be found in 'cds_preparation', 'process_reads', and 'rust' folders.
3) Final analyses were performed with Python and R scripts which can be found in 'rust_analysis' and 'analysis' folders. Common configuration files are collected under 'config' file.

1) Initial preprocessing of raw reads:

Raw reads need to be places in raw_data folder. Mapping scripts use bowtie2 and bowtie2 indexes. All indexes configured in the config files need to be available to bowtie2. Firstly a text database of read files is created by make_sample_db.py. Whole mapping pipeline is called by meta_pipeline.sh.

 /raw_data ls *fastq.gz | make_sample_db.py > ../sample_db.txt

 meta_pipeline.sh --config-file config/pipeline_conf-disome.sh --output OUTPUT_FOLDER --suffix SUFFIX
 meta_pipeline.sh --config-file config/pipeline_conf-star.sh --output OUTPUT_FOLDER --suffix SUFFIX

To get help:

 meta_pipeline.sh --help

2) Creating indexed and compressed read map files and counting:

 make_densities.sh SAMPLES

Where SAMPLES can simply created by

 cut -f 4 sample_db.txt > SAMPLES

Transcript and CDS models were prepared by:

 prepare_cds.sh TR_SAMPLES

where TR_SAMPLES is the sample names for total RNA sequencing. Counting is performed by count.sh:

 count.sh Tr TR_SAMPLES
 count.sh Mono RP_SAMPLES
 count.sh Di,all DI_SAMPLES

where RP_ and DI_SAMPLES are files containing sample names for RiboSeq and DisomeSeq, respectively.

3) Main analysis of counts can be done by the analysis.R settings.R and functions.R within the analysis folder.

Once the counts are saved under count_data folder, whole analysis, which can be configured in settings.R script can be called by

 ./analysis.R OUTPUT_DIR VERSION

where VERSION is a user-chosen version string to differentiate between different analysis runs.

RUST analysis can be performed by the rust.py script found in rust folder:

 rust.py -o OUTPUT_DIR -w 25,15 -m 30,30 --resolution codon --quantum --save-profile codon db --type Mono --cds PATH_TO_CDS_MODEL config/disome_specs.yaml config/Mmusculus.GRCm38.91.cdna.ensembl.fa analysis/density_files_db_autogen.tsv

The 'density_files_db_autogen.tsv' file is generated automatically by 'analysis.R' and saved under analysis folder. Mmusculus.GRCm38.91.cdna.ensembl.fa can be downloaded from Ensembl Database.


 
