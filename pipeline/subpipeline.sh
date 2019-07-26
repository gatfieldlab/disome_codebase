#!/usr/bin/env bash

# SUB - splits the jobs for each sample and calls MAPPER
# and performs TopHat mapping and SAM2BAM for each sample

:<<'LICENSE'
"pipeline" is a collection of shell scripts that together provide
a configurable and semi-automated pipeline to trim, filter and map
large sequence files produced by Next Generation Sequencing platforms.

Copyright (C) 2015  A. Bulak Arpat

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
LICENSE

BASE="$1"
SAMPLEDB="$3"
FILEARRAY=""
tophatproc=2
starproc=6

source $CONFIG_FILE

logfile="${__dir}/${2}/logs${SUFFIX}/${BASE}.pipeline.SUB.log"

declare -A errorlog=( ["Have files"]="not processed"
                      ["MAPPER status"]="not processed"
                      ["TOPHAT2 mapping"]="not processed"
                      ["STAR mapping"]="not processed" )

log_and_exit() {
  if [ -f ${logfile} ]; then

    logs ${SUB} $WARNING"Found a log file, most likely from previous runs.\
 Overriding it!"

    rm -f ${logfile}
  fi
  for logtype in "${!errorlog[@]}"; do
    printf "%-16s : %s\n" "$logtype" "${errorlog[$logtype]}" >> ${logfile}
  done
  exit $1
}

if [ "${SAMPLEDB}" == "skip" ]; then
  logs ${SUB} "<${BASE}> Skipping DB search by settings, calling MAPPER."
  errorlog["Have files"]="[ $(timestamp) ] OK - Skipped by --database argument"
else
  i=0
  raw_ext=${FILE_EXT["raw"]}
  while IFS= read -r file; do
      i=$(( $i+1 ))
      FILEARRAY=$(printf "$FILEARRAY\n${file:0:${#file}-${#raw_ext}}")
  done < <(grep -oE "${BASE}\S+${raw_ext}" $SAMPLEDB)
  if (( $i < 1 )); then
    logs ${SUB} $ERROR"<${BASE}> Could not find any files; terminating"
    errorlog["Have files"]="[ $(timestamp) ] ERROR - None found"
    log_and_exit 1
  fi
  logs ${SUB} "<${BASE}> Found $i files; calling MAPPER."
  errorlog["Have files"]="[ $(timestamp) ] OK - Found $i files in DB"
fi
# Spawn MAPPER processes
#-----------------------
cd $2
pipeline_bwt2_single.sh "${BASE}" "${FILEARRAY}"
if (( $? )); then

  logs ${SUB} $ERROR"<${BASE}> One or more of the MAPPER processes exited\
 with a non-zero error code; terminating"

  errorlog["MAPPER status"]="[ $(timestamp) ] ERROR -\
 One or more MAPPER processes failed"

  log_and_exit 1
else
  errorlog["MAPPER status"]="[ $(timestamp) ] OK - All finished"
fi

exit_code=0


# STAR genomic mapping
# --------------------
if (( $MAP_WITH_STAR )); then
  # Figure out what to use as input
  if (( ${#MAPTYPE[@]} == 0 )); then
    if (( ${#PRE_PROC_STEPS[@]} == 0 )); then
      FINALMAP=${INIT_PROC}
    else
      FINALMAP=${PRE_PROC_STEPS[${#PRE_PROC_STEPS[@]}-1]}
    fi
  else
    FINALMAP=${MAPTYPE[${#MAPTYPE[@]}-1]}
  fi
  star_input="fq/${BASE}${FILE_EXT[$FINALMAP]}"
  if [ ! -f ${star_input} ]; then
    exit_code=1
    errorlog["STAR mapping"]="[ $(timestamp) ] ERROR - Input file not found"
    logs ${SUB} ${ERROR}"Input file '${star_input}' could not be found for\
 [${BASE}]. STAR mapping will be skipped!"
  else
    # CPU reservation
    reserve_cpu ${starproc} "${SUB} Waiting for free CPU to map [${BASE}]\
 against genomic database using STAR"

    logs ${SUB} "<${BASE}> Entering now STAR genomic mapping"
    star_out=STAR"${SUFFIX}"/"${BASE}"/
    mkdir -p "${star_out}"
    readFileCmd=-
    use_gzip "${star_input}" && readFileCmd=zcat

    logs ${SUB} "<${BASE}> Mapping against ${BOWTIE2X['genome']}\
 genomic database using STAR..."

    STAR --runThreadN "${starproc}" --genomeDir="${STAR_INDEX['genome']}" \
    --readFilesIn "${star_input}" --outFileNamePrefix "${star_out}" \
    --readFilesCommand "${readFileCmd}" --genomeLoad LoadAndKeep \
    --outSAMtype BAM SortedByCoordinate Unsorted --alignSJDBoverhangMin 1 \
    --alignIntronMax 1000000 --outFilterType BySJout --alignSJoverhangMin 8 \
    --limitBAMsortRAM 15000000000

    cur_err=$?
    if (( $cur_err )); then
      logs ${SUB} ${ERROR}"<${BASE}> Mapping against ${BOWTIE2X['genome']}\
 genomic database with STAR failed."
      errorlog["STAR mapping"]="[ $(timestamp) ] ERROR - Mapping failed"
    else
      logs ${SUB} "<${BASE}> Mapping against ${BOWTIE2X['genome']} genomic\
 database with STAR finished OK"
      echo STAR --genomeDir="${STAR_INDEX['genome']}" \
            --genomeLoad Remove > "${CLEANUP_HOOKS}"
      logs ${SUB} "Indexing [${BASE}] genomic map"
      samtools index "${star_out}"Aligned.sortedByCoord.out.bam

      logs ${SUB} "Creating bedGraph files for [${BASE}]"

      STAR --runMode inputAlignmentsFromBAM --outWigType bedGraph \
      --inputBAMfile "${star_out}"Aligned.sortedByCoord.out.bam \
      --outWigStrand Unstranded --runThreadN "${starproc}" \
      --outFileNamePrefix "${star_out}"

      cur_err=$?
      if (( $cur_err )); then
        logs ${SUB} ${ERROR}"<${BASE}> bedGraph files could not be created"
        errorlog["STAR mapping"]="[ $(timestamp) ] ERROR - bedGraph failed"
      fi
    fi
    release_cpu ${starproc}
    if (( $cur_err )); then
      exit_code=1
    else
      errorlog["STAR mapping"]="[ $(timestamp) ] OK - All finished"
    fi
  fi 
else
  errorlog["STAR mapping"]="[ $(timestamp) ] OK - Skipped by configuration"
fi


# TOPHAT2 genomic mapping
# -----------------------
if (( $MAP_WITH_TOPHAT2 )); then
  # Figure out what to use as input
  if (( ${#MAPTYPE[@]} == 0 )); then
    if (( ${#PRE_PROC_STEPS[@]} == 0 )); then
      FINALMAP=${INIT_PROC}
    else
      FINALMAP=${PRE_PROC_STEPS[${#PRE_PROC_STEPS[@]}-1]}
    fi
  else
    FINALMAP=${MAPTYPE[${#MAPTYPE[@]}-1]}
  fi
  tophat_input="fq/${BASE}${FILE_EXT[$FINALMAP]}"
  if [ ! -f ${tophat_input} ]; then
    exit_code=1
    errorlog["TOPHAT2 mapping"]="[ $(timestamp) ] ERROR - Input file not found"
    logs ${SUB} ${ERROR}"<${BASE}> Input file '${tophat_input}' could not be\
 found; tophat2 mapping will be skipped!"

  else
    # CPU reservation

    reserve_cpu ${tophatproc} "${SUB} Waiting for free CPU to map [${BASE}]\
 against genomic database using TopHat2"

    logs ${SUB} "<${BASE}> Entering now TopHat2 genomic mapping"
    TRINDEX="${BOWTIE2X['genome']}_data/known"
    while [ ! -f $TRINDEX".rev.1.bt2" ]; do
        lockfile -r 1 -! ${TOPHAT_LOCK} > /dev/null 2>&1
        if (( $? )); then
            break
        fi
        logs ${SUB} "<${BASE}> Waiting another instance to finish\
 confirming/creating $TRINDEX ..."
        sleep 120
    done

    if [ ! -f $TRINDEX".rev.1.bt2" ]; then

        logs ${SUB} "Tophat transcriptome-index does not exist\
 for ${BOWTIE2X['genome']}."

        logs ${SUB} "Creating '${TRINDEX}' ..."

        tophat2 --GTF $GTF_DB --transcriptome-index="${TRINDEX}"\
   -p ${tophatproc} -o "${BOWTIE2X['genome']}"_data.log "${BOWTIE2X['genome']}"

        if (( $? )); then

          logs ${SUB} ${WARNING}"Tophat index '${TRINDEX}' could not be\
 created. This will probably cause errors downstream!"

        fi
    fi

    logs ${SUB} "<${BASE}> Mapping against ${BOWTIE2X['genome']}\
 genomic database using TopHat2..."

    rm -f ${TOPHAT_LOCK} &&
    tophat2 --transcriptome-index="${TRINDEX}" -p "${tophatproc}"\
   -o tophat"${SUFFIX}"/"${BASE}" "${BOWTIE2X['genome']}" "${tophat_input}"

    cur_err=$?
    release_cpu ${tophatproc}
    if (( $cur_err )); then
      exit_code=1
      errorlog["TOPHAT2 mapping"]="[ $(timestamp) ] ERROR - Something went wrong"

      logs ${SUB} ${ERROR}"<${BASE}> Mapping against ${BOWTIE2X['genome']}\
 genomic database with TopHat2 failed."

    else
      errorlog["TOPHAT2 mapping"]="[ $(timestamp) ] OK - All finished"

      logs ${SUB} "<${BASE}> Mapping against ${BOWTIE2X['genome']} genomic\
 database finished OK"

    fi
  fi
else
  errorlog["TOPHAT2 mapping"]="[ $(timestamp) ] OK - Skipped by configuration"
fi

cd ..

log_and_exit ${exit_code}
