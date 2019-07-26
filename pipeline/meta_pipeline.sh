#!/usr/bin/env bash

# META - reads arguments, environment and
# calls subpipeline for each sample
# -----------------------------------------

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

# meta_pid=$$
# gpid=$(ps opgid= "${meta_pid}")
# trap "kill -- -${gpid} 2> /dev/null" INT

trap "echo 'Forcing all processes to terminate!' && exit" INT TERM
trap "kill 0" EXIT

source pipeline_common.sh

usage(){
  cat <<EOF

${USELINE}

Usage: $0 [options] SAMPLE_FILE

-o|--output          name of the output folder.
                     default "mapping_data"
-c|--concurrent      number of samples that will be
                     processed concurrently. default 4
-d|--database        name of the sample database file
                     default "sample_db.txt". Special
                     keyword "skip" will force the pipeline
                     to skip the database search for input
                     file(s) detection and input will be
                     created specifically based on 'raw'
                     settings in config-file
-f|--config-file     configuration file to use. default
                     is 'pipeline_conf.sh'
-p|--processes       max number of processes that will
                     run in parallel. default and max
                     value is the number of CPUs
                     returned by 'nproc'
-s|--suffix          suffix to append to output folder
                     names. If set, it will be used for
                     logs, tophat and STAR folder.
                     default is ''

${LICENSE}

EOF
	exit 1
}

outputdir="mapping_data"
concurrent=4
sampledb="sample_db.txt"
NPROC=$(nproc)
proc=${NPROC}
config_file="pipeline_conf.sh"
suffix=""
cleanup_hooks=".cleanup_hooks.sh"

while [[ $# > 1 ]]; do
  key="$1"

  case $key in
    -o|--output)
    outputdir="$2"
    shift # past argument
    ;;
    -c|--concurrent)
    concurrent="$2"
    shift
    ;;
    -d|--database)
    sampledb="$2"
    shift
    ;;
    -f|--config-file)
    config_file="$2"
    shift
    ;;
    -p|--processes)
    proc=$((${NPROC}<$2?${NPROC}:$2))
    shift
    ;;
    -s|--suffix)
    suffix="$2"
    shift
    ;;
    *)
    # unknown option
    logs ${META} ${WARNING}"Skipping unknown parameter $key"
    ;;
  esac
  shift   # past argument or value
done

[ -z $1 ] && usage

samples=$1

if [ ! -f ${samples} ]; then
  logs ${META} ${ERROR}"Could not find sample file ${samples}, quitting."
  exit 1
fi

lower_sampledb=$( get_lowercase ${sampledb} )
if [ "$lower_sampledb" == "skip" ]; then
  sampledb="skip"
elif [ ! -f ${sampledb} ]; then
  logs ${META} ${ERROR}"Could not find sample database file ${sampledb},\
 quitting."
  exit 1
else
  sampledb="${__dir}/${sampledb}"
fi

if (( ${concurrent} >= $(( ${proc} / 2 )) )); then
  ${concurrent}=$(( ${proc} / 2 ))
  logs ${META} ${WARNING}"Max value allowed for ${concurrent} is half of\
 number of PROCESSES."
  logs ${META} ${WARNING}"Therefore, ${concurrent} was set to $${concurrent}"
fi

declare -a SUBFOLDER=("fq" "logs${suffix}" "sam"
                      "tophat${suffix}" "STAR${suffix}")

if [ ! -d "${TRIMMED_DIR}" ]; then
  mkdir ${TRIMMED_DIR}
  logs ${META} ${WARNING}"Trimmed folder '${TRIMMED_DIR}' did not exist,\
 created"
fi

if [ ! -d "${FILTERED_DIR}" ]; then
  mkdir ${FILTERED_DIR}
  logs ${META} ${WARNING}"Filtered folder '${FILTERED_DIR}' did not exist,\
 created"
fi

if [ ! -d "${outputdir}" ]; then
    mkdir ${outputdir}
fi

for fld in "${SUBFOLDER[@]}"
do
    if [ ! -d ${outputdir}/"$fld" ]; then
        mkdir ${outputdir}/"$fld"
    fi
done

cp ${samples} ${outputdir}/

rm -f "${CPU_LOCK}"
rm -f "${TOPHAT_LOCK}"
export NPROC=${proc}
export CONFIG_FILE="${__dir}/${config_file}"
export CLEANUP_HOOKS="${__dir}/${cleanup_hooks}"
export SUFFIX=${suffix}
echo "${proc}" > $CPU_FILE

cat ${samples} | xargs -n 1 -P ${concurrent} -I {}\
 bash -c "logs '${META} <{}> Started processing ...' && "\
"subpipeline.sh {} ${outputdir} ${sampledb} && "\
"logs '${META} <{}> Finished processing.'"

logs ${META} "Running final cleanup hooks..."
if [ ! -f ${CLEANUP_HOOKS} ]; then
  logs ${META} ${WARNING} "Could not find cleanup file ${CLEANUP_HOOKS}.\
 This might create a problem if certain data (STAR genome) is still attached\
 to shared memory!"
else
  source ${CLEANUP_HOOKS}
  rm ${CLEANUP_HOOKS}
fi

logs ${META} "Finished everything. Bye."
exit 0
