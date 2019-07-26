#!/bin/bash

# MAPPER - loops BWT2 over different DBs
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

source ${CONFIG_FILE}

# Get sample BASE name from first argument or die
if [ $# -eq 0 ] ; then
  logs ${MAPPER} ${ERROR}"Sample BASE name is missing! MAPPER terminates."
  exit 1
else
  BASE=$1
  shift
  files=($@)
fi
if (( ${#files[@]} == 0 )); then
  files=( ${BASE} )
fi

# Type is now set by conf
: ${DEFAULT_TYPE:=$(GET_TYPE ${files[0]})}
TY=$DEFAULT_TYPE

declare -a failed_steps

if command -v bowtie2 >/dev/null 2>&1; then
  BOW2VER=$(bowtie2 --version 2>&1)
else
  logs ${MAPPER} ${ERROR}"Could not find 'bowtie2'. MAPPER terminates."
  failed_steps+=("bowtie2")
fi

# Pre-processing setup & checks
# Input
declare -a INPUT
input_failed=0
all_steps=("${PRE_PROC_STEPS[@]}" "${MAPTYPE[@]}")
for file in "${files[@]}"; do
  cur_input="../${DIR_NAME[$INIT_PROC]}/${file}${FILE_EXT[$INIT_PROC]}"
  if [ ! -f ${cur_input} ]; then
    logs ${MAPPER} ${WARNING}"Could not find the INPUT file [${cur_input}] for\
 [$BASE]. MAPPER might terminate."
    if (( ${#all_steps[@]} > 0 )); then
      input_failed=1
    fi
  else
    INPUT+=(${cur_input})
  fi
done
if (( input_failed )); then
  failed_steps+=("input")
  pre_proc_pipe=""
else
  if (( ${#all_steps[@]} == 0 )); then
    pre_proc_pipe=""
  else
    use_gzip ${INPUT[0]} &&
      pre_proc_pipe="gzip -dc ${INPUT[@]}" ||
      pre_proc_pipe="cat ${INPUT[@]}"
  fi
fi
declare -a pipe_annot=("input")

# Pipe setup
for proc_step in "${PRE_PROC_STEPS[@]}"; do
  cur_log_file="../${DIR_NAME[$proc_step]}/${BASE}${SUFFIX}${LOG_EXT[$proc_step]}"
  cur_fq="../${DIR_NAME[$proc_step]}/${BASE}${FILE_EXT[$proc_step]}"
  cur_annot=(${proc_step})
  if [ "${KEEPFASTQ[$proc_step]}" = "yes" ]; then
    use_gzip ${FILE_EXT[$proc_step]} &&
      cur_fork=" | tee >(gzip > ${cur_fq})" ||
      cur_fork=" | tee ${cur_fq}"
    cur_annot+=("${proc_step}:keep_fastq")
  else
    cur_fork=""
  fi
  case "${proc_step}" in
    trim)
      case_command="cutadapt"
      case_version="--version"
      case_opts="-a ${ADAPTER} ${CUTADAPT_OPTS}"
      if (( LOGGING )); then
        case_opts="${case_opts} -"
        case_log="2> ${cur_log_file}"
      else
        case_opts="${case_opts} --quiet -"
        case_log=""
      fi
      ;;
    s-filter)
      case_command="consume"
      case_version="-v"
      case_opts="${filter_low[$TY]} ${filter_high[$TY]}"
      if (( LOGGING )); then
        case_log="2> ${cur_log_file}"
      else
        case_log="2> /dev/null"
      fi
      ;;
    q-filter)
      case_command="fastq_quality_filter"
      case_version="-h | grep -o 'FASTX Toolkit [0-9.]\{1,\}'"
      case_opts="${FASTQ_QUALITY_FILTER_OPT}"
      if (( LOGGING )); then
        case_opts="$case_opts -v"
        case_log="2> ${cur_log_file}"
      else
        case_log="2> /dev/null"
      fi
      ;;
  esac
  if command -v ${case_command} >/dev/null 2>&1; then
    version_str="$( ${case_command} ${case_version} 2>&1 )"
    if [ -z "${version_str}" ]; then
      version_str="UNDEF"
      logs ${MAPPER} ${ERROR}"Could not extract version information for\
 '${case_command}'. Setting to '${version_str}'."
    fi
    pre_proc_pipe="${pre_proc_pipe} |
${case_command} ${case_opts} ${case_log} ${cur_fork}"
    pipe_annot=("${pipe_annot[@]}" "${cur_annot[@]}")
  else
    logs ${MAPPER} ${ERROR}"Could not find command '${case_command}' for\
 the [${proc_step}] of pre-processing. MAPPER terminates."
    failed_steps+=("${proc_step}")
  fi
done

# Are we sane?
[ ${#failed_steps[@]} -eq 0 ] && sanity_passed=true || sanity_passed=false

# Logging - init
LOGFILE=logs"${SUFFIX}"/"${BASE}".pipeline.MAPPER.log
sanity_msg="$( timestamp ) Sanity checks"
if ${sanity_passed}; then
  sanity_msg="${sanity_msg} passed."
else
  sanity_msg="${sanity_msg} NOT passed!
These steps failed prior to processing:
$(IFS=, ; echo "${failed_steps[*]}")"
fi
if (( LOGGING )); then
  cat >$LOGFILE << EOL
==== '${CONF_TYPE}' pipeline MAPPER log file ====
${sanity_msg}
EOL
fi
if ! ${sanity_passed}; then
  exit 1
fi
if (( LOGGING )); then
  bowtie_params="> "
  for bkey in "${!BOWTIE_PARAMS[@]}"; do
    param_items=$( printf "%-20s : %s \n> " "$bkey" "${BOWTIE_PARAMS[$bkey]}" )
    bowtie_params="${bowtie_params}${param_items}"
  done
  cat >>$LOGFILE << EOL

== Variables ==
sample: ${BASE}
read_type: ${TY}
number of CPU cores: $( nproc )
number of parallel processes: ${NPROC}
input found: ${INPUT[@]}
  -- Pre-processing --
steps: ${PRE_PROC_STEPS[@]}
adapter: ${ADAPTER}
cutadapt params: ${CUTADAPT_OPTS}
size filter low: ${filter_low[$TY]}
size filter high: ${filter_high[$TY]}
fastq_quality_filter params: ${FASTQ_QUALITY_FILTER_OPT}
  -- Mapping --
map order: ${MAPTYPE[@]}
map tophat2: ${MAP_WITH_TOPHAT2}
tophat -gtf: ${GTF_DB}
map STAR: ${MAP_WITH_STAR}
STAR path: ${STAR_INDEX}

== Bowtie2 ==
${BOW2VER}

== Bowtie2 params ==
${bowtie_params}

EOL
fi

# Mapping checks and setup
pipe_log_file="logs${SUFFIX}/${BASE}.pipeline.PIPE.log"
pipe_exitcode="logs${SUFFIX}/${BASE}.pipeline.EXITCODES.log"
rm ${pipe_log_file} 2> /dev/null
rm ${pipe_exitcode} 2> /dev/null
mkfifo ${pipe_log_file}
pipelogs() {
    gawk -v pref="$1" -v tstr="$(timestamp)"\
    '{print tstr, pref, $0}' > "${pipe_log_file}"
}
map_pipe_can_work=true
if command -v "samtools" >/dev/null 2>&1; then
  version_str=( $( samtools --version 2>&1 ) )
  major_version=$( get_major_version ${version_str[1]} )
  if (( $major_version < 1 )); then
    logs ${SAM2BAM} ${ERROR}"'samtools' version needs to be at least 1.x"
    map_pipe_can_work=false
  fi
else
  logs ${SAM2BAM} ${ERROR}"'samtools' could not be found."
  map_pipe_can_work=false
fi

if (( ${#all_steps[@]} == 0 )); then
  logs ${MAPPER} ${WARNING}"<${BASE}> No actions were defined in configuration.\
 Skipping this and likely all others..."
  map_pipe_can_work=false
fi
map_pipe=""
i=0
for maptype in "${MAPTYPE[@]}"; do
  i=$((i+1))
  cur_fq="fq/${BASE}${FILE_EXT[${maptype}]}"
  cur_log_file="logs${SUFFIX}/${BASE}.${LIBEXTS[${maptype}]}.bw2.log"
  cur_faidx="${BOWTIE2_INDEXES}${BOWTIE2X[${maptype}]}.fa.fai"
  cur_bam="sam/$BASE.${LIBEXTS[${maptype}]}"
  cur_keep_bam="sam/${BASE}.${LIBEXTS[${maptype}]}.kept_filtered.bam"
  cur_filter_cmds=($(echo ${TOMERGE[${maptype}]} | tr "|;\-," " "))
  # setup the filter_sam command + option
  if [ "${#cur_filter_cmds[@]}" -eq 0 ]; then
    cur_filter='no'
  else
    cur_filter="${cur_filter_cmds[0]}"
    cur_keep_flag=""
  fi
  if [ "${#cur_filter_cmds[@]}" -ge 2 ]; then
    if [[ " k keep " =~ " ${cur_filter_cmds[1]} " ]]; then
      cur_keep_flag="-k -f >(samtools view -Sb -t ${cur_faidx} - > \
 ${cur_keep_bam})"
    fi
  fi
  # get map specific bowtie2 params
  if [ ${BOWTIE_PARAMS[$maptype]+isset} ]; then
    cur_bowtie="${maptype}"
  else
    cur_bowtie="default"
  fi
  # configure SAM output with filtering, sorting and converting to BAM
  if [ "${SORTSAM[${maptype}]}" = "yes" ]; then
    if [ "${KEEPUNSORTED[${maptype}]}" = "yes" ]; then
      cur_sortsam="samtools view -Shu -t ${cur_faidx} - |\
 tee >( samtools view -b - > ${cur_bam}.bam ) |\
 samtools sort -@ 2 -m 16G -O bam -T ${cur_bam} -\
 2> >( pipelogs "[samtools-sort]" ) > ${cur_bam}.sorted.bam"
    else
      cur_sortsam="samtools view -Shu -t ${cur_faidx} - |\
 samtools sort -@ 2 -m 16G -O bam -T ${cur_bam} -\
 2> >( pipelogs "[samtools-sort]" ) > ${cur_bam}.sorted.bam"
    fi
  else
    cur_sortsam="samtools view -Sb -t ${cur_faidx} -\
 2> >( pipelogs "[samtools-view]" ) > ${cur_bam}.bam"
  fi
  need_to_check_faidx=false
  if [ "${KEEPSAM[${maptype}]}" = "yes" ]; then
    need_to_check_faidx=true
    case ${cur_filter} in
      rrna)
        cur_sam=">( filter_sam -c rrna -s ${cur_faidx} ${cur_keep_flag} -\
 2> >( pipelogs "[${BASE}][${maptype}]" )\
 1> >( ${cur_sortsam} ) )"
      ;;
      filter|passtru|basic)
        cur_sam=">( filter_sam -c ${cur_filter} ${cur_keep_flag} -\
 2> >( pipelogs "[${BASE}][${maptype}]" )\
 1> >( ${cur_sortsam} ) )"
      ;;
      no)
        cur_sam="${cur_bam}.sam"
        need_to_check_faidx=false
      ;;
      *)
        logs ${SAM2BAM} ${ERROR}"Supplied command(s) for sam_filter are not\
 valid: '${TOMERGE[${maptype}]}' - could not parse/interpret!"
        map_pipe_can_work=false
      ;;
    esac
  elif [ "${KEEPSAM[$maptype]}" = "no" ]; then
    cur_sam="/dev/null"
  else
    logs ${MAPPER} ${WARNING}"<${BASE}> KEEPSAM info for [$maptype] is not\
 'yes' or 'no'! Interpreting as 'no' and discarding the SAM!"
    cur_sam="/dev/null"
  fi
  if ${need_to_check_faidx} && [ ! -f ${cur_faidx} ]; then
    logs ${SAM2BAM} ${ERROR}"<${BASE}> Required index file '${cur_faidx}' could\
 not be found."
    map_pipe_can_work=false
  fi
  # configure unmapped FASTQ output
  if [ "${KEEPFASTQ[${maptype}]}" = "yes" ]; then
    if [ "$i" -eq "${#MAPTYPE[@]}" ]; then
      use_gzip ${cur_fq} &&
        cur_fork=" >(gzip > ${cur_fq})" ||
        cur_fork=" ${cur_fq}"
    else
      use_gzip ${cur_fq} &&
        cur_fork=" >(tee >(gzip > ${cur_fq}))" ||
        cur_fork=" >(tee ${cur_fq})"
    fi
  else
    if [ "$i" -eq "${#MAPTYPE[@]}" ]; then
      cur_fork=" /dev/null"
    else
      cur_fork=" >(tee /dev/null)"
    fi
  fi
  # create the pipe
	map_pipe=${map_pipe}" |
bowtie2 ${BOWTIE_PARAMS[${cur_bowtie}]} -x ${BOWTIE2X[${maptype}]} --no-sq \
--un ${cur_fork} -S ${cur_sam} - 2>${cur_log_file}"
  pipe_annot+=(${maptype})
done

# log and check pipe status
if (( LOGGING )); then
  echo "== Pipe command [ $(timestamp) ] ==
${pre_proc_pipe} ${map_pipe}" >>$LOGFILE
fi
if ! ${map_pipe_can_work}; then
  rm ${pipe_log_file} 2> /dev/null
  rm ${pipe_exitcode} 2> /dev/null
  if (( ${#all_steps[@]} == 0 )); then
    if (( LOGGING )); then
      echo "
== Pipe exit status [ $(timestamp) ] ==
Empty pipe by configuration - execution skipped!" >>$LOGFILE
    fi
  else
    logs ${MAPPER} ${ERROR}"<${BASE}> Mapping pipe will not work, quitting."
    if (( LOGGING )); then
      echo "
== Pipe exit status [ $(timestamp) ] ==
Pipe was not executed!" >>$LOGFILE
    fi
    cur_err=1
  fi
else
  logs ${MAPPER} "<${BASE}> Now running the processing/mapping pipe"
  # cpu reservation
  PROC=$( echo ${BOWTIE_PARAMS[${cur_bowtie}]} | grep -woP '(?<=-p )[0-9]{1,}' )
  PROC=$((PROC*${#MAPTYPE[@]}*8/10))  # 80% of cpu demand is assigned
  reserve_cpu $PROC "${MAPPER} Waiting for free CPU to map $BASE"
  # pipe log reader in background, will be killed when pipe finishes running
  exec 3>&1  # this is necessary otherwise reader's output in background never
             # reaches to the stdout of the rest of the process
  $( while true; do
    if read -a logline<${pipe_log_file}; then
      opr=${logline[2]}
      tms=${logline[*]:0:2}
      msg=$( echo ${logline[*]:3:${#logline[@]}} | sed -e "s/-/${opr}/" )
      logs ${SAM2BAM} ${tms} ${msg} >&3
    fi
  done ) &
  reader_pid=$!
  trap "kill ${reader_pid} 2> /dev/null" INT
  # execute the pipe
  $( eval "${pre_proc_pipe} ${map_pipe};
    echo \${PIPESTATUS[@]} > ${pipe_exitcode}" ) &
  pipe_pid=$!
  trap "kill ${reader_pid},${pipe_pid} 2> /dev/null" INT
  while kill -0 ${pipe_pid} 2> /dev/null; do
    sleep 1
  done
  trap - INT
  kill ${reader_pid}
  wait ${reader_pid} 2>/dev/null
  if (( LOGGING )); then
    echo "
== Pipe exit status [ $(timestamp) ] ==" >>$LOGFILE
  fi
  # release cpu
  release_cpu $PROC
  # check errors and log
  read -r -a pipe_exitcode < ${pipe_exitcode}
  rm ${pipe_log_file} 2> /dev/null
  rm ${pipe_exitcode} 2> /dev/null
  cur_err=0
  if [ "${#pipe_exitcode[@]}" -ne "${#pipe_annot[@]}" ]; then
    logs ${MAPPER} ${WARNING}"<${BASE}> Pipe exit code array length does not\
 match the pipe annotation array length!"
    if (( LOGGING )); then
      echo "Pipe exit code and pipe annotation are not of same size:
  ${pipe_annot[@]}
  ${pipe_exitcode[@]}" >>$LOGFILE
    fi
    cur_err=1
  else
    for i in ${!pipe_exitcode[*]}; do
      cur_err=$((cur_err+${pipe_exitcode[$i]}))
      if (( LOGGING )); then
        echo "${pipe_annot[$i]}: ${pipe_exitcode[$i]}" >>$LOGFILE
      fi
    done
  fi
fi
final_msg="
Processing/mapping pipe"
if (( $cur_err )); then
  if (( LOGGING )); then
    echo "${final_msg}" failed, see above. Premature abortion. >>${LOGFILE}
  fi
  logs ${MAPPER} ${ERROR}"<${BASE}> Processing/mapping pipe failed!!\
 MAPPER terminates."
  exit 1
else
  if (( LOGGING )); then
    echo "${final_msg}" finished successfully. Bye! >>${LOGFILE}
  fi
  logs ${MAPPER} "<${BASE}> Finished successfully."
  exit 0
fi
