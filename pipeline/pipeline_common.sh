#!/usr/bin/env bash

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

export PROGRAM='pipeline'
export VERSION="v0.1"
export __dir="$(pwd)"
export USELINE="$PROGRAM $VERSION (C) 2015  A. Bulak Arpat"
read -r -d '' LICENSE <<'EOM'
This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute it
under certain conditions; see `LICENSE.txt' for details.
EOM
export LICENSE
export WARNING="\e[0;33m"
export ERROR="\e[0;31m"
export CPU_LOCK="${__dir}/cpu.lock"
export CPU_FILE="${__dir}/cpu.txt"
export TOPHAT_LOCK="${__dir}/tophat.lock"

# color codes for logging
export META="\e[0m\e[37m[META]\e[0m"
export SUB="\e[0m\e[35m[SUB]\e[0m"
export MAPPER="\e[0m\e[34m[MAPPER]\e[0m"
export SAM2BAM="\e[0m\e[36m[SAM2BAM]\e[0m"

# default directory structure
export RAW_DIR="raw_data"
export TRIMMED_DIR="trimmed_data"
export FILTERED_DIR="filtered_data"

logs() {
  echo -e $1 "${*:2} \e[0m"
}

timestamp() {
  date +"%D %T"
}

get_read_type_from_first_chrs() {
  echo "${2:0:$1}"
}

get_read_type_from_first_split() {
  local IFS=$1
  read -a arr <<< "$2"
  echo "${arr[0]}"
}

reserve_cpu() {
  CURPROC=$1
  MESSAGE=${*:2}
  while true; do
      lockfile $CPU_LOCK
      AVAIL_CPU=$(<$CPU_FILE)
      if (( $AVAIL_CPU < $CURPROC )); then
          rm -f $CPU_LOCK
          logs $MESSAGE
          sleep 120
      else
          REMAIN_CPU=$(( $AVAIL_CPU - $CURPROC ))
          echo $REMAIN_CPU > $CPU_FILE
          rm -f $CPU_LOCK
          break
      fi
  done
}

release_cpu () {
  CURPROC=$1
  lockfile $CPU_LOCK
  AVAIL_CPU=$(<$CPU_FILE)
  CUR_CPU=$(( $AVAIL_CPU + $CURPROC ))
  MAX_CPU=$(($CUR_CPU<$NPROC?$CUR_CPU:$NPROC))
  echo $MAX_CPU > $CPU_FILE
  rm -f $CPU_LOCK
}

use_gzip () {
  local FEXT=$1
  if [ ${FEXT:${#FEXT}-3:${#FEXT}} = '.gz' ]; then
    return 0
  else
    return 1
  fi
}

get_major_version () {
    local IFS=.
    local ver=($1)
    echo "${ver[0]}"
}

get_lowercase () {
  echo $( echo "$1" | tr '[:upper:]' '[:lower:]' )
}


export -f get_read_type_from_first_chrs
export -f get_read_type_from_first_split
export -f logs
export -f reserve_cpu
export -f release_cpu
export -f timestamp
export -f use_gzip
export -f get_major_version
export -f get_lowercase
