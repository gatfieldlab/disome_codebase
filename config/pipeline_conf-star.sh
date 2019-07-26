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

# Type of experiment (not too long, human readable, only for logs etc)
CONF_TYPE="Standalone genome mapping of 'already filtered' reads"

# Sample configuration
# ====================

DEFAULT_TYPE=""
#function GET_TYPE() { get_read_type_from_first_split "0" "$@"; }
function GET_TYPE() { get_read_type_from_first_chrs 2 "$@"; }

# Should we be logging? 1 is yes, 0 is no. No other possibilities. 0 not tested
LOGGING=1

# Pre-processing configuration
# ============================

# Where to start, we start with already trimmed, size and quality trimmed
# reads that are NOT rRNA or tRNA ('trna' from initally mapping)
INIT_PROC="trna"

# NO pre-processing steps
declare -a PRE_PROC_STEPS=( )

# File extension of tRNA - beware of gzip or not !!
declare -A FILE_EXT=( ["trna"]="_non_tRNA.fastq.gz" )

# Directory structure - ugly as we can't export hash arrays in bash
declare -A DIR_NAME=( ["trna"]="mapping_data/fq" )

# Mapping configuration
# =====================

# NO bowtie2 mapping
declare -a MAPTYPE=( )

# follow the sequential mapping with STAR against 'genome'?
# 0 is no, 1 is yes. No other possibilities.
MAP_WITH_STAR=1

# STAR index(es) for genome (other mappings are not supported yet)
# has to be FULL PATH !
declare -A STAR_INDEX=( ["genome"]="/local/databases/mouse/star/Mmusculus.GRCm38.91" )
