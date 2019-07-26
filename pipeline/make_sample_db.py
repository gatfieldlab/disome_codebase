#!/usr/bin/env python3

# utility to extract a tab-delimited flat text file to
# serve as a simple database of samples, files in the
# GNS output
# -----------------------------------------

'''
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
'''

import sys

for line in sys.stdin:
    filename = line.strip()
    base, ext = filename.split('.', 1)
    try:
        sample, run_id, lane_id, file_id = base.rsplit('_', 3)
    except ValueError:
        sys.stderr.write('Failed parsing: {}\n'.format(base))
        continue
    annot = sample.split('_')
    sample_run_id = '_'.join([sample, run_id, lane_id])
    sys.stdout.write('\t'.join([run_id, lane_id, file_id,
                                sample_run_id, base, filename] + annot) + '\n')
