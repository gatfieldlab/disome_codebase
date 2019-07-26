#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Extracts PhastCons or PhyloP score for a given UID(s) and region(s)
"""

import sys
import get_pcons
import density_interface as di


__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2019, Bulak Arpat"
__license__ = "GPLv3"
__version__ = "0.1.0"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"


class ScoresExtractor(object):
    def __init__(self, cds_file, input_file, score_file,
                 uid_pos=0, region_pos=1, seq_pos=2, den_pos=3):
        self.score_file = get_pcons.PconsFile(score_file)
        self.input_file = open(input_file)
        self.uid_pos = uid_pos
        self.reg_pos = region_pos
        self.den_pos = den_pos
        self.seq_pos = seq_pos
        self.cds = di.process_cds(cds_file, use_flag=True, include_incomplete=True)
    def iterator(self):
        for line in self.input_file:
            parsed = line.strip().split('\t')
            uid, pos, seq, den = (parsed[self.uid_pos], int(parsed[self.reg_pos]),
                                  parsed[self.seq_pos], float(parsed[self.den_pos]))
            gid, trid = uid.split("|", 1)
            tr_info = self.cds[trid]
            start = tr_info[0] + (pos - 10) * 3
            end = tr_info[0] + (pos + 10) * 3
            if (start < tr_info[0] or end > tr_info[1]):
                sys.stderr.write(
                    "Skipping {}: tr > {}:{} region > {}:{}\n".format(
                        uid, tr_info[0], tr_info[1], start, end))
                continue
            _, extract, error = self.score_file.get_pcons(trid)
            if error is not None:
                raise error
            else:
                scores = [l.split('\t')[1] for l in extract.split('\n') if l]
                yield uid, pos, den, seq, scores[start:end]

def main(argv=None):
    if argv is None:
        argv = sys.argv
    score_extractor = ScoresExtractor(argv[1], argv[2], argv[3])
    for uid, pos, den, seq, scores in score_extractor.iterator():
        sys.stdout.write("{}\n".format('\t'.join([uid, str(pos),
                                                  str(den), seq] + scores)))

if __name__ == "__main__":
    sys.exit(main())
