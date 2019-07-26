#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""gquad_predict.py
Identifies potential stretches of sequences for formation of
G-quadruplex structures. It uses the formula from:

Genes Dev. 2015 Mar 15;29(6):617-29. doi: 10.1101/gad.254631.114.
The RNA helicase MOV10L1 binds piRNA precursors to initiate
piRNA processing.

From the supplementary material of this publication:

"The genome wide G-quadruplex prediction was performed using a sensitive general
pattern: G{L1} - N{L2} - G{L3} -N{L4} - G{L5} - N{L6} - G{L7}, where G
corresponds to guanine and N to any nucleotide including G (adapted from (Ryvkin
et al., 2010; Todd et al., 2005)). L1 to L7 correspond to the length of the
stretch of G/N and are allowed to vary independently between a minimum and
maximum value. For the loose prediction L1, L3, L5 and L7 vary from 2 to 4, while
L2, L4 and L6 vary from 1 to 7. For the strict prediction, L1, L3, L5 and L7 can
be 3 or 4, while L2, L4 and L6 vary from 1 to 7."
"""


import sys
import re


__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2018,2019 Bulak Arpat"
__license__ = "GPLv3"
__version__ = "0.0.1"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"

PRED_REGEX = {
    "loose": re.compile(
        r"G{2,4}.{1,7}G{2,4}.{1,7}G{2,4}.{1,7}G{2,4}"),
    "strict": re.compile(
        r"G{3,4}.{1,7}G{3,4}.{1,7}G{3,4}.{1,7}G{3,4}"),
    "strictA": re.compile(
        r"A{3,4}.{1,7}A{3,4}.{1,7}A{3,4}.{1,7}A{3,4}"),
    "strictT": re.compile(
        r"T{3,4}.{1,7}T{3,4}.{1,7}T{3,4}.{1,7}T{3,4}"),
    "strictC": re.compile(
        r"C{3,4}.{1,7}C{3,4}.{1,7}C{3,4}.{1,7}C{3,4}")
}
DEF_PRED_TYPE = "strict"
OUTPUTS = ("seq_like", "position")


class GquadException(Exception):
    """A simple custom Exception class for gquad_predict.py
    """

    def __init__(self, m):
        message = "[gquad_predict] {}".format(m)
        super(GquadException, self).__init__(message)


class GquadPredictor(object):
    """Main class
    """

    def __init__(self, pred_type, output_type):
        if not (isinstance(pred_type, str) and pred_type in PRED_REGEX):
            raise GquadException(
                "pred_type has to be a 'str' and one of: {}".format(
                    ", ".join(PRED_REGEX.keys())))
        if not (isinstance(output_type, str) and output_type in OUTPUTS):
            raise GquadException(
                "output_type has to be a 'str' and one of: {}".format(
                    ", ".join(OUTPUTS)))
        self.pred_type = pred_type
        self.pred_regex = PRED_REGEX[pred_type]
        self.predictor = getattr(self, "_{}".format(output_type))

    def _predict(self, seq):
        predictions = self.pred_regex.finditer(seq)
        return predictions

    def _seq_like(self, seq):
        res = ""
        pos = 0
        predictions = self._predict(seq)
        for gquad in predictions:
            end = gquad.end()
            if end <= pos:
                continue
            start = max(pos, gquad.start())
            res += "N" * (start - pos)
            res += "G" * (end - start)
            pos = end
        res += "N" * (len(seq) - pos)
        return res

    def _position(self):
        raise NotImplementedError()

    def predict(self, seq):
        """Main function for predicting
        """
        if not (isinstance(seq, str) and str):
            raise GquadException("seq has to be a non-empty 'str'")
        return self.predictor(seq)


def print_fasta(handle, seq_id, seq, line_len=80):
    """Helper function to print some strings in a FASTA like format
    """
    handle.write(">{}\n".format(seq_id))
    for i in range(0, len(seq), line_len):
        handle.write("{}\n".format(seq[i:i + line_len]))


def process_fasta(fasta_file):
    """
    Reads the sequences or similar data from a fasta(-like) file
    """
    tr_seq = {}
    with open(fasta_file) as seq_handle:
        for line in seq_handle:
            line = line.strip()
            if line[0] == '>':
                try:
                    tr_split = line.split(' ')
                    tr_id = tr_split[0][1:]
                    tr_seq[tr_id] = ''
                except IOError:
                    sys.stderr.write(
                        "Parsing of fasta file '{}' failed at: {}\n".format(
                            fasta_file, line))
            else:
                tr_seq[tr_id] += line
    return tr_seq


def main():
    if len(sys.argv) > 2:
        pred_type = sys.argv[2]
    else:
        pred_type = DEF_PRED_TYPE
    gquad = GquadPredictor(pred_type=pred_type,
                           output_type="seq_like")
    tr_seq = process_fasta(sys.argv[1])
    handle = sys.stdout
    for tr_id, seq in tr_seq.items():
        print_fasta(handle, tr_id, gquad.predict(seq), line_len=80)
    return 0

if __name__ == "__main__":
    sys.exit(main())
