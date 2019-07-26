#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""extract_pcons.py
Extracts values from a BigWig file for a given GTF file and keyword(s)
Originally written for efficient extraction of phastcons or phylipcons
from large compressed BigWig files from UCSC

Output options are:
  seq_like: fasta-like (discrete)
  position: tsv (continous)
"""


import sys
import re
from collections import defaultdict
import math
# import operator as op
import pyBigWig
from accessories import utils
from genomics import conservation
import density_interface as di

__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2019 Bulak Arpat"
__license__ = "GPLv3"
__version__ = "0.0.1"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"


OUTPUTS = ("seq_like", "position")
GID_REGEX = re.compile(r'gene_id "(.+?)";')
TRID_REGEX = re.compile(r'transcript_id "(.+?)";')
CDS_START, CDS_END, TR_LEN = range(3)


class PconsException(Exception):
    """A simple custom Exception class for extract_pcons.py
    """

    def __init__(self, m):
        message = "[extract_pcons] {}".format(m)
        super(PconsException, self).__init__(message)


class PconsExtractor(object):
    """Main class
    """

    mol_key = {'cds': 'CDS',
               'transcript': 'exon'}
    def __init__(self, bigwig_file, gtf_file,
                 cds_file=None, mol_type='transcript', seq_type=None,
                 output_type='position', additional_keywords = None):
        if not (isinstance(output_type, str) and output_type in OUTPUTS):
            raise PconsException(
                "output_type has to be a 'str' and one of: {}".format(
                    ", ".join(OUTPUTS)))
        try:
            self.bw = pyBigWig.open(bigwig_file)
        except RuntimeError as err:
            raise PconsException(err)
        try:
            keywords = 'feature:{}'.format(self.mol_key[mol_type])
        except KeyError as err:
            raise PconsException("{} is not a valid mol_type".format(err))
        if additional_keywords is not None:
            keywords = ';'.join([keywords, additional_keywords])
        try:
            self.gtf = FilterGTF(gtf_file, keywords = keywords)
        except IOError as err:
            print(err)
            raise PconsException("Can not read '{}'".format(gtf_file))
        if output_type == 'seq_like':
            if mol_type == 'transcript' and seq_type == 'codon' and cds_file is None:
                raise PconsException(
                    "When mol_type is 'transcript' and seq_type is 'codon', a CDS file is required")
            self.tr_db = di.process_cds(cds_file, for_db=True, use_flag=True,
                                        include_incomplete=True)
            self.seq_formatter = {'seq_type': seq_type, 'mol_type': mol_type}
        else:
            self.tr_db = None
            self.seq_formatter = {}
        self.formatter = getattr(self, "_{}".format(output_type))
        self.values = defaultdict(list)

    def _seq_like(self):
        def _nuc_converter(x):
            return bytearray(conservation.phylop_2_chr(x[0]), 'utf-8')
        def _codon_converter(x):
            return bytearray(conservation.phylop_pos(x, 0), 'utf-8')
        def _subcodon_converter(x):
            return bytearray(conservation.phylop_subcodon(x), 'utf-8')
        if self.seq_formatter['seq_type'] == 'nuc':
            converter = _nuc_converter
            stepper = 1
        else:
            #converter = _codon_converter
            converter = _subcodon_converter
            stepper = 3
        if (self.seq_formatter['mol_type'] == 'transcript' and
            self.seq_formatter['seq_type'] == 'codon'):
            get_coord = lambda x, y : (x[CDS_START], x[CDS_END])
            init_seq = b"."
        else:
            get_coord = lambda x, y : (0, y)
            init_seq = b" "
        for uid in self.values:
            tr_id = uid.split("|")[1]
            if not tr_id in self.tr_db:
                continue
            tr_info = self.tr_db[tr_id]
            value_len = len(self.values[uid])
            start, end = get_coord(tr_info, value_len)
            print(uid, tr_info, start, end)
            # for i in range(0, value_len, 30):
            #     print(self.values[uid][i:i+30])
            seq = bytearray(init_seq * value_len)
            for i in range(start, end, stepper):
                seq[i:i+stepper] = converter(self.values[uid][i:i+stepper])
            res = get_fasta(uid, seq.decode())
            # print(res)
            # move = input("Move?")
            yield (uid, res)

    def _position(self):
        for uid in self.values:
            res = ""
            for val in self.values[uid]:
                res += "{}\t{}\n".format(uid, val)
            yield uid, res
                
    def extract(self):
        """Main function for extraction
        """
        for gtf_entry in self.gtf.filter_gtf():
            chrom = utils.ensembl_2_mm_chroms[gtf_entry[0]]
            start = int(gtf_entry[3]) - 1
            end = int(gtf_entry[4])
            reverse = gtf_entry[6] == "-"
            try:
                attrs = {k: v.strip('"') for k, v in
                         [attr.strip().split(' ', 1) for attr in
                          gtf_entry[8].split(';') if attr]}
            except ValueError:
                print("Problem: {}".format(gtf_entry[8]))
            gid = attrs['gene_id']
            trid = attrs['transcript_id']
            uid = "{}|{}".format(gid, trid)
            values = self.bw.values(chrom, start, end)
            if reverse:
                values.reverse()
            self.values[uid].extend(values)


def get_fasta(seq_id, seq, line_len=80):
    """Helper function to print some strings in a FASTA like format
    """
    res = ">{}\n".format(seq_id)
    for i in range(0, len(seq), line_len):
        res += "{}\n".format(seq[i:i + line_len])
    return res


class FilterGTF(object):
    """
    Reads a GTF-formatted file using keywords and returns an iterator
    """
    gtf_headers = ['seqname', 'source', 'feature', 'start', 'end',
                   'score', 'strand', 'frame', 'attribute']
    def __init__(self, gtf_file, keywords = ''):
        self.gtf_file = open(gtf_file)
        self.qualifiers = []
        self.values = []
        self.indx = []
        for keyval in keywords.split(";"):
            keyword, val = keyval.split(":", 1)
            if keyword in self.gtf_headers:
                ind = self.gtf_headers.index(keyword)
                self.values.append(val)
                self.indx.append(ind)
                if ind == 8:
                    self.qualifiers.append(lambda x, y: y in x)
                else:
                    self.qualifiers.append(lambda x, y : x == y)
        self.rind = range(len(self.indx))
    def _checker(self, gtf_split):
        return all([self.qualifiers[i](gtf_split[self.indx[i]], self.values[i])
                    for i in self.rind])
    def filter_gtf(self):
        for line in self.gtf_file:
            line = line.strip()
            if line[0] == '#':
                continue
            gtf_split = line.split('\t')
            if self._checker(gtf_split):
                yield gtf_split

def main():
    pcons = PconsExtractor(sys.argv[1], sys.argv[2],
                           cds_file=sys.argv[3], mol_type='transcript',
                           seq_type='codon', output_type = 'seq_like')
    pcons.extract()
    handle = open('conscores_codonpos_1_2.fa', 'w')
    for uid, res in pcons.formatter():
        print(uid)
        handle.write(res)
    return 0

if __name__ == "__main__":
    sys.exit(main())
