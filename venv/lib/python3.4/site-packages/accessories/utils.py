# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 16:12:04 2013

@author: barpat
"""

from contextlib import contextmanager
import time
import sys


mm_2_ensembl_chroms = {'chr1':'1', 'chr1_GL456210_random':'GL456210.1', 'chr1_GL456211_random':'GL456211.1',
                       'chr1_GL456212_random': 'GL456212.1', 'chr1_GL456213_random': 'GL456213.1',
                       'chr1_GL456221_random': 'GL456221.1', 'chr2': '2', 'chr3': '3', 'chr4': '4',
                       'chr4_GL456216_random': 'GL456216.1', 'chr4_GL456350_random': 'GL456350.1',
                       'chr4_JH584292_random': 'JH584292.1', 'chr4_JH584293_random': 'JH584293.1',
                       'chr4_JH584294_random': 'JH584294.1', 'chr4_JH584295_random': 'JH584295.1',
                       'chr5': '5', 'chr5_GL456354_random': 'GL456354.1',
                       'chr5_JH584296_random': 'JH584296.1', 'chr5_JH584297_random': 'JH584297.1',
                       'chr5_JH584298_random': 'JH584298.1', 'chr5_JH584299_random': 'JH584299.1',
                       'chr6': '6', 'chr7': '7', 'chr7_GL456219_random': 'GL456219.1',
                       'chr8': '8', 'chr9': '9', 'chr10': '10', 'chr11': '11',
                       'chr12': '12', 'chr13': '13', 'chr14': '14', 'chr15': '15',
                       'chr16': '16', 'chr17': '17', 'chr18': '18', 'chr19': '19',
                       'chrX': 'X', 'chrX_GL456233_random': 'GL456233.1', 'chrY': 'Y',
                       'chrY_JH584300_random': 'JH584300.1', 'chrY_JH584301_random': 'JH584301.1',
                       'chrY_JH584302_random': 'JH584302.1', 'chrY_JH584303_random': 'JH584303.1',
                       'chrUn_JH584304': 'JH584304.1', 'chrM': 'MT'}

ensembl_2_mm_chroms = {v:k for k,v in mm_2_ensembl_chroms.items()}


class ChrError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

#
# Context manager for timing
#
@contextmanager
def measureTime(title, handle=sys.stderr, unit='second'):
    t1 = time.clock()
    yield
    t2 = time.clock()
    delta = t2 - t1
    if unit == 'minute':
        delta = delta / 60
    elif unit == 'hour':
        delta = delta / 3600
    else:
        unit = 'second'
    handle.write('%s: %0.2f %s(s) elapsed\n' % (title, delta, unit))

#
# Utility function to extract a info dict from HTSeq.GenomicFeature object
#
def getInfo(feature, inclusive=True):
    d = int(inclusive)
    feature_info = {'start': feature.iv.start, 'end': feature.iv.end-d,
                    'strand': feature.iv.strand, 'source': feature.source,
                    'gene_id': feature.attr['gene_id'], 'attr': feature.attr}
    return feature_info

#
# Utility function to convert between mm and ensembl chrom names
#
def convert_chr(chrom):
    if chrom in mm_2_ensembl_chroms:
        return mm_2_ensembl_chroms[chrom]
    elif chrom in ensembl_2_mm_chroms:
        return ensembl_2_mm_chroms[chrom]
    else:
        raise ChrError("Chromosome name '{}' was not found either in mm or Ensembl\n".format(chrom))

#
# Utility function to convert various true/false strings into boolean variables
#
def check_true(s):
    s = s.lower()
    return s in ['true', '1', 't', 'y', 'yes', 'ok']
