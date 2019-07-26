# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 16:49:41 2013

@author: barpat
"""

import quicksect
import operator
from .utils import getInfo


#
# A generic counter class
#
class Counter:
    def __init__(self):
        self.count = {}
    def init_key(self, key):
        self.count[key] = 0
    def fast_tik(self, key):
        self.count[key] += 1
    def tik(self, key):
        if not self.count.has_key(key):
            self.init_key(key)
        self.count[key] += 1
    def get(self, key):
        return self.count[key]
    def destroy(self, key):
        del(self.count[key])
    def report(self, file_handler, i=1, reverse=True):
        for k,v in sorted(self.count.iteritems(), key=operator.itemgetter(i), reverse=reverse):
            file_handler.write("%s\t%i\n" % (k,v))
    def get_all(self):
        return self.count

#
# This function creates a quicksect
#

def create_quicksect_patch(loci, feature_types):
    """
    Create quicksect dictionary for looking up only exons from each single locus
    """
    quicksect_obj = {feature_type: quicksect.IntervalTree() for feature_type in feature_types}
    for locus in loci:
        for feature_type in feature_types:
            for feature in loci[locus]['features'][feature_type]:
                quicksect_obj[feature_type].insert(quicksect.Feature(feature['start'], feature['end']-1,
                                                      chr = feature['chrom'],
                                                      info = feature
                                                      ))
    return quicksect_obj
def create_quicksect(loci, feature_types):
    """
    Create quicksect dictionary for looking up only exons from each single locus
    """
    quicksect_obj = {feature_type: quicksect.IntervalTree() for feature_type in feature_types}
    for locus in loci:
        for feature_type in feature_types:
            for feature in loci[locus].features[feature_type]:
                quicksect_obj[feature_type].insert(quicksect.Feature(feature.iv.start, feature.iv.end-1,
                                                      chr = feature.iv.chrom,
                                                      info = getInfo(feature)
                                                      ))
    return quicksect_obj
#
# Generic function for querying a quicksect object for a genomic interval
#
def get_features_within_interval(contig, beg, end, quicksect_obj):
    """
    Query a quicksect object for features overlapping
    the given interval on a chromosome
    """
    all_overlaps = []
    for feature_type, feature_quicksect in quicksect_obj.iteritems():
        overlaps = feature_quicksect.find(quicksect.Feature(beg, end, chr = contig))
        if not overlaps:
            continue
        for e in overlaps:
            all_overlaps.append((feature_type, e.info["gene_id"],
                                 str(e.info["start"]), str(e.info["end"]),
                                 e.info["source"],
                                 e.info["strand"],
                                 e.info["attr"]))
    return all_overlaps

def get_features_within_interval_patch(contig, beg, end, quicksect_obj, fname_corr=None):
    """
    Query a quicksect object for features overlapping
    the given interval on a chromosome
    """
    all_overlaps = []
    for feature_type, feature_quicksect in quicksect_obj.iteritems():
        if fname_corr:
            fname = fname_corr[feature_type]
        else:
            fname = feature_type
        overlaps = feature_quicksect.find(quicksect.Feature(beg, end, chr = contig))
        if not overlaps:
            continue
        for e in overlaps:
            all_overlaps.append((fname, e.info["gene_id"],
                                 str(e.info["start"]), str(e.info["end"]),
                                 e.info[fname+"_number"],
                                 e.info["source"],
                                 e.info["strand"],
                                 e.info))
    return all_overlaps
