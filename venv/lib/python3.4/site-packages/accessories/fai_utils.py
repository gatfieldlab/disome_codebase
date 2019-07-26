# -*- coding: utf-8 -*-
"""
Created on Tue Aug  5 11:16:15 2014

@author: barpat
"""
import math



class IndexedFasta:
    """
    A class to search, parse information from fasta index files (i.e .fai), and
    retrieve related sequences from accompanying fasta files (i.e .fa)
    """
    GeneID, Length, Offset, Linebases, Linewidth = range(5)
    faidx_types = (str, int, int, int, int)
    
    def __init__(self, fasta_handle, fai_handle=None):
        self.fasta_handle = fasta_handle
        self.fai_handle = fai_handle
        
    def grep_fai(self, search_str):
        """
        A simple 'grep' like method to search index files to find lines
        with 'exact' match of the search string. Returns a list the lines
        from the index file without further parsing.
        """
        found = []
        self.fai_handle.seek(0)
        for line in self.fai_handle:
            if search_str in line:
                found.append(line)
        return found
        
    def parse_fai_line(self, fai_line):
        """
        A simple method to parse a line from the index file into a tuple
        """
        parsed = fai_line.strip().split('\t')
        return tuple([self.faidx_types[i](parsed[i]) for i in range(5)])
    
    def get_fai_info(self, search_str):
        """
        A simple wrapper to parse all the index lines that are returned
        by 'grep_fai' with the supplied the search string. Returns a list of
        parsed index lines
        """
        fai_lines = self.grep_fai(search_str)
        fai_info = []
        for line in fai_lines:
            fai_info.append(self.parse_fai_line(line))
        return fai_info
        
    def get_seq(self, fai_info):
        """
        Method to retrieve sequence from the indexed fasta file for the given
        parsed index line. Does not accept a raw line from the index, it needs
        to be parsed by parse_fai_line first! Returns the sequence string
        """
        num_lines = int(math.ceil(fai_info[self.Length] / fai_info[self.Linebases]))
        self.fasta_handle.seek(fai_info[self.Offset])
        seq = ''
        line_rf = fai_info[self.Linewidth]-fai_info[self.Linebases]
        for i in range(num_lines):
            seq += self.fasta_handle.read(fai_info[self.Linebases])
            self.fasta_handle.seek(line_rf, 1)
        seq += self.fasta_handle.read(fai_info[self.Length] - num_lines * fai_info[self.Linebases])
        return seq
    
    def iter_all_seq(self, search_str):
        """
        A simple generator method to return tuples of GeneID, Seq for the
        search string.
        """
        fai_lines = self.grep_fai(search_str)
        for line in fai_lines:
            fai_info = self.parse_fai_line(line)
            seq = self.get_seq(fai_info)
            yield (fai_info[self.GeneID], seq)
        