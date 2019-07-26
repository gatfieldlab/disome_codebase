# -*- coding: utf-8 -*-
"""
acessories/align_utils

This module hosts several utility functions
to help processing sequence alignmnets,
particularly from a SAM, i.e CIGAR and MD-tags.

Created on Thu Jul 18 15:07:42 2013
@author: Bulak Arpat, bulak.arpat@unil.ch

"""

from re import compile

# Compiled regex
MD_REGEX = compile('([0-9]+|[A-Z]+|\^[A-Z]+)')
CIGAR_REGEX = compile('[0-9]+[MIDNSHP=X]')

# Function to replace/insert str in str
def str_mod(string, pos, insert, replace=False):
    """
    This function facilitates insertion/replacement
    of a string into another one at a certain position
    """
    next_pos = pos + len(insert) * replace
    return string[:pos] + insert + string[next_pos:]

# Function to parse a CIGAR
def parse_cigar(cigar):
    return CIGAR_REGEX.findall(cigar)

# Function to convert CIGAR to a list of operations
def cigar_to_list(cigar):
    cigar_ops = []
    for op in parse_cigar(cigar):
        cigar_ops.append((op[-1], int(op[:-1])))
    return cigar_ops

# Function to calculate the length of SEQ based on cigar
def seq_len_from_cigar(cigar):
    parsed_cigar = parse_cigar(cigar)
    seq_len = 0
    for op in parsed_cigar:
        if op[-1] in ('M','I','S','=', 'X'):
            seq_len += int(op[:-1])
    return seq_len

def aln_len_from_cigar(cigar):
    aln_len = 0
    for op in parse_cigar(cigar):
        if op[-1] in ('M', 'D', 'N', 'P'):
            aln_len += int(op[:-1])
    return aln_len

# Function to parse an MD tag (not fully tested)
def parse_md(md_tag):
    return MD_REGEX.findall(md_tag)

# Function to confirm a seq, CIGAR, MD tag triple
def assert_aln(read_seq, cigar, md_tag):
    is_ok = True
    parsed_md = parse_md(md_tag)
    parsed_cigar = parse_cigar(cigar)
    md_len = 0
    md_dels = []
    for op in parsed_md:
        if op.isdigit():
            md_len += int(op)
        elif op[0] == '^':
            md_dels.append((md_len, len(op)-1))
        else:
            md_len += len(op)
    cigar_len = 0
    cigar_len2 = 0
    mod_cigar=[]
    for op in parsed_cigar:
        if op[-1] == 'D':
            cigar_del = (cigar_len2, int(op[:-1]))
            if not cigar_del in md_dels:
                is_ok = False
            else:
                md_dels.remove(cigar_del)
        elif op[-1] in ('I','S'):
            cigar_len += int(op[:-1])
        elif op[-1] not in ('H'):
            cigar_len += int(op[:-1])
            cigar_len2 += int(op[:-1])
        if op[-1] in ('I'):
            mod_cigar.append((cigar_len2, op[-1], int(op[:-1])))
    if md_dels:
        is_ok = False
    if not (cigar_len2 == md_len and cigar_len == len(read_seq)):
        is_ok = False
    if is_ok:
        return (parsed_cigar, parsed_md, mod_cigar)
    else:
        return (None, None, None)

def edit_dist_range(cigar, md_tag):
    md_by_cigar = []
    md_by_pos  = []
    parsed_cigar = parse_cigar(cigar)
    parsed_md = parse_md(md_tag)
    cur_edit = 0
    for op in parsed_md:
        if op.isdigit():
            md_by_pos += [cur_edit] * int(op)
        elif op[0] == '^':
            md_by_pos += range(cur_edit+1, cur_edit+len(op))
            cur_edit = md_by_pos[-1]
        else:
            md_by_pos += range(cur_edit+1, cur_edit+len(op)+1)
            cur_edit = md_by_pos[-1]
    cur_edit = 0
    for op in parsed_cigar:
        if op[-1] in ('S', 'H', 'P'):
            pass
        elif op[-1] in ('M', 'D'):
            md_by_cigar += [cur_edit] * int(op[:-1])
        elif op[-1] == 'I':
            cur_edit += int(op[:-1])
        elif op[-1] == 'N':
            md_by_pos = md_by_pos[:len(md_by_cigar)] + [md_by_pos[len(md_by_cigar)-1]] * int(op[:-1]) + md_by_pos[len(md_by_cigar):]
            md_by_cigar += [cur_edit] * int(op[:-1])
    if not len(md_by_cigar) ==  len(md_by_pos):
        raise("Length of edit distance arrays by cigar and md do not match")
    edit_range = [md_by_cigar[i] + md_by_pos[i] for i in range(len(md_by_cigar))]
    return edit_range


# Function to convert a CIGAR and MD tag to a gapped alignment
def sam2gapped(read_seq, cigar, md_tag):
    """
    This function takes the read sequence, CIGAR and MD tag
    as inputs and produces a gapped alignment of the read sequence
    against the reference sequence without the need of reference
    sequence as an input:

    ref     AGTGCCTTGGGTGTTCA-----ATCCCCATGCAACAACC
    aln     || ||||    ||||||     | |||||||||||||||
    read    AGAGCCTC---TGTTCACATAGACCCCCATGCAACAACC

    """
    # Assert and parse the CIGAR and MD tag
    cigar_ops, md_ops, mod_cigar = assert_aln(read_seq, cigar, md_tag)
    if not cigar_ops:
        raise Exception('Problem: {},{},{}'.format(read_seq, cigar, md_tag))

    # Initialize the ref seq and the pairwise aligment
    ref_seq = read_seq
    pairwise = '|' * len(read_seq)

    # First loop over the MD tag
    # Skip soft clipped bps (S) at the 5'
    soft_pos = 0
    for i in (0,min(1, len(cigar_ops)-1)):
        if cigar_ops[i][-1] not in ('H','S'):
            break
        elif cigar_ops[i][-1] == 'S':
            soft_pos += int(cigar_ops[i][:-1])
    cur_pos = soft_pos
    for md_op in md_ops:
        #print cur_pos, md_op
        while mod_cigar and mod_cigar[0][0] <= cur_pos-soft_pos:
            insertion = mod_cigar.pop(0)
            cur_pos += insertion[2]
        # a match
        if md_op.isdigit():
            op_len = int(md_op)
        # a deletion
        elif md_op[0] == '^':
            op_len = len(md_op)-1
            read_seq = str_mod(read_seq, cur_pos, '-'*op_len)
            ref_seq = str_mod(ref_seq, cur_pos, md_op[1:])
            pairwise = str_mod(pairwise, cur_pos, ' '*op_len)
        # a mismatch
        else:
            op_len = len(md_op)
            ref_seq = str_mod(ref_seq, cur_pos, md_op, replace=True)
            pairwise = str_mod(pairwise, cur_pos, '.'*op_len, replace=True)
        cur_pos += op_len

    # Then loop over CIGAR
    cur_pos = 0
    for cigar_op in cigar_ops:
        op_len = int(cigar_op[:-1])
        op_com = cigar_op[-1]
        if  op_com in ('I','S'):
            ref_seq = str_mod(ref_seq, cur_pos, '-'*op_len, replace=True)
            pairwise = str_mod(pairwise, cur_pos, ' '*op_len, replace=True)
        if not op_com in ('N','H','P'):
            cur_pos += op_len
    return '\n'.join([ref_seq, pairwise, read_seq])
