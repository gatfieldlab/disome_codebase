#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Simple methods to pretty and functional printing various data types"""


__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2017, Bulak Arpat"
__license__ = "GPLv3"
__version__ = "0.1.0"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"


def pretty_list(lst, max_len=6):
    """Prettify the string representation of a long list by presenting few
    values from the beginning and ending and putting ... in the middle without
    commas:
    [1, 2, 3, 4, 5, 6, 7] ->  '[ 1 2 3 ... 5 6 7]'

    Args:
        l: :obj:`list`
        max_len: :obj:`int` defining the max length of list after which it will
        be prettified
    Returns: :obj:`str`
    """
    if not isinstance(lst, list):
        raise TypeError('pretty_list arg 1 has to be a list')
    if not isinstance(max_len, int):
        raise TypeError('pretty_list arg 2 has to be an integer')
    if not max_len > 1:
        raise ValueError('pretty_list arg 2 has to be > 1')
    if len(lst) > max_len:
        right = max_len // 2
        left = max_len - right
        lst = lst[:left] + ['...'] + lst[-right:]
    return "[ " + " ".join([str(item) for item in lst]) + "]"


def pretty_dict(dic):
    """Prettify the string representation of a dictionary:
    {'c': 39, 'ded': 2, 'aa': 21} -> '{ aa: 21, c: 39, ded: 2,}'

    Arg: :obj:`dict` to be prettified
    """
    if not isinstance(dic, dict):
        raise TypeError('pretty_dict arg has to a dictionary')
    d_str = "{"
    for k in sorted(dic.keys()):
        d_str += " {}: {},".format(k, dic[k])
    d_str += "}"
    return d_str
