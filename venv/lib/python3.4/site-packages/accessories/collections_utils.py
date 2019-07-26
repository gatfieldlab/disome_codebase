# -*- coding: utf-8 -*-
"""
Utilities for collections items

@author: barpat
"""
import collections

def fix_depth_recursive_defaultdict(depth=3, t=int):
    """Generates a recursive defaultdict object with a given depth. The
    innermost values of the object will the of type 't'"""
    if depth==1:
        return collections.defaultdict(t)
    else:
        return collections.defaultdict(
            lambda: fix_depth_recursive_defaultdict(depth-1, t))


recursive_defaultdict = lambda: collections.defaultdict(recursive_defaultdict)
recursive_defaultdict.__doc__ = "Generates an endless recursive defaultdict object"

def int_factory():
    return collections.defaultdict(int)

def dict_factory():
    return collections.defaultdict(dict)
