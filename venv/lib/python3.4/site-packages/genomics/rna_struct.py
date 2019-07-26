#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This module provides tables for RNA secondary
structures and tools
"""

# import itertools

__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2017, Bulak Arpat"
__license__ = "GPLv3"
__version__ = "0.2.0"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"


# G-quadruplex translation table
GQUAD_LETTERS = ["N", "G"]
# GQUAD_TABLE = {codon: GQUAD_LETTERS["G" in codon] for codon in
#                ["".join(word) for word in itertools.product(
#                   GQUAD_LETTERS, repeat=3)]}

# NGN is intentially undefined. We should never come across this key
# if we do it must be an error
GQUAD_TABLE = {
    "GGG": "G",
    "GGN": "G",
    "GNN": "G",
    "NNN": "N",
    "NNG": "G",
    "NGG": "G",
    "GNG": "G"
}
