#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""This script aims at annotating EnsemblID based text outputs by searching
and appending gene and transcript names from an annotation file.
It is very specific to the annotation file and not intended for general use
"""

import sys
import re
import argparse
import fileinput

__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2016, 2017, Bulak Arpat"
__license__ = "GPLv3"
__version__ = "0.1.0"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"


ID_REGEX = re.compile(r'(ENSMUS[GT]\d{11})')


def get_argparser():
    """Generates the argument parser"""
    parser = argparse.ArgumentParser(prog="pretty_names.py",
                                     description=__version__)
    parser.add_argument(
        "-a", "--annotation", type=str, required=True, help="Annotation file")
    parser.add_argument(
        "-d", "--delimiter", type=str, default="\t", help="delimiter in input")
    parser.add_argument(
        "-m", "--mode", type=str, choices=["append", "replace"],
        default="append", help='mode of action for new annotations')
    parser.add_argument(
        "input", nargs="*", help='input file(s)', default="-")
    return parser


def main():
    """Main function that reads the annotation file and appends the user-
    friendly annotations to the IDs used in the input files"""
    parser = get_argparser()
    args = parser.parse_args()
    id_annotation = {}
    try:
        with open(args.annotation) as annot_file:
            for line in annot_file:
                parsed = line.strip().split()
                id_annotation[parsed[0]] = parsed[1]
                id_annotation[parsed[3]] = parsed[4]
    except OSError as e:
        sys.stderr.write(
            "Error '{}' occured tying to read from '{}': {}\n".format(
                e.errno, args.annotation, e.strerror))
        return 1
    except IndexError:
        sys.stderr.write(
            "Could not properly parse the annotation file: {}\n".format(
                args.annotation))
        return 1

    for line in fileinput.input(args.input):
        ids_found = ID_REGEX.findall(line)
        if ids_found:
            unique_ids = list(set(ids_found))
            unique_ids.sort(key=lambda x: ids_found.index(x))
            annotations = [id_annotation[iden] for iden in unique_ids]
            id_col = args.delimiter.join(annotations) + args.delimiter
        else:
            id_col = ""
        sys.stdout.write("{}{}".format(id_col, line))


if __name__ == "__main__":
    sys.exit(main())
