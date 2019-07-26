#!/bin/env python3

import sys

with open(sys.argv[1]) as infile:
    if sys.argv[3] == "5":
        mod_pos = 10
    elif sys.argv[3] == "3":
        mod_pos = 9
    reg_len = int(sys.argv[2]) * 3
    min_len = reg_len * 2
    for line in infile:
        parsed = line.strip().split("\t")
        if int(parsed[8]) < min_len:
            continue
        parsed[mod_pos] = str(int(parsed[9]) + reg_len)
        sys.stdout.write("{}\n".format("\t".join(parsed)))
