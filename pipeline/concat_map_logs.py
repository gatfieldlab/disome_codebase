#!/usr/bin/env python3

import os
import re
import sys



reads = re.compile('^(\d+) reads')
align = re.compile('(\d+) .+ aligned (.+)$')
log_path = sys.argv[2]
star_log_path = sys.argv[3]

def get_log_files(log_path, base_name, extension="bw2.log"):
    file_regex = re.compile('^{}\.(.+)\.{}'.format(base_name, extension))
    file_list = os.listdir(log_path)
    log_files = []
    for file in file_list:
        file_match = file_regex.match(file)
        if file_match:
            log_type = file_match.group(1)
            log_files.append((log_type, os.path.join(log_path, file)))
    return log_files

def get_STAR_log(log_path, base_name):
    log_file = os.path.join(log_path, base_name, "Log.final.out")
    with open(log_file) as f:
        for line in f:
            try:
                key, value = tuple([val.strip() for val in line.strip().split("|")])
            except ValueError:
                continue
            if key == "Uniquely mapped reads number":
                uniq_map = int(value)
            elif key == "Number of reads mapped to multiple loci":
                multi_map = int(value)
            elif key == "Number of reads mapped to too many loci":
                super_map = int(value)
            elif key == "Number of input reads":
                all_reads = int(value)
        try:
            total_mapped = uniq_map + multi_map + super_map
            unmapped = all_reads - total_mapped
        except NameError as err:
            print("Error while parsing {}".format(log_file))
            print(err)
            sys.exit(1)
        return (all_reads, total_mapped, unmapped)


with open(sys.argv[1]) as samples:
    for sample in samples:
        sample = sample.strip()
        star_log = get_STAR_log(star_log_path, sample)
        log_files = get_log_files(log_path, sample)
        for log_type, log_file in log_files:
            align_stat = {}
            with open(log_file) as f:
                for line in f:
                    read_match = reads.match(line)
                    if read_match:
                        tot_reads = read_match.group(1)
                    align_match = align.search(line)
                    if align_match:
                        align_stat[align_match.group(2)] = align_match.group(1)
                align_stat['total_aligned'] = str(int(tot_reads) - int(align_stat['0 times']))
                keys = ['total_aligned']
                for t in keys:
                    print("\t".join([sample, log_type, t, align_stat[t]]))
                if log_type == 'mouse_cDNA':
                    #print("\t".join([sample, "unmapped", "total_aligned", align_stat['0 times']]))
                    assert int(align_stat['0 times']) == star_log[0]
        print("\t".join([sample, "genomic", "total_aligned", str(star_log[1])]))
        print("\t".join([sample, "unmapped", "total_aligned", str(star_log[2])]))

        #     except IOError as e:
        #         sys.stderr.write("Error: problem accessing data file 'mapping_data/logs/{0}.{1}.bw2.log'\n".format(file_base, log_type))
        #         sys.stderr.write(e.message)
        #             #sys.exit(1)
        # if sample_base in tophat_processed:
        #     continue
        # tophat_processed.append(sample_base)
        # with open("mapping_data/tophat/"+sample_base+"/align_summary.txt") as tophat_alns:
        #     tophat_align = {}
        #     for line in tophat_alns:
        #         align_match = tophat.search(line)
        #         if align_match:
        #             tophat_align[align_match.group(1)] = int(align_match.group(2))
        # try:
        #     align_stat['genomic'] = {'0 times': tophat_align['Input'] - tophat_align['Mapped'],
        #                              'exactly 1 time': tophat_align['Mapped'] - tophat_align['of these'],
        #                              '>1 times': tophat_align['of these']}
        # except:
        #     print "Error: could not extract tophat stats properly for: {}".format(sample_base)
        #     print "Unexpected error:", sys.exc_info()[0]
        #     sys.exit(1)
        # for t in ('0 times', 'exactly 1 time', '>1 times'):
        #     print "\t".join([sample_name, 'all', 'genomic', t, str(align_stat['genomic'][t])])
