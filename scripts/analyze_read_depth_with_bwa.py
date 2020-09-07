import argparse
import pickle
import os
import numpy as np
from utils import get_reverse_complement
from filter_corrected_alleles import parse_perfect_sam, parse_fasta
from parse_contig_realign import parse_CIGAR
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fa', '--fn_alleles',
        help = 'input target allele fasta file'
    )
    parser.add_argument(
        '-fs', '--fn_sam',
        help = 'input short reads to allele sam file'
    )
    parser.add_argument(
        '-t', '--thrsd', type=int,
        help = 'mininum coverage requirement for reads in read depth analysis'
    )
    
    parser.add_argument(
        '-foc', '--fo_calling_report',
        help = 'output calling report file'
    )
    parser.add_argument(
        '-fop', '--fo_grouping_pickle',
        help = 'output allele-read group pickle file'
    )
    parser.add_argument(
        '-fv', '--fn_verify_annotation',
        help = 'input annotation file for verification'
    )
    args = parser.parse_args()
    return args

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def parse_perfect_sam_with_S(fn_sam):
    list_perfect_fields = []
    list_mismatch_fields = []
    with open(fn_sam, 'r') as f_s:
        for line in f_s:
            if line[0] != '@':
                fields = line.strip().split()
                if "NM:i:0" in line:
                    list_perfect_fields.append(fields)
                else:
                    list_mismatch_fields.append(fields)
    return list_perfect_fields, list_mismatch_fields


def histogram_read_depth(dict_allele_histogram, list_perfect_fields, dict_allele_reads, thrsd=100):
    for fields in list_perfect_fields:
        allele_name = fields[2]
        read_name = fields[0]
        start_pos = int(fields[3])
        cigar = fields[5]
        ref_end_pos = len(dict_allele_histogram[allele_name]) + 1
        min_coverage_thrsd = min(thrsd, ref_end_pos-1)
        number, operate = parse_CIGAR(cigar)
        if len(operate) > 3:
            eprint("Error on CIGAR!", operate)
            continue
        M_idx = None
        for idx, ele in enumerate(operate):
            if ele == 'M':
                M_idx = idx
        if number[M_idx] < min_coverage_thrsd:
            continue
        end_pos = number[M_idx] + start_pos
        if start_pos == 1: # read map to the left side
            if operate[-1] == 'M' or end_pos == ref_end_pos: # right side is open or till the end
                dict_allele_histogram[allele_name][start_pos-1:end_pos-1] += 1
        elif end_pos == ref_end_pos: # read map to the left side
            if operate[0] == 'M': # left side is open
                dict_allele_histogram[allele_name][start_pos-1:end_pos-1] += 1
        elif operate == ['M']: # read is including in the allele
            dict_allele_histogram[allele_name][start_pos-1:end_pos-1] += 1

        dict_allele_reads[allele_name].add(read_name)



if __name__ == '__main__':
    args = parse_args()
    fn_alleles = args.fn_alleles
    fn_sam = args.fn_sam
    thrsd = args.thrsd
    fo_calling_report  = args.fo_calling_report
    fo_grouping_pickle = args.fo_grouping_pickle
    fn_verify_annotation = args.fn_verify_annotation

    dict_allele = parse_fasta(fn_alleles)
    dict_allele_histogram = { name.split()[0]:np.zeros(len(SEQ)) for name, SEQ in dict_allele.items() }
    dict_allele_reads = { name.split()[0]:set() for name in dict_allele.keys() }

    list_perfect_fields, list_mismatch_fields = parse_perfect_sam_with_S(fn_sam)
    histogram_read_depth(dict_allele_histogram, list_perfect_fields, dict_allele_reads, thrsd)

    list_min_depth = [ [min(histogram), None, name.split('|')[1]] for name, histogram in  dict_allele_histogram.items() ]

    set_annotation = set()
    if fn_verify_annotation:
        with open(fn_verify_annotation, "r") as f_v:
            for line in f_v:
                set_annotation.add(line.split(',')[0])

    for idx, ele in enumerate(list_min_depth):
        min_depth, tag, name = ele
        if "corrected" in name:
            list_min_depth[idx][1] = 4
            continue
        if name in set_annotation:
            if min_depth > 0:
                list_min_depth[idx][1] = 3
            else:
                list_min_depth[idx][1] = 1
        else:
            if min_depth > 0:
                list_min_depth[idx][1] = 2
            else:
                list_min_depth[idx][1] = 0

    f_oc = open(fo_calling_report,'w')
    for element in sorted(list_min_depth, reverse=True):
        min_depth, tag, name = element
        if fn_verify_annotation:
            if tag==0:
                continue
            f_oc.write(name + '\t' + str(min_depth) + '\t' + str(tag))
        else:
            f_oc.write(name + '\t' + str(min_depth))
        f_oc.write('\n')
    f_oc.close()

    f_op = open(fo_grouping_pickle, 'wb')
    pickle.dump(dict_allele_reads, f_op)

    

