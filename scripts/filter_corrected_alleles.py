import argparse
import pickle
import os
import numpy as np
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fs', '--fn_sam',
        help = 'input sam file where alleles.fasta align to corrected_alleles.fasta'
    )
    parser.add_argument(
        '-fca', '--fn_corrected_alleles',
        help = 'input corrected allele file'
    )
    
    parser.add_argument(
        '-fof', '--fo_filtered_alleles',
        help = 'output filtered corrected allele file'
    )
    args = parser.parse_args()
    return args

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def parse_fasta(fn_fasta):
    '''parse the fasta file into a dictionary'''
    # dict_name_SEQ {}
    #  - keys: seq_name
    #  - values: seq_SEQ
    dict_name_SEQ = {}
    with open(fn_fasta, 'r') as f_f:
        seq_name = ""
        seq_SEQ = ""
        for line in f_f:
            if line[0] == '>':
                if seq_name != "":
                    if dict_name_SEQ.get(seq_name):
                        print("WARNING! Duplicate sequence name:", seq_name)
                    dict_name_SEQ[seq_name] = seq_SEQ
                seq_name = line.strip()[1:]
                seq_SEQ = ""
            else:
                seq_SEQ += line.strip()
        if dict_name_SEQ.get(seq_name):
            print("WARNING! Duplicate sequence name:", seq_name)
        dict_name_SEQ[seq_name] = seq_SEQ
    return dict_name_SEQ


def parse_perfect_sam(fn_sam):
    list_perfect_fields = []
    list_mismatch_fields = []
    with open(fn_sam, 'r') as f_s:
        for line in f_s:
            if line[0] != '@':
                if "NM:i:0" in line:
                    fields = line.strip().split()
                    list_perfect_fields.append(fields)
                else:
                    fields = line.strip().split()
                    list_mismatch_fields.append(fields)
    return list_perfect_fields, list_mismatch_fields


if __name__ == "__main__":
    args = parse_args()
    fn_sam = args.fn_sam
    fn_corrected_alleles = args.fn_corrected_alleles
    fo_filtered_alleles = args.fo_filtered_alleles

    dict_allele_SEQ = parse_fasta(fn_corrected_alleles)
    list_perfect_fields, list_mismatch_fields = parse_perfect_sam(fn_sam)
    set_duplicate_allele = set()
    for fields in list_perfect_fields:
        # get the ref (allele) info
        set_duplicate_allele.add(fields[2])

    for ele in set_duplicate_allele:
        dict_allele_SEQ.pop(ele)

    f_o = open(fo_filtered_alleles, 'w')
    for allele_name in sorted(dict_allele_SEQ.keys()):
        f_o.write(">" + allele_name + '\n')
        f_o.write(dict_allele_SEQ[allele_name] + '\n')
    f_o.close()


