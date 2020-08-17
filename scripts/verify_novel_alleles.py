import argparse
import pickle
import os
import numpy as np
from filter_corrected_alleles import parse_perfect_sam
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fs1', '--fn_sam_H1',
        help = 'input sam file where alleles.fasta align to corrected_alleles.fasta'
    )
    parser.add_argument(
        '-fs2', '--fn_sam_H2',
        help = 'input sam file where alleles.fasta align to corrected_alleles.fasta'
    )
    
    parser.add_argument(
        '-for', '--fo_verification_report',
        help = 'output filtered corrected allele file'
    )
    args = parser.parse_args()
    return args

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)



if __name__ == "__main__":
    args = parse_args()
    fn_sam_H1 = args.fn_sam_H1
    fn_sam_H2 = args.fn_sam_H2
    fo_verification_report = args.fo_verification_report

    list_perfect_fields_1, list_mismatch_fields_1 = parse_perfect_sam(fn_sam_H1)
    list_perfect_fields_2, list_mismatch_fields_2 = parse_perfect_sam(fn_sam_H2)
    list_perfect_fields_1.sort(key=lambda fields:(fields[2], fields[3]), reverse=True)
    list_perfect_fields_2.sort(key=lambda fields:(fields[2], fields[3]), reverse=True)
    print("There are", len(list_perfect_fields_1), "true novel alleles in H1:")
    for fields in list_perfect_fields_1:
        print("\t", fields[0] + '\t\t' + fields[2] + '\t' + fields[3])
    print("There are", len(list_perfect_fields_2), "true novel alleles in H1:")
    for fields in list_perfect_fields_2:
        print("\t", fields[0] + '\t\t' + fields[2] + '\t' + fields[3])

    set_mismatch_fields_1 = {fields[0] for fields in list_mismatch_fields_1}
    set_mismatch_fields_2 = {fields[0] for fields in list_mismatch_fields_2}
    print("Useless corrected alleles are:\t", sorted(set_mismatch_fields_1.intersection(set_mismatch_fields_2)))


