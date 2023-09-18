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
    f_or = open(fo_verification_report, 'w')
    print("There are ", len(list_perfect_fields_1), "true novel alleles in H1:")
    f_or.write("There are" + str(len(list_perfect_fields_1)) + " true novel alleles in H1:\n")
    for fields in list_perfect_fields_1:
        print("\t", fields[0] + '\t\t' + fields[2] + '\t' + fields[3])
        f_or.write("\t" + fields[0] + '\t\t' + fields[2] + '\t' + fields[3] + '\n')
    print("There are ", len(list_perfect_fields_2), "true novel alleles in H2:")
    f_or.write("There are" + str(len(list_perfect_fields_2)) + " true novel alleles in H2:\n")
    for fields in list_perfect_fields_2:
        print("\t", fields[0] + '\t\t' + fields[2] + '\t' + fields[3])
        f_or.write("\t" + fields[0] + '\t\t' + fields[2] + '\t' + fields[3] + '\n')

    dict_mismatch_fields_1 = {}
    for fields in list_mismatch_fields_1:
        try:
            dict_mismatch_fields_1[fields[0]] = (fields[5], fields[11])
        except:
            dict_mismatch_fields_1[fields[0]] = ('*', '*')
    dict_mismatch_fields_2 = {}
    for fields in list_mismatch_fields_2:
        try:
            dict_mismatch_fields_2[fields[0]] = (fields[5], fields[11])
        except:
            dict_mismatch_fields_2[fields[0]] = ('*', '*')
    
    set_mismatch_fields_1 = set(dict_mismatch_fields_1.keys())
    set_mismatch_fields_2 = set(dict_mismatch_fields_2.keys())
    set_common_mismatch_fields = set_mismatch_fields_1.intersection(set_mismatch_fields_2)
    
    if len(set_common_mismatch_fields) == 0:
        print("No common mismatched alleles between H1 and H2.")
        f_or.write("No common mismatched alleles between H1 and H2.")
    else:
        print("Common mismatched alleles / CIGAR in H1 / NM:i in H1 / CIGAR in H2 / NM:i in H2")
        f_or.write("Common mismatched alleles / CIGAR in H1 / NM:i in H1 / CIGAR in H2 / NM:i in H2\n")
        for allele_name in sorted(set_common_mismatch_fields):
            pair_H1 = dict_mismatch_fields_1[allele_name]
            pair_H2 = dict_mismatch_fields_2[allele_name]
            print("\t", allele_name, ' / ', pair_H1[0], ' / ', pair_H1[1], ' / ', pair_H2[0], ' / ', pair_H2[1])
            f_or.write("\t" + allele_name + ' / ' + pair_H1[0] + ' / ' + pair_H1[1] + ' / ' + pair_H2[0] + ' / ' + pair_H2[1] + '\n')
    f_or.close()


