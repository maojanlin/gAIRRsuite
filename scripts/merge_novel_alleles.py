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
        help = 'input original allele fasta file'
    )
    parser.add_argument(
        '-fn', '--fn_novel',
        help = 'input fasta file with novel alleles'
    )
    
    parser.add_argument(
        '-fom', '--fo_merged_fasta',
        help = 'output merged fasta file'
    )
    args = parser.parse_args()
    return args

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)



if __name__ == '__main__':
    args = parse_args()
    fn_alleles = args.fn_alleles
    fn_novel   = args.fn_novel
    fo_merged_fasta = args.fo_merged_fasta

    dict_original_allele = parse_fasta(fn_alleles)
    dict_novel_allele    = parse_fasta(fn_novel)

    f_om = open(fo_merged_fasta, 'w')
    for name in sorted(dict_original_allele.keys()):
        f_om.write(">" + name.split()[0] + "\n")
        f_om.write(dict_original_allele[name] + "\n")
    for name in sorted(dict_novel_allele.keys()):
        if ("extend" in name) == False:
            f_om.write(">|" + name + "|\n")
            f_om.write(dict_novel_allele[name] + "\n")
    f_om.close()

