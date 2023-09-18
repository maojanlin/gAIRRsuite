import argparse
import pickle
import os
import numpy as np
from parse_sam_haplotyping import parse_CIGAR, warp_SEQ

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fa', '--target_allele', 
        help = 'interested allele name to call from sam'
    )
    parser.add_argument(
        '-fs', '--fn_sam_file', 
        help = 'sam file including the target allele'
    )
    args = parser.parse_args()
    return args


def fetch_sam(fn_sam_file, target_allele):
    dict_operate = {'M': 0, 'I': 1, 'D': -1, 'S':1, 'H':1}
    
    with open(fn_sam_file, 'r') as f_r:
        match_flag = False
        for line in f_r:
            if line[0] != '@': # neglect the headers
                fields = line.split()
                allele_name = fields[0].split('|')[1]
                if target_allele == allele_name: # find the match
                    cigar = fields[5]
                    allele_SEQ = fields[9]
                    allele_SEQ = warp_SEQ(allele_SEQ, cigar)
                    print(fields[2] + '\t' + fields[3] + '-' + str(int(fields[3]) + len(allele_SEQ)) + '\t' + allele_SEQ )

if __name__ == '__main__':
    args = parse_args()
    fn_sam_file = args.fn_sam_file
    target_allele = args.target_allele

    fetch_sam(fn_sam_file, target_allele)

