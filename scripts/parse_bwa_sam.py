import argparse
import pickle
import os
import numpy as np
from parse_sam_haplotyping import parse_CIGAR

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fs', '--fn_sam',
        help = 'input sam file'
    )
    parser.add_argument(
        '-td', '--thrsd', type = int,
        help = 'edit distance threshold'
    )
    parser.add_argument(
        '-fo', '--fn_output_file',
        help = 'output report file'
    )
    args = parser.parse_args()
    return args

def parse_edit_distance(fn_sam, thrsd=0):
    with open(fn_sam, 'r') as f_o:
        for line in f_o:
            if line[0] != '@': # real alignment information
                fields = line.split()
                #print(fields[11])
                eDist = int(fields[11].split(':')[2])
                cigar = fields[5]
                if 'S' in cigar:
                    continue
                if fields[2] != '*' and eDist <= thrsd:
                    print_word = fields[0] + ' ' + fields[2].split('_')[5] + ' ' + fields[2]
                    #print_word = fields[0] + '\t' + fields[2] + '\t' + fields[11]
                    print(print_word)


if __name__ == '__main__':
    args = parse_args()
    fn_sam = args.fn_sam
    thrsd = args.thrsd 

    parse_edit_distance(fn_sam, thrsd);
    

