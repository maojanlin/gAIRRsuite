'''
This script reads a set of genomic sequences and generate the 
position weight matrix (PWM) of the read set.
'''
import argparse
import numpy as np
from utils import read_from_file

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', '--fn_input',
        help='input file, can be .txt or .fasta'
    )
    parser.add_argument(
        '-o', '--fn_out',
        help='output file (names)'
    )
    args = parser.parse_args()
    return args


def build_pwm(list_seqs):
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    mat_count = np.zeros((4, len(list_seqs[0])))
    for seq in list_seqs:
        for i, nuc in enumerate(seq):
            mat_count[mapping[nuc], i] += 1
    print (mat_count)
    print (mat_count / len(list_seqs))
    return

if __name__ == '__main__':
    args = parse_args()
    fn_input = args.fn_input

    list_seqs = read_from_file(fn_input)
    pwm = build_pwm(list_seqs)