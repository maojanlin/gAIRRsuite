'''
Process fasta files for TCR/BCR alleles, and return a dictionary
that stores the length of the sequences.
'''

import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', '--fn_input', action = 'append', nargs = '+',
        help = 'FASTA file we want to include'
    )
    parser.add_argument(
        '-o', '--fn_output',
        help = 'output file including only alleles [None]'
    )
    args = parser.parse_args()
    return args

def read_fasta(fn, dc):
    num_before = len(dc)
    with open(fn, 'r') as f:
        seq = ''
        for line in f:
            if line[0] == '>':
                if len(seq) > 0:
                    dc[name] = len(seq)
                    seq = ''
                name = line.split('|')[1]
            else:
                seq += line.rstrip()
        if len(seq) > 0:
            dc[name] = len(seq)
    num_after = len(dc)
    sys.stderr.write('{} seqs added from {}\n'.format(num_after-num_before, fn))
    return dc

if __name__ == '__main__':
    args = parse_args()

    sys.stderr.write(str(args.fn_input) + '\n')

    dc = {}
    for fn in args.fn_input:
        read_fasta(fn[0], dc)
    
    print ('ALLELE\tLENGTH')
    for k in dc.keys():
        print ('{}\t{}'.format(k, dc[k]))
