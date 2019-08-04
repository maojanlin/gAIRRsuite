'''
This script reads a set of genomic sequences and generate the 
position weight matrix (PWM) of the read set.
'''
import argparse
import numpy as np

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

def read_from_file(fn_input, ext):
    list_seqs = []
    assert ext in ['txt', 'fasta', 'fa']
    with open(fn_input, 'r') as f:
        for line in f:
            if ext in ['fasta', 'fa']:
                #: name
                if line[0] == '>':
                    continue
            line = line.rstrip().upper()
            
            #: doesn't allow nuc other than ACGT
            count_A = line.count('A')
            count_C = line.count('C')
            count_G = line.count('G')
            count_T = line.count('T')
            assert count_A + count_C + count_G + count_T == len(line)
            
            list_seqs.append(line)

    #: check if all the seqs are equal in length
    for s in list_seqs:
        assert len(s) == len(list_seqs[0])
    
    return list_seqs

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

    ext_input = fn_input[fn_input.rfind('.') + 1:]
    list_seqs = read_from_file(fn_input, ext_input)
    pwm = build_pwm(list_seqs)