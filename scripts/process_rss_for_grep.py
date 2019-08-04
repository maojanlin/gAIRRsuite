import argparse
from build_pwm import read_from_file

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', '--fn_input',
        help='input file, can be .txt or .fasta'
    )
    args = parser.parse_args()
    return args

def get_reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement

def process_seqs_for_grep(list_seqs):
    list_rc = []
    for seq in list_seqs:
        list_rc.append(get_reverse_complement(seq))
    set_all = set(list_seqs).union(set(list_rc))
    return set_all

if __name__ == '__main__':
    args = parse_args()
    fn_input = args.fn_input

    ext_input = fn_input[fn_input.rfind('.') + 1:]
    list_seqs = read_from_file(fn_input, ext_input)

    print ('Length of list {0} = {1}'.format(fn_input, len(list_seqs)))
    set_all = process_seqs_for_grep(list_seqs)
    print ('Size of set for forward/reverse seqs = {}'.format(len(set_all)))
    print ('\|'.join(set_all))