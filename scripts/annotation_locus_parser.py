import argparse
import pickle
import os

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fa', '--fn_annotated',
        help = 'input annotation file'
    )
    parser.add_argument(
        '-fo', '--fn_output',
        help = 'output file with the contig_name and interval information'
    )
    args = parser.parse_args()
    return args


def append_locus(fn_annotated, fn_output):
    f_o = open(fn_output, 'a')
    contig_name = fn_annotated.split('/')[-1]
    contig_name = contig_name.split('.')[0]
    with open(fn_annotated) as f_a:
        for line in f_a:
            if line[0:4] != "Name":
                f_o.write(contig_name + ',')
                f_o.write(line)
    f_o.close()


if __name__ == '__main__':
    args = parse_args()
    fn_annotated = args.fn_annotated
    fn_output = args.fn_output

    append_locus(fn_annotated, fn_output)
