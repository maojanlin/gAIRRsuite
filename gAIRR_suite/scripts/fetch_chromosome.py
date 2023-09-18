import argparse
import pickle
import os
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fg', '--fn_genome',
        help = 'genome fasta file'
    )
    parser.add_argument(
        '-fc', '--chromosome_name', 
        help = 'target chromosome_name'
    )
    parser.add_argument(
        '--interval',
        help = 'interested interval of the chromosome'
    )
    parser.add_argument(
        '-fo', '--fn_output_fasta_file',
        default = "chromosome.fasta",
        help = 'output target fasta file'
    )
    args = parser.parse_args()
    return args


def fetch_chromosome_reads(chromosome_name, fn_genome, fn_output_fasta_file):
    f_o = open(fn_output_fasta_file, 'w')
    chromosome_SEQ = ""

    with open(fn_genome, 'r') as f_r:
        match_flag = False
        for line in f_r:
            if line[0] == '>':
                if (line[1:-1]) == chromosome_name:
                    f_o.write(line)
                    match_flag = True
                elif match_flag:
                    break
            elif match_flag == True:
                f_o.write(line)
                chromosome_SEQ += line.strip()

    f_o.close()
    return chromosome_SEQ

if __name__ == '__main__':
    args = parse_args()
    fn_genome = args.fn_genome
    chromosome_name = args.chromosome_name
    fn_output_fasta_file = args.fn_output_fasta_file
    interval = args.interval

    chromosome_SEQ = ""
    # load chromosome file if it's already there
    if os.path.exists(fn_output_fasta_file):
        print ('Dumped chromosome file {} has existed, load for it instead of dumping'.format(fn_output_fasta_file))
        with open(fn_output_fasta_file, 'r') as f_r:
            for line in f_r:
                if line[0] != '>':
                    chromosome_SEQ += line.strip()
    else:
        chromosome_SEQ = fetch_chromosome_reads(chromosome_name, fn_genome, fn_output_fasta_file)

    if args.interval:
        positions = interval.split('-')
        print(chromosome_SEQ[int(positions[0])-1:int(positions[1])-1])


