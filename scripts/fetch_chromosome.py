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
        '-fo', '--fn_output_fasta_file',
        default = "chromosome.fasta",
        help = 'output target fasta file'
    )
    args = parser.parse_args()
    return args


def fetch_chromosome_reads(chromosome_name, fn_genome, fn_output_fasta_file):
    f_o = open(fn_output_fasta_file, 'w')

    with open(fn_genome, 'r') as f_r:
        match_flag = False
        for line in f_r:
            if line[0] == '>':
                if (line[1:-1]) == chromosome_name:
                    f_o.write(line)
                    
                    read_name = line
                    read_SEQs = ""
                    match_flag = True
                elif match_flag:
                    break
            elif match_flag == True:
                f_o.write(line)

    f_o.close()

if __name__ == '__main__':
    args = parse_args()
    fn_genome = args.fn_genome
    chromosome_name = args.chromosome_name
    fn_output_fasta_file = args.fn_output_fasta_file

    fetch_chromosome_reads(chromosome_name, fn_genome, fn_output_fasta_file)

