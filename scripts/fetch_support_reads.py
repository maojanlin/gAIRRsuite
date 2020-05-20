import argparse
import pickle
import os
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fr', '--fn_read',
        help = 'reads fasta file'
    )
    parser.add_argument(
        '-fsup', '--fn_support_read_pickle', 
        help = 'input pickle file indicating support reads of alleles'
    )
    parser.add_argument(
        '-fa', '--target_allele', 
        help = 'interested allele name to call support reads'
    )
    parser.add_argument(
        '-fo', '--fn_output_fasta_file',
        default = "support_reads.fasta",
        help = 'output fasta file composed of support reads'
    )
    args = parser.parse_args()
    return args


def read_pickle(pickle_file):
    f = open(pickle_file, 'rb')
    x = pickle.load(f)
    return x


def fetch_support_reads(target_allele, fn_read, dict_support_read, fn_support_read_pickle):
    f_o = open(fn_output_fasta_file, 'w')
    total_sup_reads = 0

    set_target_allele = dict_support_read[target_allele]
    with open(fn_read, 'r') as f_r:
        read_name = ""
        read_SEQs = ""
        match_flag = False
        for line in f_r:
            if line[0] == '>':
                if (line[1:-1]) in dict_support_read[target_allele]:
                    f_o.write(read_name)
                    f_o.write(read_SEQs)
                    
                    read_name = line
                    read_SEQs = ""
                    match_flag = True
                    total_sup_reads += 1
                else:
                    match_flag = False
            elif match_flag == True:
                read_SEQs = read_SEQs + line
        f_o.write(read_name)
        f_o.write(read_SEQs)

    f_o.close()
    print("Total " + str(total_sup_reads) + " support reads!")
    return total_sup_reads



if __name__ == '__main__':
    args = parse_args()
    fn_read = args.fn_read
    fn_support_read_pickle = args.fn_support_read_pickle
    target_allele = args.target_allele
    fn_output_fasta_file = args.fn_output_fasta_file

    dict_support_read = read_pickle(fn_support_read_pickle)
    support_reads = fetch_support_reads(target_allele, fn_read, dict_support_read, fn_support_read_pickle)

