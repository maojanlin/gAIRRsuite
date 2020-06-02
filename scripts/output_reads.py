import argparse
import pickle
import os
import numpy as np
from utils import get_hamming_dist
from utils import get_reverse_complement
import sys
from coverage_analysis import fetch_reads_alleles, add_reads, add_alleles


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fr', '--fn_read',
        help = 'reads fasta file'
    )
    parser.add_argument(
        '-fa', '--fn_allele',
        help = 'allele fasta file'
    )
    parser.add_argument(
        '-fnp', '--fn_pickle_file', 
        help = 'input pickle file indicating reads and alleles in the clusters'
    )
    parser.add_argument(
        '-fo', '--fn_output_read_fasta',
        help = 'output reads in a cluster [None]'
    )
    args = parser.parse_args()
    return args

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def read_pickle(pickle_file):
    f = open(pickle_file, 'rb')
    x = pickle.load(f)

def add_diploid_reads(dict_cluster_info):
    dict_diploid = {'1':'_2', '2':'_1'}
    dict_cluster_diploid_info = dict(dict_cluster_info)
    for cluster_id in dict_cluster_info.keys():
        set_reads_name = dict_cluster_info[cluster_id][1]
        set_diploid_reads = set()
        for name in set_reads_name:
            read_name = name.split('_')[0]
            diploidity = name.split('_')[1]
            diploid_name = read_name + dict_diploid[diploidity[0]] + diploidity[1:]
            set_diploid_reads.add(diploid_name)
        dict_cluster_diploid_info[cluster_id][1] |= set_diploid_reads
    return dict_cluster_diploid_info


def output_reads(dict_cluster_info, fn_read, fn_allele, fn_output_read_fasta):
    dict_cluster_diploid_info = add_diploid_reads(dict_cluster_info)
    dict_read_allele_clusters = fetch_reads_alleles(fn_read, fn_allele, dict_cluster_diploid_info)
    
    for cluster_id in dict_read_allele_clusters.keys():
        H1_name = fn_output_read_fasta.split('.')[0] + '_' + cluster_id + '_H1.' + fn_output_read_fasta.split('.')[1]
        H2_name = fn_output_read_fasta.split('.')[0] + '_' + cluster_id + '_H2.' + fn_output_read_fasta.split('.')[1]
        allele_file = fn_output_read_fasta.split('.')[0] + '_' + cluster_id + '_allele.' + fn_output_read_fasta.split('.')[1]
        
        f1_o = open(H1_name, 'w')
        f2_o = open(H2_name, 'w')
        f3_o = open(allele_file, 'w')
        
        dict_alleles = dict_read_allele_clusters[cluster_id][0]
        dict_reads = dict_read_allele_clusters[cluster_id][1]
        sorted_reads_name = sorted(dict_reads)
        for key in sorted_reads_name:
            if key.split('_')[1][0] == '1':
                f1_o.write('>' + key.split('_')[0] + ' ')
                f1_o.write(key.split('_')[1] + '\n')
                f1_o.write(dict_reads[key] + '\n')
            else:
                f2_o.write('>' + key.split('_')[0] + ' ')
                f2_o.write(key.split('_')[1] + '\n')
                f2_o.write(dict_reads[key] + '\n')
        f1_o.close()
        f2_o.close()

        for key in dict_alleles.keys():
            f3_o.write('>' + key + '\n')
            f3_o.write(dict_alleles[key] + '\n')
        f3_o.close()

    print("Dumping files finished!")



if __name__ == '__main__':
    args = parse_args()
    fn_read = args.fn_read
    fn_allele = args.fn_allele
    fn_pickle_file = args.fn_pickle_file
    fn_output_read_fasta = args.fn_output_read_fasta

    # load pickle if it's already there
    if os.path.exists(fn_pickle_file):
        print ('Pickle file {} has existed, load it...'.format(fn_pickle_file))
        f = open(fn_pickle_file, 'rb')
        dict_cluster_info = pickle.load(f)
        f.close()

        output_reads(dict_cluster_info, fn_read, fn_allele, fn_output_read_fasta)
    
