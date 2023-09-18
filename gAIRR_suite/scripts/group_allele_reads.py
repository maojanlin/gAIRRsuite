import argparse
import pickle
import os
import sys
import pyfastx


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fr1', '--fn_read_1',
        help = 'fasta file with all R1 reads'
    )
    parser.add_argument(
        '-fr2', '--fn_read_2',
        help = 'fasta file with all R2 reads'
    )
    parser.add_argument(
        '-fa', '--fn_allele',
        help = 'fasta file with all target alleles'
    )
    parser.add_argument(
        '-fp', '--fn_pickle_file',
        help = 'input pickle file indicating supporting reads of the alleles'
    )
    parser.add_argument(
        '-name', '--target_name',
        help = 'the target allele names'
    )
    parser.add_argument(
        '-fod', '--fo_output_dir',
        help = 'target output directory containing all alleles and reads'
    )
    args = parser.parse_args()
    return args

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def read_pickle(pickle_file):
    f = open(pickle_file, 'rb')
    x = pickle.load(f)
    return x


def group_allele_reads(dict_allele_support_reads, dict_allele, dict_read_1, dict_read_2, target_name, fo_output_dir):
    idx = 0
    f_omap = open((fo_output_dir + target_name + '_group_allele_map.txt'), 'w')
    for allele_name, list_read_name in sorted(dict_allele_support_reads.items()):
        if len(list_read_name) == 0:
            continue
        f_oa  = open((fo_output_dir + target_name + '_' + str(idx) + '_allele.fasta'), 'w')
        f_oa.write('>' + allele_name + '\n')
        f_oa.write(dict_allele[allele_name] + '\n')
        f_oa.close()
        
        if dict_read_2:
            f_or1 = open((fo_output_dir + target_name + '_' + str(idx) + '_read_R1.fasta'), 'w')
            f_or2 = open((fo_output_dir + target_name + '_' + str(idx) + '_read_R2.fasta'), 'w')
            for read_name in sorted(list_read_name):
                f_or1.write('>' + read_name + ' 1\n')
                f_or1.write(dict_read_1[read_name] + '\n')
                f_or2.write('>' + read_name + ' 2\n')
                f_or2.write(dict_read_2[read_name] + '\n')
            f_or1.close()
            f_or2.close()
        else:
            f_or1 = open((fo_output_dir + target_name + '_' + str(idx) + '_read.fasta'), 'w')
            for read_name in sorted(list_read_name):
                f_or1.write('>' + read_name + '\n')
                f_or1.write(dict_read_1[read_name] + '\n')
            f_or1.close()
        
        f_omap.write(str(idx) + ',' + allele_name + '\n')
        idx += 1
    f_omap.close()
    return idx
    

def parse_fastx(fn_fastx):
    "parser for both fasta and fastq file"
    dict_read = {}
    f_fastx = pyfastx.Fastx(fn_fastx)
    for seq_info in f_fastx:
        dict_read[seq_info[0]] = seq_info[1]
    return dict_read



if __name__ == '__main__':
    args = parse_args()
    fn_read_1 = args.fn_read_1
    fn_read_2 = args.fn_read_2
    fn_allele = args.fn_allele
    fn_pickle_file = args.fn_pickle_file
    target_name = args.target_name
    fo_output_dir = args.fo_output_dir

    # load the pickle, allele and R1 R2 reads
    dict_allele_support_reads = read_pickle(fn_pickle_file)
    try:
        dict_allele_support_reads = {name.split('|')[1]:list_reads for name, list_reads in dict_allele_support_reads.items()}
    except:
        dict_allele_support_reads = {name.split()[0]:list_reads for name, list_reads in dict_allele_support_reads.items()}
    #print(dict_allele_support_reads.keys())
    dict_allele = parse_fastx(fn_allele)
    try:
        dict_allele = {name.split('|')[1]:SEQ for name, SEQ in dict_allele.items()}
    except:
        dict_allele = {name.split()[0]:SEQ for name, SEQ in dict_allele.items()}
    #print(dict_allele.keys())
    dict_read_1 = parse_fastx(fn_read_1)
    group_num = 0
    if fn_read_2:
        dict_read_2 = parse_fastx(fn_read_2)
        group_num = group_allele_reads(dict_allele_support_reads, dict_allele, dict_read_1, dict_read_2, target_name, fo_output_dir)
    else:
        group_num = group_allele_reads(dict_allele_support_reads, dict_allele, dict_read_1, None, target_name, fo_output_dir)
    print(group_num)
