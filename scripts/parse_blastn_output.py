'''
# BLASTN 2.10.0+
# Query: M01931:79:000000000-CHMNT:1:1101:21945:1007 1:N:0:46
# Database: ../../../IMGT_DBSEQ/TRBV.fasta.txt
# Fields: query id, subject id, evalue, bit score
# 11 hits found
M01931:79:000000000-CHMNT:1:1101:21945:1007    X57607|TRBV7-7*02|Homo    0.45    24.7
M01931:79:000000000-CHMNT:1:1101:21945:1007    L36092|TRBV7-7*01|Homo    0.45    24.7
M01931:79:000000000-CHMNT:1:1101:21945:1007    X58806|TRBV7-6*02|Homo    0.45    24.7
M01931:79:000000000-CHMNT:1:1101:21945:1007    L36092|TRBV7-6*01|Homo    0.45    24.7
M01931:79:000000000-CHMNT:1:1101:21945:1007    L13762|TRBV7-4*02|Homo    1.6    22.9
'''

import argparse
# import operator
# import collections

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', '--fn_input',
        help = 'input file (blastn log file)'
    )
    parser.add_argument(
        '-n', '--top_n',
        help = 'number of top alleles printed (10)',
        default = 10,
        type = int
    )
    parser.add_argument(
        '-o', '--fn_output',
        help = 'output file'
    )
    args = parser.parse_args()
    return args

def process_blastn_log(fn_input):
    f = open(fn_input, 'r')
    new_record_flag = True
    dict_allele_count = {}
    for line in f:
        if line[0] == '#':
            new_record_flag = True
            continue
        elif new_record_flag:
            fields = line.split()
            if dict_allele_count.get(fields[1]):
                dict_allele_count[fields[1]] += 1
            else:
                dict_allele_count[fields[1]] = 1
            new_record_flag = False
    return dict_allele_count

def print_dict_allele_count(dict_allele_count, top_n):
    sum_value = 0
    for i, (key, value) in enumerate(dict_allele_count.items()):
        # print (key, value)
        sum_value += value
    print ("Number of reads processed:", sum_value)
    sorted_dict = sorted(dict_allele_count.items(), key=lambda x: x[1], reverse=True)
    print (sorted_dict[: top_n])

if __name__ == '__main__':
    args = parse_args()
    fn_input = args.fn_input
    fn_output = args.fn_output
    top_n = args.top_n

    dict_allele_count = process_blastn_log(fn_input)
    print_dict_allele_count(dict_allele_count, top_n)
