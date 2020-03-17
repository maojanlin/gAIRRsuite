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
        '-s', '--scoring',
        help = 'scoring of a hit, ["bit_score", "count"] ("bit_score")',
        default = "bit_score"
    )
    parser.add_argument(
        '--only-take-top',
        help = 'enabled to only use the top-1 result, otherwise split the score/count between tied hits (False)',
        action = 'store_true'
    )
    parser.add_argument(
        '-o', '--fn_output',
        help = 'output file'
    )
    args = parser.parse_args()
    return args

def process_blastn_log_chunk(chunk_records, scoring, only_take_top):
    size = len(chunk_records)
    # scoring
    if scoring == 'count':
        score = 1 / size
    elif scoring == 'bit_score':
        score = float(chunk_records[0][3]) / size
    
    if only_take_top:
        if dict_allele_count.get(chunk_records[0][1]):
            dict_allele_count[chunk_records[0][1]] += score
        else:
            dict_allele_count[chunk_records[0][1]] = score
    else:
        for i, _ in enumerate(chunk_records):
            if dict_allele_count.get(chunk_records[i][1]):
                dict_allele_count[chunk_records[i][1]] += score
            else:
                dict_allele_count[chunk_records[i][1]] = score
    
    return dict_allele_count

def process_blastn_log(fn_input, scoring, only_take_top):
    f = open(fn_input, 'r')
    new_record_flag = True
    dict_allele_count = {}
    chunk_records = []
    for line in f:
        if line[0] == '#':
            if len(chunk_records) > 0:
                dict_allele_count = process_blastn_log_chunk(chunk_records, scoring, only_take_top)
            # initialize chunk
            new_record_flag = True
            chunk_records = []
            continue
        elif new_record_flag:
            fields = line.split()
            # the first hit
            if len(chunk_records) == 0:
                chunk_records.append(fields)
            # if tied bit score
            elif fields[3] == chunk_records[-1][3]:
                chunk_records.append(fields)
            else:
                new_record_flag = False
    # process the last record
    if len(chunk_records) > 0:
        dict_allele_count = process_blastn_log_chunk(chunk_records, scoring, only_take_top)

    return dict_allele_count

def print_dict_allele_count(dict_allele_count, top_n):
    sum_value = 0
    for i, (key, value) in enumerate(dict_allele_count.items()):
        # print (key, value)
        sum_value += value
    print ('Number of reads processed:', sum_value)
    sorted_dict = sorted(dict_allele_count.items(), key=lambda x: x[1], reverse=True)
    print (sorted_dict[: top_n])

if __name__ == '__main__':
    args = parse_args()
    fn_input = args.fn_input
    fn_output = args.fn_output
    top_n = args.top_n
    scoring = args.scoring
    assert scoring in ['bit_score', 'count']
    only_take_top = args.only_take_top

    dict_allele_count = process_blastn_log(fn_input, scoring, only_take_top)
    print_dict_allele_count(dict_allele_count, top_n)
