'''
An example of input file:

old
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

new
# BLASTN 2.10.0+
# Query: M01931:79:000000000-CHMNT:1:1101:21945:1007 1:N:0:46
# Database: ../../20200317_filtered_V_alleles_for_probe_design/all_probs.fasta
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 1 hits found
M01931:79:000000000-CHMNT:1:1101:21945:1007	AE000661|TRAV37*01|Homo	100.000	13	0	0	140	152	144	156	3.0	24.7
# BLASTN 2.10.0+
# Query: M01931:79:000000000-CHMNT:1:1101:21666:1012 1:N:0:46
# Database: ../../20200317_filtered_V_alleles_for_probe_design/all_probs.fasta
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 2 hits found
M01931:79:000000000-CHMNT:1:1101:21666:1012	IMGT000035|IGHV(III)-2-1*02|Homo	95.000	20	1	0	231	250	119	138	0.020	32.8
M01931:79:000000000-CHMNT:1:1101:21666:1012	AB019441|IGHV(III)-2-1*01|Homo	95.000	20	1	0	231	250	119	138	0.020	32.8
# BLASTN 2.10.0+
# Query: M01931:79:000000000-CHMNT:1:1101:16553:1013 1:N:0:46
# Database: ../../20200317_filtered_V_alleles_for_probe_design/all_probs.fasta
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 122 hits found
M01931:79:000000000-CHMNT:1:1101:16553:1013	D86996|IGLV11-55*01|Homo	100.000	98	0	0	2	99	98	1	2.15e-46	178
M01931:79:000000000-CHMNT:1:1101:16553:1013	KM455555|IGLV11-55*02|Homo	98.980	98	1	0	2	99	98	1	9.15e-45	173
'''

import argparse
from utils import read_allele_cluster, get_allele_cluster
import pickle

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', '--fn_input',
        help = 'input file (blastn log file)'
    )
    parser.add_argument(
        '-c', '--fn_cluster',
        help = 'allele cluster report (from clustal-omega)'
    )
    parser.add_argument(
        '-n', '--top_n',
        help = 'number of top alleles printed (10)',
        default = 10,
        type = int
    )
    parser.add_argument(
        '--identity_thrsd',
        help = 'minimum identity level that we take into account ([0,100]) [0]',
        default = 0,
        type = int
    )
    parser.add_argument(
        '-s', '--scoring',
        help = 'scoring of a hit, ["bit_score", "count"] ("bit_score")',
        default = "bit_score"
    )
    parser.add_argument(
        '--only_take_top',
        help = 'enabled to only use the top-1 result, otherwise split the score/count between tied hits (False)',
        action = 'store_true'
    )
    parser.add_argument(
        '--fn_allele_len',
        help = 'TSV file specifying allele lengths'
    )
    parser.add_argument(
        '-o', '--fn_output',
        help = 'output file including only alleles [None]'
    )
    parser.add_argument(
        '-ocp', '--fn_output_cluster_pickle',
        help = 'output pickle file that stores cluster_info dict [None]'
    )
    args = parser.parse_args()
    return args

def update_cluster_info(dict_cluster_info, allele_name, read_name, dict_allele_cluster):
    '''
    `dict_cluster_info` is mutable
    '''
    cluster_id = dict_allele_cluster.get(allele_name)
    if not cluster_id:
        # allele_name not in cluster
        return
    if dict_cluster_info.get(cluster_id):
        dict_cluster_info[cluster_id][0].add(allele_name)
        dict_cluster_info[cluster_id][1].add(read_name)
    else:
        dict_cluster_info[cluster_id] = [{allele_name}, {read_name}]

def process_blastn_log_chunk(dict_allele_count, chunk_records, scoring, only_take_top, dict_allele_len, dict_allele_cluster, dict_cluster_info):
    size = len(chunk_records)
    # size = 1

    # scoring
    if scoring == 'count':
        score = 1 / size
    elif scoring == 'bit_score':
        score = float(chunk_records[0][11]) / size
    
    read_name = chunk_records[0][0]
    if only_take_top:
        name = chunk_records[0][1].split('|')[1]
        if dict_allele_len:
            allele_len = dict_allele_len[name]
        else:
            allele_len = 1

        if dict_allele_count.get(name):
            dict_allele_count[name] += float(score / allele_len)
        else:
            dict_allele_count[name] = float(score / allele_len)

        update_cluster_info(dict_cluster_info, name, read_name, dict_allele_cluster)
    else:
        for i, _ in enumerate(chunk_records):
            name = chunk_records[i][1].split('|')[1]
            if dict_allele_len:
                allele_len = dict_allele_len[name]
            else:
                allele_len = 1

            if dict_allele_count.get(name):
                dict_allele_count[name] += float(score / allele_len)
            else:
                dict_allele_count[name] = float(score / allele_len)

            update_cluster_info(dict_cluster_info, name, read_name, dict_allele_cluster)
    
    return dict_allele_count

def process_blastn_log(fn_input, scoring, only_take_top, dict_allele_len, dict_allele_cluster, identity_thrsd=0):
    f = open(fn_input, 'r')
    new_record_flag = True
    dict_allele_count = {}
    chunk_records = []
    num_records = 0
    dict_cluster_info = {}
    # dict_cluster_info
    #   - keys: cluster_id
    #   - values: [{a set of allele names}, {a set of read names}]

    # one chunk associates with one read
    for line in f:
        if line[0] == '#':
            if len(chunk_records) > 0:
                num_records += 1
                dict_allele_count = process_blastn_log_chunk(dict_allele_count, chunk_records, scoring, only_take_top, dict_allele_len, dict_allele_cluster, dict_cluster_info)
            # initialize chunk
            new_record_flag = True
            chunk_records = []
            continue
        elif new_record_flag:
            fields = line.split()
            # the first hiti
            if len(chunk_records) == 0:
                # check if the identity >= identity_thrsd
                if (float(fields[2]) >= identity_thrsd):
                    chunk_records.append(fields)
                else:
                    new_record_flag = True
            # if tied bit score
            elif fields[11] == chunk_records[-1][11]:
                if (float(fields[2]) >= identity_thrsd):
                    chunk_records.append(fields)
                else:
                    new_record_flag = True
            else:
                new_record_flag = False
    # process the last record
    if len(chunk_records) > 0:
        num_records += 1
        dict_allele_count = process_blastn_log_chunk(dict_allele_count, chunk_records, scoring, only_take_top, dict_allele_cluster, dict_cluster_info)
    
    print ('Number of records processed:', num_records)

    return dict_allele_count, dict_cluster_info

def print_dict_allele_count(dict_allele_count, top_n, fn_output):
    sorted_dict = sorted(dict_allele_count.items(), key=lambda x: x[1], reverse=True)
    if top_n > len(sorted_dict):
        print ('Warning: top_n({}) exceeds the number of alleles found, only print {} alleles'.format(top_n, len(sorted_dict)))
        top_n = len(sorted_dict)
    for i in range(top_n):
        print (sorted_dict[i][0], ':', int(sorted_dict[i][1]))
    # print (sorted_dict[: top_n])
    
    if fn_output:
        with open(fn_output, 'w') as f_o:
            for i in range(top_n):
                f_o.write(sorted_dict[i][0] + '\t' + str(sorted_dict[i][1]) + '\n')
                # f_o.write(sorted_dict[i][0].split('|')[1] + '\t' + str(sorted_dict[i][1]) + '\n')

def get_allele_len(fn_allele_len):
    '''
    Read allele lengths TSV file
    '''
    dc = {}
    with open(fn_allele_len, 'r') as f:
        cnt = 0
        for line in f:
            if cnt == 0:
                cnt += 1
            else:
                dc[line.split('\t')[0]] = int(line.split('\t')[1])
    return dc

if __name__ == '__main__':
    args = parse_args()
    fn_input = args.fn_input
    fn_cluster = args.fn_cluster
    fn_output = args.fn_output
    top_n = args.top_n
    scoring = args.scoring
    identity_thrsd = args.identity_thrsd
    fn_allele_len = args.fn_allele_len
    assert scoring in ['bit_score', 'count']
    only_take_top = args.only_take_top
    fn_output_cluster_pickle = args.fn_output_cluster_pickle

    if fn_allele_len:
        dict_allele_len = get_allele_len(fn_allele_len)
    else:
        dict_allele_len = None

    if fn_cluster:
        dict_allele_cluster = read_allele_cluster(fn_cluster)

    dict_allele_count, dict_cluster_info = process_blastn_log(fn_input, scoring, only_take_top, dict_allele_len, dict_allele_cluster, identity_thrsd)
    print_dict_allele_count(dict_allele_count, top_n, fn_output)

    if fn_output_cluster_pickle:
        f = open(fn_output_cluster_pickle, 'wb')
        pickle.dump(dict_cluster_info, f)
        f.close()
