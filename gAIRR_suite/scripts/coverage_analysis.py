import argparse
import pickle
import os
import numpy as np
from scripts.utils import get_hamming_dist
from scripts.utils import get_reverse_complement
import sys

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
        '-ans', '--fn_annotation',
        help = 'answer file, for development'
    )
    parser.add_argument(
        '-md', '--min_depth',
        default=0,
        help = 'minimun required read-depth for hc allele'
    )
    parser.add_argument(
        '-fnp', '--fn_pickle_file', 
        help = 'input pickle file indicating reads and alleles in the clusters'
    )
    parser.add_argument(
        '-fop', '--fn_output_cluster_pickle',
        help = 'output file of the cluster with reads and alleles [None]'
    )
    parser.add_argument(
        '-fsp', '--fn_output_support_read_pickle',
        help = 'output file of indicating reads that support an alleles [None]'
    )
    args = parser.parse_args()
    return args

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def read_pickle(pickle_file):
    f = open(pickle_file, 'rb')
    x = pickle.load(f)
    return x

def add_alleles(allele_name, allele_SEQs, dict_cluster_info, dict_read_allele_clusters):
    for cluster_id in dict_cluster_info:
        if allele_name in dict_cluster_info[cluster_id][0]:
            if dict_read_allele_clusters.get(cluster_id):
                dict_read_allele_clusters[cluster_id][0][allele_name] = allele_SEQs
            else:
                dict_read_allele_clusters[cluster_id] = [{allele_name: allele_SEQs}, {}]
    return dict_read_allele_clusters

def add_reads(read_name, read_SEQs, dict_cluster_info, dict_read_allele_clusters):
    for cluster_id in dict_cluster_info:
        if read_name in dict_cluster_info[cluster_id][1]:
            if dict_read_allele_clusters.get(cluster_id):
                dict_read_allele_clusters[cluster_id][1][read_name] = read_SEQs
            else:
                print("Warning: cluster_id errors!")
    return dict_read_allele_clusters

def fetch_reads_alleles(fn_read, fn_allele, dict_cluster_info):
    dict_read_allele_clusters = {}
    # dict_read_allele_clusters
    #   - keys: cluster_id
    #   - values: [{a dict with {allele names: ALLELE_SEQs}, {read names: READ_SEQs}]

    with open(fn_allele, 'r') as f_a:
        allele_name = ""
        allele_SEQs = ""
        for line in f_a:
            if line[0] == '>':
                dict_read_allele_clusters = add_alleles(allele_name, allele_SEQs, dict_cluster_info, dict_read_allele_clusters)

                #set new allele name and reset SEQs
                allele_name = line.split('|')[1]
                allele_SEQs = ""
            else:
                allele_SEQs = allele_SEQs + line.strip()
        dict_read_allele_clusters = add_alleles(allele_name, allele_SEQs, dict_cluster_info, dict_read_allele_clusters)

    with open(fn_read, 'r') as f_r:
        read_name = ""
        read_SEQs = ""
        for line in f_r:
            if line[0] == '>':
                dict_read_allele_clusters = add_reads(read_name, read_SEQs, dict_cluster_info, dict_read_allele_clusters)

                #set new read name and reset SEQs
                read_name = line.split()[0]
                read_name = read_name[1:]
                read_SEQs = ""
            else:
                read_SEQs = read_SEQs + line.strip()
        dict_read_allele_clusters = add_reads(read_name, read_SEQs, dict_cluster_info, dict_read_allele_clusters)
                
    return dict_read_allele_clusters

def hamming_traverse(
    seq_allele,
    read,
    min_coverage = 100,
    min_hamming  = 1
):
    '''
    get the allele and read, traverse all possible relative position until match found.
    return triplet (True/False, start position, end position)
    does not concern the reverse_complement
    '''
    if len(read) < min_coverage:
        return (False,0,0)
    if len(seq_allele) > len(read):
        for i in range(len(seq_allele) - len(read) + 1):
            if get_hamming_dist(
                str_a=seq_allele[i: i+len(read)],
                str_b=read,
                thrsd=int(min_hamming),
                try_both_orient=False
            ):
                return (True, i, i+len(read))
        
        for i in range(min_coverage, len(read)):
            if get_hamming_dist(
                str_a=seq_allele[0:i],
                str_b=read[len(read)-i: len(read)],
                thrsd=int(min_hamming),
                try_both_orient=False
            ):
                return (True, 0, i)

        for i in range(len(seq_allele)-len(read), len(seq_allele)-min_coverage):
            if get_hamming_dist(
                str_a=seq_allele[i:],
                str_b=read[: len(seq_allele)-i],
                thrsd=int(min_hamming),
                try_both_orient=False
            ):
                return (True, i, len(seq_allele))
    else: # len(seq_allele) <= len(read)
        for i in range(len(read) - len(seq_allele) + 1):
            if get_hamming_dist(
                str_a=seq_allele,
                str_b=read[i: i+len(seq_allele)],
                thrsd=int(min_hamming),
                try_both_orient=False
            ):
                return (True, 0, len(seq_allele))
        
        for i in range(min_coverage, len(seq_allele)):
            if get_hamming_dist(
                str_a=seq_allele[0:i],
                str_b=read[len(read)-i: len(read)],
                thrsd=int(min_hamming),
                try_both_orient=False
            ):
                return (True, 0, i)

        for i in range(0, len(seq_allele)-min_coverage):
            if get_hamming_dist(
                str_a=seq_allele[i:],
                str_b=read[: len(seq_allele)-i],
                thrsd=int(min_hamming),
                try_both_orient=False
            ):
                return (True, i, len(seq_allele))
    
    return (False, 0, 0)
    

def coverage_analysis(
    dict_read_allele_clusters,
    fn_annotation,
    required_min_depth = 0,
    required_single_coverage = 50,
    required_single_identity = 1,
):

    dict_hc_calls = {}
    dict_sup_reads = {}

    # tmp: for dev
    list_annotated = []
    #f_tmp = open('./NA12878_annotated_all.txt', 'r')
    f_tmp = open(fn_annotation, 'r')
    for line in f_tmp:
        list_annotated.append(line.rstrip())

    # tmp
    list_answer = []

    # for each cluster
    #for cluster_id in dict_read_allele_clusters.keys():
    for cluster_id in range(55,len(dict_read_allele_clusters.keys()),50):
        print("Cluster: " + str(cluster_id))
        eprint("============= Cluster: " + str(cluster_id) + " ==============")
        cluster = dict_read_allele_clusters[str(cluster_id)]
        dict_allele = cluster[0]
        dict_read = cluster[1]

        # for each allele in a cluster
        for allele in dict_allele.keys():
            print(allele)
            seq_allele = dict_allele[allele]
            seq_coverage = np.zeros( len(seq_allele) )
            dict_sup_reads[allele] = set()

            # tmp
            if allele in list_annotated:
                list_answer.append(allele)
            
            # need simplify
            for read in dict_read:
                seq_read = dict_read[read]
                # ignore the reads that are too short
                if len(seq_read) < required_single_coverage:
                    continue

                traverse_result = hamming_traverse(seq_allele, seq_read, required_single_coverage, required_single_identity)
                if traverse_result[0]:
                    seq_coverage[traverse_result[1]:traverse_result[2]] += 1
                    dict_sup_reads[allele].add(read)
                else:
                    r_seq_read = get_reverse_complement(seq_read)
                    traverse_result = hamming_traverse(seq_allele, r_seq_read, required_single_coverage, required_single_identity)
                    if traverse_result[0]:
                        seq_coverage[traverse_result[1]:traverse_result[2]] += 1
                        dict_sup_reads[allele].add(read)
            
            if min(seq_coverage) > required_min_depth:
                if dict_hc_calls.get(allele):
                    print("Warning! evaluate two times")
                    dict_hc_calls[allele] += 1
                else:
                    dict_hc_calls[allele] = min(seq_coverage)
                print("OOO: " + str(min(seq_coverage)) + ' ' + str(sum(seq_coverage)/len(seq_coverage)) + ' ' + str(max(seq_coverage)))
            else:
                print("XXX: " + str(min(seq_coverage)) + ' ' + str(sum(seq_coverage)/len(seq_coverage)) + ' ' + str(max(seq_coverage)))
            print(seq_coverage)


    print (dict_hc_calls)
    print (list_answer)

    # tmp
    print ('Num. high-confidence calls')
    print (len(set(dict_hc_calls)))
    print ('Num. answer')
    print (len(set(list_answer)))
    print ('Num. intersection')
    print (len(set(list_answer).intersection(set(dict_hc_calls))))

    #print ("Support reads of alleles:")
    #print (dict_sup_reads)

    return dict_sup_reads
    




if __name__ == '__main__':
    args = parse_args()
    fn_read = args.fn_read
    fn_allele = args.fn_allele
    fn_annotation = args.fn_annotation
    min_depth = int(args.min_depth)
    fn_pickle_file = args.fn_pickle_file
    fn_output_cluster_pickle = args.fn_output_cluster_pickle
    fn_output_support_read_pickle = args.fn_output_support_read_pickle

    # load pickle if it's already there
    if os.path.exists(fn_output_cluster_pickle):
        print ('Pickle file {} has existed, load for it instead of re-calculating'.format(fn_output_cluster_pickle))
        f = open(fn_output_cluster_pickle, 'rb')
        dict_read_allele_clusters = pickle.load(f)
        f.close()
    else:
        dict_cluster_info = read_pickle(fn_pickle_file)
        dict_read_allele_clusters = fetch_reads_alleles(fn_read, fn_allele, dict_cluster_info)

        if fn_output_cluster_pickle:
            f = open(fn_output_cluster_pickle, 'wb')
            pickle.dump(dict_read_allele_clusters, f)
            f.close()

    dict_sup_reads = coverage_analysis(dict_read_allele_clusters, fn_annotation, min_depth)
    if fn_output_support_read_pickle:
        f = open(fn_output_support_read_pickle, 'wb')
        pickle.dump(dict_sup_reads, f)
        f.close()

