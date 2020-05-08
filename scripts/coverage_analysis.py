import argparse
import pickle
import os
import numpy as np
from utils import get_hamming_dist
from utils import get_reverse_complement

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
        '-fop', '--fn_output_cluster_pickle',
        help = 'output file of the cluster with reads and alleles [None]'
    )
    args = parser.parse_args()
    return args

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
                dict_read_allele_clusters[cluster_id] = [{allele_name: allele_SEQs}, set()]
    return dict_read_allele_clusters

def add_reads(read_name, read_SEQs, dict_cluster_info, dict_read_allele_clusters):
    for cluster_id in dict_cluster_info:
        if read_name in dict_cluster_info[cluster_id][1]:
            if dict_read_allele_clusters.get(cluster_id):
                dict_read_allele_clusters[cluster_id][1].add(read_SEQs)
            else:
                print("Warning: cluster_id errors!")
    return dict_read_allele_clusters


def fetch_reads_alleles(fn_read, fn_allele, dict_cluster_info):
    dict_read_allele_clusters = {}
    # dict_read_allele_clusters
    #   - keys: cluster_id
    #   - values: [{a dict with {allele names: ALLELE_SEQs}, {a set of read SEQs}]

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
                
    return dict_read_allele_clusters

def coverage_analysis(
    dict_read_allele_clusters,
    required_single_coverage = 100,
    required_total_coverage  = 1,
    required_single_identity = 1,
):

    dict_hcover_calls = {}

    # tmp: for dev
    list_annotated = []
    f_tmp = open('./NA12878_annotated_all.txt', 'r')
    for line in f_tmp:
        list_annotated.append(line.rstrip())

    # tmp
    list_answer = []

    # for each cluster
    for cluster_id in dict_read_allele_clusters.keys():
    #for cluster_id in range(9,25,15):
        print(cluster_id)
        cluster = dict_read_allele_clusters[str(cluster_id)]
        dict_allele = cluster[0]
        set_read = cluster[1]

        # for each allele in a cluster
        for allele in dict_allele.keys():
            seq_allele = dict_allele[allele]
            seq_coverage = np.zeros( len(seq_allele) )

            # tmp
            if allele in list_annotated:
                list_answer.append(allele)
            
            # need simplify
            for read in set_read:
                if len(read) < required_single_coverage:
                    continue
                if len(seq_allele) > len(read):
                    read_found = False
                    for i in range(len(seq_allele) - len(read) + 1):
                        # if dist <= threshold
                        if get_hamming_dist(
                            str_b=read, # reverse complement is on read
                            str_a=seq_allele[i: i+len(read)],
                            thrsd=int(required_single_identity)
                        ):
                            seq_coverage[i:i+len(read)] = 1
                            read_found = True
                            break
                    # the read position is already found
                    if read_found:
                        continue
                    
                    r_read = get_reverse_complement(read)
                    for i in range(required_single_coverage, len(read)):
                        # if dist <= threshold
                        if get_hamming_dist(
                            str_a=seq_allele[0:i],
                            str_b=read[len(read)-i: len(read)],
                            thrsd=int(required_single_identity),
                            try_both_orient=False
                        ):
                            seq_coverage[0:i] = 1
                            read_found = 1
                            break
                        if get_hamming_dist(
                            str_a=seq_allele[0:i],
                            str_b=r_read[len(read)-i: len(read)],
                            thrsd=int(required_single_identity),
                            try_both_orient=False
                        ):
                            seq_coverage[0:i] = 1
                            read_found = 1
                            break
                    # the read position is already found
                    if read_found:
                        continue

                    for i in range(len(seq_allele)-len(read), len(seq_allele)-required_single_coverage):
                        if get_hamming_dist(
                            str_a=seq_allele[i:],
                            str_b=read[: len(seq_allele)-i],
                            thrsd=int(required_single_identity),
                            try_both_orient=False
                        ):
                            seq_coverage[i:] = 1
                            break
                        if get_hamming_dist(
                            str_a=seq_allele[i:],
                            str_b=r_read[: len(seq_allele)-i],
                            thrsd=int(required_single_identity),
                            try_both_orient=False
                        ):
                            seq_coverage[i:] = 1
                            break
                else:
                    read_found = False
                    for i in range(len(read) - len(seq_allele) + 1):
                        # if dist <= threshold
                        if get_hamming_dist(
                            str_a=seq_allele,
                            str_b=read[i: i+len(seq_allele)],
                            thrsd=int(required_single_identity)
                        ):
                            read_found = True
                            if dict_hcover_calls.get(allele):
                                dict_hcover_calls[allele] += 1
                            else:
                                dict_hcover_calls[allele] = 1
                            break
                    # read_found means the allele is totally covered
                    if read_found:
                        break
                    
                    r_read = get_reverse_complement(read)
                    for i in range(required_single_coverage, len(seq_allele)):
                        # if dist <= threshold
                        if get_hamming_dist(
                            str_a=seq_allele[0:i],
                            str_b=read[len(read)-i: len(read)],
                            thrsd=int(required_single_identity),
                            try_both_orient=False
                        ):
                            seq_coverage[0:i] = 1
                            read_found = 1
                            break
                        if get_hamming_dist(
                            str_a=seq_allele[0:i],
                            str_b=r_read[len(read)-i: len(read)],
                            thrsd=int(required_single_identity),
                            try_both_orient=False
                        ):
                            seq_coverage[0:i] = 1
                            read_found = 1
                            break
                    # the read position is already found
                    if read_found:
                        continue

                    for i in range(0, len(seq_allele)-required_single_coverage):
                        if get_hamming_dist(
                            str_a=seq_allele[i:],
                            str_b=read[: len(seq_allele)-i],
                            thrsd=int(required_single_identity),
                            try_both_orient=False
                        ):
                            seq_coverage[i:] = 1
                            break
                        if get_hamming_dist(
                            str_a=seq_allele[i:],
                            str_b=r_read[: len(seq_allele)-i],
                            thrsd=int(required_single_identity),
                            try_both_orient=False
                        ):
                            seq_coverage[i:] = 1
                            break

            
            if sum(seq_coverage)/len(seq_coverage) >= required_total_coverage:
                if dict_hcover_calls.get(allele):
                    dict_hcover_calls[allele] += 1
                else:
                    dict_hcover_calls[allele] = sum(seq_coverage)/len(seq_coverage)


    print (dict_hcover_calls)
    print (list_answer)

    # tmp
    print ('Num. high-confidence calls')
    print (len(set(dict_hcover_calls)))
    print ('Num. answer')
    print (len(set(list_answer)))
    print ('Num. intersection')
    print (len(set(list_answer).intersection(set(dict_hcover_calls))))
    


def get_high_confidence_calls(
    dict_read_allele_clusters,
    required_coverage = 280,
    required_identity = 0.99
):
    '''
    Get high-confidence allele calls.
    High-confidence means high identity and high allele coverage
    '''

    dict_hc_calls = {}

    # tmp: for dev
    list_annotated = []
    f_tmp = open('./NA12878_annotated_all.txt', 'r')
    for line in f_tmp:
        list_annotated.append(line.rstrip())

    # tmp
    list_answer = []

    # for each cluster
    for cluster_id in dict_read_allele_clusters.keys():
        cluster = dict_read_allele_clusters[cluster_id]
        dict_allele = cluster[0]
        set_read = cluster[1]

        # for each allele in a cluster
        for allele in dict_allele.keys():
            seq_allele = dict_allele[allele]
            # tmp
            if allele in list_annotated:
                list_answer.append(allele)
            # tmp
            if len(seq_allele) > required_coverage:
                continue
            else:
                # print (allele)
                for read in set_read:
                    for i in range(len(read) - len(seq_allele) + 1):
                        # if dist <= threshold
                        if get_hamming_dist(
                            str_a=seq_allele,
                            str_b=read[i: i+len(seq_allele)],
                            thrsd=int(len(seq_allele)*(1-required_identity))
                        ):
                            if dict_hc_calls.get(allele):
                                dict_hc_calls[allele] += 1
                            else:
                                dict_hc_calls[allele] = 1
                            break
        
    print (dict_hc_calls)
    print (list_answer)

    # tmp
    print ('Num. high-confidence calls')
    print (len(set(dict_hc_calls)))
    print ('Num. answer')
    print (len(set(list_answer)))
    print ('Num. intersection')
    print (len(set(list_answer).intersection(set(dict_hc_calls))))



if __name__ == '__main__':
    args = parse_args()
    fn_read = args.fn_read
    fn_allele = args.fn_allele
    fn_pickle_file = args.fn_pickle_file
    fn_output_cluster_pickle = args.fn_output_cluster_pickle

    # load pickle if it's already there
    if os.path.exists(fn_output_cluster_pickle):
        print ('Pickle file {} has existed, load for it instread re-calculating'.format(fn_output_cluster_pickle))
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

    coverage_analysis(dict_read_allele_clusters)

