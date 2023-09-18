import argparse
import pickle
import os
import numpy as np
from parse_sam_haplotyping import parse_CIGAR
from coverage_analysis import read_pickle

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fna', '--fn_allele_contig',
        help = 'allele_contig file show both contig and allele information'
    )
    parser.add_argument(
        '-fpa', '--fn_allele_pickle',
        help = 'pickle file store the allele locus information'
    )
    parser.add_argument(
        '-fpf', '--fn_flanking_pickle',
        help = 'pickle file store the flanking coverage information'
    )
    args = parser.parse_args()
    return args


def flanking_coverage(dict_allele, dict_flanking):
    list_contig_name = sorted(dict_allele)
    # dict from dict_flanking
    dict_FP = {}
    dict_TP = {}
    # dict from dict_allele
    dict_FN = {}
    for contig_name in list_contig_name:
        if dict_flanking.get(contig_name):
            set_allele = dict_allele[contig_name]
            set_flanking = dict_flanking[contig_name]
            for allele_region in set_allele:
                cover_flag = False
                for flanking_region in set_flanking:
                    if flanking_region[0] <= allele_region[0] and flanking_region[1] >= allele_region[1]:
                        cover_flag = True
                        if dict_TP.get(contig_name):
                            dict_TP[contig_name].add(flanking_region)
                        else:
                            dict_TP[contig_name] = {flanking_region}
                        #break
                if cover_flag == False:
                    if dict_FN.get(contig_name):
                        dict_FN[contig_name].add(allele_region)
                    else:
                        dict_FN[contig_name] = {allele_region}
        else:
            dict_FN[contig_name] = dict_allele[contig_name]
    
    num_TP = 0
    num_FN = 0
    num_annotated_locus = 0
    num_flanking_locus = 0
    for contig_name in sorted(dict_TP):
        print(contig_name)
        print(sorted(dict_TP[contig_name]))
        num_TP += len(dict_TP[contig_name])
    print("===============================")
    for contig_name in sorted(dict_FN):
        print(contig_name)
        print(sorted(dict_FN[contig_name]))
        num_FN += len (dict_FN[contig_name])
    for contig_name in dict_allele.keys():
        #print(contig_name + '\t' + str(len(dict_allele[contig_name])))
        num_annotated_locus += len(dict_allele[contig_name])
    for contig_name in dict_flanking.keys():
        num_flanking_locus += len(dict_flanking[contig_name])
    print("number of annotated locus:" + str(num_annotated_locus))
    print("number of TP:" + str(num_TP))
    print("number of FN:" + str(num_FN))
    print("number of flanking locus:" + str(num_flanking_locus))


def read_annotation_contigs(fn_allele_contig):
    dict_locus = {}
    list_locus = []
    with open (fn_allele_contig) as f_a:
        for line in f_a:
            list_locus.append(line.strip())

    for locus in list_locus:
        fields = locus.split(',')
        contig_name = fields[0]
        if dict_locus.get(contig_name):
            dict_locus[contig_name].add((int(fields[2]),int(fields[3]), fields[1]))
        else:
            dict_locus[contig_name] = {(int(fields[2]),int(fields[3]), fields[1])}
    
    return dict_locus


if __name__ == '__main__':
    args = parse_args()
    fn_allele_contig = args.fn_allele_contig
    fn_allele_pickle = args.fn_allele_pickle
    fn_flanking_pickle = args.fn_flanking_pickle

    dict_allele = {}
    # load pickle if it's already there
    if os.path.exists(fn_allele_pickle):
        print ('Pickle file {} has existed, load for it instead of re-calculating'.format(fn_allele_pickle))
        f = open(fn_allele_pickle, 'rb')
        dict_allele = pickle.load(f)
        f.close()
    else:
        dict_allele = read_annotation_contigs(fn_allele_contig)
        if fn_allele_pickle:
            f = open(fn_allele_pickle, 'wb')
            pickle.dump(dict_allele, f)
            f.close()

    dict_flanking = read_pickle(fn_flanking_pickle)
    
    flanking_coverage(dict_allele, dict_flanking)
