import argparse
import pickle
import os
import numpy as np
import sys

# make sure the package modules is in the path
sys.path.append(os.path.dirname(__file__))
sys.path.append(os.path.dirname(__file__)+'/..')

from utils import get_reverse_complement, eprint
from filter_corrected_alleles import parse_fasta
from parse_cluster_realign import cluster_separate, variant_link_graph, haplotyping_link_graph, haplotyping_link_graph, mark_edit_region, output_contig_correction
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fs', '--fn_sam',
        help = 'input target shrinking file'
    )
    parser.add_argument(
        '-fc', '--fn_contig',
        help = 'input full length contig fasta file'
    )
    parser.add_argument(
        '-fr', '--fn_allele_region',
        help = 'input allele region of contigs'
    )
    parser.add_argument(
        '-ext', '--len_extend', type=int,
        help = 'the flanking sequence length cropped from the contigs'
    )
    
    parser.add_argument(
        '-foh', '--fo_flanking_haplotype',
        help = 'output haplotyped flanking sequence fasta file'
    )
    parser.add_argument(
        '-foc', '--fo_shrinked_contigs',
        help = 'output shrinked contigs fasta file'
    )
    args = parser.parse_args()
    return args

def parse_allele_region(region_report):
    dict_contig_region = {}
    f_r = open(region_report, 'r')
    for line in f_r:
        region, contig_name = line.strip().split(',')
        dict_contig_region[contig_name] = (int(region.split('-')[0]), int(region.split('-')[1]))
    return dict_contig_region



if __name__ == '__main__':
    args = parse_args()
    fn_sam = args.fn_sam
    fn_contig = args.fn_contig
    fn_allele_region = args.fn_allele_region
    len_extend = args.len_extend
    fo_flanking_haplotype = args.fo_flanking_haplotype
    fo_shrinked_contigs = args.fo_shrinked_contigs

    dict_contig_region = parse_allele_region(fn_allele_region)
    dict_cluster_contig = parse_fasta(fn_contig)
    f_oc = open(fo_shrinked_contigs, 'w')
    for contig_name, SEQ in sorted(dict_cluster_contig.items()):
        region = dict_contig_region[contig_name]
        cut_start = max(0,region[0]-len_extend)
        cut_end = min(len(SEQ), region[1]+len_extend)
        dict_contig_region[contig_name] = (cut_start, cut_end)
        shrinked_SEQ = SEQ[cut_start:cut_end]
        f_oc.write('>' + contig_name + '\n')
        f_oc.write(shrinked_SEQ + '\n')
    f_oc.close()

    dict_contig = cluster_separate(fn_contig, fn_sam)

    for contig_name, contig_info in sorted(dict_contig.items()):
        #parse the sam file and generate
        edit_histogram, cov_histogram, list_read_info = mark_edit_region(contig_name, contig_info, True)
        
        #determine the region contains alternative flanking region
        edit_region = []
        interest_region = dict_contig_region[contig_name]
        for idx, ele in enumerate(edit_histogram):
            print(str(idx) + ':\t' + str(cov_histogram[idx])  + '\t' + str(ele))
            if ele > cov_histogram[idx]/4.2 and (interest_region[0] < idx < interest_region[1]):
                edit_region.append(idx)

        eprint(contig_name, interest_region)
        print(contig_name, edit_region)
        
        contig_SEQ = dict_contig[contig_name][2]
        interest_region_str = str(interest_region[0]) + '-' + str(interest_region[1])
        interest_edit_region = edit_region
        if interest_edit_region != [] and min(cov_histogram[interest_region[0]:interest_region[1]]) > 20:
            print("=========== allele correction ==============")
            eprint("CORRECT", contig_name.split('_')[0], min(cov_histogram[interest_region[0]:interest_region[1]]), interest_edit_region)
            dict_link_graph, dict_var_weight, dict_link_outward, dict_link_inward = variant_link_graph(interest_edit_region, list_read_info)
            haplotype_0, haplotype_1 = haplotyping_link_graph(dict_link_graph, dict_var_weight, dict_link_outward, dict_link_inward, interest_region_str)
            #output_contig_correction(contig_SEQ, region_st, region_ed, haplotype_0, haplotype_1, contig_name, corrected_contig_output_file)
            output_contig_correction(contig_SEQ.lower(), interest_region[0], interest_region[1], haplotype_0, haplotype_1, contig_name.split('_')[0], fo_flanking_haplotype)
        #elif interest_edit_region != []:
        elif min(cov_histogram[interest_region[0]:interest_region[1]]) < 20:
            eprint("Deficient", contig_name.split('_')[0], min(cov_histogram[interest_region[0]:interest_region[1]]), interest_edit_region)
            print("=== cov not efficient:", min(cov_histogram[1:]), "=======")
        else:
            eprint("No variant", contig_name.split('_')[0], min(cov_histogram[interest_region[0]:interest_region[1]]), interest_edit_region)
            f_c = open(fo_flanking_haplotype, 'a')
            f_c.write(">" + contig_name.split('_')[0] + "/original\n")
            f_c.write(contig_SEQ.lower()[interest_region[0]:interest_region[1]] + "\n")
            f_c.close()
            print("============ No novel allele ===============")
     



