import argparse
import pickle
import os
import numpy as np
from utils import get_reverse_complement, eprint
from filter_corrected_alleles import parse_fasta
from compare_annotation import parse_annotation
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fasm1', '--fn_asm_1',
        help = 'input assembly fasta file'
    )
    parser.add_argument(
        '-fasm2', '--fn_asm_2',
        help = 'input assembly fasta file'
    )
    parser.add_argument(
        '-fa', '--fn_annotation',
        help = 'input annotation txt file indicating the position of each alleles'
    )
    parser.add_argument(
        '-ext', '--len_extend', type=int,
        help = 'the flanking sequence length cropped from the contigs'
    )
    
    parser.add_argument(
        '-fof', '--fo_asm_flanking',
        help = 'output flanking sequence fasta file from asm'
    )
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    fn_asm_1 = args.fn_asm_1
    fn_asm_2 = args.fn_asm_2
    fn_annotation = args.fn_annotation
    len_extend = args.len_extend
    fo_asm_flanking = args.fo_asm_flanking

    dict_contig_H1 = parse_fasta(fn_asm_1)
    if fn_asm_2:
        dict_contig_H2 = parse_fasta(fn_asm_2)
    list_annotation = parse_annotation(fn_annotation)
    
    # dict_flank_SEQ {}
    # - keys: allele_name (novel)
    # - values: SEQ set {}
    dict_flank_SEQ = {}
    for annotation_info in list_annotation:
        allele_name = annotation_info[0]
        contig_name = annotation_info[1]
        start_pos = int(annotation_info[2])
        end_pos   = start_pos + int(annotation_info[3])
        contig_SEQ = ""
        if fn_asm_2:
            if dict_contig_H1.get(contig_name):
                contig_SEQ = dict_contig_H1[contig_name]
            elif dict_contig_H2.get(contig_name):
                contig_SEQ = dict_contig_H2[contig_name]
            else:
                eprint("Fatal Error! contig name " + contig_name + " not found!")
        else:
            contig_SEQ = dict_contig_H1[contig_name]
        left_flank  = max(0,start_pos-len_extend-1) 
        right_flank = min(len(contig_SEQ), end_pos+len_extend-1)
        flanking_SEQ = contig_SEQ[left_flank:right_flank]
        if len(annotation_info) > 4: # with mismatch
            allele_name = allele_name + "/novel"
        if dict_flank_SEQ.get(allele_name):
            if flanking_SEQ in dict_flank_SEQ[allele_name] or get_reverse_complement(flanking_SEQ) in dict_flank_SEQ[allele_name]:
                pass
            else:
                dict_flank_SEQ[allele_name].add(flanking_SEQ)
        else:
            dict_flank_SEQ[allele_name] = {flanking_SEQ}

    f_of = open(fo_asm_flanking, 'w')
    for allele_name, set_allele_SEQ in sorted(dict_flank_SEQ.items()):
        for idx, allele_SEQ in enumerate(sorted(set_allele_SEQ)):
            f_of.write('>' + allele_name + '-' + str(idx) + '\n')
            f_of.write(allele_SEQ.lower() + '\n')
    f_of.close()

