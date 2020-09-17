import argparse
import pickle
import os
import numpy as np
from utils import get_reverse_complement
from filter_corrected_alleles import parse_perfect_sam, parse_fasta
from parse_contig_realign import parse_CIGAR
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fa', '--fn_list_alleles',
        help = 'input allele fasta file list'
    )
    
    parser.add_argument(
        '-fof', '--fo_merged_fasta',
        help = 'output merged fasta file'
    )
    parser.add_argument(
        '-for', '--fo_merged_report',
        help = 'output merging report'
    )
    args = parser.parse_args()
    return args

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)



if __name__ == '__main__':
    args = parse_args()
    fn_list_alleles  = args.fn_list_alleles
    fo_merged_fasta  = args.fo_merged_fasta
    fo_merged_report = args.fo_merged_report
    print(fn_list_alleles)

    # dict_database {}
    # - keys: allele_name
    # - values: dict_SEQ {}
    #           - keys: SEQ
    #           - values: [person_name_0, person_name_1, ...]
    dict_database = {}
    for file_name in fn_list_alleles.split():
        dict_allele = parse_fasta(file_name)
        name_st = max(file_name.find("novel_") + 6, file_name.find("flanking_") + 9)
        name_ed = file_name.find(".fasta")
        person_name = file_name[name_st:name_ed]
        for allele_name, SEQ in dict_allele.items():
            allele_name = allele_name[:allele_name.rfind('-')]
            SEQ = SEQ.lower()
            if dict_database.get(allele_name):
                dict_SEQ = dict_database[allele_name]
                if dict_SEQ.get(SEQ):
                    dict_SEQ[SEQ].append(person_name)
                elif dict_SEQ.get(get_reverse_complement(SEQ)):
                    dict_SEQ[get_reverse_complement(SEQ)].append(person_name)
                else:
                    dict_SEQ[SEQ] = [person_name]
            else:
                dict_database[allele_name] = {SEQ:[person_name]}

    f_of = open(fo_merged_fasta,  'w')
    f_or = open(fo_merged_report, 'w')
    for allele_name, dict_SEQ in sorted(dict_database.items()):
        if allele_name == "":
            continue
        for idx, (SEQ, list_person) in enumerate(dict_SEQ.items()):
            f_of.write(">" + allele_name + '-' + str(idx) + '\n')
            f_of.write(SEQ + '\n')

            f_or.write(allele_name + '-' + str(idx) + ',' + str(len(list_person)) + ',' + ','.join(list_person) + '\n' )
    f_of.close()
    f_or.close()


