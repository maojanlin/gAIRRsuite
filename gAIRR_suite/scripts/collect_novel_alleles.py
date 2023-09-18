import argparse
import pickle
import os
import numpy as np
from .utils import get_reverse_complement
from .filter_corrected_alleles import parse_perfect_sam, parse_fasta
from .parse_contig_realign import parse_CIGAR
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fa', '--fn_list_alleles',
        help = 'input allele fasta file list'
    )
    parser.add_argument(
        '-type', '--merge_type',
        help = '\"simple\": without reference, \"reference\": with novel reference'
    )
    parser.add_argument(
        '-fnr', '--fn_novel_reference',
        help = 'input novel allele reference fasta file'
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
    merge_type = args.merge_type.lower()
    dict_novel_serial = {}
    list_novel_SEQ = None
    if merge_type == "reference":
        # dict_novel_serial {}
        # - keys: SEQ
        # - values: merged novel serial
        dict_temp = parse_fasta(args.fn_novel_reference)
        dict_novel_serial = {SEQ:name for (name,SEQ) in dict_temp.items()}
        list_novel_SEQ = sorted(dict_novel_serial.keys())
    elif merge_type != "simple":
        print("WARNING! Incorrect Merge Type:", merge_type)

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
        # get the personal name from file name
        name_st = max(file_name.find("novel_") + 6, file_name.find("flanking_") + 9)
        name_ed = file_name.find(".fasta")
        person_name = file_name[name_st:name_ed]
        for allele_name, SEQ in dict_allele.items():
            if allele_name == "":
                continue
            novel_allele_name = "" 
            if 'novel' in file_name:
                novel_allele_name = allele_name[:allele_name.rfind('/')]
                novel_allele_name += '/n'
            elif 'flanking' in file_name:
                if (merge_type == "reference") and ("novel" in allele_name):
                    #Find the same name in novel allele reference database
                    for ref_SEQ in list_novel_SEQ:
                        if ref_SEQ in SEQ:
                            novel_allele_name = dict_novel_serial[ref_SEQ]
                            novel_allele_name += '/f'
                else:
                    novel_allele_name = allele_name[:allele_name.rfind('-')]
                    novel_allele_name += '/f'
            else:
                print("WARNING! Incorrect naming in file", person_name)
            SEQ = SEQ.lower()
            if dict_database.get(novel_allele_name):
                dict_SEQ = dict_database[novel_allele_name]
                if dict_SEQ.get(SEQ):
                    dict_SEQ[SEQ].append(person_name)
                elif dict_SEQ.get(get_reverse_complement(SEQ)):
                    dict_SEQ[get_reverse_complement(SEQ)].append(person_name)
                else: # add the SEQ into dict_SEQ
                    dict_SEQ[SEQ] = [person_name]
            else:
                dict_database[novel_allele_name] = {SEQ:[person_name]}

    f_of = open(fo_merged_fasta,  'w')
    f_or = open(fo_merged_report, 'w')
    f_or.write('allele_name\tnumber_of_found_in_database\tsamples_possessing_the_allele\n')
    for allele_name, dict_SEQ in sorted(dict_database.items()):
        for idx, (SEQ, list_person) in enumerate(sorted(dict_SEQ.items(), key = lambda pair : len(pair[1]), reverse=True)):
            f_of.write(">" + allele_name + str(idx+1).zfill(2) + '\n')
            f_of.write(SEQ + '\n')

            f_or.write(allele_name + str(idx+1).zfill(2) + '\t' + str(len(list_person)) + '\t' + ','.join(list_person) + '\n' )
    f_of.close()
    f_or.close()


