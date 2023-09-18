import argparse
import pickle
import os
import numpy as np
from .utils import eprint, get_reverse_complement
from .filter_corrected_alleles import parse_fasta

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fab', '--fn_annotation_base',
        help = 'base reliable annotation'
    )
    parser.add_argument(
        '-fan', '--fn_annotation_new',
        help = 'new called annotation'
    )
    parser.add_argument(
        '-fh', '--fn_haplotypes',
        help = 'new sequences fasta file'
    )

    parser.add_argument(
        '-for', '--fo_report',
        help = 'comparison report'
    )
    parser.add_argument(
        '-fos', '--fo_summary',
        help = 'comparison summary'
    )
    args = parser.parse_args()
    return args


def parse_annotation(fn_annotation):
    list_annotation = []
    f_a = open(fn_annotation, 'r')
    for line in f_a:
        list_annotation.append(line.strip().split(','))
    return list_annotation


if __name__ == "__main__":
    args = parse_args()
    fn_annotation_base = args.fn_annotation_base
    fn_annotation_new  = args.fn_annotation_new
    fn_haplotypes = args.fn_haplotypes
    fo_report  = args.fo_report
    fo_summary = args.fo_summary

    list_annotation_base = parse_annotation(fn_annotation_base)
    list_annotation_new  = parse_annotation(fn_annotation_new)

    list_annotate_only_info = []
    list_captured_only_info = []
    
    idx_new = 0
    used_seqs = set()
    f_or = open(fo_report, 'w')
    for idx_base in range(len(list_annotation_base)):
        fields = list_annotation_base[idx_base]
        if len(fields) < 2:
            continue
        allele_name = fields[0]
        contig_name = fields[1]
        start_pos = int(fields[2])
        end_pos   = int(fields[3]) + start_pos
        allele_len = end_pos - start_pos
        miss_flag = True
        for idx_search in range(idx_new, len(list_annotation_new)):
            new_fields = list_annotation_new[idx_search]
            if len(new_fields) < 2:
                continue
            search_name = new_fields[0]
            used_seqs.add(search_name)
            search_contig = new_fields[1]
            new_start_pos = int(new_fields[2])
            new_end_pos   = int(new_fields[3]) + new_start_pos
            if contig_name == search_contig:
                if new_start_pos <= start_pos < end_pos <= new_end_pos:
                    miss_flag = False
                    l_flank_len = start_pos - new_start_pos
                    r_flank_len = new_end_pos - end_pos
                    mismatch_tag = ""
                    if len(new_fields) == 6:
                        mismatch_tag = ',' + new_fields[4] + ',' + new_fields[5]
                    elif len(new_fields) == 5:
                        mismatch_tag = ',' + new_fields[4]
                    
                    f_or.write(contig_name + ',' + str(start_pos) + ',' + str(allele_len) + \
                            ',' + str(l_flank_len) + ',' + str(r_flank_len) + ',' \
                            + allele_name + ',' + search_name + mismatch_tag + '\n')
                    
                    for idx_more in range(idx_new, idx_search):
                        more_fields = list_annotation_new[idx_more]
                        more_name = more_fields[0]
                        more_contig = more_fields[1]
                        more_start_pos = int(more_fields[2])
                        more_flank_len = int(more_fields[3])
                        mismatch_tag   = ""
                        if len(more_fields) == 6:
                            mismatch_tag = ',' + more_fields[4] + ',' + more_fields[5]
                        elif len(more_fields) == 5:
                            mismatch_tag = ',' + more_fields[4]
                        f_or.write(more_contig + ',*,*,' + str(more_start_pos) + ',' + str(more_flank_len) + ',' + more_name + ',' + mismatch_tag + '\n')
                        list_captured_only_info.append(more_contig + ',' + more_name + ',' + str(more_start_pos) + ',' + str(more_flank_len) + ',' + mismatch_tag)
                    idx_new = idx_search + 1
                    break
        if miss_flag:
            f_or.write(contig_name + ',' + str(start_pos) + ',' + str(allele_len) + ',*,*,' + allele_name + ',' + '*' + '\n')
            list_annotate_only_info.append(contig_name + ',' + allele_name + ',' + str(start_pos) + ',' + str(allele_len))
    f_or.close()

    f_os = open(fo_summary, 'w')
    if len(list_annotate_only_info) > 0:
        f_os.write("Annotation only position:\n")
        f_os.write("contig,allele_name,position,len,mismatch_tag\n")
        print("Annotation only position:")
        for info in list_annotate_only_info:
            f_os.write(info + '\n')
            print(info)
        f_os.write("\n")
    else:
        f_os.write("All annotated places are covered by AIRRCall.\n\n")
        print("All annotated place are covered by AIRRCall.\n")
    if len(list_captured_only_info) > 0:
        f_os.write("AIRRCall only position:\n")
        f_os.write("contig,allele_name,position,len,mismatch_tag\n")
        print("AIRRCall only position:")
        for info in list_captured_only_info:
            f_os.write(info + '\n')
            print(info)
        f_os.write("\n")
    else:
        f_os.write("No AIRRCall only places.\n\n")
        print("No AIRRCall only places.\n")

    if fn_haplotypes:
        dict_haplotypes = parse_fasta(fn_haplotypes)
        pool_seqs = set(dict_haplotypes.keys())
        f_os.write("Redundant AIRRCall seqs:\n")
        print("Redundant AIRRCall seqs:")
        for seq_name in sorted(pool_seqs - used_seqs):
            f_os.write(seq_name + '\n')
            print("\t" + seq_name)
    f_os.close()

    
    
