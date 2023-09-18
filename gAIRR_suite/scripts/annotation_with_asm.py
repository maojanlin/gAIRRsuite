import argparse
import pickle
import os
import numpy as np
from scripts.filter_corrected_alleles import parse_perfect_sam, parse_fasta
from scripts.parse_contig_realign import parse_CIGAR
from scripts.utils import get_reverse_complement
import sys

def parse_args(arguments):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fs1', '--fn_sam_H1',
        help = 'input sam file where alleles.fasta align to corrected_alleles.fasta'
    )
    parser.add_argument(
        '-fs2', '--fn_sam_H2',
        help = 'input sam file where alleles.fasta align to corrected_alleles.fasta'
    )
    parser.add_argument(
        '-fasm1', '--fn_asm_H1',
        help = 'input assembly fasta file'
    )
    parser.add_argument(
        '-fasm2', '--fn_asm_H2',
        help = 'input assembly fasta file'
    )
    parser.add_argument(
        '-fp1', '--fn_pickle_H1',
        help = 'pickle summary of the assembly fasta file'
    )
    parser.add_argument(
        '-fp2', '--fn_pickle_H2',
        help = 'pickle summary of the assembly fasta file'
    )
    parser.add_argument(
        '-ext', '--len_extend', type=int, default=200,
        help = 'the flanking sequence length cropped from the contigs'
    )
    
    parser.add_argument(
        '-fos', '--fo_annotation_summary',
        help = 'output allele with perfect match relative to whole genome assembly'
    )
    parser.add_argument(
        '-foa', '--fo_perfect_annotation_report',
        help = 'output allele with perfect match relative to whole genome assembly'
    )
    parser.add_argument(
        '-foma', '--fo_mismatched_annotation_report',
        help = 'output allele with perfect match relative to whole genome assembly'
    )
    parser.add_argument(
        '-fom', '--fo_mismatched_fasta',
        help = 'output allele with perfect match relative to whole genome assembly'
    )
    parser.add_argument(
        '-fof', '--fo_flanking_fasta',
        help = 'output flanking sequence fasta file from asm'
    )
    args = parser.parse_args(arguments)
    return args

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def check_occupied(dict_occupied_place, list_fields, dict_contig):
    for fields in list_fields:
        allele_name = fields[0]
        allele_tag = int(fields[1])
        contig_name = fields[2]
        contig_len = len(dict_contig[contig_name])
        pos_start = int(fields[3])
        cigar = fields[5]
        if (contig_name =='*'):# or ('S' in cigar) or ('H' in cigar):
            continue
        # figure out how long the allele is on the reference
        align_length = 0
        clip_length  = 0
        number, operate = parse_CIGAR(cigar)
        for idx, op in enumerate(operate):
            if (op == 'M') or (op == 'D'):
                align_length += number[idx]
            elif (op == 'H') or (op == 'S'):
                clip_length += number[idx]
       
        # if the clip portion is too large, ignore the alignment
        if clip_length > 0.1*align_length: 
            continue
        
        if (operate[0] == 'S') or (operate[0] == 'H'):
            pos_start = pos_start - number[0]
        pos_end = pos_start + align_length + clip_length
        # trim the occupied region if they surpass the contig boundary
        if pos_start < 1:
            clip_length += pos_start
            pos_start = 1
        if pos_end > contig_len+1:
            clip_length = clip_length + contig_len+1 - pos_end
            pos_end = contig_len+1

        mismatch = int(fields[11][5:]) + clip_length # get NM:i: number
        if mismatch == 0:
            if dict_occupied_place.get(contig_name):
                dict_occupied_place[contig_name].append( (pos_start, pos_end, mismatch, allele_name, allele_tag) )
            else:
                dict_occupied_place[contig_name] = [ (pos_start, pos_end, mismatch, allele_name, allele_tag) ]
        else: # contain mismatches!
            if dict_occupied_place.get(contig_name):
                occupied_flag = False
                for idx, element in enumerate(dict_occupied_place[contig_name]): # check if the place is already occupied
                    e_start = element[0]
                    e_end = element[1]
                    #if (pos_end < e_start) == False and (e_end < pos_start) == False:
                    if min((pos_end - e_start), (e_end - pos_start)) > 0.8*min((pos_end-pos_start), (e_end-e_start)):
                        e_mismatch = element[2]
                        occupied_flag = True
                        if mismatch < e_mismatch:
                            MD_tag = fields[12][5:]
                            dict_occupied_place[contig_name][idx] = (pos_start, pos_end, mismatch, allele_name, allele_tag, cigar, MD_tag, clip_length)
                        break
                if occupied_flag == False:
                    MD_tag = fields[12][5:]
                    dict_occupied_place[contig_name].append( (pos_start, pos_end, mismatch, allele_name, allele_tag, cigar, MD_tag, clip_length) )
            else:
                MD_tag = fields[12][5:]
                dict_occupied_place[contig_name] = [ (pos_start, pos_end, mismatch, allele_name, allele_tag, cigar, MD_tag, clip_length) ]


def output_dict_occupied(dict_occupied_place, f_o):
    contig_num = 0
    allele_num = 0
    novel_num = 0
    for contig_name in sorted(dict_occupied_place.keys(), reverse=True):
        contig_num += 1
        for pairs in sorted(dict_occupied_place[contig_name], reverse=True):
            allele_num += 1
            mismatch_tag = ""
            if pairs[2] > 0:
                novel_num += 1
                mismatch_tag = ',NM:i:' + str(pairs[2]-pairs[7])
                if pairs[7]:
                    mismatch_tag += ',HS:i:' + str(pairs[7])
            
            try:
                f_o.write(pairs[3].split('|')[1] + ',' + contig_name + ',' + str(pairs[0]) + ',' + str(pairs[1]-pairs[0]) + mismatch_tag)
            except:
                f_o.write(pairs[3].split()[0] + ',' + contig_name + ',' + str(pairs[0]) + ',' + str(pairs[1]-pairs[0]) + mismatch_tag)
            f_o.write('\n')
    return (novel_num, allele_num, contig_num)


def occupied_annotation(dict_occupied_place_1, dict_occupied_place_2, list_fields_1, list_fields_2, fo_annotation_report, dict_contig_H1, dict_contig_H2, fo_annotation_summary):
    # dict_occupied_place {}
    #  - key: contig_name
    #  - value: list_occupied_pair[ (pos_start, pos_end, NM:i), ... ]
    
    check_occupied(dict_occupied_place_1, list_fields_1, dict_contig_H1)
    check_occupied(dict_occupied_place_2, list_fields_2, dict_contig_H2)
    
    f_o = open(fo_annotation_report, 'w')
    f_os = open(fo_annotation_summary, 'w')
    (novel_num, allele_num, contig_num) = output_dict_occupied(dict_occupied_place_1, f_o)
    f_o.write("There are " + str(novel_num) + " novel alleles in total " + str(allele_num) + " alleles from " + str(contig_num) + " contigs in H1.\n")
    f_os.write("There are " + str(novel_num) + " novel alleles in total " + str(allele_num) + " alleles from " + str(contig_num) + " contigs in H1.\n")
    (novel_num, allele_num, contig_num) = output_dict_occupied(dict_occupied_place_2, f_o)
    if dict_contig_H2 != {}:
        f_o.write("There are " + str(novel_num) + " novel alleles in total " + str(allele_num) + " alleles from " + str(contig_num) + " contigs in H2.\n")
        f_os.write("There are " + str(novel_num) + " novel alleles in total " + str(allele_num) + " alleles from " + str(contig_num) + " contigs in H2.\n")
    f_o.close()
    f_os.close()


def reference_recovery(allele_SEQ, cigar, MD_tag):
    if 'I' in cigar:
        number, operate = parse_CIGAR(cigar)
        allele_cursor = 0
        for idx, op in enumerate(operate):
            if op == 'M':
                allele_cursor += number[idx]
            elif op == 'H':
                allele_SEQ = allele_SEQ[:allele_cursor] + allele_SEQ[allele_cursor+number[idx]:]
                allele_cursor += number[idx]
    allele_SEQ = allele_SEQ.lower()

    MD_num = 0
    MD_D_SEQ = ""
    MD_D_flag = False
    MD_cursor = 0
    for char in MD_tag:
        if char.isdigit():
            MD_num = MD_num*10 + int(char)
            if MD_D_flag:
                allele_SEQ = allele_SEQ[:MD_cursor] + MD_D_SEQ + allele_SEQ[MD_cursor:]
                MD_D_SEQ  = ""
                MD_D_flag = False
        elif MD_D_flag == True:
            MD_D_SEQ += char
        elif char == '^':
            MD_cursor += MD_num
            MD_num = 0
            MD_D_flag = True
        else: # common substitution
            MD_cursor += MD_num
            MD_num = 0
            allele_SEQ = allele_SEQ[:MD_cursor] + char + allele_SEQ[MD_cursor+1:]
            MD_cursor += 1

    return allele_SEQ


def correct_allele(dict_occupied_place, dict_SEQ, dict_corrected_alleles, dict_flanking_alleles, dict_contig, len_extend):
    # dict_corrected_alleles {}
    #  - keys: allele_name
    #  - values: corrected_SEQ_set {corrected_SEQ_1, corrected_SEQ_2}
    
    # dict_flanking_alleles {}
    #  - keys: allele_name
    #  - values: corrected_SEQ_set {flanking_SEQ_1, flanking_SEQ_2}
    for contig_name, list_contig in dict_occupied_place.items():
        contig_SEQ = dict_contig[contig_name]
        contig_len = len(contig_SEQ)
        for pairs in list_contig:
            pos_start = pairs[0]
            pos_end   = pairs[1]
            allele_name = pairs[3]
            try:
                allele_name = allele_name.split('|')[1]
            except:
                allele_name = allele_name.split()[0]
            flag = pairs[4]
            flanking_SEQ = contig_SEQ[max(0, pos_start-len_extend-1):min(contig_len, pos_end+len_extend-1)].lower()
            if pairs[2] != 0: # mismatched alleles but remain in contig
                corrected_SEQ = contig_SEQ[pos_start-1:pos_end-1].lower()
                if flag % 32 >= 16:
                    corrected_SEQ = get_reverse_complement(corrected_SEQ)
                    flanking_SEQ  = get_reverse_complement(flanking_SEQ)

                if dict_corrected_alleles.get(allele_name):
                    dict_corrected_alleles[allele_name].add(corrected_SEQ)
                else:
                    dict_corrected_alleles[allele_name] = {corrected_SEQ}
                
                flanking_name = allele_name + "/novel"
                if dict_flanking_alleles.get(flanking_name):
                    dict_flanking_alleles[flanking_name].add(flanking_SEQ)
                else:
                    dict_flanking_alleles[flanking_name] = {flanking_SEQ}
            else:
                if flag % 32 >= 16:
                    flanking_SEQ  = get_reverse_complement(flanking_SEQ)
                flanking_name = allele_name
                if dict_flanking_alleles.get(flanking_name):
                    dict_flanking_alleles[flanking_name].add(flanking_SEQ)
                else:
                    dict_flanking_alleles[flanking_name] = {flanking_SEQ}

    return dict_corrected_alleles, dict_flanking_alleles


def get_SEQ_from_sam_list(list_fields, dict_SEQ):
    for fields in list_fields:
        if fields[9] != '*':
            name = fields[0]
            flag = int(fields[1])
            if dict_SEQ.get(name):
                continue
            else:
                if flag % 32 >= 16: # SEQ is reverse_complemented
                    dict_SEQ[name] = get_reverse_complement(fields[9])
                else: # SEQ is original one
                    dict_SEQ[name] = fields[9]



def main(arguments=None):
    args = parse_args(arguments)
    # input information
    fn_sam_H1 = args.fn_sam_H1
    fn_sam_H2 = args.fn_sam_H2
    fn_asm_H1 = args.fn_asm_H1
    fn_asm_H2 = args.fn_asm_H2
    fn_pickle_H1 = args.fn_pickle_H1
    fn_pickle_H2 = args.fn_pickle_H2
    len_extend = args.len_extend

    # output files
    fo_annotation_summary = args.fo_annotation_summary
    fo_perfect_annotation_report = args.fo_perfect_annotation_report
    fo_mismatched_annotation_report = args.fo_mismatched_annotation_report
    fo_mismatched_fasta = args.fo_mismatched_fasta
    fo_flanking_fasta = args.fo_flanking_fasta

    if fn_sam_H2: # if there are two genome H1, H2 to analyze
        dict_contig_H1 = None
        dict_contig_H2 = None
        if os.path.exists(fn_pickle_H1):
            print ('Pickle file', fn_pickle_H1, 'has existed, load for it instead of re-fetching')
            f = open(fn_pickle_H1, 'rb')
            dict_contig_H1 = pickle.load(f)
            f.close()
        else:
            dict_contig_H1 = parse_fasta(fn_asm_H1)
            if fn_pickle_H1:
                f = open(fn_pickle_H1, 'wb')
                pickle.dump(dict_contig_H1, f)
                f.close()
        if os.path.exists(fn_pickle_H2):
            print ('Pickle file', fn_pickle_H2, 'has existed, load for it instead of re-fetching')
            f = open(fn_pickle_H2, 'rb')
            dict_contig_H2 = pickle.load(f)
            f.close()
        else:
            dict_contig_H2 = parse_fasta(fn_asm_H2)
            if fn_pickle_H2:
                f = open(fn_pickle_H2, 'wb')
                pickle.dump(dict_contig_H2, f)
                f.close()
        list_perfect_fields_1, list_mismatch_fields_1 = parse_perfect_sam(fn_sam_H1)
        list_perfect_fields_2, list_mismatch_fields_2 = parse_perfect_sam(fn_sam_H2)
        set_perfect_allele_1 = set(fields[0] for fields in list_perfect_fields_1)
        set_perfect_allele_2 = set(fields[0] for fields in list_perfect_fields_2)
        print("========== Annotation of IMGT Alleles ==========")
        print("There are", len(list_perfect_fields_1), "allele sites and", len(set_perfect_allele_1), "alleles in H1.")
        print("There are", len(list_perfect_fields_2), "allele sites and", len(set_perfect_allele_2), "alleles in H2.")
        print("There are", len(set_perfect_allele_1.intersection(set_perfect_allele_2)), "common alleles in H1 and H2.") 
    
        # output the perfect annotations
        dict_occupied_place_1 = {}
        dict_occupied_place_2 = {}
        occupied_annotation(dict_occupied_place_1, dict_occupied_place_2, list_perfect_fields_1, list_perfect_fields_2, \
                            fo_perfect_annotation_report, dict_contig_H1, dict_contig_H2, fo_annotation_summary)
        
        # output the mismatched annotations
        if fo_mismatched_annotation_report:
            occupied_annotation(dict_occupied_place_1, dict_occupied_place_2, list_mismatch_fields_1, list_mismatch_fields_2, \
                                fo_mismatched_annotation_report, dict_contig_H1, dict_contig_H2, fo_annotation_summary)
            print("========== Annotation of Imperfect Matches ==========")
            allele_num_H1 = sum([len(list_contig) for list_contig in dict_occupied_place_1.values()])
            print("There are", allele_num_H1, "potential alleles in H1 among", len(dict_occupied_place_1), "contigs.")
            allele_num_H2 = sum([len(list_contig) for list_contig in dict_occupied_place_2.values()])
            print("There are", allele_num_H2, "potential alleles in H2 among", len(dict_occupied_place_2), "contigs.")
    
            if fo_mismatched_fasta:
                dict_SEQ = {}
                get_SEQ_from_sam_list(list_perfect_fields_1, dict_SEQ)
                get_SEQ_from_sam_list(list_perfect_fields_2, dict_SEQ)
                get_SEQ_from_sam_list(list_mismatch_fields_1, dict_SEQ)
                get_SEQ_from_sam_list(list_mismatch_fields_2, dict_SEQ)
    
                dict_corrected_alleles = {}
                dict_flanking_alleles  = {}
                correct_allele(dict_occupied_place_1, dict_SEQ, dict_corrected_alleles, dict_flanking_alleles, dict_contig_H1, len_extend)
                correct_allele(dict_occupied_place_2, dict_SEQ, dict_corrected_alleles, dict_flanking_alleles, dict_contig_H2, len_extend)
                f_n = open(fo_mismatched_fasta, 'w')
                for allele_name in sorted(dict_corrected_alleles.keys()):
                    for idx, SEQ in enumerate(sorted(dict_corrected_alleles[allele_name])):
                        f_n.write('>' + allele_name + '/novel-' + str(idx) + '\n')
                        f_n.write(SEQ + '\n')
                f_n.close()
                print("Output novel alleles.")
                f_f = open(fo_flanking_fasta, 'w')
                for allele_name in sorted(dict_flanking_alleles.keys()):
                    for idx, SEQ in enumerate(sorted(dict_flanking_alleles[allele_name])):
                        f_f.write('>' + allele_name + '-' + str(idx) + '\n')
                        f_f.write(SEQ + '\n')
                f_f.close()
                print("Output flanking sequences")
            else:
                print("Corrected mismatched files not specified.")
    else: # if there is only one genome to analyze
        dict_contig_H1 = None
        if os.path.exists(fn_pickle_H1):
            print ('Pickle file', fn_pickle_H1, 'has existed, load for it instead of re-fetching')
            f = open(fn_pickle_H1, 'rb')
            dict_contig_H1 = pickle.load(f)
            f.close()
        else:
            dict_contig_H1 = parse_fasta(fn_asm_H1)
            if fn_pickle_H1:
                f = open(fn_pickle_H1, 'wb')
                pickle.dump(dict_contig_H1, f)
                f.close()
        list_perfect_fields_1, list_mismatch_fields_1 = parse_perfect_sam(fn_sam_H1)
        set_perfect_allele_1 = set(fields[0] for fields in list_perfect_fields_1)
        print("========== Annotation of IMGT Alleles ==========")
        print("There are", len(list_perfect_fields_1), "allele sites and", len(set_perfect_allele_1), "alleles in genome.")
    
        # output the perfect annotations
        dict_occupied_place_1 = {}
        occupied_annotation(dict_occupied_place_1, {}, list_perfect_fields_1, [], fo_perfect_annotation_report, dict_contig_H1, {}, fo_annotation_summary)
        
        # output the mismatched annotations
        if fo_mismatched_annotation_report:
            occupied_annotation(dict_occupied_place_1, {}, list_mismatch_fields_1, [], fo_mismatched_annotation_report, dict_contig_H1, {}, fo_annotation_summary)
            print("========== Annotation of Imperfect Matches ==========")
            allele_num_H1 = sum([len(list_contig) for list_contig in dict_occupied_place_1.values()])
            print("There are", allele_num_H1, "potential alleles in the genome among", len(dict_occupied_place_1), "contigs.")
    
            if fo_mismatched_fasta:
                dict_SEQ = {}
                get_SEQ_from_sam_list(list_perfect_fields_1, dict_SEQ)
                get_SEQ_from_sam_list(list_mismatch_fields_1, dict_SEQ)
    
                dict_corrected_alleles = {}
                dict_flanking_alleles  = {}
                correct_allele(dict_occupied_place_1, dict_SEQ, dict_corrected_alleles, dict_flanking_alleles, dict_contig_H1, len_extend)
                f_n = open(fo_mismatched_fasta, 'w')
                for allele_name in sorted(dict_corrected_alleles.keys()):
                    for idx, SEQ in enumerate(sorted(dict_corrected_alleles[allele_name])):
                        f_n.write('>' + allele_name + '/novel-' + str(idx) + '\n')
                        f_n.write(SEQ + '\n')
                f_n.close()
                print("Output novel alleles.")
                f_f = open(fo_flanking_fasta, 'w')
                for allele_name in sorted(dict_flanking_alleles.keys()):
                    for idx, SEQ in enumerate(sorted(dict_flanking_alleles[allele_name])):
                        f_f.write('>' + allele_name + '-' + str(idx) + '\n')
                        f_f.write(SEQ + '\n')
                f_f.close()
                print("Output flanking sequences")
            else:
                print("Corrected mismatched files not specified.")
            


if __name__ == "__main__":
    main()
