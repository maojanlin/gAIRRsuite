import argparse
import pickle
import os
import numpy as np
from parse_sam_haplotyping import parse_CIGAR

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fs1', '--fn_sam_H1',
        help = 'input H1 sam file'
    )
    parser.add_argument(
        '-fs2', '--fn_sam_H2',
        help = 'input H2 sam file'
    )
    parser.add_argument(
        '-td', '--thrsd', type = int,
        help = 'edit distance threshold'
    )
    parser.add_argument(
        '-fo1', '--fn_output_file_H1',
        help = 'output H1 pickle file'
    )
    parser.add_argument(
        '-fo2', '--fn_output_file_H2',
        help = 'output H2 pickle file H2'
    )
    args = parser.parse_args()
    return args

def position_allele_contig(fn_sam):
    list_allele_contig = []

    dict_operate = {'M': 1, 'I': 0, 'D': 1, 'S':0, 'H':0, '*':0}
    with open(fn_sam, 'r') as f_o:
        for line in f_o:
            if line[0] != '@': # real alignment information
                fields = line.split()
                cigar = fields[5]
                number, operate = parse_CIGAR(cigar)
                op_val = [dict_operate[op] for op in operate]

                start_pos = int(fields[3]) 
                end_pos   = start_pos + np.dot(op_val, number)
                
                edit_num  = int(fields[11].split(':')[2])

                contig_name = fields[2]
                allele_name = fields[0].split('|')
                if len(allele_name) > 1:
                    allele_name = allele_name[1]
                else:
                    allele_name = allele_name[0]
                
                list_allele_contig.append((allele_name, contig_name, start_pos, end_pos, edit_num))

    return sorted(list_allele_contig)


def add_dict(dict_locus, items):
    # dict_locus {}
    #   - keys: "contig_names"
    #   - values: {a set of region pairs(1300,1600)}
    contig_name = items[1]
    start_pos = items[2]
    end_pos   = items[3]
    edit_num  = items[4]
    if dict_locus.get(contig_name):
        dict_locus[contig_name].add((start_pos,end_pos,edit_num))
    else:
        dict_locus[contig_name] = {(start_pos,end_pos,edit_num)}
    return dict_locus


def locus_thred(list_allele_contig_H1, list_allele_contig_H2, thrsd=50):
    # dict_locus {}
    #   - keys: "contig_names"
    #   - values: {a set of region pairs(1300,1600)}
    dict_locus_H1 = {}
    dict_locus_H2 = {}

    for idx, items_1 in enumerate(list_allele_contig_H1):
        items_2 = list_allele_contig_H2[idx]
        determine = ""
        #print (items_1[0] + '\tH1: ' + str(items_1[4]) + '\tH2:' + str(items_2[4]))
        if items_2[4] - items_1[4] > thrsd:
            dict_locus_H1 = add_dict(dict_locus_H1, items_1)
            determine = 'keep H1, discard ' + items_2[1] + '\t' + str(items_2[2]) + ':' + str(items_2[3])
        elif items_1[4] - items_2[4] > thrsd:
            dict_locus_H2 = add_dict(dict_locus_H2, items_2)
            determine = 'keep H2, discard ' + items_1[1] + '\t' + str(items_1[2]) + ':' + str(items_1[3])
        else:
            dict_locus_H1 = add_dict(dict_locus_H1, items_1)
            dict_locus_H2 = add_dict(dict_locus_H2, items_2)
            determine = 'keep BOTH'
        print (items_1[0] + '\t' +  str(items_1[4]) + '\t' + str(items_2[4]) + '\t' + determine)

    return dict_locus_H1, dict_locus_H2


def print_list_allele_contig(list_allele_contig):
    list_contig = sorted(list_allele_contig, key=lambda element: element[1])
    pre_contig_name = ""
    pre_interval = ""
    for idx, element in enumerate(list_contig):
        allele_name = element[0]
        contig_name = element[1]
        interval = str(element[2]) + ':' + str(element[3])
        edit_num = str(element[4])
        if contig_name == pre_contig_name:
            contig_name = " "*len(contig_name)
        else:
            pre_contig_name = contig_name
        if interval == pre_interval:
            interval = " "*len(interval)
        else:
            pre_interval = interval

        print (contig_name + '\t' + interval + '\t' + allele_name + '\t' + edit_num)


if __name__ == '__main__':
    args = parse_args()
    fn_sam_H1 = args.fn_sam_H1
    fn_sam_H2 = args.fn_sam_H2
    thrsd = args.thrsd
    fn_output_file_H1 = args.fn_output_file_H1
    fn_output_file_H2 = args.fn_output_file_H2

    list_allele_contig_H1 = position_allele_contig(fn_sam_H1);
    list_allele_contig_H2 = position_allele_contig(fn_sam_H2);
    
    dict_locus_H1, dict_locus_H2 = locus_thred(list_allele_contig_H1, list_allele_contig_H2, thrsd)

    #print('H1:')
    #print_list_allele_contig(list_allele_contig_H1)
    #print('H2:')
    #print_list_allele_contig(list_allele_contig_H2)

    if args.fn_output_file_H1:
        f = open(fn_output_file_H1, 'wb')
        pickle.dump(dict_locus_H1, f)
        f.close()
    if args.fn_output_file_H2:
        f = open(fn_output_file_H2, 'wb')
        pickle.dump(dict_locus_H2, f)
        f.close()
    

