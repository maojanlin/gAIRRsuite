import argparse
import pickle
import os
import numpy as np
from parse_sam_haplotyping import parse_CIGAR

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fs', '--fn_sam',
        help = 'input sam file'
    )
    parser.add_argument(
        '-td', '--thrsd', type = int,
        help = 'edit distance threshold'
    )
    parser.add_argument(
        '-fo', '--fn_output_file',
        help = 'output report file'
    )
    args = parser.parse_args()
    return args

def position_allele_contig(fn_sam, thrsd=0):
    list_allele_contig = []
    # dict_locus {}
    #   - keys: "contig_names"
    #   - values: {a set of region pairs(1300,1600)}
    dict_locus = {}
    dict_operate = {'M': 1, 'I': 0, 'D': 1, 'S':0, 'H':0, '*':0}
    with open(fn_sam, 'r') as f_o:
        for line in f_o:
            if line[0] != '@': # real alignment information
                fields = line.split()
                cigar = fields[5]
                number, operate = parse_CIGAR(cigar)
                op_val = [dict_operate[op] for op in operate]

                start_pos = int(fields[3]) 
                end_pos = start_pos + np.dot(op_val, number)

                contig_name = fields[2]
                allele_name = fields[0].split('|')
                if len(allele_name) > 1:
                    allele_name = allele_name[1]
                else:
                    allele_name = allele_name[0]
                
                list_allele_contig.append((contig_name, str(start_pos)+":"+str(end_pos), allele_name))
                if dict_locus.get(contig_name):
                    dict_locus[contig_name].add((start_pos,end_pos))
                else:
                    dict_locus[contig_name] = {(start_pos,end_pos)}

    return sorted(list_allele_contig), dict_locus



if __name__ == '__main__':
    args = parse_args()
    fn_sam = args.fn_sam
    thrsd = args.thrsd
    fn_output_file = args.fn_output_file

    list_allele_contig, dict_locus = position_allele_contig(fn_sam, thrsd);
    pre_contig_name = ""
    pre_interval = ""
    for idx, element in enumerate(list_allele_contig):
        contig_name, interval, allele_name = element
        if contig_name == pre_contig_name:
            contig_name = " "*len(contig_name)
        else:
            pre_contig_name = contig_name
        if interval == pre_interval:
            interval = " "*len(interval)
        else:
            pre_interval = interval

        print (contig_name + '\t' + interval + '\t' + allele_name)

    if args.fn_output_file:
        f = open(fn_output_file, 'wb')
        pickle.dump(dict_locus, f)
        f.close()
    

