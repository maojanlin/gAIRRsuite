'''
input format is like this:

M01931:79:000000000-CHMNT:1:1104:10725:24346	163	NA12878-S1544-H1-000005F	20885	60	301M	=	21219	622	AAGAAATGTCATGCATTGATGGTGGGAATTCAAAACGTTACATGCACAAAATGAGATTTTTTGGCATTTTTTAAATAGAGATAAAAGTAGAGTTAAAATGTGAACTTGTGCTTGTGTTCCAAAGTATTTACAACGTTGATTCAGAAATTGATGTTTACAAAGATTGATTCAGAGGAAGTTCTGTATCAGATTTGTTAATGTGATTCATTCTACAATCCCTGAAATTTGCTTGCAGAATAAATGTTGTATGAAAAATATCTCAAATAACTAAAATCCTGTCCACTCAAGCCCTTTTCCAGGG	CCCCCGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFFGFFFGGGGGGGGGBFFGGGGGGGGDGGGGGFGGGGGGGGGGFFGGGGGGFGGGGGGGGGGGGCEGGGGC;EGGGGGGGGFGGGGGGGGGGCGGCEG?DD@;EFGGGA==FFFG8FG,@=FDGCCFFGGF=:DFGGGFGFFGG?=FFDDFGD7D?9AFF88B7CFGFFFCFFC;88AFF+53>AFFFFEFFEEEBEB24*9@58>;)55>EFDB>).7@EFF@	NM:i:2	MD:Z:120G172G7	MC:Z:252M13I36M	AS:i:291	XS:i:231
M01931:79:000000000-CHMNT:1:2107:19461:4454	163	NA12878-S1544-H1-000005F	20885	60	301M	=	21168	584	AAGAAATGTCATGCATTGATGGTGGGAATTCAAAACGTTACATGCACAAAATGAGATTTTTTGGCATTTTTTAAATAGAGATAAAAGTAGAGTTAAAATGTGAACTTGTGCTTGTGTTCCGAAGTATTTACAACGTTGATTCAGAAATTGATGTTTACAAAGATTGATTCAGAGGAAGTTCTGTATCAGATTTGTTAATGTGATTCATTCTACAATCCCTGAAATTTGCTTGCAGAATAAATGTTGTATGAAAAATATCCCAAATAACTAAAATCCTGTCCACTCAAGCCCTTGTCCAGGG	CCCCCGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGDEGFG9FGGFGFFGGFFFFGACFGGEFGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGEEB;FDFGFFFGGFGGGGGGGGGGGGFGGGGGGGGGFGGAEA;9EGGGGFAFGGGGDGGGECFFGGGGGGGFFDFGGFFF9C7FF?CCCFF?8@F?7F7+96;+?@FEGGF6@+6AFF4EFA+16A3499+*4@6===07;EE2:@F(9=>?FFF	NM:i:1	MD:Z:259T41	MC:Z:301M	AS:i:296	XS:i:236
'''

import argparse
import pickle
import os
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fs', '--fn_sam',
        help = 'input sam file'
    )
    parser.add_argument(
        '-fc', '--chromosome_file', 
        help = 'dumped target chromosome file'
    )
    parser.add_argument(
        '-fas', '--allele_SEQ', 
        help = 'target allele_SEQ'
    )
    parser.add_argument(
        '--interval', required = True,
        help = '\"target_allele_start_position\"-\"target_allele_end_position\"'
    )
    parser.add_argument(
        '-fp', '--target_position', nargs='+', 
        help = 'target haplotyping position'
    )
    parser.add_argument(
        '-fo', '--fn_output_file',
        help = 'output report file'
    )
    args = parser.parse_args()
    return args

def parse_CIGAR(cigar):
    number = []
    operate = []
    tmp_num = 0
    for char in cigar:
        if char.isdigit():
            # since cigar will always starts in number
            tmp_num = tmp_num*10 + int(char)
        else:
            # since cigar will always ends in operator
            number.append(tmp_num)
            tmp_num = 0
            operate.append(char)
    
    return (number, operate)


def warp_SEQ(query_SEQ, cigar):
    # number of CIGAR operation and operator
    dict_operate = {'M': 0, 'I': 1, 'D': -1, 'S':1, 'H':1}
    number, operate = parse_CIGAR(cigar)
    
    list_split = [sum(number[:i+1]) for i in range(len(number))]
    list_split.insert(0,0)
    align_SEQ = query_SEQ
    for idx, op in enumerate(reversed(operate)):
        try:
            op_value = dict_operate[op]
            if op_value == 0:
                continue
            elif op_value == 1:
                op_num = number[-idx-1]
                split_pos = list_split[-idx-1]
                align_SEQ = align_SEQ[:split_pos] + align_SEQ[split_pos+op_num:]
            else:
                op_num = number[-idx-1]
                split_pos = list_split[-idx-1]
                align_SEQ = align_SEQ[:split_pos] + '_'*op_num + align_SEQ[split_pos:]
        except KeyError:
            print("Warning! strange cigar code " + op)

    return align_SEQ



def haplotyping(fn_sam, target_position, allele_start_position, allele_end_position):
    dict_bp = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    histogram = {}
    haplo_map = np.zeros(4*(allele_end_position-allele_start_position))
    with open(fn_sam, 'r') as f_r:
        match_flag = False
        for line in f_r:
            fields = line.split('\t')
            cigar = fields[5] #cigar code
            
            query_SEQ = fields[9]
            # warping the query_SEQ according to CIGAR code
            align_SEQ = warp_SEQ(query_SEQ, cigar)
            start_pos = int(fields[3])
            
            ##### check the base on the specific positions #####
            haplotype = ""
            try:
                for position in target_position:
                    if position - start_pos < 0:
                        haplotype += '_'
                    else:
                        haplotype += (align_SEQ[position - start_pos])
            except IndexError:
                haplotype += '_'
            
            if histogram.get(haplotype):
                histogram[haplotype] += 1
            else:
                histogram[haplotype] = 1
            
            ##### check the whole allele region #####
            for idx in range(allele_end_position-allele_start_position):
                try:
                    relative_idx = idx + allele_start_position - start_pos
                    if relative_idx > 0:
                        base = align_SEQ[relative_idx]
                        haplo_map[idx*4 + dict_bp[base]] += 1
                except KeyError:
                    pass
                except IndexError:
                    pass

    return(histogram, haplo_map)


def print_suspect_SNPs(interval_length, haplo_map, chromosome_SEQ, allele_SEQ):
    dict_r_bp = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    for idx in range(interval_length):
        calling = haplo_map[4*idx:4*idx+4]
        max_base = np.argmax(calling)
        print_flag = False
        print_chromosome = ""
        print_allele = ""
        r_base = dict_r_bp[max_base]
        if calling[max_base]/sum(calling) < 0.90:
            print_flag = True
        if chromosome_SEQ:
            c_base = chromosome_SEQ[idx]
            print_chromosome = ", contig base: " + c_base
            if c_base != r_base:
                print_flag = True
        if allele_SEQ:
            a_base = allele_SEQ[idx]
            print_allele = ", allele base: " + a_base
            if a_base != r_base:
                print_flag = True
        if print_flag:
            str_calling = r_base + "  |  "
            for idy, depth in enumerate(calling):
                str_calling += dict_r_bp[idy] + ": " + str(int(depth)).rjust(5) + " (" + str(int(depth/sum(calling)*100)) + "%)\t"
            print("position: " + str(allele_start_position + idx) + print_chromosome + print_allele + ", reads: " + str_calling)



if __name__ == '__main__':
    args = parse_args()
    fn_sam = args.fn_sam
    interval = args.interval.split('-')
    allele_start_position = int(interval[0])
    allele_end_position = int(interval[1])
    
    allele_SEQ = args.allele_SEQ
    chromosome_SEQ = None
    if args.chromosome_file:
        whole_chromosome_SEQ = ""
        with open(args.chromosome_file, 'r') as f_r:
            for line in f_r:
                if line[0] != '>':
                    whole_chromosome_SEQ += line.strip()
        chromosome_SEQ = whole_chromosome_SEQ[allele_start_position-1:allele_end_position-1]
    if args.target_position:
        target_position = list(map(int, args.target_position))
    else:
        target_position = []
    fn_output_file = args.fn_output_file

    histogram, haplo_map = haplotyping(fn_sam, target_position, allele_start_position, allele_end_position)
    if fn_output_file:
        with open(fn_output_file,'w') as f_o:
            list_haplotype = sorted(histogram)
            for haplotype in list_haplotype:
                f_o.write(haplotype + '\t' + str(histogram[haplotype]) + '\n')
    elif args.target_position:
        list_haplotype = sorted(histogram)
        #list_haplotype = np.sort(histogram.keys())
        for haplotype in list_haplotype:
            print(haplotype + '\t' + str(histogram[haplotype]))

    print_suspect_SNPs(allele_end_position - allele_start_position, haplo_map, chromosome_SEQ, allele_SEQ)
    

