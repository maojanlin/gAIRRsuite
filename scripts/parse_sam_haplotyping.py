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
        '-fp', '--target_position', nargs='+', 
        help = 'target haplotyping position'
    )
    parser.add_argument(
        '-fo', '--fn_output_file',
        help = 'output report file'
    )
    args = parser.parse_args()
    return args


def haplotyping(fn_sam, target_position):
    dict_operate = {'M': 0, 'I': 1, 'D': -1, 'S':1, 'H':1}
    histogram = {}
    with open(fn_sam, 'r') as f_r:
        match_flag = False
        for line in f_r:
            fields = line.split('\t')
            cigar = fields[5] #cigar code
            number = []
            operate = []
            tmp_num = 0
            
            #parse cigar code
            for char in cigar:
                if char.isdigit():
                    # since cigar will always starts in number
                    tmp_num = tmp_num*10 + int(char)
                else:
                    # since cigar will always ends in operator
                    number.append(tmp_num)
                    tmp_num = 0
                    operate.append(char)

            #print(number)
            #print(operate)
            #print(cigar)

            list_split = [sum(number[:i+1]) for i in range(len(number))]
            list_split.insert(0,0)
            query_SEQ = fields[9]
            start_pos = int(fields[3])
            align_SEQ = fields[9]
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
                        align_SEQ = align_SEQ[:split_pos] + 'X'*op_num + align_SEQ[split_pos:]
                except KeyError:
                    print("Warning! starange cigar code " + op)
            
            haplotype = ""
            try:
                for position in target_position:
                    haplotype += (align_SEQ[position - start_pos])
            except IndexError:
                pass
            
            if histogram.get(haplotype):
                histogram[haplotype] += 1
            else:
                histogram[haplotype] = 1

    return(histogram)


if __name__ == '__main__':
    args = parse_args()
    fn_sam = args.fn_sam
    target_position = list(map(int, args.target_position))
    fn_output_file = args.fn_output_file

    histogram = haplotyping(fn_sam, target_position)
    if fn_output_file:
        with open(fn_output_file,'w') as f_o:
            list_haplotype = sorted(histogram)
            for haplotype in list_haplotype:
                f_o.write(haplotype + '\t' + str(histogram[haplotype]) + '\n')
    else:
        list_haplotype = sorted(histogram)
        #list_haplotype = np.sort(histogram.keys())
        for haplotype in list_haplotype:
            print(haplotype + '\t' + str(histogram[haplotype]))

