'''
parse_contig_realign.py deal with both the CIGAR and the MD tag to see where is the edit region
Sample input sam:

@SQ	SN:NODE_1_length_1979_cov_43.519438	LN:1979
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1198-dirty	CL:../../../bwa/bwa mem ../asm_contigs/TCRV_contigs_225.fasta TCRV_remain_225_P1.fasta
M01931:79:000000000-CHMNT:1:1101:10805:2092	0	NODE_1_length_1979_cov_43.519438	875	60	301M	*	0	0	GGTATCGACAAGACCCAGGCATGGGGCTGAGGCTGATTCATTACTCAGTTGGTGAGGGTACAACTGCCAAAGGAGAGGTCCCTGATGGCTACAATGTCTCCAGATTAAAAAAACAGAATTTCCTGCTGGGGTTGGAGTCGGCTGCTCCCTCCCAAACATCTGTGTACTTCTGTGCCAGCAGTTACTCCACAGTGCTGCACGGCTGTCTCCTCTCTGCACAGAAAGGCAAGGGAAGGTGCTGCCCTCCCCCGCAGCACAGATTCAGCGATGCCCTTGGTCCTAGCACCGGAAACTTTGGAGC	*	NM:i:2	MD:Z:247T40A12	AS:i:291	XS:i:0
M01931:79:000000000-CHMNT:1:1101:13808:15579	0	NODE_1_length_1979_cov_43.519438	955	60	301M	*	0	0	CCTGATGGCTACAATGTCTCCAGATTAAAAAAACAGAATTTCCTGCTGGGGTTGGAGTCGGCTGCTCCCTCCCAAACATCTGTGTACTTCTGTGCCAGCAGTTACTCCACAGTGCTGCACGGCTGTCTCCTCTCTGCACAGAAAGGCAAGGGAAGGTGCTGCCCTCCTCCGCAGCACAGATTCAGCGATGCCCTTGGTCCTAGCACCGAAAACTTTGGAGCCCCAATGGGCCCGGGCAGTGCGAGCCTTCATCTGTGCCAGGTGCCTCTGCAGTCGGTCTCGGCCAGGCCTGGATCGGTCC	*	NM:i:0	MD:Z:301	AS:i:301	XS:i:0
M01931:79:000000000-CHMNT:1:1101:14223:20821	0	NODE_1_length_1979_cov_43.519438	902	60	301M	*	0	0	TGAGGCTGATTCATTACTCAGTTGGTGAGGGTACAACTGCCAAAGGAGAGGTCCCTGATGGCTACAATGTCTCCAGATTAAAAAAACAGAATTTCCTGCTGGGGTTGGAGTCGGCTGCTCCCTCCCAAACATCTGTGTACTTCTGTGCCAGCAGTTACTCCACAGTGCTGCACGGCTGTCTCCTCTCTGCACAGAAAGGCAAGGGAAGGTTCTGCCCTCCTCCGCAGCACAGATTCAGCGATGCCCTTGGTCCTAGCACCGAAAACTTTGGAGCCCCAATGGGCCCGGGCAGTGCGAACCT	*	NM:i:2	MD:Z:210G86G3	AS:i:292	XS:i:0
'''
import argparse
import pickle
import os
import numpy as np
from parse_sam_haplotyping import parse_CIGAR
np.set_printoptions(threshold=5000)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fs', '--fn_sam',
        help = 'sam file of reads realign to contig'
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


def parse_MD(MD_tag):
    tmp_num = 0
    ref_idx = 1 # the sam files are started with 1
    mis_region = []
    for char in MD_tag:
        if char.isdigit():
            tmp_num = tmp_num*10 + int(char)
        else:
            ref_idx += tmp_num
            tmp_num = 0
            if char != '^':
                mis_region.append(ref_idx)
                ref_idx += 1
    return(mis_region)


def mark_edit_region(fn_sam, fn_output_file):
    edit_histogram = None
    cov_histogram  = None
    list_match = []
    list_edit  = []
    contig_len = 0
    dict_reads = {}
    even_odd_flag = 1
    with open(fn_sam, 'r') as f_s:
        for line in f_s:
            if line[0] == '@': # header, information of the contig
                if line.find('LN:') != -1:
                    contig_len = int(line[line.find('LN:')+3:-1]) + 1 # the number system start with 1
                    edit_histogram = np.zeros(contig_len)
                    cov_histogram  = np.zeros(contig_len)
            else: # real alignment information
                fields    = line.split()
                read_name = fields[0]
                read_SEQ  = fields[9]
                edit_dist = int(fields[11].split(':')[2])
                cigar     = fields[5]
                MD_tag    = fields[12].split(':')[2]
                start_pos = int(fields[3])
                
                number, operate = parse_CIGAR(cigar)
                mis_region_MD = parse_MD(MD_tag)
                if operate[0] == 'S':
                    mis_region_MD = [ele + number[0] + start_pos - 1 for ele in mis_region_MD]
                else:
                    mis_region_MD = [ele + start_pos - 1 for ele in mis_region_MD]

                mis_region_I = []   # insertion boundary region
                diff_len = 0        # len contribution of D and I
                if 'I' in operate or 'D' in operate:
                    idx_I = start_pos - 1
                    for idx, op in enumerate(operate):
                        if op == 'I':
                            diff_len -= number[idx]
                            mis_region_I.append(idx_I)
                            mis_region_I.append(idx_I+1)
                        else:
                            idx_I += number[idx]
                            if op == 'D':
                                diff_len += number[idx]
                
                edit_histogram[mis_region_MD] += 1
                edit_histogram[mis_region_I]  += 1
                
                end_pos   = int(fields[3]) + len(fields[9]) + diff_len
                cov_histogram[start_pos:end_pos] += 1
                
                # record the reads information
                dict_reads[(read_name, even_odd_flag)] = read_SEQ
                if edit_dist == 0:
                    list_match.append((start_pos, end_pos, read_name, even_odd_flag, read_SEQ))
                else:
                    list_edit.append((start_pos, end_pos, read_name, even_odd_flag, read_SEQ))
                if even_odd_flag == 1:
                    even_odd_flag = 2
                else:
                    even_odd_flag = 1

    return edit_histogram, cov_histogram, list_match, list_edit, dict_reads


def pop_perfect_reads(edit_region, list_match, dict_reads, fn_output_file):
    for read_info in list_match:
        start_pos, end_pos, read_name, even_odd_flag, read_SEQ = read_info
        for spot in edit_region:
            # if perfect match cover the edit spot
            if start_pos <= spot and end_pos >= spot:
                # if the read and its pair is in dict_read
                if dict_reads.get((read_name, even_odd_flag)):
                    # pop the reads
                    dict_reads.pop((read_name, 1))
                    dict_reads.pop((read_name, 2))
                break

    # output the pair end reads
    list_read_pair = sorted(dict_reads)
    if fn_output_file:
        fn_o_name_1 = ""
        fn_o_name_2 = ""
        if fn_output_file.find('.') == -1:
            fn_o_name_1 = fn_output_file + "_P1.fasta"
            fn_o_name_2 = fn_output_file + "_P2.fasta"
        else:
            fn_o_name_1 = fn_output_file.split('.')[0] + "_P1.fasta"
            fn_o_name_2 = fn_output_file.split('.')[0] + "_P2.fasta"
        
        fn_o_1 = open(fn_o_name_1, 'w')
        fn_o_2 = open(fn_o_name_2, 'w')
        for read_pair in list_read_pair:
            if read_pair[1] == 1:
                fn_o_1.write('>' + read_pair[0] + ' 1\n')
                fn_o_1.write(dict_reads[read_pair]+'\n')
            elif read_pair[1] == 2:
                fn_o_2.write('>' + read_pair[0] + ' 2\n')
                fn_o_2.write(dict_reads[read_pair]+'\n')
            else:
                print("Warning!! read pair doesn't exist!!")
        fn_o_1.close()
        fn_o_2.close()



if __name__ == '__main__':
    args = parse_args()
    fn_sam = args.fn_sam
    thrsd = args.thrsd
    fn_output_file = args.fn_output_file

    #parse the sam file and generate
    edit_histogram, cov_histogram, list_match, list_edit, dict_reads = mark_edit_region(fn_sam, fn_output_file)
    
    #determine the region contains alternative flanking region
    edit_region = []
    for idx, ele in enumerate(edit_histogram):
        print(str(idx) + ':\t' + str(cov_histogram[idx])  + '\t' + str(ele))
        if ele > cov_histogram[idx]/4:
            edit_region.append(idx)

    #print("edit_region")
    #print(edit_region)
    #print(list_match)
    #print(list_edit)

    pop_perfect_reads(edit_region, list_match, dict_reads, fn_output_file)

    

