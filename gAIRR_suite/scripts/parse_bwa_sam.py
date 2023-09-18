'''
Sample input sam:

TRAV38-1*02	0	NODE_8_length_2093_cov_53.250254	1044	60	280M	*	0	0	GCCCAGACAGTCACTCAGTCTCAACCAGAGATGTCTGTGCAGGAGGCAGAGACTGTGACCCTGAGTTGCACATATGACACCAGTGAGAATGATTATTATTTGTTCTGGTACAAGCAGCCTCCCAGCAGGCAGATGATTCTCGTTATTCGCCAAGAAGCTTATAAGCAACAGAATGCAACGGAGAATCGTTTCTCTGTGAACTTCCAGAAAGCAGCCAAATCCTTCAGTCTCAAGATCTCAGACTCACAGCTGGGGGACACTGCGATGTATTTCTGTGCTT	*	NM:i:1	MD:Z:90A189	AS:i:275	XS:i:237	XA:Z:NODE_21_length_1799_cov_58.263158,-772,280M,9;
TRAV14/DV4*01	16	NODE_141_length_563_cov_68.215596	108	25	290M	*	0	0	CCCTCTCTCATTGCACAGAAGTACATTGCTGAGTCCCCCAGTTGTGAAGCGGAGATGACAAGGTTGGCGGATTTTCTTGCCTTCTGGAAATTCAATGAGTAGCGACCTTCTGTTGCATTTTGCTGGTCATAAGACCCCTGATAAATAAGAAAAATCATTTCCCCACTGCTGGGCTGCTTGTACCAGAATAGACCATAACTTGGATCACTGGTGTCATATGTGCAGTCCAGAGTCACAGCCTCCTTTTCCTGCACGAACATTCCTGGTTGGGTTTGAGTTATCTTCTGGGC	*	NM:i:0	MD:Z:290	AS:i:290	XS:i:275	XA:Z:NODE_4_length_2191_cov_35.952035,+1061,290M,3;
TRAV14/DV4*03	16	NODE_141_length_563_cov_68.215596	116	5	282M	*	0	0	CATTGCACAGAAATACATTGCTGAGTCCCCCAGTTGTGAAGCGGAGATGACAAGGTTGGCGGATTTTCTTGCCTTCTGGAAATTCAATGAGTAGCGACCTTCTGTTGCATTTTGCTGGTCATAAGACCCCTGATAAATAAGAAAAATCATTTCCCCACTGCTGGGCTGCTTGTACCAGAATAGACCATAACTTGGATCACTGGTGTCATATGTGCAGTCCAGAGTCACAGCCTCCTTTTCCTGCACGAACATTCCTGGTTGGGTTTGAGTTATCTTCTGGGC	*	NM:i:1	MD:Z:12G269	AS:i:277	XS:i:272	XA:Z:NODE_4_length_2191_cov_35.952035,+1061,282M,2;
'''
import argparse
import pickle
import os
import numpy as np
import sys

# make sure the package modules is in the path
sys.path.append(os.path.dirname(__file__))

from parse_sam_haplotyping import parse_CIGAR
from utils import eprint, get_reverse_complement

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fs', '--fn_sam',
        help = 'input sam file'
    )
    parser.add_argument(
        '-fc', '--fn_contig',
        help = 'input contig file to call'
    )
    parser.add_argument(
        '-td', '--thrsd', type = int,
        help = 'edit distance threshold'
    )
    parser.add_argument(
        '--cluster_id',
        help = 'suffix after Nodes'
    )
    parser.add_argument(
        '-fo', '--fn_output_file',
        help = 'output report file'
    )
    parser.add_argument(
        '-fr', '--fn_output_flanking_region',
        help = 'output flanking region fasta file'
    )
    parser.add_argument(
        '-fsize', '--flanking_size', type = int,
        help = 'user-defined extending flanking size'
    )
    parser.add_argument(
        '-frs', '--fn_output_flanking_size',
        help = 'output flanking region fasta file extending only fsize'
    )
    args = parser.parse_args()
    return args

def fetch_contig(fn_contig):
    dict_contig = {}
    contig_name = ""
    contig_SEQ  = ""
    try:
        with open(fn_contig) as f_o:
            for line in f_o:
                if line[0] == '>':
                    if dict_contig.get(contig_name):
                        print("Warning! Duplicate contig name: " + contig_name)
                    else:
                        dict_contig[contig_name] = contig_SEQ
                    contig_name = line[1:-1]
                    contig_SEQ = ""
                else:
                    contig_SEQ += line.strip()
            dict_contig[contig_name] = contig_SEQ
    except:
        pass
    return dict_contig


def parse_edit_distance(fn_sam, fn_output_file, fn_output_flanking_region, fn_output_flanking_size, dict_contig, cluster_id, thrsd=0, flanking_size=100):
    f_report = open(fn_output_file, 'a')
    f_flank = open(fn_output_flanking_region, 'a')
    f_flank_size = open(fn_output_flanking_size, 'a')
    with open(fn_sam, 'r') as f_o:
        for line in f_o:
            if line[0] != '@': # real alignment information
                fields = line.split()
                #print(fields[11])
                eDist = int(fields[11].split(':')[2])
                cigar = fields[5]
                if 'S' in cigar:
                    continue
                contig_name = fields[2]
                if contig_name != '*' and eDist <= thrsd:
                    print_word = fields[0] + ' ' + contig_name.split('_')[2] + ' ' + contig_name + ' ' + str(cluster_id) + '\n'
                    #print_word = fields[0] + '\t' + fields[2] + '\t' + fields[11]
                    f_report.write(print_word)
                    
                    if dict_contig.get(contig_name):
                        contig_SEQ = dict_contig[contig_name]
                        allele_name = fields[0]
                        allele_print = allele_name + '_cluster_' + cluster_id
                        f_flank.write('>' + allele_print + '\n')
                        if (int(fields[1]) % 32) >= 16:
                            f_flank.write(get_reverse_complement(contig_SEQ) + '\n')
                            print(str(len(contig_SEQ) - int(fields[3]) - len(fields[9]) +1) + '-' + str(len(contig_SEQ) - int(fields[3]) +1) + ',' + allele_print)
                        else:
                            f_flank.write(contig_SEQ + '\n')
                            print(str(int(fields[3]) -1) + '-' + str(int(fields[3]) -1 + len(fields[9])) + ',' + allele_print)
                        
                        start_pos = int(fields[3]) -1 - flanking_size
                        end_pos   = int(fields[3]) -1 + len(fields[9]) + flanking_size
                        #print(str(start_pos) + '-' + str(end_pos))
                        if start_pos < 0:
                            start_pos = 0
                        if end_pos > len(contig_SEQ):
                            end_pos = len(contig_SEQ)
                        f_flank_size.write('>' + allele_print + '\n')
                        f_flank_size.write(contig_SEQ[start_pos:end_pos] + '\n')
                    else:
                        eprint("Warning! Contig name does not exist! " + contig_name)
    f_report.close()
    f_flank.close()
    f_flank_size.close()


if __name__ == '__main__':
    args = parse_args()
    fn_sam = args.fn_sam
    fn_contig = args.fn_contig
    thrsd = args.thrsd
    cluster_id = args.cluster_id
    fn_output_file = args.fn_output_file
    fn_output_flanking_region = args.fn_output_flanking_region
    flanking_size = args.flanking_size
    fn_output_flanking_size = args.fn_output_flanking_size

    dict_contig = fetch_contig(fn_contig)
    if len(dict_contig) > 0:
        parse_edit_distance(fn_sam, fn_output_file, fn_output_flanking_region, fn_output_flanking_size, dict_contig, cluster_id, thrsd, flanking_size)
    else:
        eprint("No contig to process.")


