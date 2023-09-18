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

from scripts.parse_sam_haplotyping import parse_CIGAR
from scripts.utils import get_reverse_complement
from scripts.coverage_analysis import eprint
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
        '-it_region', '--interest_region',
        help = '\"start-end\" the region interested, only variants inside are considered'
    )
    parser.add_argument(
        '-fo', '--fn_output_file',
        help = 'output report file'
    )
    
    parser.add_argument(
        '--contig_file',
        help = 'fasta file showing the contig information'
    )
    parser.add_argument(
        '--allele_file',
        help = 'fasta file containing the allele name information'
    )
    parser.add_argument(
        '--corrected_contig_output_file',
        help = 'output corrected contig as flanking region fasta file'
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


def mark_edit_region(fn_sam, fn_output_file, contig_file):
    edit_histogram = None
    cov_histogram  = None
    #list_read_info: [ (start_pos, end_pos, read_name, even_odd_flag, mis_region) ]
    list_read_info = []
    contig_len = 0
    contig_name = ""
    # dict_reads{}
    #  - key: (read_name, pair_number)
    #  - values: read_SEQ
    dict_reads = {}
    even_odd_flag = 1
    with open(fn_sam, 'r') as f_s:
        for line in f_s:
            if line[0] == '@': # header, information of the contig
                if line.find('LN:') != -1:
                    # sometimes SPAdes would produce more than 1 contig, but the short one are not very useful
                    # so we discard the short contigs and reads align to them
                    if contig_len == 0: 
                        contig_len = int(line[line.find('LN:')+3:-1]) + 1 # the number system start with 1
                        contig_name = line.split(':')[1][:-3] 
                        edit_histogram = np.zeros(contig_len)
                        cov_histogram  = np.zeros(contig_len)
            else: # real alignment information
                fields    = line.split()
                # if the read align to shorter contigs, pass
                if contig_name != fields[2]: 
                    dict_reads[(read_name, even_odd_flag)] = read_SEQ
                    list_read_info.append((0, 0, read_name, even_odd_flag, [], "", read_SEQ))
                    if even_odd_flag == 1:
                        even_odd_flag = 2
                    else:
                        even_odd_flag = 1
                    continue
                read_name = fields[0]
                read_SEQ  = fields[9]
                cigar     = fields[5]
                sam_flag  = int(fields[1])
                # if the alignment is a supplementary alignment, pass
                # read BWA manual "Supplementary Alignment" for more information
                if sam_flag > 1024:
                    continue
                # if cigar == '*', means alignment is bad, pass
                if cigar == '*':
                    dict_reads[(read_name, even_odd_flag)] = read_SEQ
                    #list_read_info.append((start_pos, end_pos, read_name, even_odd_flag, mis_region))
                    list_read_info.append((0, 0, read_name, even_odd_flag, [], "", read_SEQ))
                    if even_odd_flag == 1:
                        even_odd_flag = 2
                    else:
                        even_odd_flag = 1
                    continue

                edit_dist = int(fields[11].split(':')[2])
                MD_tag    = fields[12].split(':')[2]
                start_pos = int(fields[3])
                
                number, operate = parse_CIGAR(cigar)
                mis_region_MD = parse_MD(MD_tag)
                #if operate[0] == 'S':
                #    mis_region_MD = [ele + number[0] + start_pos - 1 for ele in mis_region_MD]
                #else:
                mis_region_MD = [ele + start_pos - 1 for ele in mis_region_MD]

                mis_region_I = []   # insertion boundary region
                diff_len = 0        # len contribution of D and I
                if 'I' in operate or 'D' in operate:
                    idx_I = start_pos - 1 # index in reference
                    for idx, op in enumerate(operate):
                        if op == 'I':
                            diff_len -= number[idx]
                            mis_region_I.append(idx_I)
                            mis_region_I.append(idx_I+1)
                        else:
                            if op == 'S':
                                diff_len -= number[idx]
                            else:
                                idx_I += number[idx]
                                if op == 'D':
                                    diff_len += number[idx]
                
                #print(fields[0])
                #print(mis_region_MD)
                #print(mis_region_I)
                #print(mis_region)
                mis_region = mis_region_MD + mis_region_I
                mis_region.sort()
                

                edit_histogram[mis_region] += 1
                
                end_pos   = start_pos + len(fields[9]) + diff_len
                cov_histogram[start_pos:end_pos] += 1
                
                # record the reads information
                if int(sam_flag/16)%2 == 1:
                    dict_reads[(read_name, even_odd_flag)] = get_reverse_complement(read_SEQ.upper())
                else:
                    dict_reads[(read_name, even_odd_flag)] = read_SEQ
                list_read_info.append((start_pos, end_pos, read_name, even_odd_flag, mis_region, cigar, read_SEQ))
                if even_odd_flag == 1:
                    even_odd_flag = 2
                else:
                    even_odd_flag = 1

    contig_SEQ = ""
    with open(contig_file, 'r') as f_c:
        contig_flag = False
        for line in f_c:
            if line[0] == '>':
                tmp_name = line[1:].strip()
                if tmp_name == contig_name:
                    contig_flag = True
                else:
                    contig_flag = False
            elif contig_flag:
                contig_SEQ += line.strip()

    return edit_histogram, cov_histogram, list_read_info, dict_reads, contig_SEQ


def pop_perfect_reads(edit_region, list_read_info, dict_reads, fn_output_file):
    # pop perfect matches in the edit_region
    for read_info in list_read_info:
        start_pos, end_pos, read_name, even_odd_flag, mis_region, cigar, read_SEQ = read_info
        for spot in edit_region:
            # if match cover the edit spot
            if start_pos <= spot and end_pos >= spot:
                support_flag = False
                for ele in mis_region:
                    if ele in edit_region: # the read support some edit regions
                        support_flag = True
                        break
                # if the read and its pair is in dict_read
                if support_flag:
                    # the read support some edit regions
                    break
                else: # pop the reads
                    #if dict_reads.get((read_name, even_odd_flag)):
                    if dict_reads.get((read_name, 1)):
                        dict_reads.pop((read_name, 1))
                    if dict_reads.get((read_name, 2)):
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


def report_variant_base(start_pos, cigar, mis_region, read_SEQ):
    list_var_pair = []
    if len(mis_region) < 1:
        return []
    mis_idx = 0
    ref_cursor = start_pos
    read_cursor = 0

    number, operate = parse_CIGAR(cigar)
    for idx, op in enumerate(operate):
        if mis_idx >= len(mis_region): # if the CIGAR is out of mis_region, skip
            break
        # Not operate now ---> We just ignore the soft-clip and hard-clip reads for precision purpose
        if op == 'M':# or op == 'S':
            tmp_ref_cursor = ref_cursor + number[idx]
            while (mis_idx < len(mis_region) ) and (mis_region[mis_idx] < tmp_ref_cursor):
                dist = mis_region[mis_idx] - ref_cursor
                list_var_pair.append( (mis_region[mis_idx], read_SEQ[read_cursor+dist]) )
                mis_idx += 1
            ref_cursor  += number[idx]
            read_cursor += number[idx]
        elif op == 'D':
            num_D = number[idx]
            #list_var_pair += [ (D_idx,'-') for D_idx in range(mis_region[mis_idx], mis_region[mis_idx]+num_D) ]
            if (mis_idx < len(mis_region)) and (mis_region[mis_idx] <= ref_cursor):
                list_var_pair.append( (mis_region[mis_idx], 'D'+str(num_D)) )
                mis_idx  += num_D
            ref_cursor += num_D
        elif op == 'I':
            num_I = number[idx]
            if (mis_idx < len(mis_region)) and (mis_region[mis_idx] <= ref_cursor):
                list_var_pair.append( (mis_region[mis_idx], 'I'+read_SEQ[read_cursor:read_cursor+num_I]) )
                mis_idx + 1
            read_cursor += num_I
        elif op == 'S':
            num_S = number[idx]
            read_cursor += num_S
    return list_var_pair


def variant_link_graph(edit_region, list_read_info):
    # the function group the reads that cover at least two variant and link the variants for haplotyping
    #
    # dict_link_graph = {}      # record all the links on the interested position
    #  - keys: (position, base)
    #  - values: dict_link {}
    #             - keys: (node_idx, ((position_0, base_0), (position_1, base_1), ... (position_n, base_n)))
    #                     (position, base) is the (position_node_idx, base_node_idx)
    #             - values: weight
    #
    # dict_var_weight = {}      # record all the bases (not concerning the links) on the interested position
    #  - keys: position
    #  - values: dict_base_weight {}
    #            - keys: base
    #            - values: weight
    #
    # dict_link_outward = {}    # record all right side links with one shift relative to (position, base)
    #  - keys: (position, base)
    #  - values: dict_1_link {}
    #             - keys: ((position, base), (position_0, base_0))
    #             - values: weight
    #
    # dict_link_inward = {}     # record all left side links with one shift relative to (position, base)
    #  - keys: (position, base)
    #  - values: dict_1_link {}
    #             - keys: ((position_inward, base_inward), (position, base))
    #             - values: weight

    dict_link_graph = {}
    dict_var_weight = { pos:{} for pos in edit_region }
    dict_link_outward = {}
    dict_link_inward  = {}
   
    list_idx = -2
    len_list_read_info = len(list_read_info)-3
    #for list_idx in range(0,len(list_read_info),2):
    while list_idx < len_list_read_info:
        list_idx += 2
        start_pos_0, end_pos_0, read_name_0, even_odd_flag_0, mis_region_0, cigar_0, read_SEQ_0 = list_read_info[list_idx]
        start_pos_1, end_pos_1, read_name_1, even_odd_flag_1, mis_region_1, cigar_1, read_SEQ_1 = list_read_info[list_idx+1]
        if read_name_0 != read_name_1:
            eprint("Error with read name:", read_name_0, "and", read_name_1, list_idx)
            list_idx -= 1
            continue
        # we suppose that one of the read with soft-clip means that the read pairs should not be here
        '''
        if 'S' in cigar_0 or 'S' in cigar_1:
            if start_pos_0 <= 583 <= end_pos_0 or start_pos_1 <= 583 <= end_pos_1:
                print("Soft Clip!", read_name_0, start_pos_0, end_pos_0, "----", start_pos_1, end_pos_1)
            continue
        else:
            if start_pos_0 <= 583 <= end_pos_0:
                print("pos_0", start_pos_0, end_pos_0)
            if start_pos_1 <= 583 <= end_pos_1:
                print("pos_1", start_pos_1, end_pos_1)
        '''
        covered_edit_region_0 = []
        covered_edit_region_1 = []
        for spot in edit_region:
            # if match cover the edit spot
            if start_pos_0 <= spot <= end_pos_0:
                covered_edit_region_0.append(spot)
            if start_pos_1 <= spot <= end_pos_1:
                covered_edit_region_1.append(spot)

        # if the read pairs cover any interested spot, record the base information in dict_var_weight
        if len(covered_edit_region_0) + len(covered_edit_region_1) > 0:
            list_var_pair_0 = report_variant_base(start_pos_0, cigar_0, covered_edit_region_0, read_SEQ_0)
            list_var_pair_1 = report_variant_base(start_pos_1, cigar_1, covered_edit_region_1, read_SEQ_1)
            haplotype_frag = tuple(sorted(list_var_pair_0 + list_var_pair_1))
            # record all the variant proportions at every variant position
            for var_pair in haplotype_frag:
                if dict_var_weight[var_pair[0]].get(var_pair[1]):
                    dict_var_weight[var_pair[0]][var_pair[1]] += 1
                else:
                    dict_var_weight[var_pair[0]][var_pair[1]] = 1
            
            # if there are long range coverage between mate pairs, then record the links in dict_link_graph
            if len(haplotype_frag) > 1:
                for node_idx, var_pair in enumerate(haplotype_frag):
                    # dict_link_graph
                    if dict_link_graph.get(var_pair):
                        if dict_link_graph[var_pair].get((node_idx, haplotype_frag)):
                            dict_link_graph[var_pair][(node_idx, haplotype_frag)] += 1
                        else:
                            dict_link_graph[var_pair][(node_idx, haplotype_frag)] = 1
                    else:
                        dict_link_graph[var_pair] = {(node_idx, haplotype_frag):1}
                    # dict_link_outward
                    if node_idx < len(haplotype_frag)-1:
                        if dict_link_outward.get(var_pair):
                            if dict_link_outward[var_pair].get((var_pair, haplotype_frag[node_idx+1])):
                                dict_link_outward[var_pair][(var_pair, haplotype_frag[node_idx+1])] += 1
                            else:
                                dict_link_outward[var_pair][(var_pair, haplotype_frag[node_idx+1])] = 1
                        else:
                            dict_link_outward[var_pair] = {(var_pair, haplotype_frag[node_idx+1]):1}
                    # dict_link_inward
                    if node_idx > 0:
                        if dict_link_inward.get(var_pair):
                            if dict_link_inward[var_pair].get((haplotype_frag[node_idx-1], var_pair)):
                                dict_link_inward[var_pair][(haplotype_frag[node_idx-1], var_pair)] += 1
                            else:
                                dict_link_inward[var_pair][(haplotype_frag[node_idx-1], var_pair)] = 1
                        else:
                            dict_link_inward[var_pair] = {(haplotype_frag[node_idx-1], var_pair):1}

    return dict_link_graph, dict_var_weight, dict_link_outward, dict_link_inward


def find_double_pos(pos_start_idx, list_pos_weight, haplotype_0, haplotype_1, hap_cursor_0, hap_cursor_1):
    while pos_start_idx < len(list_pos_weight):
        if len(list_pos_weight[pos_start_idx][1]) == 0:
            print("There is no reads covered the position", list_pos_weight[0][0])
            pos_start_idx += 1
        elif len(list_pos_weight[pos_start_idx][1]) == 1:
            print("There is no variant detected in position", list_pos_weight[0][0])
            haplotype_0.append((list_pos_weight[pos_start_idx][0], list_pos_weight[pos_start_idx][1][0][0]))
            haplotype_1.append((list_pos_weight[pos_start_idx][0], list_pos_weight[pos_start_idx][1][0][0]))
            pos_start_idx += 1
            hap_cursor_0 += 1
            hap_cursor_1 += 1
        else:
            haplotype_0.append((list_pos_weight[pos_start_idx][0], list_pos_weight[pos_start_idx][1][0][0]))
            haplotype_1.append((list_pos_weight[pos_start_idx][0], list_pos_weight[pos_start_idx][1][1][0]))
            pos_start_idx += 1
            return pos_start_idx, haplotype_0, haplotype_1, hap_cursor_0, hap_cursor_1
    return pos_start_idx, haplotype_0, haplotype_1, hap_cursor_0, hap_cursor_1


def trim_dict(dict_target, thrsd=4):
    sorted_items = sorted(dict_target.items(), key=lambda pair:pair[1], reverse=True)
    for pair_id in range(1, len(sorted_items)):
        if sorted_items[ pair_id-1 ][1] > sorted_items[ pair_id ][1]*thrsd:
            for i in range(pair_id, len(sorted_items)):
                dict_target.pop(sorted_items[i][0])
            break


def get_farthest_ext(dict_link, hap_base_info):
    list_ext = [ext[1][0] for ext in dict_link.keys()]
    position = hap_base_info[0]
    subst_base = hap_base_info[1]
    if subst_base[0] == 'D':
        min_ext = position + int(subst_base[1:])
        list_ext.append(min_ext)
    elif subst_base[0] == 'H':
        min_ext = position + len(subst_base[1:])
        list_ext.append(min_ext)
    else:
        list_ext.append(position)
    return max(list_ext)


def haplotyping_link_graph(dict_link_graph, dict_var_weight, dict_link_outward, dict_link_inward, edit_region):
    # sort the potential variants on the interested site, can only use these variants bases
    list_pos_weight = []
    print("Trimming the significant bases at interested site:")
    print("Original site-base dict", dict_var_weight)
    '''
    for key in sorted(dict_var_weight.keys()):
        list_pos_base = list(dict_var_weight[key].items()) #[(pos_key,pos_value) for pos_key, pos_value in dict_var_weight[key].items()]
        list_pos_base.sort(key=lambda pair:pair[1], reverse=True)
        for idx in range(1, len(list_pos_base)):
            if list_pos_base[idx-1][1] > list_pos_base[idx][1]*10:
                list_pos_base = list_pos_base[:idx]
                break
        list_pos_weight.append((key,list_pos_base))
    '''
    for key in sorted(dict_var_weight.keys()):
        dict_part = dict_var_weight[key]
        trim_dict(dict_part, 10)
        list_pos_weight.append((key, sorted(dict_part.items(), key=lambda pair:pair[1], reverse=True)))

    print("Final site-base list:", list_pos_weight)
    if list_pos_weight == []:
        print("There is no variant detected!")
        return [], []
    
    print("+++++++++++++++++++", "dict_link_graph", "+++++++++++++++++++")
    for key in sorted(dict_link_graph.keys()):
        print(key, dict_link_graph[key])
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    
    # initializing the haplotype list, the cursor, and the last_ext
    haplotype_0 = []        # record the (position, base) pair of the haplotype
    hap_cursor_0 = 0        # record the position got the linking information (still useless in this version)
    break_flag_0 = False    # the flag indicating of the haplotype is breaked
    haplotype_1 = []
    hap_cursor_1 = 0
    break_flag_1 = False
    pos_start_idx = 0
   
    # find the first variant site with two variants
    pos_start_idx, haplotype_0, haplotype_1, hap_cursor_0, hap_cursor_1 = find_double_pos(pos_start_idx, list_pos_weight, haplotype_0, haplotype_1, hap_cursor_0, hap_cursor_1)

    # haplotyping from list_pos_weight:
    for pos_idx in range(pos_start_idx, len(list_pos_weight)):
        pos_weight = list_pos_weight[pos_idx]
        position = pos_weight[0]
        list_pos_base = pos_weight[1]
        print("XXXXXXXXXXXXXX", position, "XXXXXXXXXXXXXXXX")

        # deal with haplotype_0's outward link
        dict_outward_0 = {}
        if dict_link_outward.get(haplotype_0[hap_cursor_0]):
            dict_outward_0 = dict_link_outward[haplotype_0[hap_cursor_0]]
        trim_dict(dict_outward_0)
        if position > get_farthest_ext(dict_outward_0, haplotype_0[hap_cursor_0]):
            break_flag_0 = True
            eprint("Haplotype 0 has a break at", haplotype_0[hap_cursor_0], "to", position)
        print(dict_outward_0)
        # deal with haplotype_1's outward link
        print("--------------------")
        dict_outward_1 = {}
        if dict_link_outward.get(haplotype_1[hap_cursor_1]):
            dict_outward_1 = dict_link_outward[haplotype_1[hap_cursor_1]]
        trim_dict(dict_outward_1)
        if position > get_farthest_ext(dict_outward_1, haplotype_1[hap_cursor_1]):
            break_flag_1 = True
            eprint("Haplotype 1 has a break at", haplotype_1[hap_cursor_1], "to", position)
        print(dict_outward_1)
        # deal with position's inward link
        print("--------------------")
        dict_inward_0 = {}
        if dict_link_inward.get((position, list_pos_base[0][0])):
            dict_inward_0 = dict_link_inward[(position, list_pos_base[0][0])]
        trim_dict(dict_inward_0)
        print(dict_inward_0)
        if len(list_pos_base) > 1:
            print("--------------------")
            dict_inward_1 = {}
            if dict_link_inward.get((position, list_pos_base[1][0])):
                dict_inward_1 = dict_link_inward[(position, list_pos_base[1][0])]
            trim_dict(dict_inward_1)
            print(dict_inward_1)

        connect_info_0 = None
        connect_info_1 = None
        # There must be at least one kind of base in the position
        for (outward_key, weight) in sorted(dict_outward_0.items(), key=lambda pair:pair[1], reverse=True):
            if dict_inward_0.get(outward_key):
                print("Potential Connect: ", outward_key, 0, 0)
                connect_info_0 = (dict_outward_0[outward_key], (position, outward_key[1][1]))
                break
        for (outward_key, weight) in sorted(dict_outward_1.items(), key=lambda pair:pair[1], reverse=True):
            if dict_inward_0.get(outward_key):
                print("Potential Connect: ", outward_key, 1, 0)
                connect_info_1 = (dict_outward_1[outward_key], (position, outward_key[1][1]))
                break
        
        # if there are two variants in the position
        if len(list_pos_base) > 1:
            for (outward_key, weight) in sorted(dict_outward_0.items(), key=lambda pair:pair[1], reverse=True):
                if dict_inward_1.get(outward_key):
                    print("Potential Connect: ", outward_key, 0, 1)
                    if connect_info_0 == None or connect_info_0[0] < weight:
                        connect_info_0 = (dict_outward_0[outward_key], (position, outward_key[1][1]))
                        break
            for (outward_key, weight) in sorted(dict_outward_1.items(), key=lambda pair:pair[1], reverse=True):
                if dict_inward_1.get(outward_key):
                    print("Potential Connect: ", outward_key, 1, 1)
                    if connect_info_1 == None or connect_info_1[0] < weight:
                        connect_info_1 = (dict_outward_1[outward_key], (position, outward_key[1][1]))
                        break
        
        # the case that two haplotypes may collapse into one
        if connect_info_0 and connect_info_1:
            if connect_info_0[1] == connect_info_1[1]: # two haplotypes are collapsed
                for redouble_idx in range(pos_idx, len(list_pos_weight)): # make it double again
                    rd_pos_weight = list_pos_weight[redouble_idx]
                    rd_position = rd_pos_weight[0]
                    rd_list_pos_base = rd_pos_weight[1]
                    if(len(rd_list_pos_base)) >= 2: # if there are two variants at the site
                        # call the potential connections
                        last_info_0 = haplotype_0[hap_cursor_0]
                        last_info_1 = haplotype_1[hap_cursor_1]
                        dict_info_0 = dict_link_graph[last_info_0]
                        dict_info_1 = dict_link_graph[last_info_1]
                        # connect them
                        rd_info_0 = None
                        rd_info_1 = None
                        for rd_link_info, rd_weight in sorted(dict_info_0.items(), key=lambda pair:pair[1], reverse=True):
                            variant_flag = False
                            for info_pair in rd_link_info[1]:
                                tmp_rd_info = []
                                if info_pair == connect_info_0[1]:
                                    variant_flag = True
                                    tmp_rd_info.append(info_pair)
                                if variant_flag:
                                    tmp_rd_info.append(info_pair)
                                    if info_pair[0] == rd_position:
                                        rd_info_0 = tmp_rd_info
                                        break
                            if rd_info_0:
                                break
                            
                        for rd_link_info, rd_weight in sorted(dict_info_1.items(), key=lambda pair:pair[1], reverse=True):
                            for info_pair in rd_link_info[1]:
                                tmp_rd_info = []
                                if info_pair == connect_info_1[1]:
                                    variant_flag = True
                                    tmp_rd_info.append(info_pair)
                                if variant_flag:
                                    tmp_rd_info.append(info_pair)
                                    if info_pair[0] == rd_position:
                                        rd_info_1 = tmp_rd_info
                                        break
                            if rd_info_1:
                                break

                        if rd_info_0 != rd_info_1:
                            if rd_info_0:
                                haplotype_0 += rd_info_0
                                hap_cursor_0 += len(rd_info_0)
                            else:
                                break_flag_0 = True
                            if rd_info_1:
                                haplotype_1 += rd_info_1
                                hap_cursor_1 += len(rd_info_1)
                            else:
                                break_flag_1 = True
                            break
                         
                print("Crossing the single base variant site...")
                continue



        # update the nodes if the connection is found
        if connect_info_0:    
            haplotype_0.append(connect_info_0[1])
            hap_cursor_0 += 1
            if break_flag_1 and len(list_pos_base) >1:
                for idx in range(2):
                    potential_base = list_pos_base[idx][0] 
                    if potential_base != connect_info_0[1][1]:
                        eprint("Link rebuilt on Haplotype 1 at", haplotype_1[hap_cursor_1] , "to", position)
                        haplotype_1.append((position, potential_base))
                        hap_cursor_1 += 1
                        break_flag_1 = False
                        break
        if connect_info_1:    
            haplotype_1.append(connect_info_1[1])
            hap_cursor_1 += 1
            if break_flag_0 and len(list_pos_base) >1:
                for idx in range(2):
                    potential_base = list_pos_base[idx][0] 
                    if potential_base != connect_info_1[1][1]:
                        eprint("Link rebuilt on Haplotype 0 at", haplotype_0[hap_cursor_0] , "to", position)
                        haplotype_0.append((position, potential_base))
                        hap_cursor_0 += 1
                        break_flag_0 = False
                        break

        if break_flag_0 and break_flag_1:
            eprint("BREAKING LINKS FOR BOTH HAPLOTYPE AT", position, "!!!!")
            eprint("Breaking links cannot resloved, we guess...")
            haplotype_0.append((position, list_pos_base[0][0]))
            hap_cursor_0 += 1
            if len(list_pos_base) > 1:
                haplotype_1.append((position, list_pos_base[1][0]))
                hap_cursor_1 += 1
    
    print(haplotype_0)
    print(haplotype_1)
    return haplotype_0, haplotype_1


def sequence_substitution(sequence, list_correct_pairs):
    hap_offset = 0
    for pair in list_correct_pairs:
        position = pair[0]
        operation = pair[1]
        if operation[0] == 'D':
            D_len = int(operation[1:])
            sequence = sequence[:position+hap_offset] + sequence[position+hap_offset+D_len:]
            hap_offset -= D_len
        elif operation[0] == 'I':
            I_len = len(operation[1:])
            print("I situation", sequence[:position+hap_offset-1], operation, sequence[position+hap_offset-1:])
            sequence = sequence[:position+hap_offset-1] + operation[1:] + sequence[position+hap_offset-1:]
            hap_offset += I_len
        else:
            sequence = sequence[:position+hap_offset-1] + operation + sequence[position+hap_offset:]
    return sequence


def get_allele_name(allele_file):
    with open(allele_file, 'r') as f_a:
        for line in f_a:
            if line[0] == '>':
                allele_name = line[1:].strip()
                return allele_name


def output_contig_correction(contig_SEQ, region_st, region_ed, haplotype_0, haplotype_1, allele_name, corrected_contig_output_file, suffix="/haplotype-"):
    corrected_contig_SEQ_0 = sequence_substitution(contig_SEQ, haplotype_0)
    corrected_contig_SEQ_0 = corrected_contig_SEQ_0[region_st:region_ed]

    corrected_contig_SEQ_1 = sequence_substitution(contig_SEQ, haplotype_1)
    corrected_contig_SEQ_1 = corrected_contig_SEQ_1[region_st:region_ed]

    f_c = open(corrected_contig_output_file, 'a')
    f_c.write(">" + allele_name + suffix + "0\n")
    f_c.write(corrected_contig_SEQ_0 + "\n")
    f_c.write(">" + allele_name + suffix + "1\n")
    f_c.write(corrected_contig_SEQ_1 + "\n")
    f_c.close()
    return (corrected_contig_SEQ_0, corrected_contig_SEQ_1)


def output_original_contig(contig_SEQ, region_st, region_ed, allele_file, corrected_contig_output_file):
    allele_name = ""
    with open(allele_file, 'r') as f_a:
        for line in f_a:
            if line[0] == '>':
                allele_name = line[1:].strip()
                break

    f_c = open(corrected_contig_output_file, 'a')
    f_c.write(">" + allele_name + "/original\n")
    f_c.write(contig_SEQ[region_st:region_ed] + "\n")
    f_c.close()
    


if __name__ == '__main__':
    args = parse_args()
    fn_sam = args.fn_sam
    thrsd = args.thrsd
    interest_region = args.interest_region
    fn_output_file = args.fn_output_file
    # parameter for corrected contig
    contig_file = args.contig_file
    allele_file = args.allele_file
    corrected_contig_output_file = args.corrected_contig_output_file

    #parse the sam file and generate
    edit_histogram, cov_histogram, list_read_info, dict_reads, contig_SEQ = mark_edit_region(fn_sam, fn_output_file, contig_file)
    
    #determine the region contains alternative flanking region
    edit_region = []
    for idx, ele in enumerate(edit_histogram):
        print(str(idx) + ':\t' + str(cov_histogram[idx])  + '\t' + str(ele))
        if ele > cov_histogram[idx]/4:
            edit_region.append(idx)

    if interest_region:
        print("interested_region:")
        print(interest_region)
        print("edit_region:")
        print(edit_region)
        region_st = int(interest_region.split('-')[0])
        region_ed = int(interest_region.split('-')[1])
        break_flag = False
        interest_edit_region = []
        for ele in edit_region:
            if region_st <= ele <= region_ed:
                interest_edit_region.append(ele)
        if len(interest_edit_region) > 0:
            pop_perfect_reads(edit_region, list_read_info, dict_reads, fn_output_file)
                
            print("=========== start haplotyping ==============")
            dict_link_graph, dict_var_weight, dict_link_outward, dict_link_inward = variant_link_graph(interest_edit_region, list_read_info)
            haplotype_0, haplotype_1 = haplotyping_link_graph(dict_link_graph, dict_var_weight, dict_link_outward, dict_link_inward, interest_region)
            output_contig_correction(contig_SEQ, region_st, region_ed, haplotype_0, haplotype_1, get_allele_name(allele_file), corrected_contig_output_file)
        else:
            output_original_contig(contig_SEQ, region_st, region_ed, allele_file, corrected_contig_output_file)
            print("No variant detected in the interested region, output original contig...")
    else:
        print("edit_region:")
        print(edit_region)
        if len(edit_region) > 0:
            pop_perfect_reads(edit_region, list_read_info, dict_reads, fn_output_file)
                
            print("=========== start haplotyping ==============")
            dict_link_graph, dict_var_weight, dict_link_outward, dict_link_inward = variant_link_graph(edit_region, list_read_info)
            haplotyping_link_graph(dict_link_graph, dict_var_weight, dict_link_outward, dict_link_inward, edit_region)
            output_contig_correction(contig_SEQ, region_st, region_ed, haplotype_0, haplotype_1, get_allele_name(allele_file), corrected_contig_output_file)
        else:
            print("No variant detected!")

    

