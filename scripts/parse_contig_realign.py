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
from utils import get_reverse_complement
from coverage_analysis import eprint
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
                    list_read_info.append((0, 0, read_name, even_odd_flag, []))
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
                                diff_len += number[idx]
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
                
                '''
                mis_base = []
                mis_idx = 0
                ref_cursor = start_pos
                read_cursor = 0
                for idx, op in enumerate(operate):
                    print(mis_base)
                    print(ref_cursor, '---', read_cursor)
                    if op == 'M' or op == 'S':
                        tmp_ref_cursor = ref_cursor + number[idx]
                        while (mis_idx < len(mis_region) ) and (mis_region[mis_idx] < tmp_ref_cursor):
                            dist = mis_region[mis_idx] - ref_cursor
                            mis_base.append(read_SEQ[read_cursor+dist])
                            mis_idx += 1
                        ref_cursor  += number[idx]
                        read_cursor += number[idx]
                    elif op == 'D':
                        num_D = number[idx]
                        mis_idx  += num_D
                        mis_base += ['-']*num_D
                        ref_cursor += num_D
                    elif op == 'I':
                        num_I = number[idx]
                        mis_idx + 1
                        mis_base.append('I')
                        mis_base.append(read_SEQ[read_cursor:read_cursor+num_I])
                        read_cursor += num_I
                print(mis_base)
                '''


                edit_histogram[mis_region] += 1
                #edit_histogram[mis_region_MD] += 1
                #edit_histogram[mis_region_I]  += 1
                
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

    return edit_histogram, cov_histogram, list_read_info, dict_reads


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
                elif dict_reads.get((read_name, even_odd_flag)):
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


def report_variant_base(start_pos, cigar, mis_region, read_SEQ):
    #list_var_pair = []
    if len(mis_region) < 1:
        return []
    mis_base = []
    mis_idx = 0
    ref_cursor = start_pos
    read_cursor = 0

    number, operate = parse_CIGAR(cigar)
    for idx, op in enumerate(operate):
        if mis_idx >= len(mis_region):
            break
        if op == 'M':# or op == 'S' or op == 'H':
            tmp_ref_cursor = ref_cursor + number[idx]
            while (mis_idx < len(mis_region) ) and (mis_region[mis_idx] < tmp_ref_cursor):
                dist = mis_region[mis_idx] - ref_cursor
                mis_base.append( (mis_region[mis_idx], read_SEQ[read_cursor+dist]) )
                mis_idx += 1
            ref_cursor  += number[idx]
            read_cursor += number[idx]
        elif op == 'D':
            num_D = number[idx]
            #mis_base += [ (D_idx,'-') for D_idx in range(mis_region[mis_idx], mis_region[mis_idx]+num_D) ]
            mis_base.append( (mis_region[mis_idx], 'D'+str(num_D)) )
            mis_idx  += num_D
            ref_cursor += num_D
        elif op == 'I':
            num_I = number[idx]
            mis_base.append( (mis_region[mis_idx], 'I'+read_SEQ[read_cursor:read_cursor+num_I]) )
            mis_idx + 1
            read_cursor += num_I
    return mis_base
    #return list_var_pair


def variant_link_graph(edit_region, list_read_info):
    # the function group the reads that cover at least two variant and link the variants for haplotyping
    # dict_link_graph = {}
    #  - keys: (position, base)
    #  - values: dict_link {}
    #             - keys: ((position_0, base_0), (position_1, base_1), ... (position_n, base_n))
    #                       position_0 is not necessary position, but position is one of position_x
    #             - values: weight
    #
    # dict_var_weight = {}
    #  - keys: position
    #  - values: dict_base_weight {}
    #            - keys: base
    #            - values: weight

    dict_link_graph = {}
    dict_var_weight = { pos:{} for pos in edit_region }
    
    for list_idx in range(0,len(list_read_info),2):
        read_info_0 = list_read_info[list_idx]
        read_info_1 = list_read_info[list_idx+1]
        start_pos_0, end_pos_0, read_name_0, even_odd_flag_0, mis_region_0, cigar_0, read_SEQ_0 = read_info_0
        start_pos_1, end_pos_1, read_name_1, even_odd_flag_1, mis_region_1, cigar_1, read_SEQ_1 = read_info_1
        if read_name_0 != read_name_1:
            eprint("Error with read name:", read_name_0, "and", read_name_1)
            continue
        if 'S' in cigar_0 or 'S' in cigar_1:
            continue

        covered_edit_region_0 = []
        covered_edit_region_1 = []
        for spot in edit_region:
            # if match cover the edit spot
            if start_pos_0 <= spot and end_pos_0 >= spot:
                covered_edit_region_0.append(spot)
            elif start_pos_1 <= spot and end_pos_1 >= spot:
                covered_edit_region_1.append(spot)

        #if len(covered_edit_region_0) > 0 and len(covered_edit_region_1) > 0:
        if len(covered_edit_region_0) + len(covered_edit_region_1) > 1:
            # if there are long range coverage between mate pairs
            list_var_pair_0 = report_variant_base(start_pos_0, cigar_0, covered_edit_region_0, read_SEQ_0)
            list_var_pair_1 = report_variant_base(start_pos_1, cigar_1, covered_edit_region_1, read_SEQ_1)
            #print("================================")
            #print(read_name_0)
            #print(covered_edit_region_0)
            #print(list_var_pair_0)
            #print(covered_edit_region_1)
            #print(list_var_pair_1)
            haplotype_frag = tuple(sorted(list_var_pair_0 + list_var_pair_1))
            # record all the variant proportions at every variant position
            for var_pair in haplotype_frag:
                if dict_var_weight[var_pair[0]].get(var_pair[1]):
                    dict_var_weight[var_pair[0]][var_pair[1]] += 1
                else:
                    dict_var_weight[var_pair[0]][var_pair[1]] = 1
            
            # record all the links
            for node_idx, var_pair in enumerate(haplotype_frag):
                if dict_link_graph.get(var_pair):
                    if dict_link_graph[var_pair].get((node_idx, haplotype_frag)):
                        dict_link_graph[var_pair][(node_idx, haplotype_frag)] += 1
                    else:
                        dict_link_graph[var_pair][(node_idx, haplotype_frag)] = 1
                else:
                    dict_link_graph[var_pair] = {(node_idx, haplotype_frag):1}
    #for key in sorted(dict_link_graph.keys()):
    #    print (dict_link_graph[key])
    return dict_link_graph, dict_var_weight


def trim_dict(dict_target):
    sorted_items = sorted(dict_target.items(), key=lambda pair:pair[1], reverse=True)
    for pair_id in range(1, len(sorted_items)):
        if sorted_items[ pair_id-1 ] > sorted_items[ pair_id ]*4:
            for i in range(pair_id, len(sorted_items)):
                dict_target.pop(sorted_items[i][0])
            break


def haplotyping_link_graph(dict_link_graph, dict_var_weight, edit_region):
    list_pos_weight = []
    for key in sorted(dict_var_weight.keys()):
        #print("+++++++++++++++++++")
        list_pos_base = [(pos_key,pos_value) for pos_key, pos_value in dict_var_weight[key].items()]
        list_pos_base.sort(key=lambda pair:pair[1], reverse=True)
        print(list_pos_base)
        for idx in range(len(list_pos_base)-1):
            if list_pos_base[idx][1] > list_pos_base[idx+1][1]*10:
                list_pos_base = list_pos_base[:idx+1]
                break
        list_pos_weight.append((key,list_pos_base))
        print(list_pos_base)
    print(list_pos_weight)
    print("+++++++++++++++++++")

    # initializing the haplotype list and the cursor
    haplotype_0 = []
    hap_cursor_0 = 0
    haplotype_1 = []
    hap_cursor_1 = 0
    if len(list_pos_weight[0]) < 2:
        eprint("ERROR ON START OF THE HAPLOTYPING!")
    else: # len(list_pos_weight[0]) == 2:
        haplotype_0.append((list_pos_weight[0][0], list_pos_weight[0][1][0][0]))
        haplotype_1.append((list_pos_weight[0][0], list_pos_weight[0][1][1][0]))

    # haplotyping from list_pos_weight:
    for pos_idx in range(1, len(list_pos_weight)):
        pos_weight = list_pos_weight[pos_idx]
        position = pos_weight[0]
        list_pos_base = pos_weight[1]
        print("XXXXXXXXXXXXXX", position, "XXXXXXXXXXXXXXXX")
        # deal with haplotype_0's outward link
        dict_connection = dict_link_graph[haplotype_0[hap_cursor_0]]
        dict_outward_0  = {}
        for key in dict_connection.keys():
            node_idx = key[0]
            path     = key[1]
            if len(path) > node_idx+1:
                if dict_outward_0.get((path[node_idx], path[node_idx+1])):
                    dict_outward_0[(path[node_idx], path[node_idx+1])] += dict_connection[key]
                else:
                    dict_outward_0[(path[node_idx], path[node_idx+1])] = dict_connection[key]
                #print(path[node_idx], path[node_idx+1], dict_connection[key])
        
        trim_dict(dict_outward_0)
        print(dict_outward_0)
        # deal with haplotype_1's outward link
        print("----------------")
        dict_connection = dict_link_graph[haplotype_1[hap_cursor_1]]
        dict_outward_1  = {}
        for key in dict_connection.keys():
            node_idx = key[0]
            path     = key[1]
            if len(path) > node_idx+1:
                if dict_outward_1.get((path[node_idx], path[node_idx+1])):
                    dict_outward_1[(path[node_idx], path[node_idx+1])] += dict_connection[key]
                else:
                    dict_outward_1[(path[node_idx], path[node_idx+1])] = dict_connection[key]
                #print(path[node_idx], path[node_idx+1], dict_connection[key])
        
        trim_dict(dict_outward_1)
        print(dict_outward_1)
        # deal with position's inward link
        print("----------------")
        dict_connection = dict_link_graph[(position, list_pos_base[0][0])]
        dict_inward_0   = {}
        for key in dict_connection.keys():
            node_idx = key[0]
            path     = key[1]
            if node_idx > 0:
                if dict_inward_0.get((path[node_idx-1], path[node_idx])):
                    dict_inward_0[(path[node_idx-1], path[node_idx])] += dict_connection[key]
                else:
                    dict_inward_0[(path[node_idx-1], path[node_idx])] = dict_connection[key]
                #print(path[node_idx-1], path[node_idx], dict_connection[key])

        print(dict_inward_0)
        trim_dict(dict_inward_0)
        print(dict_inward_0)
        #print(dict_link_graph[(position, list_pos_base[1][0])])
        if len(list_pos_base) > 1:
            print("----------------")
            dict_connection = dict_link_graph[(position, list_pos_base[1][0])]
            dict_inward_1   = {}
            for key in dict_connection.keys():
                node_idx = key[0]
                path     = key[1]
                if node_idx > 0:
                    if dict_inward_1.get((path[node_idx-1], path[node_idx])):
                        dict_inward_1[(path[node_idx-1], path[node_idx])] += dict_connection[key]
                    else:
                        dict_inward_1[(path[node_idx-1], path[node_idx])] = dict_connection[key]
                    #print(path[node_idx-1], path[node_idx], dict_connection[key])
            trim_dict(dict_inward_1)
            print(dict_inward_1)

        for outward_key in sorted(dict_outward_0.keys()):
            for inward_key in sorted(dict_inward_0.keys()):
                if outward_key == inward_key:
                    print("Connect: ", outward_key, 0, 0)
                    haplotype_0.append((position, outward_key[1][1]))
                    hap_cursor_0 += 1
                    break
        for outward_key in sorted(dict_outward_1.keys()):
            for inward_key in sorted(dict_inward_0.keys()):
                if outward_key == inward_key:
                    print("Connect: ", outward_key, 1, 0)
                    haplotype_1.append((position, outward_key[1][1]))
                    hap_cursor_1 += 1
                    break
        
        if len(list_pos_base) > 1:
            for outward_key in sorted(dict_outward_0.keys()):
                for inward_key in sorted(dict_inward_1.keys()):
                    if outward_key == inward_key:
                        print("Connect: ", outward_key, 0, 1)
                        haplotype_0.append((position, outward_key[1][1]))
                        hap_cursor_0 += 1
                        break
            for outward_key in sorted(dict_outward_1.keys()):
                for inward_key in sorted(dict_inward_1.keys()):
                    if outward_key == inward_key:
                        print("Connect: ", outward_key, 1, 1)
                        haplotype_1.append((position, outward_key[1][1]))
                        hap_cursor_1 += 1
                        break
    
    print(haplotype_0)
    print(haplotype_1)

    '''
    sorted_key_list = sorted(dict_link_graph.keys())
    list_var_weight = []
    ref_cursor = sorted_key_list[0][0]
    for key in sorted_key_list:
        print(key, dict_link_graph[key])
        
        if key[0] != ref_cursor:
            # update ref_cursor
            ref_cursor = key[0]
            print ("---------------------")
            if len(list_var_weight) == 1:
                print ("only one variant")
                print (list_var_weight[0])
            elif len(list_var_weight) == 2:
                print ("Normal Case")
                print (list_var_weight[0])
                print (list_var_weight[1])
            else: # multiple variants, take the largest two
                list_var_weight.sort(reverse=True)
                # evaluate which two? variant are the most common
                print (list_var_weight)
                print (list_var_weight[0])
                print (list_var_weight[1])
            
            list_var_weight = []
            print ("---------------------")
        '''
        #list_var_weight.append( (sum(dict_link_graph[key].values()), key[1]) )
        #print( ref_cursor, key[1], sum(dict_link_graph[key].values()) )
        #print( key, sum(dict_link_graph[key].values()), dict_link_graph[key] )



if __name__ == '__main__':
    args = parse_args()
    fn_sam = args.fn_sam
    thrsd = args.thrsd
    interest_region = args.interest_region
    fn_output_file = args.fn_output_file

    #parse the sam file and generate
    edit_histogram, cov_histogram, list_read_info, dict_reads = mark_edit_region(fn_sam, fn_output_file)
    
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
            dict_link_graph, dict_var_weight = variant_link_graph(interest_edit_region, list_read_info)
            haplotyping_link_graph(dict_link_graph, dict_var_weight, interest_region)
        else:
            print("No variant detected in the interested region!")
    else:
        print("edit_region:")
        print(edit_region)
        if len(edit_region) > 0:
            pop_perfect_reads(edit_region, list_read_info, dict_reads, fn_output_file)
            dict_link_graph, dict_var_weight = variant_link_graph(edit_region, list_read_info)
            haplotyping_link_graph(dict_link_graph, dict_var_weight, edit_region)
        else:
            print("No variant detected!")

    

